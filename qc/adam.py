"""Reference ADaM derivations (ADSL, ADAE, ADTTE-PFS) from the SDTM-like tables, one named rule per step.

alt: a set of known alternative interpretations; switching one on reproduces a specific, common derivation mistake.
  trtsdt_zero        TRTSDT = first EX record even if the dose was 0 (spec: first dose > 0)
  saffl_any          SAFFL = Y for any EX record (spec: received any dose > 0)
  agegr_gt           AGEGR1 '>=65' only for AGE > 65 (spec: AGE >= 65)
  impute_first       partial AE start always imputed to the 1st of month / Jan 1 (spec: TRTSDT if same month/year)
  day0               ASTDY = ASTDT - TRTSDT (spec: no day 0, +1 on/after TRTSDT)
  teae_window_other  TEAE window off by 2 days from the spec
  ignore_newtx       PFS ignores new anticancer therapy (spec: censor at last adequate assessment before it)
  ignore_missed      PFS ignores >= 2 missed assessments before the event (spec: censor at last adequate assessment)
  aval_noplus1       AVAL = ADT - STARTDT (spec: + 1)
"""
import datetime as dt
from .sdtm import D, iso

ALTS = ["trtsdt_zero", "saffl_any", "agegr_gt", "impute_first", "day0", "teae_window_other", "ignore_newtx", "ignore_missed",
        "aval_noplus1"]


def by_subj(rows):
    out = {}
    for x in rows:
        out.setdefault(x["USUBJID"], []).append(x)
    return out


def derive(T, spec, alt=frozenset(), extra=frozenset()):
    """spec: dict(teae_window, visit_weeks, window). extra: bug kinds outside the diagnosis library (held-out only)."""
    EX, RS, DS, CM, AE = (by_subj(T[k]) for k in ("EX", "RS", "DS", "CM", "AE"))
    adsl, adae, adtte = {}, {}, {}
    for dm in T["DM"]:
        u = dm["USUBJID"]
        ex = sorted(EX.get(u, []), key=lambda e: e["EXSTDTC"])
        dosed = [e for e in ex if e["EXDOSE"] > 0]
        trts = (ex[0] if "trtsdt_zero" in alt else dosed[0])["EXSTDTC"] if dosed else None
        trte = (ex[-1] if "trtedt_any" in extra else dosed[-1])["EXSTDTC"] if dosed else None
        ag = dm["AGE"] > 65 if "agegr_gt" in alt else dm["AGE"] >= 65
        dth = DS[u][0]["DSSTDTC"] if u in DS else None
        adsl[u] = dict(TRTSDT=trts, TRTEDT=trte, SAFFL="Y" if (ex if "saffl_any" in alt else dosed) else "N",
                       AGEGR1=">=65" if ag else "<65", RANDDT=dm["RANDDT"], DTHDT=dth)
        # ADAE
        win = spec["teae_window"] + (2 if "teae_window_other" in alt else 0)
        for ae in AE.get(u, []):
            s = ae["AESTDTC"]
            if len(s) == 10:
                ast = s
            elif len(s) == 7:
                ast = trts if (trts and trts[:7] == s and "impute_first" not in alt) else s + "-01"
            else:
                ast = trts if (trts and trts[:4] == s and "impute_first" not in alt) else s + "-01-01"
            dy = None
            if trts:
                diff = (D(ast) - D(trts)).days
                dy = diff + (1 if diff >= 0 and "day0" not in alt else 0)
            te = "Y" if trts and D(trts) <= D(ast) <= D(trte) + dt.timedelta(days=win) else None
            adae[(u, ae["AESEQ"])] = dict(ASTDT=ast, ASTDY=dy, TRTEMFL=te)
        # ADTTE: PFS from randomisation
        rand = D(dm["RANDDT"])
        post = sorted((x for x in RS.get(u, []) if x["RSORRES"] != "BASELINE" and D(x["RSDTC"]) > rand), key=lambda x: x["RSDTC"])
        adequate = [x for x in post if x["RSORRES"] != "NE"]
        pd = next((D(x["RSDTC"]) for x in adequate if x["RSORRES"] == "PD"), None)
        death = D(dth) if dth else None
        if death and "death_censored" in extra:
            death = None
        newtx = D(CM[u][0]["CMSTDTC"]) if u in CM else None
        gap2 = dt.timedelta(days=2 * 7 * spec["visit_weeks"] + spec["window"])
        event_dt = min(x for x in (pd, death) if x) if (pd or death) else None
        before = lambda lim: [D(x["RSDTC"]) for x in adequate if D(x["RSDTC"]) <= lim and x["RSORRES"] != "PD"]
        if not adequate and not (death and death - rand <= gap2):
            adt, cnsr, why = rand, 1, "NO POST-BASELINE ASSESSMENT"
        elif newtx and (event_dt is None or newtx < event_dt) and "ignore_newtx" not in alt:
            last = before(newtx)
            adt, cnsr, why = (last[-1] if last else rand), 1, "NEW ANTICANCER THERAPY"
        elif event_dt:
            last = before(event_dt)
            prev = last[-1] if last else rand
            if event_dt - prev > gap2 and "ignore_missed" not in alt:
                adt, cnsr, why = prev, 1, ">=2 MISSED ASSESSMENTS BEFORE EVENT"
            else:
                adt, cnsr, why = event_dt, 0, "PROGRESSION" if pd and event_dt == pd else "DEATH"
        else:
            adt, cnsr, why = (D(adequate[-1]["RSDTC"]) if adequate else rand), 1, "LAST ADEQUATE ASSESSMENT"
        aval = (adt - rand).days + (0 if "aval_noplus1" in alt else 1)
        adtte[u] = dict(PARAMCD="PFS", STARTDT=iso(rand), ADT=iso(adt), CNSR=cnsr, AVAL=aval, EVNTDESC=why)
    return dict(ADSL=adsl, ADAE=adae, ADTTE=adtte)
