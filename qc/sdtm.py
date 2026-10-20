"""Synthetic oncology study in SDTM-like tables (DM, EX, RS, DS, CM, AE). Dates are ISO strings; AE dates may be partial."""
import datetime as dt
import random

DEV = dict(n=120, start="2025-01-06", visit_weeks=8, window=7, miss_p=0.08, ne_p=0.03, newtx_p=0.15, death_p=0.35,
           zero_first_p=0.12, ae_partial_p=0.25, cutoff_days=760, cycle_days=21)


def D(s):
    return dt.date.fromisoformat(s)


def iso(d):
    return d.isoformat()


def make(seed, cfg=None):
    c = {**DEV, **(cfg or {})}
    r = random.Random(seed)
    start, cutoff = D(c["start"]), D(c["start"]) + dt.timedelta(days=c["cutoff_days"])
    T = dict(DM=[], EX=[], RS=[], DS=[], CM=[], AE=[], cutoff=iso(cutoff), cfg=c)
    for i in range(c["n"]):
        u = f"AST-{i + 1:04d}"
        rand = start + dt.timedelta(days=r.randrange(0, 200))
        age = r.choice([r.randrange(28, 85), 65, 64, 66])
        first = rand + dt.timedelta(days=r.randrange(0, 8))
        T["DM"].append(dict(USUBJID=u, AGE=age, SEX=r.choice("MF"), ARM=r.choice(["DRUG", "SOC"]), RANDDT=iso(rand),
                            RFSTDTC=iso(first)))
        pd_day = int(r.expovariate(1 / 260)) + 20
        pd_true = rand + dt.timedelta(days=pd_day)
        newtx = rand + dt.timedelta(days=r.randrange(30, max(31, pd_day))) if r.random() < c["newtx_p"] else None
        death = None
        if r.random() < c["death_p"]:
            death = (pd_true if r.random() < 0.7 else rand) + dt.timedelta(days=r.randrange(10, 260))
        end_tx = min(x for x in (pd_true, newtx or cutoff, death or cutoff, cutoff))
        # dosing: 21-day cycles; sometimes a zero-dose first record (dose held at cycle 1)
        d, k = first, 0
        while d < end_tx and d < cutoff:
            dose = 0 if (k == 0 and r.random() < c["zero_first_p"]) else 200
            T["EX"].append(dict(USUBJID=u, EXSTDTC=iso(d), EXDOSE=dose))
            d += dt.timedelta(days=c["cycle_days"]); k += 1
        if r.random() < 0.04:                          # randomised, dose held at cycle 1, never treated
            T["EX"] = [e for e in T["EX"] if e["USUBJID"] != u] + [dict(USUBJID=u, EXSTDTC=iso(first), EXDOSE=0)]
        elif not any(e["USUBJID"] == u and e["EXDOSE"] > 0 for e in T["EX"]):
            T["EX"].append(dict(USUBJID=u, EXSTDTC=iso(first + dt.timedelta(days=c["cycle_days"])), EXDOSE=200))
        # tumour assessments: baseline, then every visit_weeks +- window until PD / death / cutoff
        T["RS"].append(dict(USUBJID=u, RSDTC=iso(rand - dt.timedelta(days=r.randrange(1, 14))), RSORRES="BASELINE"))
        j = 1
        while True:
            a = rand + dt.timedelta(days=7 * c["visit_weeks"] * j + r.randrange(-c["window"], c["window"] + 1))
            j += 1
            if a > cutoff or (death and a >= death):
                break
            if r.random() < c["miss_p"]:
                continue
            if a >= pd_true:
                T["RS"].append(dict(USUBJID=u, RSDTC=iso(a), RSORRES="PD"))
                break
            T["RS"].append(dict(USUBJID=u, RSDTC=iso(a), RSORRES="NE" if r.random() < c["ne_p"] else r.choice(["SD", "SD", "PR", "CR"])))
        if death and death <= cutoff:
            T["DS"].append(dict(USUBJID=u, DSDECOD="DEATH", DSSTDTC=iso(death)))
        if newtx and newtx <= cutoff:
            T["CM"].append(dict(USUBJID=u, CMCAT="NEW ANTICANCER THERAPY", CMSTDTC=iso(newtx)))
        for _ in range(r.randrange(0, 5)):
            ad = first + dt.timedelta(days=r.randrange(-25, 240))
            if ad > cutoff:
                continue
            p = r.random()
            s = iso(ad)[:7] if p < c["ae_partial_p"] else (iso(ad)[:4] if p < c["ae_partial_p"] * 1.15 else iso(ad))
            T["AE"].append(dict(USUBJID=u, AESEQ=len([x for x in T["AE"] if x["USUBJID"] == u]) + 1, AESTDTC=s,
                                AEDECOD=r.choice(["NAUSEA", "FATIGUE", "ANAEMIA", "RASH", "DIARRHOEA", "NEUTROPENIA"])))
    return T
