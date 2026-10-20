"""Second-programmer QC: compare an agent's ADaM to an independent derivation, then say which SAP rule explains the gap."""
from .adam import ALTS, derive

DOC = {"trtsdt_zero": "TRTSDT taken from a zero-dose EX record (spec: first dose > 0)",
       "saffl_any": "SAFFL = Y for subjects who never received a dose > 0",
       "agegr_gt": "AGE = 65 put in '<65' (spec: >= 65)",
       "impute_first": "partial AE dates imputed to the 1st even in the first-dose month/year (spec: impute to TRTSDT)",
       "day0": "ASTDY has a day 0 (spec/SDTMIG: no day 0, +1 on or after TRTSDT)",
       "teae_window_other": "treatment-emergent window differs from the SAP",
       "ignore_newtx": "PFS not censored at the last adequate assessment before new anticancer therapy (FDA 2018 PFS conventions)",
       "ignore_missed": "PFS event counted after >= 2 missed assessments (FDA 2018: censor at last adequate assessment)",
       "aval_noplus1": "AVAL = ADT - STARTDT without +1"}


def diff(agent, ref):
    out = []
    for ds in ref:
        for key in sorted(set(ref[ds]) | set(agent.get(ds, {})), key=str):
            a, r = agent.get(ds, {}).get(key), ref[ds].get(key)
            if a is None or r is None:
                out.append((ds, key, "*", a, r))
                continue
            for v in r:
                if a.get(v) != r[v]:
                    out.append((ds, key, v, a.get(v), r[v]))
    return out


def qc(T, spec, agent):
    """-> dict(decision PASS / FIX / REVIEW, rule, mismatches, subjects, example)."""
    ref = derive(T, spec)
    d = diff(agent, ref)
    if not d:
        return dict(decision="PASS", rule=None, mismatches=0)
    subj = sorted({k[0] if isinstance(k, tuple) else k for _, k, *_ in d})
    for alt in ALTS:
        if derive(T, spec, frozenset([alt])) == agent:
            return dict(decision="FIX", rule=alt, why=DOC[alt], mismatches=len(d), subjects=len(subj), example=d[0])
    return dict(decision="REVIEW", rule=None, why="differs from the independent derivation; no known interpretation explains it",
                mismatches=len(d), subjects=len(subj), example=d[0])


def p21_like(agent):
    """Structural/conformance-style checks only (what a validator sees without the SAP)."""
    issues = []
    for u, x in agent["ADSL"].items():
        if x["AGEGR1"] not in ("<65", ">=65") or x["SAFFL"] not in ("Y", "N"):
            issues.append(("ADSL", u, "CT"))
    for k, x in agent["ADAE"].items():
        if x["ASTDY"] == 0:
            issues.append(("ADAE", k, "study day 0"))
        if x["TRTEMFL"] not in ("Y", None):
            issues.append(("ADAE", k, "CT"))
    for u, x in agent["ADTTE"].items():
        if x["CNSR"] not in (0, 1) or x["ADT"] < x["STARTDT"]:
            issues.append(("ADTTE", u, "CNSR/ADT"))
    return issues
