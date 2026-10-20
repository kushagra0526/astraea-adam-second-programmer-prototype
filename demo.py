"""An agent-built ADaM package for a synthetic oncology study, checked by an independent second program."""
from qc import sdtm, adam, compare

SPEC = dict(teae_window=30, visit_weeks=8, window=7)
T = sdtm.make(1)
ref = adam.derive(T, SPEC)
for bug in ("ignore_newtx", "impute_first", "trtsdt_zero"):
    agent = adam.derive(T, SPEC, frozenset([bug]))           # what an agent's program produced
    r = compare.qc(T, SPEC, agent)
    ds, key, var, got, want = r["example"]
    print(f"agent package with '{bug}':")
    print(f"  conformance-style check: {len(compare.p21_like(agent))} issues")
    print(f"  second programmer: {r['decision']} - {r['mismatches']} cells, {r['subjects']} subjects differ")
    print(f"  rule: {r['why']}")
    print(f"  e.g. {ds} {key} {var}: agent {got!r}, spec {want!r}\n")
agent = adam.derive(T, SPEC, extra=frozenset(["death_censored"]))
r = compare.qc(T, SPEC, agent)
print(f"agent package with a mistake outside the library: {r['decision']} - {r['why']} ({r['mismatches']} cells)")
