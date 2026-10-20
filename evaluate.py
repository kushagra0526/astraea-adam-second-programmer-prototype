"""Planted ADaM derivation bugs: found and named vs a conformance-style check.  python evaluate.py [dev|heldout] [out.json]"""
import collections, json, random, sys
from qc import sdtm, adam, compare


def studies(which):
    if which == "heldout":
        h = json.load(open("heldout/config.json"))
        seeds, cfg, spec, extra = h["seeds"], h["cfg"], h["spec"], h["unknown_bugs"]
    else:
        seeds, cfg, spec, extra = list(range(1, 61)), None, dict(teae_window=30, visit_weeks=8, window=7), []
    for s in seeds:
        r = random.Random(s * 31)
        T = sdtm.make(s, cfg)
        kinds = adam.ALTS + extra
        bug = r.choice(kinds) if r.random() < 0.8 else None
        alt = frozenset([bug]) if bug in adam.ALTS else frozenset()
        ext = frozenset([bug]) if bug in extra else frozenset()
        yield s, T, spec, bug, bug in extra, adam.derive(T, spec, alt, ext)


def run(which):
    R = collections.Counter()
    per = collections.defaultdict(collections.Counter)
    rows = []
    for s, T, spec, bug, unknown, agent in studies(which):
        res = compare.qc(T, spec, agent)
        p21 = compare.p21_like(agent)
        visible = bug is not None and adam.derive(T, spec) != agent       # a bug that changes no value in this study is invisible
        key = "unknown" if unknown else ("clean" if bug is None else "known")
        R[key] += 1
        if bug and not visible:
            R["invisible"] += 1
        if visible:
            R[key + "_visible"] += 1
            R[key + "_flagged"] += res["decision"] != "PASS"
            R[key + "_p21"] += bool(p21)
            R[key + "_named_right"] += res.get("rule") == bug
            R[key + "_named_wrong"] += res.get("rule") not in (None, bug)
            per[bug]["n"] += 1; per[bug]["named"] += res.get("rule") == bug; per[bug]["p21"] += bool(p21)
            per[bug]["review"] += res["decision"] == "REVIEW"
        if bug is None:
            R["clean_flagged"] += res["decision"] != "PASS"; R["clean_p21"] += bool(p21)
        rows.append(dict(seed=s, bug=bug, decision=res["decision"], rule=res.get("rule"), mismatches=res["mismatches"], p21=len(p21)))
    return dict(set=which, summary=dict(R), per_bug={k: dict(v) for k, v in per.items()}, studies=rows)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "dev"
    r = run(which)
    S = r["summary"]
    print(which, dict(sorted(S.items())))
    for k, v in sorted(r["per_bug"].items()):
        print(f"  {k:18s} studies {v['n']:2d}  named right {v['named']:2d}  REVIEW {v.get('review', 0):2d}  conformance check flagged {v['p21']:2d}")
    if len(sys.argv) > 2:
        json.dump(r, open(sys.argv[2], "w"), indent=1, default=str)
