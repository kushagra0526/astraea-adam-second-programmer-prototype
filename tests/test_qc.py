from qc import sdtm, adam, compare

SPEC = dict(teae_window=30, visit_weeks=8, window=7)


def tiny():
    return dict(DM=[dict(USUBJID="S1", AGE=65, SEX="F", ARM="DRUG", RANDDT="2025-01-06", RFSTDTC="2025-01-08")],
                EX=[dict(USUBJID="S1", EXSTDTC="2025-01-08", EXDOSE=0), dict(USUBJID="S1", EXSTDTC="2025-01-29", EXDOSE=200)],
                RS=[dict(USUBJID="S1", RSDTC="2025-03-03", RSORRES="SD"), dict(USUBJID="S1", RSDTC="2025-04-28", RSORRES="SD"),
                    dict(USUBJID="S1", RSDTC="2025-06-23", RSORRES="PD")],
                DS=[], CM=[dict(USUBJID="S1", CMCAT="NEW ANTICANCER THERAPY", CMSTDTC="2025-05-15")],
                AE=[dict(USUBJID="S1", AESEQ=1, AESTDTC="2025-01", AEDECOD="NAUSEA"), dict(USUBJID="S1", AESEQ=2, AESTDTC="2025-01-29", AEDECOD="RASH")])


def test_reference_rules():
    a = adam.derive(tiny(), SPEC)
    assert a["ADSL"]["S1"]["TRTSDT"] == "2025-01-29"                     # first dose > 0
    assert a["ADSL"]["S1"]["AGEGR1"] == ">=65"
    assert a["ADAE"][("S1", 1)]["ASTDT"] == "2025-01-29"                 # partial date in first-dose month -> TRTSDT
    assert a["ADAE"][("S1", 2)]["ASTDY"] == 1                            # no day 0
    t = a["ADTTE"]["S1"]
    assert t["CNSR"] == 1 and t["ADT"] == "2025-04-28" and t["EVNTDESC"] == "NEW ANTICANCER THERAPY"
    assert t["AVAL"] == 113


def test_each_known_mistake_is_named():
    T = sdtm.make(1)
    for alt in adam.ALTS:
        agent = adam.derive(T, SPEC, frozenset([alt]))
        r = compare.qc(T, SPEC, agent)
        assert r["decision"] == "PASS" or r["rule"] == alt


def test_clean_passes():
    T = sdtm.make(2)
    assert compare.qc(T, SPEC, adam.derive(T, SPEC))["decision"] == "PASS"


def test_unknown_mistake_goes_to_review_not_a_wrong_name():
    T = sdtm.make(3)
    agent = adam.derive(T, SPEC, extra=frozenset(["death_censored"]))
    assert compare.qc(T, SPEC, agent)["decision"] == "REVIEW"


def test_conformance_check_misses_semantic_bugs():
    T = sdtm.make(1)
    assert compare.p21_like(adam.derive(T, SPEC, frozenset(["ignore_newtx"]))) == []


def test_diff_reports_cells():
    T = sdtm.make(1)
    d = compare.diff(adam.derive(T, SPEC, frozenset(["aval_noplus1"])), adam.derive(T, SPEC))
    assert d and all(v == "AVAL" for _, _, v, _, _ in d)
