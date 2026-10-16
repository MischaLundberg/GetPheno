
import subprocess
import hashlib
import warnings
import pandas as pd
import numpy as np
import pytest
from datetime import datetime, date

import get_pheno as gp

# Ensure globals are set to safe defaults for tests
gp.only_ICD10 = False
gp.only_ICD9 = False
gp.only_ICD8 = False
gp.cluster_run = "Default"
gp.verbose = False
gp.DayFirst = False
gp.DateFormat = "%Y-%m-%d"
#gp.logger = setup_logger("test_get_pheno", to_console=False)

def test_match_codes_exact_and_prefix():
    s = pd.Series(["ICD10:F32", "ICD10:F321", "ATC:N06A", "OTHER"])
    # exact match
    assert gp.match_codes(s, ["ICD10:F32"], exact=True).tolist() == [True, False, False, False]
    # prefix match
    assert gp.match_codes(s, ["ICD10:F32"], exact=False).tolist() == [True, True, False, False]
    # ATC prefix matching
    assert gp.match_codes(pd.Series(["ATC:N06A", "ATC:N06AB"]), ["ATC:N06A"], exact=False).tolist() == [True, True]


def test_wildcards_are_allowed_with_and_without_exact_match():
    df = pd.DataFrame({
        "iid": ["a", "b", "c", "d", "e", "f"],
        "diagnosis": ["F250", "F339", "F32.2", "F32.20", "F321", "G40"],
    })
    requested = ["F25*", "F33*", "F32.2"]

    prefix_hits = gp.map_cases(requested, exact_match=False, df1=df, diagcol="diagnosis")
    exact_hits = gp.map_cases(requested, exact_match=True, df1=df, diagcol="diagnosis")

    assert set(prefix_hits["iid"]) == {"a", "b", "c", "d"}
    assert set(exact_hits["iid"]) == {"a", "b", "c"}


def test_expanded_ranges_follow_exact_match_mode():
    df = pd.DataFrame({
        "iid": ["a", "b", "c", "d", "e", "f"],
        "diagnosis": ["ICD10:T36", "ICD10:T360", "ICD10:T37", "ICD10:T38", "ICD10:T380", "ICD10:T39"],
    })
    requested = gp.expand_ranges("ICD10:T36-ICD10:T38")

    prefix_hits = gp.map_cases(requested, exact_match=False, df1=df, diagcol="diagnosis")
    exact_hits = gp.map_cases(requested, exact_match=True, df1=df, diagcol="diagnosis")

    assert requested == ["ICD10:T36", "ICD10:T37", "ICD10:T38"]
    assert set(prefix_hits["iid"]) == {"a", "b", "c", "d", "e"}
    assert set(exact_hits["iid"]) == {"a", "c", "d"}


def test_remove_leading_icd_and_format_numeric():
    assert gp.remove_leading_icd("ICD10:F32") == "F32"
    assert gp.remove_leading_icd("ICD8:123.4") == "123.4"
    # non-string is returned as-is
    assert gp.remove_leading_icd(123) == 123
    # format_numeric with DK cluster mode pads integer and decimal parts
    # Use a mode that is in DK_clusters
    res = gp.format_numeric(12.3, "CHB_DBDS")
    # Expect integer padded to 3 digits and decimal padded to 2 digits -> "012.30"
    assert res == "012.30"


def test_normalize_atc_codes_preserves_prefix_and_wildcards():
    assert gp.normalize_atc_codes(["atc:n06a*", " n05 "], expect_prefix=None) == ["ATC:N06A*", "N05"]
    assert gp.normalize_atc_codes(["n06a"], expect_prefix=True) == ["ATC:N06A"]
    assert gp.normalize_atc_codes(["ATC:N06A"], expect_prefix=False) == ["N06A"]


def test__to_datetime_series_and_convert_if_not_datetime(tmp_path):
    s = pd.Series(["20220101", "2022-01-02", "03/04/2022", None, "2022-13-01"])
    out = gp._to_datetime_series(s, fmt=None, dayfirst=False)
    # first parsed as 2022-01-01 via YYYYMMDD
    assert out.iloc[0].date() == date(2022, 1, 1)
    # second parsed directly
    assert out.iloc[1].date() == date(2022, 1, 2)
    # third will parse with fallback - dayfirst False -> mm/dd/YYYY fails -> coerce or attempt variants
    assert isinstance(out.iloc[2], (pd.Timestamp, type(pd.NaT)))  # either parsed or NaT
    # convert_if_not_datetime
    assert gp.convert_if_not_datetime(None) is pd.NaT
    dt = datetime(2000, 1, 2)
    assert gp.convert_if_not_datetime(dt) == dt
    assert gp.convert_if_not_datetime("2001-01-01").date() == date(2001, 1, 1)


def test__as_list_and__to_dt_list():
    assert gp._as_list(5) == [5]
    assert gp._as_list([1, 2, ""]) == [1, 2]
    dt_list = gp._to_dt_list(["2020-01-01", pd.NaT, "1999-12-31"])
    assert all(isinstance(x, pd.Timestamp) for x in dt_list)
    assert len(dt_list) == 2


def test_normalize_iid_series_and_auto():
    s = pd.Series(["001", "002.0", "ABC", "003"])
    # target str keeps leading zeros, removes .0
    out = gp.normalize_iid_series(s, target="str")
    assert out.tolist()[0] == "001"
    assert out.tolist()[1] == "002"
    # target int falls back to string if non-numeric present
    out2 = gp.normalize_iid_series(s, target="int")
    assert out2.dtype.name == "string" or out2.dtype.name == "object"
    # auto chooses Int64 if all numeric
    s2 = pd.Series(["10", "20.0", "30"])
    auto = gp.normalize_iid_series(s2, target="int")
    assert getattr(auto.dtype, "name", "") in ("Int64", "int64") or pd.api.types.is_integer_dtype(auto)


def test_convert_if_not_datetime_various():
    assert gp.convert_if_not_datetime(pd.NaT) is pd.NaT
    assert gp.convert_if_not_datetime(date(2000, 1, 1)) == date(2000, 1, 1)
    assert gp.convert_if_not_datetime("2000-01-02").date() == date(2000, 1, 2)


def test_ipsych_cohort_restrictions_do_not_filter_on_parental_birthplace():
    df = pd.DataFrame({
        "iid": ["keep", "early", "end", "foreign", "no_mother", "not_age_one", "parents_foreign", "twin_a", "twin_b"],
        "birthdate": [
            "1981-05-01", "1981-04-30", "2009-01-01", "1990-01-01", "1990-01-01",
            "1990-01-01", "1990-01-01", "1991-02-03", "1991-02-03",
        ],
        "fkode": [100, 100, 100, 6000, 100, 100, 100, 100, 100],
        "pnr_m": ["m1", "m2", "m3", "m4", pd.NA, "m6", "m7", "twins", "twins"],
        "stat": ["1", "1", "1", "1", "1", "2", "1", "1", "1"],
        "statd": [
            "1981-05-02", "1981-05-02", "2009-01-02", "1990-01-02", "1990-01-02",
            "1990-04-01", "1990-01-02", "1991-02-04", "1991-02-04",
        ],
        "DKDK": [True, True, True, True, True, True, False, True, True],
    })

    result = gp.apply_ipsych_cohort_restrictions(
        df,
        iidcol="iid",
        options={"enabled": True, "dk_parents_col": "DKDK"},
    )

    assert result["iid"].tolist() == ["keep", "parents_foreign"]
    assert result.set_index("iid").loc["parents_foreign", "both_parents_DK"] == False


def test_ipsych_cohort_warns_when_dkdk_indicator_is_unavailable(capsys):
    df = pd.DataFrame({
        "iid": ["keep"],
        "birthdate": ["1990-01-01"],
        "fkode": [100],
        "pnr_m": ["m1"],
        "stat": ["1"],
        "statd": ["1990-01-02"],
    })

    result = gp.apply_ipsych_cohort_restrictions(
        df,
        iidcol="iid",
        options={"enabled": True},
    )

    assert result["iid"].tolist() == ["keep"]
    assert "Danish-born-parents (DKDK) indicator could not be reported" in capsys.readouterr().out


def test_ipsych_restrictions_use_status_columns_from_additional_input(capsys):
    stam = pd.DataFrame({
        "iid": ["keep", "exclude"],
        "birthdate": ["1990-01-01", "1990-01-01"],
        "fkode": [100, 100],
        "pnr_m": ["m1", "m2"],
        "DKDK": [True, True],
    })
    additional = pd.DataFrame({
        "iid": ["keep", "exclude"],
        "status_code": ["1", "2"],
        "status_date": ["1990-01-02", "1990-04-01"],
    })

    filtered_stam, filtered_additional = gp.restrict_ipsych_iid_inputs(
        stam,
        additional,
        iidcol="iid",
        birthdatecol="birthdate",
        options={
            "enabled": True,
            "dk_parents_col": "DKDK",
            "status_col": "status_code",
            "status_date_col": "status_date",
        },
    )

    assert filtered_stam["iid"].tolist() == ["keep"]
    assert filtered_additional["iid"].tolist() == ["keep"]
    assert "neither status column" not in capsys.readouterr().out


def test_death_filters_support_age_range_living_and_deceased_modes():
    metadata = pd.DataFrame({
        "iid": ["young_death", "older_death", "alive"],
        "birthdate": ["1980-01-01"] * 3,
        "D_DODSDATO": ["2000-01-01", "2020-01-01", ""],
    })

    age_only = gp.apply_death_filters(
        metadata,
        deathdatecol="D_DODSDATO",
        min_max_age_death="20,30",
    )
    age_and_alive = gp.apply_death_filters(
        metadata,
        deathdatecol="D_DODSDATO",
        min_max_age_death="20,30",
        keep_alive=True,
    )
    older_dead = gp.apply_death_filters(
        metadata,
        deathdatecol="D_DODSDATO",
        min_max_age_death="30,50",
        only_dead=True,
    )

    assert age_only["iid"].tolist() == ["young_death"]
    assert age_and_alive["iid"].tolist() == ["young_death", "alive"]
    assert older_dead["iid"].tolist() == ["older_death"]
    assert age_and_alive.set_index("iid").loc["young_death", "AgeAtDeath"] == 20
    assert pd.isna(age_and_alive.set_index("iid").loc["alive", "AgeAtDeath"])


def test_death_filters_reject_conflicting_status_modes():
    with pytest.raises(ValueError, match="cannot be enabled together"):
        gp.apply_death_filters(
            pd.DataFrame({
                "birthdate": ["1980-01-01"],
                "D_DODSDATO": ["2020-01-01"],
            }),
            deathdatecol="D_DODSDATO",
            keep_alive=True,
            only_dead=True,
        )


def test_ipsych_death_date_enforces_survival_to_first_birthday():
    df = pd.DataFrame({
        "iid": ["alive", "before_one", "on_first_birthday"],
        "birthdate": ["1990-01-01"] * 3,
        "D_DODSDATO": ["", "1990-12-31", "1991-01-01"],
        "fkode": [100] * 3,
        "pnr_m": ["m1", "m2", "m3"],
        "stat": ["1", "1", "1"],
        "DKDK": [True] * 3,
    })

    result = gp.apply_ipsych_cohort_restrictions(
        df,
        iidcol="iid",
        options={
            "enabled": True,
            "dk_parents_col": "DKDK",
            "deathdate_col": "D_DODSDATO",
        },
    )

    assert result["iid"].tolist() == ["alive", "on_first_birthday"]


def test_ipsych_status_only_does_not_require_death_date():
    df = pd.DataFrame({
        "iid": ["keep"],
        "birthdate": ["1990-01-01"],
        "fkode": [100],
        "pnr_m": ["m1"],
        "stat": ["1"],
        "DKDK": [True],
    })

    result = gp.apply_ipsych_cohort_restrictions(
        df,
        iidcol="iid",
        options={"enabled": True, "dk_parents_col": "DKDK"},
    )

    assert result["iid"].tolist() == ["keep"]


def test_ipsych_diagnosis_window_uses_date_in_and_date_out_overlap():
    df = pd.DataFrame({
        "iid": ["starts_on_boundary", "crosses_end", "starts_at_end", "before_window", "undated"],
        "date_in": ["2009-12-01", "2015-12-31", "2016-01-01", "2008-01-01", None],
        "date_out": ["2010-01-01", "2016-01-02", "2016-01-03", "2009-12-31", None],
    })

    result = gp.filter_records_to_ipsych_date_window(
        df,
        input_date_in_name="date_in",
        input_date_out_name="date_out",
        options={
            "enabled": True,
            "diagnosis_start": "2010-01-01",
            "diagnosis_end": "2016-01-01",
        },
    )

    assert result["iid"].tolist() == ["starts_on_boundary", "crosses_end"]


def test_ipsych_diagnosis_window_requires_both_date_columns():
    with pytest.raises(ValueError, match="requires both date-in and date-out"):
        gp.filter_records_to_ipsych_date_window(
            pd.DataFrame({"date_in": ["2015-01-01"]}),
            input_date_in_name="date_in",
            input_date_out_name="date_out",
            options={"enabled": True},
        )


def test_remove_duplicates_preserve_order():
    lst = [1, 2, 2, 3, 1, 4]
    assert gp.remove_duplicates_preserve_order(lst) == [1, 2, 3, 4]


def test_generate_cpr_enc_and_random_date():
    s = gp.generate_cpr_enc()
    assert isinstance(s, str) and len(s) >= 40
    d1 = date(2000, 1, 1)
    d2 = date(2000, 1, 10)
    r = gp.generate_random_date(d1, d2)
    assert isinstance(r, date)
    assert d1 <= r <= d2


def test_setup_logger_and_usage(monkeypatch, tmp_path):
    # capture stdout logs via stream handler, ensure setup_logger returns logger and handlers configured
    log_fn = tmp_path / "test.log"
    logger = gp.setup_logger("testscript", to_console=False, to_file=str(log_fn))
    assert isinstance(logger, type(gp.logger))
    # usage() calls psutil.Process, ensure it runs without error
    gp.logger = logger  # set module logger to prevent None references
    gp.usage()


def test_expand_ranges():
    assert gp.expand_ranges("T36-T38,T40-T42,T45") == ["T36", "T37", "T38", "T40", "T41", "T42", "T45"]
    assert gp.expand_ranges(["T36-T38"]) == ["T36", "T37", "T38"]
    assert gp.expand_ranges(["36-38"]) == ["36", "37", "38"]
    assert gp.expand_ranges("ICD10:T36-T38") == ["ICD10:T36", "ICD10:T37", "ICD10:T38"]
    assert gp.expand_ranges(["X20"]) == ["X20"]


def test_update_icd_coding_and_process_entry(monkeypatch):
    # ensure processing simple ICD10 entry works
    gp.cluster_run = "Default"
    out = gp.update_icd_coding(["ICD10:F32", "ATC:N06A"], eM=False, skip=False, remove_point_in_diag_request=True, ICDCM=False, no_Fill=False, noLeadingICD=True, icdprefix="")
    # When remove_point_in_diag_request True, points removed
    assert any("F32" in str(x) or "N06A" in str(x) for x in out)
    # process_entry behavior for ATC and ICD8
    assert gp.process_entry("ATC:N06A", remove_leading=True, eM=False, mode="Default", icdprefix="", remove_point=True, ICDCM=False) == "N06A"
    # numeric entry interpreted as ICD8 formatting when mode in DK_clusters; test with integer input
    res = gp.process_entry(123, remove_leading=False, eM=False, mode="CHB_DBDS", icdprefix="", remove_point=False, ICDCM=False)
    assert isinstance(res, str)


def test_icd8_non_exact_short_code_is_as_inclusive_as_explicit_subcodes(monkeypatch):
    monkeypatch.setattr(gp, "cluster_run", "CHB_DBDS")
    monkeypatch.setattr(gp, "only_ICD10", False)
    monkeypatch.setattr(gp, "only_ICD9", False)
    monkeypatch.setattr(gp, "only_ICD8", False)

    short_codes = gp.update_icd_coding(
        ["ICD8:296.2", "ICD8:298.0", "ICD8:300.4"],
        eM=False,
        remove_point_in_diag_request=True,
        noLeadingICD=True,
    )
    long_codes = gp.update_icd_coding(
        ["ICD8:296.2", "ICD8:296.29", "ICD8:298.0", "ICD8:298.09", "ICD8:300.4", "ICD8:300.49"],
        eM=False,
        remove_point_in_diag_request=True,
        noLeadingICD=True,
    )

    df = pd.DataFrame({
        "iid": ["a", "b", "c", "d", "e", "f"],
        "diagnosis": ["29620", "29629", "29800", "29809", "30040", "30049"],
    })

    short_hits = gp.map_cases(short_codes, exact_match=False, df1=df, diagcol="diagnosis")
    long_hits = gp.map_cases(long_codes, exact_match=False, df1=df, diagcol="diagnosis")

    assert short_codes == ["2962", "2980", "3004"]
    assert set(short_hits["iid"]) == set(long_hits["iid"])


def test_icd8_wildcard_remains_prefix_match_even_with_exact_match(monkeypatch):
    monkeypatch.setattr(gp, "cluster_run", "CHB_DBDS")
    monkeypatch.setattr(gp, "only_ICD10", False)
    monkeypatch.setattr(gp, "only_ICD9", False)
    monkeypatch.setattr(gp, "only_ICD8", False)

    codes = gp.update_icd_coding(
        ["ICD8:300.4*"],
        eM=True,
        remove_point_in_diag_request=True,
        noLeadingICD=True,
    )
    df = pd.DataFrame({
        "iid": ["a", "b", "c"],
        "diagnosis": ["30040", "30049", "30050"],
    })

    hits = gp.map_cases(codes, exact_match=True, df1=df, diagcol="diagnosis")

    assert codes == ["3004*"]
    assert set(hits["iid"]) == {"a", "b"}


def test_dict_update_icd_coding_basic():
    curr = pd.DataFrame({
        "Disorder": ["X"],
        "Disorder Codes": ["ICD10:F32,ICD10:F33"]
    })
    pheno, updated_list = gp.dict_update_icd_coding(curr.copy(), exact_match=False, skip_icd_update=False,
                                                   remove_point_in_diag_request=False, ICDCM=False, noLeadingICD=False, icdprefix="")
    # pheno should be a DataFrame and updated_list be a list
    assert isinstance(pheno, pd.DataFrame)
    assert isinstance(updated_list, list)


def test_load_phenotypes_normalizes_all_multi_phenotype_rows(tmp_path, monkeypatch):
    monkeypatch.setattr(gp, "cluster_run", "Default")
    request = tmp_path / "multi.tsv"
    request.write_text(
        "\n".join([
            "mdd\tF32,F33",
            "skiz\tF20",
            "skizospek\tF20,F21,F22,F23,F24,F25,F28,F29",
            "bipol\tF30,F31",
            "affek\tF30,F31,F32,F33,F34,F35,F36,F37,F38,F39",
            "autism\tF84.0,F84.1,F84.5,F84.8,F84.9",
            "adhd\tF90.0",
            "anorek\tF50.0,F50.1",
        ]),
        encoding="utf-8",
    )

    pheno, multi, normalized = gp.load_phenotypes(
        pheno_request=str(request),
        pheno_requestcol="diagnosis",
        skip_icd_update=True,
        exact_match=False,
        remove_point_in_diag_request=False,
    )

    assert multi is True
    assert set(normalized) == {
        "F20", "F21", "F22", "F23", "F24", "F25", "F28", "F29",
        "F30", "F31", "F32", "F33", "F34", "F35", "F36", "F37", "F38", "F39",
        "F50.0", "F50.1", "F84.0", "F84.1", "F84.5", "F84.8", "F84.9", "F90.0",
    }
    assert set(gp.flatten_disorder_codes(pheno)) == set(normalized)


def test_parse_pheno_rules_basic():
    # simple rule with range
    gp.cluster_run = "Default"
    r = gp.parse_pheno_rules("main=T36-T38;sub=F10,F11", exact_match=False, skip_icd_update=False, remove_point_in_diag_request=False, ICDCM=False, noLeadingICD=False, icdprefix="")
    r2 = gp.parse_pheno_rules("T36-T38", exact_match=False, skip_icd_update=False, remove_point_in_diag_request=False, ICDCM=False, noLeadingICD=False, icdprefix="")
    r3 = gp.parse_pheno_rules("main=T36-T38;sub=F10,F11", exact_match=False, skip_icd_update=False, remove_point_in_diag_request=True, ICDCM=False, noLeadingICD=False, icdprefix="D")
    r4 = gp.parse_pheno_rules("main=T36-T38;sub=F10,F11", exact_match=False, skip_icd_update=False, remove_point_in_diag_request=True, ICDCM=False, noLeadingICD=True, icdprefix="D")
    r5 = gp.parse_pheno_rules("main=F;sub=T36-T38,F10,F11", exact_match=False, skip_icd_update=False, remove_point_in_diag_request=True, ICDCM=False, noLeadingICD=False, icdprefix="D")
    gp.cluster_run = "CHB_DBDS"
    r6 = gp.parse_pheno_rules("main=T36-T38;sub=F10,F11", exact_match=False, skip_icd_update=False, remove_point_in_diag_request=True, ICDCM=False, noLeadingICD=False, icdprefix="D")
    r7 = gp.parse_pheno_rules("main=F;sub=T36-T38,F10,F11", exact_match=False, skip_icd_update=False, remove_point_in_diag_request=True, ICDCM=False, noLeadingICD=False, icdprefix="D")
    r8 = gp.parse_pheno_rules("main=F;sub=T36-T50,T52-T60", exact_match=False, skip_icd_update=False, remove_point_in_diag_request=True, ICDCM=False, noLeadingICD=False, icdprefix="D")
    assert "main" in r and "sub" in r
    assert isinstance(r["main"], list)
    assert r['main'] == ['ICD10:T36', 'ICD10:T37', 'ICD10:T38']
    assert r['sub'] == ['ICD10:F10', 'ICD10:F11']
    assert "ranges" in r
    assert isinstance(r["ranges"], list)
    assert r2['ranges'] == ['ICD10:T36', 'ICD10:T37', 'ICD10:T38']
    assert r3['main'] == ['ICD10:DT36', 'ICD10:DT37', 'ICD10:DT38']
    assert r3['sub'] == ['ICD10:DF10', 'ICD10:DF11']
    assert r4['main'] == ['DT36', 'DT37', 'DT38']
    assert r4['sub'] == ['DF10', 'DF11']
    assert r5['main'] == ['ICD10:DF']
    assert r5['sub'] == ['ICD10:DT36', 'ICD10:DT37', 'ICD10:DT38', 'ICD10:DF10', 'ICD10:DF11']
    assert r6['main'] == ['ICD10:DT36', 'ICD10:DT37', 'ICD10:DT38']
    assert r6['sub'] == ['ICD10:DF10', 'ICD10:DF11']
    assert r7['main'] == ['ICD10:DF']
    assert r7['sub'] == ['ICD10:DT36', 'ICD10:DT37', 'ICD10:DT38', 'ICD10:DF10', 'ICD10:DF11']


def test_merge_IIDs_aggregation():
    # Create a small tmp_result_df with two diagnoses, different dates
    df = pd.DataFrame({
        "cpr": ["A", "A", "B"],
        "diagnosis": ["ICD10:F32", "ICD10:F33", "ICD10:F32"],
        "date_in": [pd.Timestamp("2020-01-01"), pd.Timestamp("2021-01-01"), pd.Timestamp("2022-01-01")],
        "date_out": [pd.Timestamp("2020-01-02"), pd.Timestamp("2021-01-02"), pd.Timestamp("2022-01-02")],
    })
    merged = gp.merge_IIDs(
        tmp_result_df=df,
        diagnostic_col="diagnosis",
        birthdatecol="birthdate",
        input_date_in_name="date_in",
        input_date_out_name="date_out",
        iidcol="cpr",
        verbose=False,
        Cases=True,
        Covariates=False,
        BuildEntryExitDates=False
    )
    # Expect one row per IID
    assert merged["cpr"].nunique() == 2
    assert "diagnoses" in merged.columns
    assert "first_dx" in merged.columns or "date_in" in merged.columns


def test_select_by_iid_and_diag_optimized_and_h5(tmp_path):
    # create small dataframe and write to h5 store
    df = pd.DataFrame({
        "pnr": [1, 2, 3, 2],
        "diagnosis": ["A01", "B02", "A01", "A019"],
        "date_in": [pd.Timestamp("2020-01-01")] * 4
    })
    h5p = str(tmp_path / "teststore.h5")
    with pd.HDFStore(h5p, mode="w") as st:
        st.append("df", df, format="table", data_columns=["pnr", "diagnosis"])
    # select by diag prefix
    out = gp.select_by_iid_and_diag_optimized(h5_path=h5p, table_name="df", iidcol="pnr", iids=[2,3], diagcol="diagnosis", diags=["A01"])
    assert isinstance(out, pd.DataFrame)
    # should find pnr 3 and pnr 2 (A019 startswith A01)
    assert set(out["pnr"].tolist()).issubset({2, 3})
    # exact matching (no prefix)
    out_exact = gp.select_by_iid_and_diag_optimized(h5_path=h5p, table_name="df", iidcol="pnr", iids=None, diagcol="diagnosis", diags=["A01"], prefix_all=False)
    assert isinstance(out_exact, pd.DataFrame)


def test_get_h5_cases_integration(tmp_path):
    # create store and call get_h5_cases wrapper
    df = pd.DataFrame({"iid": [10, 20, 30], "diagnosis": ["ICD10:F32", "ICD10:F33", "OTHER"], "date_in": pd.to_datetime(["2020-01-01"]*3)})
    h5p = str(tmp_path / "test2.h5")
    with pd.HDFStore(h5p, mode="w") as st:
        st.append("df", df, format="table", data_columns=["iid", "diagnosis"])
    out = gp.get_h5_cases(h5file=h5p, iids=[10, 30], iidcol="iid", diags=["ICD10:F32"], diagcol="diagnosis", directmapping=True, table_name="df")
    assert isinstance(out, pd.DataFrame)


def test_load_data_file_and_finalize(tmp_path):
    # create a small csv
    csvp = tmp_path / "stam.csv"
    df = pd.DataFrame({"pnr": [1, 2], "birthdate": ["01/01/1980", "02/02/1990"], "diagnosis": ["ICD10:F32", "ICD10:F33"], "sex": ["M", "F"]})
    df.to_csv(csvp, index=False)
    res = gp.load_data_file(data_file=str(csvp), sep=",", birthdatecol="birthdate", diagnostic_col="diagnosis", sexcol="sex", stam_cols_to_read_as_date=["birthdate"])
    assert "birthdate" in res.columns
    # finalize_lpr_data should rename columns accordingly
    df2 = pd.DataFrame({"c_adiag": ["A"], "d_inddto": [pd.Timestamp("2020-01-01")]})
    fin = gp.finalize_lpr_data(df2.copy(), diagnostic_col="c_adiag", birthdatecol="d_inddto", ctype_col="", verbose=False)
    assert "diagnosis" in fin.columns


def test_load_additional_information_files_merges_columns_by_iid(tmp_path):
    standard = tmp_path / "person.tsv"
    extra = tmp_path / "laboratory.tsv"
    pd.DataFrame({
        "pnr": ["1", "2"],
        "birthdate": ["1980-01-01", "1990-02-02"],
        "region": ["north", None],
    }).to_csv(standard, sep="\t", index=False)
    pd.DataFrame({
        "pnr": ["2", "3"],
        "region": ["south", "west"],
        "marker": [12, 18],
    }).to_csv(extra, sep="\t", index=False)

    result = gp.load_additional_information_files(
        data_file=f"{standard},{extra}",
        sep="\t",
        birthdatecol="birthdate",
        diagnostic_col="diagnosis",
        sexcol="sex",
        stam_cols_to_read_as_date=["birthdate"],
        iidcol="pnr",
    )

    assert result.columns.tolist() == ["pnr", "birthdate", "region", "marker"]
    assert result["pnr"].tolist() == ["1", "2", "3"]
    assert result.set_index("pnr").loc["2", "region"] == "south"
    assert result.set_index("pnr").loc["3", "marker"] == 18


def test_additional_files_preserve_large_numeric_looking_iids(tmp_path):
    person = tmp_path / "t_person.tsv"
    deaths = tmp_path / "t_dodsaarsag_2.tsv"
    person.write_text(
        "cpr_enc\tD_FODDATO\tS_KON\n"
        "9007199254740993\t1980-01-01\tF\n"
        "9007199254740995\t1990-02-02\tM\n",
        encoding="utf-8",
    )
    # The empty IID makes pandas infer float without an explicit IID dtype;
    # values above 2**53 would then lose precision before the merge.
    deaths.write_text(
        "cpr_enc\tD_DODSDATO\tD_FINDEDATO\n"
        "9007199254740993\t2020-03-04\t2020-03-05\n"
        "9007199254740995\t2021-04-05\t2021-04-06\n"
        "\t2022-01-01\t2022-01-02\n",
        encoding="utf-8",
    )

    result = gp.load_additional_information_files(
        data_file=f"{person},{deaths}",
        sep="\t",
        birthdatecol="D_FODDATO",
        diagnostic_col="diagnosis",
        sexcol="S_KON",
        stam_cols_to_read_as_date=["D_FODDATO", "D_DODSDATO", "D_FINDEDATO"],
        iidcol="cpr_enc",
    ).set_index("cpr_enc")

    assert result.loc["9007199254740993", "D_DODSDATO"] == pd.Timestamp("2020-03-04")
    assert result.loc["9007199254740995", "D_FINDEDATO"] == pd.Timestamp("2021-04-06")


def test_multiple_lpr_files_map_per_file_columns_and_fill_missing(tmp_path, capsys):
    diagnoses = tmp_path / "lpr3_diagnoses.tsv"
    procedures = tmp_path / "lpr3_procedures.tsv"
    pd.DataFrame({
        "cpr_enc": ["1"],
        "source": ["lpr3_diagnoses"],
        "date_in": ["2020-01-01"],
        "date_out": ["2020-01-02"],
        "contact_type": ["inpatient"],
        "diagnosis_type": ["A"],
        "diagnosis": ["F32"],
    }).to_csv(diagnoses, sep="\t", index=False)
    pd.DataFrame({
        "cpr_enc": ["2"],
        "source": ["lpr3_procedures"],
        "procedure_date": ["2021-03-04"],
        "date_out": ["2021-03-05"],
        "code": ["KABC"],
    }).to_csv(procedures, sep="\t", index=False)

    schema = gp.build_lpr_input_schema(
        lpr_file=f"{diagnoses},{procedures}",
        diagnostic_col="diagnosis,code",
        iidcol="cpr_enc",
        input_date_in_name="date_in,procedure_date",
        input_date_out_name="date_out",
        ctype_col="diagnosis_type",
    )
    headers = [
        {"path": str(diagnoses), "columns": gp._read_input_header(str(diagnoses), "\t")},
        {"path": str(procedures), "columns": gp._read_input_header(str(procedures), "\t")},
    ]
    gp.validate_lpr_input_schema_headers(schema, headers)
    result = gp.process_lpr_data(
        lpr_file=f"{diagnoses},{procedures}",
        lpr2nd_file="",
        dta_input=False,
        sep="\t",
        lpr_cols_to_read_as_date=[],
        DateFormat="%Y-%m-%d",
        potential_lpr_cols_to_read_as_date=[],
        diagnostic_col="diagnosis",
        diagnostic2nd_col="",
        lpr_recnummer="",
        lpr2nd_recnummer="",
        iidcol="cpr_enc",
        lpr_schema=schema,
        input_date_in_name="date_in",
        input_date_out_name="date_out",
        ctype_col="diagnosis_type",
    )

    assert result["diagnosis"].tolist() == ["F32", "KABC"]
    assert result["cpr_enc"].tolist() == ["1", "2"]
    assert result.loc[1, "date_in"] == pd.Timestamp("2021-03-04")
    assert pd.isna(result.loc[1, "diagnosis_type"])
    warning = capsys.readouterr().out
    assert "--DiagTypecol" in warning
    assert str(procedures) in warning
    assert warning.count("--DiagTypecol") == 1


def test_explicit_blank_per_file_column_adds_na_without_warning(tmp_path, capsys):
    first = tmp_path / "diagnoses.tsv"
    second = tmp_path / "procedures.tsv"
    pd.DataFrame({"iid": ["1"], "diagnosis": ["F32"], "dtype": ["A"]}).to_csv(first, sep="\t", index=False)
    pd.DataFrame({"iid": ["2"], "code": ["KABC"]}).to_csv(second, sep="\t", index=False)
    schema = gp.build_lpr_input_schema(
        f"{first},{second}", "diagnosis,code", "iid", ctype_col="dtype,"
    )
    headers = [
        {"path": str(first), "columns": gp._read_input_header(str(first), "\t")},
        {"path": str(second), "columns": gp._read_input_header(str(second), "\t")},
    ]
    gp.validate_lpr_input_schema_headers(schema, headers)
    result = gp.process_lpr_data(
        f"{first},{second}", "", False, "\t", [], "%Y-%m-%d", [],
        "diagnosis", "", "", "", "iid", lpr_schema=schema, ctype_col="dtype",
    )
    assert pd.isna(result.loc[1, "dtype"])
    assert "--DiagTypecol" not in capsys.readouterr().out


def test_multiple_lpr_files_work_in_low_memory_batch(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    diagnoses = tmp_path / "diagnoses.tsv"
    procedures = tmp_path / "procedures.tsv"
    pd.DataFrame({"iid": ["1", "2"], "diagnosis": ["F32", "F33"]}).to_csv(diagnoses, sep="\t", index=False)
    pd.DataFrame({"iid": ["2", "3"], "code": ["PROC2", "PROC3"]}).to_csv(procedures, sep="\t", index=False)
    schema = gp.build_lpr_input_schema(
        f"{diagnoses},{procedures}", "diagnosis,code", "iid"
    )

    result = gp.batch_load_lprfile(
        df=pd.DataFrame(),
        lprfile=f"{diagnoses},{procedures}",
        lpr_recnummer="",
        lpr2nd_file="",
        lpr2nd_recnummer="",
        iidcol="iid",
        iid_batch=["2", "3"],
        batch_num=0,
        sep="\t",
        potential_lpr_cols_to_read_as_date=[],
        lpr_cols_to_read_as_date=[],
        verbose=False,
        dta_input=False,
        DateFormat="%Y-%m-%d",
        diagnostic_col="diagnosis",
        diagnostic2nd_col="",
        lpr_schema=schema,
    )

    assert set(result["iid"]) == {"2", "3"}
    assert set(result["diagnosis"]) == {"F33", "PROC2", "PROC3"}
    assert not list(tmp_path.glob("*.filtered_temp.csv"))


def test_lpr_column_mapping_count_must_match_file_count():
    with pytest.raises(ValueError, match="--fcol supplies 2 column names for 3 -f files"):
        gp.build_lpr_input_schema("one.tsv,two.tsv,three.tsv", "diagnosis,code", "iid")


def test_select_input_file_variant_keeps_matching_column_names():
    path, columns = gp._select_input_file_variant(
        "psyk.tsv,lpr.tsv",
        1,
        {
            "diagnostic_col": "psyk_diag,lpr_diag",
            "iidcol": "iid",
            "ctype_col": "psyk_type,",
        },
    )

    assert path == "lpr.tsv"
    assert columns == {
        "diagnostic_col": "lpr_diag",
        "iidcol": "iid",
        "ctype_col": "",
    }


def test_build_temp_file_and_load_mapping_rows(monkeypatch, tmp_path):
    # create a small file
    csvp = tmp_path / "big.csv"
    df = pd.DataFrame({"id": [f"id{i}" for i in range(10)], "v": list(range(10))})
    df.to_csv(csvp, index=False)
    # load_mapping_rows should identify row indices
    rows = gp.load_mapping_rows(str(csvp), iidcol="id", target_iids=["id1", "id3"], sep=",")
    assert isinstance(rows, list)
    # Monkeypatch subprocess.run used in build_temp_file to avoid awk on system
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: None)
    # call build_temp_file with the indices from load_mapping_rows
    tempf = tmp_path / "filtered_temp.csv"
    index_file = tmp_path / "row_indices.txt"
    gp.build_temp_file(str(csvp), [1, 3], temp_file=str(tempf), index_file=str(index_file), verbose=True)
    # Because we monkeypatched subprocess.run, file may not be created; ensure the function ran without exception


def test_map_cases_and_build_phenotype_cases_basic():
    df = pd.DataFrame({
        "iid": ["A", "A", "B"],
        "diagnosis": ["ICD10:F32", "ICD10:F33", "ICD10:F32"],
        "date_in": [pd.Timestamp("2020-01-01"), pd.Timestamp("2021-01-02"), pd.Timestamp("2022-01-01")],
        "date_out": [pd.Timestamp("2020-01-02"), pd.Timestamp("2021-01-03"), pd.Timestamp("2022-01-02")]
    })
    # map_cases basic
    out = gp.map_cases(values_to_match=["ICD10:F32"], exact_match=False, df1=df, diagcol="diagnosis")
    assert not out.empty
    # build_phenotype_cases returns aggregated structure
    merged = gp.build_phenotype_cases(df1=df, exact_match=False, values_to_match=["ICD10:F32"], diagnostic_col="diagnosis", birthdatecol="birthdate",
                                      iidcol="iid", input_date_in_name="date_in", input_date_out_name="date_out", verbose=False, Covariates=False)
    assert isinstance(merged, pd.DataFrame)


def test_process_ophold_minimal(tmp_path):
    # Build minimal stam and ophold tables
    stam = pd.DataFrame({"pnr": ["A", "B"], "fkode": [6000, 6000], "fkode_m": [6000, 6000], "fkode_f": [6000, 6000], "birthdate": [pd.Timestamp("1990-01-01"), pd.Timestamp("1980-01-01")]})
    ophold = pd.DataFrame({"pnr": ["A", "A", "B"], "stat": [10, 10, 10], "statd": [pd.Timestamp("2000-01-01"), pd.Timestamp("2001-01-01"), pd.Timestamp("2002-01-01")], "tflytd": [pd.Timestamp("2000-01-01"), pd.Timestamp("2001-01-01"), pd.Timestamp("2002-01-01")], "fflytd": [pd.Timestamp("2000-01-01"), pd.Timestamp("2001-01-01"), pd.Timestamp("2002-01-01")], "komkod": ["x", "y", "z"], "orig": ["o", "o", "o"], "opholdnr": [1,2,3]})
    res = gp.process_ophold(ophold=ophold, stam=stam, tmp_result_df="", ophold_out_file="", birthdatecol="birthdate", iidcol="pnr", verbose=False)
    assert "moved_to_dk" in res.columns


def test_update_DxDates_multi_exclusion_basic():
    # Prepare a simple DataFrame with lists in columns
    all_data = pd.DataFrame({
        "IID": ["A"],
        "diagnoses": [["ICD10:F32", "ICD10:F33"]],
        "in_dates": [[pd.Timestamp("2020-01-01"), pd.Timestamp("2021-01-01")]],
        "DUD": [["DUD1"]],
        "DUD_In_Dates": [[pd.Timestamp("2019-12-01")]],
    })
    out = gp.update_DxDates_multi_exclusion(all_data.copy(),
                                           exclusion_type="1yprior",
                                           exc_diag="DUD",
                                           exc_diag_date="DUD_In_Dates",
                                           exc_diag_inflicted_changes="DUD_Inflicted_changes",
                                           diag_excode=7,
                                           level2codes="Level2_diagnoses",
                                           level2dates="Level2_dates",
                                           level2datemodifiercodes="diagnoses_Level2_modifier",
                                           level2datemodifierdates="date_Level2_modifier",
                                           level2datemodifierDXs="disorder_Level2_modifier",
                                           date_format="%Y-%m-%d",
                                           verbose=False)
    assert isinstance(out, pd.DataFrame)



def create_test_h5(tmp_path):
    """Helper to create test HDF5 store with known data"""
    h5path = str(tmp_path / "test.h5")
    
    # Create test data with various formats
    df = pd.DataFrame({
        'iid': ['1', '2', '3', '4', '5', '1', '2'],  # some duplicates
        'diagnosis': ['ICD10:F32', 'ICD10:F33', 'ICD10:F32.1', 'OTHER', 'ICD10:F32.0', 'ICD10:G40', 'ICD10:F32'],
        'extra': ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    })
    
    # Store with different dtypes to test robustness
    store = pd.HDFStore(h5path, mode='w')
    store.append('df', df, data_columns=True, min_itemsize={'diagnosis': 30})
    store.close()
    
    return h5path, df

def test_select_by_iid_and_diag_basic(tmp_path):
    """Test basic functionality with exact matches"""
    h5path, orig_df = create_test_h5(tmp_path)
    
    # Test 1: Exact ICD code match
    result = gp.select_by_iid_and_diag_optimized(
        h5_path=h5path,
        table_name='df',
        iidcol='iid',
        iids=None,  # all IIDs
        diagcol='diagnosis',
        diags=['ICD10:F32'],
        prefix_all=False
    )
    # Should find all rows with exactly ICD10:F32
    assert len(result) == 2
    assert all(d == 'ICD10:F32' for d in result['diagnosis'])

def test_select_by_iid_and_diag_prefix(tmp_path):
    """Test prefix matching"""
    h5path, orig_df = create_test_h5(tmp_path)
    
    # Test prefix matching (should match F32, F32.0, F32.1)
    result = gp.select_by_iid_and_diag_optimized(
        h5_path=h5path,
        table_name='df',
        iidcol='iid',
        iids=None,
        diagcol='diagnosis',
        diags=['ICD10:F32'],
        prefix_all=True  # Enable prefix matching
    )
    
    # Should find all F32* codes
    assert len(result) == 4
    assert all(d.startswith('ICD10:F32') for d in result['diagnosis'])

def test_select_by_iid_and_diag_combined(tmp_path):
    """Test combined IID and diagnosis filtering"""
    h5path, orig_df = create_test_h5(tmp_path)
    
    # Test combined IID and diagnosis filtering
    result = gp.select_by_iid_and_diag_optimized(
        h5_path=h5path,
        table_name='df',
        iidcol='iid',
        iids=['1', '2'],  # Only these IIDs
        diagcol='diagnosis',
        diags=['ICD10:F32'],
        prefix_all=True
    )
    
    # Should find F32* codes only for IIDs 1 and 2
    assert len(result) == 2
    assert all(i in ['1', '2'] for i in result['iid'])
    assert all(d.startswith('ICD10:F32') for d in result['diagnosis'])

def test_select_by_iid_and_diag_numeric_iids(tmp_path):
    """Test with numeric IIDs to catch type conversion issues"""
    h5path, _ = create_test_h5(tmp_path)
    
    # Test with numeric IIDs
    result = gp.select_by_iid_and_diag_optimized(
        h5_path=h5path,
        table_name='df',
        iidcol='iid',
        iids=[1, 2],  # Numeric IIDs
        diagcol='diagnosis',
        diags=['ICD10:F32'],
        prefix_all=False
    )
    
    assert len(result) > 0
    assert all(int(i) in [1, 2] for i in result['iid'])

def test_select_by_iid_and_diag_empty_results(tmp_path):
    """Test handling of no matches"""
    h5path, _ = create_test_h5(tmp_path)
    
    # Should return empty DataFrame with correct columns
    result = gp.select_by_iid_and_diag_optimized(
        h5_path=h5path,
        table_name='df',
        iidcol='iid',
        iids=['999'],  # Non-existent IID
        diagcol='diagnosis',
        diags=['NON:EXISTENT'],
        prefix_all=False
    )
    
    assert len(result) == 0
    assert 'iid' in result.columns
    assert 'diagnosis' in result.columns


def test_preflight_input_files_reads_headers_and_reports_missing(tmp_path):
    stam = tmp_path / "stam.tsv"
    stam.write_text("iid\tbirthdate\tsex\n1\t2000-01-01\tF\n", encoding="utf-8")

    headers = gp.preflight_input_files([
        ("stam", str(stam), "\t", True),
        ("optional", "", "\t", False),
    ])

    assert headers["stam"][0]["columns"] == ["iid", "birthdate", "sex"]
    assert headers["optional"] == []

    with pytest.raises(FileNotFoundError):
        gp.preflight_input_files([("missing", str(tmp_path / "missing.tsv"), "\t", True)])


def test_read_input_header_uses_explicit_engine_for_regex_separator(tmp_path):
    input_file = tmp_path / "input.tsv"
    input_file.write_text("iid\tstatus\n1\tactive\n", encoding="utf-8")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        columns = gp._read_input_header(str(input_file), sep=r"\t")

    assert columns == ["iid", "status"]
    assert not any(isinstance(item.message, pd.errors.ParserWarning) for item in caught)


def test_load_phenotypes_does_not_emit_applymap_futurewarning(tmp_path):
    request = tmp_path / "phenotypes.tsv"
    request.write_text(" MDD \t F32,F33 \n", encoding="utf-8")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        phenotypes, _, _ = gp.load_phenotypes(
            pheno_request=str(request),
            pheno_requestcol="diagnosis",
            skip_icd_update=True,
        )

    assert phenotypes.iloc[0].tolist() == ["MDD", ["F32", "F33"]]
    assert not any(
        isinstance(item.message, FutureWarning) and "applymap" in str(item.message)
        for item in caught
    )


def test_contains_atc_codes_detects_nested_main_and_exclusion_values(tmp_path):
    exclusion_file = tmp_path / "exclusions.tsv"
    exclusion_file.write_text("Disorder\tDisorder Codes\nAUD\tICD10:F10\nMED\tmain=ATC:N06A;sub=ICD10:F32\n", encoding="utf-8")

    assert gp.contains_atc_codes(pd.Series(["ICD10:F32", "ATC:N06A"]))
    assert gp.contains_atc_codes(pd.DataFrame({"Disorder Codes": ["ICD10:F32", "main=ATC:N05;sub=F10"]}))
    assert gp.file_contains_atc_codes(str(exclusion_file))
    assert not gp.contains_atc_codes(["ICD10:F32", "ICD9:303"])


def test_batch_helpers_limit_every_frame_to_current_iid_universe():
    df3 = pd.DataFrame({"iid": ["1", "2.0", "3", "4"], "birthdate": pd.to_datetime(["2000-01-01"] * 4)})
    df1 = pd.DataFrame({"iid": ["1", "2", "2", "5"], "diagnosis": ["A", "B", "C", "D"]})

    normalized_df3, iid_list, batch_size, num_batches = gp.prepare_iid_batches(df3, "iid", 2)
    assert normalized_df3["iid"].tolist() == ["1", "2", "3", "4"]
    assert iid_list == ["1", "2", "3", "4"]
    assert batch_size == 2
    assert num_batches == 2

    batch = iid_list[:2]
    assert gp.subset_to_iid_batch(df1, "iid", batch)["iid"].tolist() == ["1", "2", "2"]

    unique_df3 = gp.enforce_iid_universe(df3, "iid", batch, frame_name="df3", unique_iids=True)
    assert unique_df3["iid"].tolist() == ["1", "2"]


def test_fill_casecontrol_status_columns_marks_missing_phenotype_values_as_controls():
    df = pd.DataFrame({
        "iid": ["1", "2", "3"],
        "MDD": ["Case", np.nan, ""],
        "OCD": [np.nan, "Case", ""],
        "MDD_Codes": ["F32", "", ""],
    })

    out = gp.fill_casecontrol_status_columns(df, ["MDD", "OCD", "BPD"])

    assert out["MDD"].tolist() == ["Case", "Control", "Control"]
    assert out["OCD"].tolist() == ["Control", "Case", "Control"]
    assert out["BPD"].tolist() == ["Control", "Control", "Control"]
    assert out["MDD_Codes"].tolist() == ["F32", "", ""]


def test_cohort_membership_mask_supports_new_chb_and_old_degen_columns():
    df = pd.DataFrame({
        "iid": ["1", "2", "3", "4", "5"],
        "chb": ["TRUE", "FALSE", "", np.nan, "0"],
        "degen": ["FALSE", "TRUE", "", "", ""],
        "degen_new": ["FALSE", "FALSE", "TRUE", "", ""],
        "degen_old": ["FALSE", "FALSE", "FALSE", "1", ""],
    })

    mask = gp.cohort_membership_mask(df, ["chb", "degen", "degen_new", "degen_old"])

    assert mask.tolist() == [True, True, True, True, False]


def test_add_chb_degen_compat_columns_aliases_new_and_old_names():
    chb_only = gp.add_chb_degen_compat_columns(pd.DataFrame({"chb": ["TRUE", "FALSE"]}))
    degen_only = gp.add_chb_degen_compat_columns(pd.DataFrame({"degen": ["TRUE", "FALSE"]}))

    assert chb_only["degen"].tolist() == ["TRUE", "FALSE"]
    assert degen_only["chb"].tolist() == ["TRUE", "FALSE"]


def test_group_cohort_columns_keeps_only_present_columns_together():
    cols = ["iid", "birthdate", "chb", "sex", "dbds", "degen", "diagnosis"]

    assert gp.group_cohort_columns(cols) == [
        "iid", "birthdate", "chb", "dbds", "degen", "sex", "diagnosis"
    ]
    assert gp.group_cohort_columns(["iid", "birthdate", "sex"]) == ["iid", "birthdate", "sex"]


def test_expected_output_schema_uses_stam_first_and_source_specific_columns():
    pheno = pd.DataFrame({
        "Disorder": ["AUD", "MED", "MIXED"],
        "Disorder Codes": ["ICD10:F10", "ATC:N06A", "ICD10:F32,ATC:N05"],
    })
    lifetime = pd.DataFrame({"Disorder": ["EXC"], "Disorder Codes": ["ATC:N07"]})

    cols = gp.build_expected_output_columns(
        iidcol="iid",
        df3_columns=["iid", "birthdate", "sex"],
        df4_columns=["extra_info"],
        in_pheno_codes=pheno,
        lifetime_exclusions=lifetime,
        oneYearPrior_exclusions=pd.DataFrame(),
        post_exclusions=pd.DataFrame(),
        covariates=pd.DataFrame(),
        cluster_run="CHB_DBDS",
        ctype_col="",
        sexcol="sex",
        main_pheno_name="MainPheno",
    )

    assert cols[:4] == ["iid", "birthdate", "sex", "extra_info"]
    assert cols.index("diagnosis") < cols.index("AUD")
    assert "AUD_admissiontype" in cols
    assert "AUD_apk" not in cols
    assert "MED_apk" in cols
    assert "MED_admissiontype" not in cols
    assert "MIXED_admissiontype" in cols
    assert "MIXED_apk" in cols
    assert "EXC_apk" in cols


def test_expected_output_schema_does_not_add_missing_cohort_aliases():
    pheno = pd.DataFrame({"Disorder": ["MDD"], "Disorder Codes": ["ICD10:F32"]})

    cols = gp.build_expected_output_columns(
        iidcol="iid",
        df3_columns=["iid", "birthdate", "chb", "sex", "dbds"],
        df4_columns=[],
        in_pheno_codes=pheno,
        lifetime_exclusions=pd.DataFrame(),
        oneYearPrior_exclusions=pd.DataFrame(),
        post_exclusions=pd.DataFrame(),
        covariates=pd.DataFrame(),
        cluster_run="CHB_DBDS",
        ctype_col="",
        sexcol="sex",
        main_pheno_name="MainPheno",
    )

    assert cols[:5] == ["iid", "birthdate", "chb", "dbds", "sex"]
    assert "degen" not in cols
    assert "degen_new" not in cols
    assert "degen_old" not in cols


def test_coalesce_suffixed_columns_prefers_new_merge_values():
    df = pd.DataFrame({
        "iid": ["1", "2"],
        "sex_x": ["old", "old"],
        "sex_y": ["F", None],
        "birthdate": [None, "1990-01-01"],
        "birthdate_y": ["2000-01-01", "1990-01-01"],
    })

    out = gp.coalesce_suffixed_columns(df, ["sex", "birthdate"])

    assert out["sex"].tolist() == ["F", None]
    assert out["birthdate"].tolist() == ["2000-01-01", "1990-01-01"]
    assert "sex_x" not in out.columns
    assert "sex_y" not in out.columns
    assert "birthdate_y" not in out.columns


def test_fillna_blank_and_sanitize_tsv_cells_handle_arrays_and_newlines():
    df = pd.DataFrame({
        "text": ["a\tb", "c\nd", None],
        "items": [["x", "y\nz"], np.array(["a\tb"]), None],
        "num": [1.0, np.nan, 3.0],
    })

    filled = gp.fillna_blank(df)
    assert filled.loc[2, "text"] == ""
    assert filled.loc[1, "num"] == ""

    sanitized = gp.sanitize_tsv_cells(df)
    assert sanitized.loc[0, "text"] == "a b"
    assert sanitized.loc[1, "text"] == "c d"
    assert "\n" not in sanitized.loc[1, "items"]
    assert "\t" not in sanitized.loc[0, "items"]


def test_checksum_file_is_nonempty_and_uses_sha256_format(tmp_path):
    output = tmp_path / "result.tsv"
    output.write_text("iid\tdiagnosis\n1\tCase\n", encoding="utf-8")
    fgwa = tmp_path / "result.tsv.fgwa.pheno"
    fgwa.write_text("FID\tIID\tCaseControl\n1\t1\t1\n", encoding="utf-8")

    checksum_file = gp.checksum_path_for_output(str(output))
    targets = gp.collect_existing_checksum_targets(str(output), write_fastGWA_format=True)
    gp.write_checksum_file(targets, checksum_file)

    lines = (tmp_path / "result.checksums.sha256").read_text(encoding="utf-8").splitlines()
    expected_digest = hashlib.sha256(output.read_bytes()).hexdigest()
    assert lines
    assert lines[0] == f"{expected_digest}  result.tsv"
    assert any(line.endswith("  result.tsv.fgwa.pheno") for line in lines)


def test_generate_readme_uses_actual_checksum_file_name():
    readme = gp.generate_readme(
        flags_used="-g: file pheno.tsv\n--iidcol: value IID\n--sexcol: value sex\n--bdcol: value birthdate\n--iidstatus: value \n--iidstatusdate: value ",
        default_args="",
        checksum_file="result.checksums.sha256",
    )

    assert "sha256sum -c ./result.checksums.sha256" in readme
    assert "sha256sum -c ./checksums.sha256" not in readme
# Note: main() and other CLI/integration-heavy functions are not invoked here.
# The tests above exercise logic-heavy code paths and ensure basic functionality.
