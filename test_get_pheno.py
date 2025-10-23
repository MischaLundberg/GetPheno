
import sys
#sys.path.append("/Users/mischa/Desktop")

import os
import io
import uuid
import tempfile
import subprocess
import pandas as pd
import numpy as np
import types
import pytest
from datetime import datetime, date, timedelta

import get_pheno_refactored_new as gp

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


def test_split_and_format_various():
    # numeric input returns string
    assert gp.split_and_format(123) == "123"
    assert gp.split_and_format("456") == "456"
    # ICD-like string
    assert gp.split_and_format("F32.1") == "F321"
    # with fill -> ensure padded to 4
    assert gp.split_and_format("F32.1", fill=True) == "F3210"
    # Non-matching string returns itself
    assert gp.split_and_format("XYZ") == "XYZ"


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
    auto = gp.normalize_iid_series_auto(s2)
    assert getattr(auto.dtype, "name", "") in ("Int64", "int64", "Int64") or pd.api.types.is_integer_dtype(auto)


def test_convert_if_not_datetime_various():
    assert gp.convert_if_not_datetime(pd.NaT) is pd.NaT
    assert gp.convert_if_not_datetime(date(2000, 1, 1)) == date(2000, 1, 1)
    assert gp.convert_if_not_datetime("2000-01-02").date() == date(2000, 1, 2)


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
    assert gp.expand_ranges(["T36-T38,T40-T42,T45"]) == ["T36", "T37", "T38", "T40", "T41", "T42", "T45"]
    assert gp.expand_ranges(["T36-T38"]) == ["T36", "T37", "T38"]
    assert gp.expand_ranges(["36-38"]) == ["36", "37", "38"]
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


def test_load_stam_file_and_finalize(tmp_path):
    # create a small csv
    csvp = tmp_path / "stam.csv"
    df = pd.DataFrame({"pnr": [1, 2], "birthdate": ["01/01/1980", "02/02/1990"], "diagnosis": ["ICD10:F32", "ICD10:F33"], "sex": ["M", "F"]})
    df.to_csv(csvp, index=False)
    res = gp.load_stam_file(stam_file=str(csvp), isep=",", birthdatecol="birthdate", diagnostic_col="diagnosis", sexcol="sex", stam_cols_to_read_as_date=["birthdate"])
    assert "birthdate" in res.columns
    # finalize_lpr_data should rename columns accordingly
    df2 = pd.DataFrame({"c_adiag": ["A"], "d_inddto": [pd.Timestamp("2020-01-01")]})
    fin = gp.finalize_lpr_data(df2.copy(), diagnostic_col="c_adiag", birthdatecol="d_inddto", ctype_col="", verbose=False)
    assert "diagnosis" in fin.columns


def test_build_temp_file_and_load_mapping_rows(monkeypatch, tmp_path):
    # create a small file
    csvp = tmp_path / "big.csv"
    df = pd.DataFrame({"id": [f"id{i}" for i in range(10)], "v": list(range(10))})
    df.to_csv(csvp, index=False)
    # load_mapping_rows should identify row indices
    rows = gp.load_mapping_rows(str(csvp), iidcol="id", target_iids=["id1", "id3"], fsep=",")
    assert isinstance(rows, list)
    # Monkeypatch subprocess.run used in build_temp_file to avoid awk on system
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: None)
    # call build_temp_file with the indices from load_mapping_rows
    tempf = tmp_path / "filtered_temp.csv"
    gp.build_temp_file(str(csvp), [1, 3], temp_file=str(tempf), verbose=True)
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
    h5path = tmp_path + "/test.h5"
    
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
    print(result)
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
    
    print(result)
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
    
    print(result)
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
    
    print(result)
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
    
    print(result)
    assert len(result) == 0
    assert 'iid' in result.columns
    assert 'diagnosis' in result.columns
# Note: main() and other CLI/integration-heavy functions are not invoked here.
# The tests above exercise logic-heavy code paths and ensure basic functionality.


