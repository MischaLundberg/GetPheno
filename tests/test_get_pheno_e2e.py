import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import get_pheno as gp

FIXTURE_ROOT = REPO_ROOT / "tests" / "fixtures" / "e2e"


@pytest.mark.parametrize("flag,expected", [("--version", f"GetPheno {gp.__version__}"),
                                          ("--help", "usage: get_pheno.py")])
def test_cli_information_flags_need_no_input_arguments(tmp_path, flag, expected):
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "get_pheno.py"), flag],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert expected in result.stdout
    assert not list(tmp_path.iterdir())


def _read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _copy_fixture_inputs(case_name, tmp_path):
    src = FIXTURE_ROOT / case_name / "input"
    dst = tmp_path / "input"
    shutil.copytree(src, dst)
    return dst


def _run_getpheno(
    input_dir,
    output_path,
    *,
    lowmem=False,
    batchsize=100000,
    ipsych_options=None,
    additional_file="",
    stam_files="stam.csv",
    qced_file="",
    build_entry_exit=False,
    status_col="",
    status_date_col="",
    diagnosis_files="lpr.csv",
    diagnostic_col="diagnosis",
    iid_col="iid",
    date_in_col="date_in",
    date_out_col="date_out",
    diagnosis_type_col="c_types",
    ini_file=None,
    default_args=None,
    deathdate_col="",
    death_age_range="0,0",
    keep_alive=False,
    only_dead=False,
    write_omop_cohort=False,
    omop_person_map="",
    omop_definitions="",
    omop_output_dir="",
    phenotype_file="pheno.tsv",
    register_col="",
    run_psyk=False,
    run_lpr=False,
    verylowmem=False,
    date_format="%Y-%m-%d",
    exact_match=True,
    opr_options=None,
    lifetime_file="",
    post_file="",
    prior_file="",
    match_fi=False,
    type_exclusions="",
    type_inclusions="",
):
    gp.main(
        ",".join(str(input_dir / path) for path in diagnosis_files.split(",")),
        str(input_dir / phenotype_file) if phenotype_file else "",
        ",".join(str(input_dir / path) for path in stam_files.split(",")),
        ",".join(str(input_dir / path) for path in additional_file.split(","))
        if additional_file else "",
        False,
        "",
        diagnostic_col,
        "diagnosis",
        iid_col,
        "birthdate",
        "sex",
        ",",
        ",",
        ",",
        "\t",
        str(output_path),
        exact_match,
        date_in_col,
        date_out_col,
        str(input_dir / qced_file) if qced_file else "",
        type_exclusions,
        type_inclusions,
        diagnosis_type_col,
        str(input_dir / lifetime_file) if lifetime_file else "",
        str(input_dir / post_file) if post_file else "",
        str(input_dir / prior_file) if prior_file else "",
        False,
        "",
        "",
        False,
        False,
        False,
        match_fi,
        True,
        date_format,
        status_col,
        status_date_col,
        "",
        False,
        1,
        "MainPheno",
        build_entry_exit,
        False,
        False,
        False,
        False,
        [],
        [],
        "0,0",
        False,
        False,
        False,
        lowmem,
        verylowmem,
        batchsize,
        False,
        "",
        "",
        "",
        "diagnosis",
        "",
        "",
        "",
        [],
        run_lpr,
        run_psyk,
        ",",
        "",
        str(ini_file or output_path.parent / "missing.ini"),
        False,
        False,
        False,
        False,
        "",
        "",
        "",
        default_args or [],
        "",
        ipsych_options,
        deathdate_col,
        death_age_range,
        keep_alive,
        only_dead,
        write_omop_cohort,
        str(omop_person_map) if omop_person_map else "",
        "iid",
        "person_id",
        str(omop_definitions) if omop_definitions else "",
        str(omop_output_dir) if omop_output_dir else "",
        "\t",
        register_col=register_col,
        opr_options=opr_options,
    )


@pytest.mark.parametrize("exact_match", [False, True])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("lowmem", [False, True])
def test_internal_wildcard_phenotype_matches_in_all_input_modes(tmp_path, monkeypatch, exact_match, indexed, lowmem):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gp.socket, "gethostname", lambda: "local-test")
    iids = ["001", "002", "003", "004", "005", "006"]
    pd.DataFrame({"iid": iids, "birthdate": ["1985-01-01"] * 6, "sex": ["F"] * 6}).to_csv(
        tmp_path / "stam.csv", index=False)
    pd.DataFrame({
        "iid": iids,
        "diagnosis": ["ICD8:297.29", "ICD8:297.89", "ICD8:297.291", "ICD8:297X29", "ICD8:298.29", "ICD8:297.9"],
        "date_in": ["2015-01-01"] * 6, "date_out": ["2015-01-02"] * 6,
    }).to_csv(tmp_path / "lpr.csv", index=False)
    (tmp_path / "pheno.tsv").write_text("MATCH\tICD8:297.*9\n")
    if indexed:
        gp.index_diag_file(str(tmp_path / "lpr.csv"), index_columns=["iid", "diagnosis"], separator=",")
    output = tmp_path / "wildcards.tsv"
    _run_getpheno(tmp_path, output, exact_match=exact_match, lowmem=lowmem, batchsize=2)
    result = _read_tsv(output).set_index("iid")
    assert result.index.tolist() == iids
    assert set(result.index[result["MATCH"] == "Case"]) == {"001", "002", "006"}


@pytest.mark.parametrize("exact_match", [False, True])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("lowmem", [False, True])
def test_wildcard_category_ranges_include_subcodes_in_all_modes(tmp_path, monkeypatch, exact_match, indexed, lowmem):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gp.socket, "gethostname", lambda: "local-test")
    iids = ["001", "002", "003", "004", "005", "006"]
    pd.DataFrame({"iid": iids, "birthdate": ["1985-01-01"] * 6, "sex": ["F"] * 6}).to_csv(
        tmp_path / "stam.csv", index=False)
    pd.DataFrame({
        "iid": iids,
        "diagnosis": ["ICD10:F20", "ICD10:F20.1", "ICD10:F29.8", "ICD10:F21.1", "ICD10:F30", "ICD10:F29"],
        "date_in": ["2015-01-01"] * 6, "date_out": ["2015-01-02"] * 6,
    }).to_csv(tmp_path / "lpr.csv", index=False)
    (tmp_path / "pheno.tsv").write_text("PLAIN\tICD10:F20-F29\nWILDCARD\tICD10:F20*-F29*\n")
    if indexed:
        gp.index_diag_file(str(tmp_path / "lpr.csv"), index_columns=["iid", "diagnosis"], separator=",")
    output = tmp_path / "ranges.tsv"
    _run_getpheno(tmp_path, output, exact_match=exact_match, lowmem=lowmem, batchsize=2)
    result = _read_tsv(output).set_index("iid")
    expected = {"001", "002", "003", "004", "006"}
    assert set(result.index[result["WILDCARD"] == "Case"]) == expected
    assert set(result.index[result["PLAIN"] == "Case"]) == ({"001", "006"} if exact_match else expected)


@pytest.mark.parametrize("exact_match", [False, True])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("lowmem", [False, True])
def test_shared_decimal_range_matches_in_all_input_modes(tmp_path, monkeypatch, exact_match, indexed, lowmem):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gp.socket, "gethostname", lambda: "local-test")
    iids = [f"{i:03d}" for i in range(1, 19)]
    pd.DataFrame({"iid": iids, "birthdate": ["1985-01-01"] * 18, "sex": ["F"] * 18}).to_csv(
        tmp_path / "stam.csv", index=False)
    pd.DataFrame({
        "iid": iids,
        "diagnosis": [f"ICD8:{i}.39" for i in range(300, 306)]
                     + ["ICD8:300.38", "ICD8:305.390", "ICD8:306.39", "ICD8:299.39",
                        "ICD8:298.29", "ICD8:298.30", "ICD8:298.99", "ICD8:298.28", "ICD8:298.999",
                        "ICD8:303.29", "ICD8:303.291", "ICD8:303.28"],
        "date_in": ["2015-01-01"] * 18, "date_out": ["2015-01-02"] * 18,
    }).to_csv(tmp_path / "lpr.csv", index=False)
    (tmp_path / "pheno.tsv").write_text(
        "PLAIN\tICD8:300.39-305.39\nWILDCARD\tICD8:300.39*-305.39*\n"
        "INTERVAL\tICD8:298.29-298.99\nINTERVAL_PREFIX\tICD8:298.29*-298.99*\nSUFFIX\tICD8:303.*9\n"
    )
    if indexed:
        gp.index_diag_file(str(tmp_path / "lpr.csv"), index_columns=["iid", "diagnosis"], separator=",")
    output = tmp_path / "decimals.tsv"
    _run_getpheno(tmp_path, output, exact_match=exact_match, lowmem=lowmem, batchsize=3)
    result = _read_tsv(output).set_index("iid")
    assert result.index.tolist() == iids
    exact_cases = set(iids[:6])
    prefix_cases = exact_cases | {"008"}
    assert set(result.index[result["PLAIN"] == "Case"]) == (exact_cases if exact_match else prefix_cases)
    assert set(result.index[result["WILDCARD"] == "Case"]) == prefix_cases
    interval_cases = {"011", "012", "013"}
    assert set(result.index[result["INTERVAL"] == "Case"]) == (interval_cases if exact_match else interval_cases | {"015"})
    assert set(result.index[result["INTERVAL_PREFIX"] == "Case"]) == interval_cases | {"015"}
    assert set(result.index[result["SUFFIX"] == "Case"]) == {"004", "016"}


def _write_opr_inputs(tmp_path):
    iids = [f"{i:03d}" for i in range(1, 7)]
    pd.DataFrame({"iid": iids, "birthdate": ["1985-01-01"] * 6, "sex": ["F"] * 6}).to_csv(
        tmp_path / "stam.csv", index=False)
    pd.DataFrame({
        "iid": ["001", "001", "002", "004", "006"],
        "diagnosis": ["ICD10:F32", "KAAA001", "ICD10:G40", "OPR:KBBB001", "ICD10:F32"],
        "date_in": ["2015-01-01"] * 5, "date_out": ["2015-01-02"] * 5,
    }).to_csv(tmp_path / "lpr.csv", index=False)
    pd.DataFrame({
        "person": ["002", "003", "004", "999"],
        "procedure": ["kaaa001", "OPR:KAAA001", "KBBB002", "KAAA001"],
        "performed": ["2016-06-01", "2017-07-01", "2018-08-01", "2015-01-01"],
    }).to_csv(tmp_path / "opr.tsv", sep="\t", index=False)
    (tmp_path / "pheno.tsv").write_text(
        "DX\tICD10:F32\nOP\tOPR:KAAA001\nOP_PREFIX\tOPR:KBBB*\nMIX\tICD10:F32,OPR:KAAA001\n")
    return {"file": str(tmp_path / "opr.tsv"), "code_col": "procedure", "date_col": "performed",
            "iid_col": "person", "sep": "\t"}


@pytest.mark.parametrize("lowmem", [False, True])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("separate_opr", [False, True])
def test_opr_conditional_source_and_fallback_across_input_modes(tmp_path, monkeypatch, lowmem, indexed, separate_opr):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gp.socket, "gethostname", lambda: "local-test")
    options = _write_opr_inputs(tmp_path)
    if indexed:
        gp.index_diag_file(str(tmp_path / "lpr.csv"), index_columns=["iid", "diagnosis"], separator=",")
    output = tmp_path / "procedures.tsv"
    _run_getpheno(tmp_path, output, opr_options=options if separate_opr else None,
                  lowmem=lowmem, batchsize=2)
    result = _read_tsv(output).set_index("iid")
    assert result.index.tolist() == [f"{i:03d}" for i in range(1, 7)]
    assert set(result.index[result["DX"] == "Case"]) == {"001", "006"}
    assert set(result.index[result["OP"] == "Case"]) == ({"002", "003"} if separate_opr else {"001"})
    assert set(result.index[result["OP_PREFIX"] == "Case"]) == {"004"}
    assert set(result.index[result["MIX"] == "Case"]) == ({"001", "002", "003", "006"} if separate_opr else {"001", "006"})
    assert result.loc["002" if separate_opr else "001", "OP_earliest_date"].startswith("2016-06-01" if separate_opr else "2015-01-01")


def test_opr_inputs_are_not_opened_without_opr_requests(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    options["file"] = str(tmp_path / "unused-missing.tsv")
    (tmp_path / "pheno.tsv").write_text("DX\tICD10:F32\n")
    monkeypatch.setattr(gp, "configure_opr_input", lambda *args: pytest.fail("Unused OPR source was opened"))
    output = tmp_path / "unused.tsv"
    _run_getpheno(tmp_path, output, opr_options=options)
    assert (_read_tsv(output)["DX"] == "Case").sum() == 2


def test_requested_opr_source_must_exist_and_have_required_columns(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    options["file"] = str(tmp_path / "missing.tsv")
    with pytest.raises(FileNotFoundError, match="opr_file"):
        _run_getpheno(tmp_path, tmp_path / "missing-output.tsv", opr_options=options)
    options["file"] = str(tmp_path / "opr.tsv")
    options["date_col"] = "missing_date"
    with pytest.raises(ValueError, match="missing_date"):
        _run_getpheno(tmp_path, tmp_path / "missing-column.tsv", opr_options=options)


def test_opr_matchfi_includes_procedure_only_iids(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    output = tmp_path / "matchfi.tsv"
    _run_getpheno(tmp_path, output, opr_options=options, match_fi=True)
    assert "003" in set(_read_tsv(output)["iid"])


def test_opr_multiple_files_can_use_different_columns(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    second = tmp_path / "opr2.tsv"
    pd.DataFrame({"other_iid": ["003"], "code": ["KBBB099"], "date": ["2019-01-01"]}).to_csv(
        second, sep="\t", index=False)
    options.update(file=options["file"] + "," + str(second), code_col="procedure,code",
                   iid_col="person,other_iid", date_col="performed,date")
    output = tmp_path / "multiple.tsv"
    _run_getpheno(tmp_path, output, opr_options=options)
    result = _read_tsv(output).set_index("iid")
    assert set(result.index[result["OP_PREFIX"] == "Case"]) == {"003", "004"}


def test_opr_file_can_be_stata(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    frame = pd.read_csv(options["file"], sep="\t", dtype=str)
    frame["performed"] = pd.to_datetime(frame["performed"])
    path = tmp_path / "procedures.dta"
    frame.to_stata(path, write_index=False)
    options["file"] = str(path)
    output = tmp_path / "stata.tsv"
    _run_getpheno(tmp_path, output, opr_options=options)
    assert set(_read_tsv(output).query("OP == 'Case'")["iid"]) == {"002", "003"}


def test_opr_in_exclusion_file_triggers_procedure_loading(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    (tmp_path / "pheno.tsv").write_text("DX\tICD10:F32\n")
    (tmp_path / "lifetime.tsv").write_text("SURG\tOPR:KAAA001\n")
    pd.DataFrame({"person": ["006"], "procedure": ["KAAA001"], "performed": ["2014-01-01"]}).to_csv(
        options["file"], sep="\t", index=False)
    output = tmp_path / "exclusion.tsv"
    _run_getpheno(tmp_path, output, opr_options=options, lifetime_file="lifetime.tsv")
    result = _read_tsv(output).set_index("iid")
    assert "OPR:KAAA001" in result.loc["006", "SURG"]
    assert result.loc["006", "SURG_earliest_date"].startswith("2014-01-01")


def test_opr_combines_with_separate_atc_input(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    (tmp_path / "pheno.tsv").write_text("MULTI\tICD10:F32,ATC:N06A*,OPR:KAAA001\n")
    pd.DataFrame({"iid": ["005"], "drug": ["ATC:N06AB"], "prescribed": ["2016-01-01"]}).to_csv(
        tmp_path / "atc.csv", index=False)
    ini = tmp_path / "medication.ini"
    ini.write_text(f"[Default]\natc_file = {tmp_path / 'atc.csv'}\natc_diag_col = drug\natc_date_col = prescribed\n")
    output = tmp_path / "mixed.tsv"
    _run_getpheno(tmp_path, output, opr_options=options, ini_file=ini,
                  default_args=["--atc", "--atccol", "--atcdatecol"])
    assert set(_read_tsv(output).query("MULTI == 'Case'")["iid"]) == {"001", "002", "003", "005", "006"}


@pytest.mark.parametrize('apply_filters', [False, True])
def test_opr_cli_ini_and_explicit_column_override(tmp_path, apply_filters):
    options = _write_opr_inputs(tmp_path)
    ini = tmp_path / "procedure.ini"
    ini.write_text(f"[Default]\nopr_file = {options['file']}\nopr_diag_col = wrong\n"
                   "opr_iid_col = person\nopr_date_col = performed\nopr_sep = \\t\n")
    output = tmp_path / "cli.tsv"
    if apply_filters:
        ini.write_text(ini.read_text() + 'opr_filters = person=003\n')
    filter_args = ['--oprFilter', 'person=002', '--oprFilter', 'procedure=KAAA*',
                   '--oprExclude', 'person=003'] if apply_filters else []
    completed = subprocess.run([
        sys.executable, str(REPO_ROOT / "get_pheno.py"), "--ini", str(ini),
        "-i", str(tmp_path / "stam.csv"), "-f", str(tmp_path / "lpr.csv"),
        "-g", str(tmp_path / "pheno.tsv"), "-o", str(output), "--iidcol", "iid",
        "--fcol", "diagnosis", "--din", "date_in", "--don", "date_out",
        "--DateFormat", "%Y-%m-%d", "--oprcol", "procedure", "--eM",
    ] + filter_args, cwd=tmp_path, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr + completed.stdout[-2000:]
    assert set(_read_tsv(output).query("OP == 'Case'")["iid"]) == ({'002'} if apply_filters else {'002', '003'})


@pytest.mark.parametrize("indexed", [False, True])
def test_empty_opr_input_produces_controls(tmp_path, monkeypatch, indexed):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    (tmp_path / "pheno.tsv").write_text("OP\tOPR:KAAA001\n")
    pd.DataFrame(columns=["person", "procedure", "performed"]).to_csv(options["file"], sep="\t", index=False)
    if indexed:
        gp.index_diag_file(str(tmp_path / "lpr.csv"), index_columns=["iid", "diagnosis"], separator=",")
    output = tmp_path / "empty.tsv"
    _run_getpheno(tmp_path, output, opr_options=options)
    result = _read_tsv(output)
    assert result["OP"].tolist() == ["Control"] * 6


def test_opr_records_respect_qc_iid_universe(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    (tmp_path / "qced.txt").write_text("002\n")
    output = tmp_path / "qc.tsv"
    _run_getpheno(tmp_path, output, opr_options=options, qced_file="qced.txt")
    result = _read_tsv(output)
    assert result["iid"].tolist() == ["002"]
    assert result["OP"].tolist() == ["Case"]


def test_separate_opr_events_respect_comparison_date_window(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    frame = pd.read_csv(options["file"], sep="\t", dtype=str)
    frame.loc[frame["person"] == "002", "performed"] = "2015-12-31"
    frame.to_csv(options["file"], sep="\t", index=False)
    output = tmp_path / "window.tsv"
    _run_getpheno(tmp_path, output, opr_options=options, ipsych_options={"enabled": True})
    assert set(_read_tsv(output).query("OP == 'Case'")["iid"]) == {"002"}


def test_opr_compound_rule_can_exclude_diagnosis_cases(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    frame = pd.read_csv(tmp_path / "lpr.csv", dtype=str)
    frame.loc[frame["iid"] == "002", "diagnosis"] = "ICD10:F32"
    frame["c_types"] = "A"
    frame.to_csv(tmp_path / "lpr.csv", index=False)
    (tmp_path / "pheno.tsv").write_text("OP\tmain=OPR:KAAA001;rule_out=ICD10:F32\n")
    output = tmp_path / "compound.tsv"
    _run_getpheno(tmp_path, output, opr_options=options)
    assert set(_read_tsv(output).query("OP == 'Case'")["iid"]) == {"003"}


@pytest.mark.parametrize('separate', [False, True])
@pytest.mark.parametrize('indexed,lowmem', [(False, False), (True, True)])
@pytest.mark.parametrize('explicit_filters', [False, True])
@pytest.mark.parametrize('source_flag', ['PSYK', 'LPR'])
def test_opr_filter_policy_is_identical_in_both_input_routes(
        tmp_path, monkeypatch, separate, indexed, lowmem, explicit_filters, source_flag):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gp.socket, 'gethostname', lambda: 'local-test')
    monkeypatch.setattr(gp, 'DK_clusters', (*gp.DK_clusters, 'Default'))
    options = _write_opr_inputs(tmp_path)
    procedures = pd.read_csv(options['file'], sep='\t', dtype=str).iloc[:2].copy()
    procedures['procedure_origin'] = ['somatic', 'lpr_procedures']
    procedures['procedure_kind'] = ['M', 'B']
    procedures.to_csv(options['file'], sep='\t', index=False)
    primary = pd.read_csv(tmp_path / 'lpr.csv', dtype=str)
    primary = primary[~primary['diagnosis'].str.contains('KAAA001', na=False)].copy()
    primary['register'] = 'psyc_diag'
    primary['c_types'] = 'A'
    fallback = pd.DataFrame({
        'iid': procedures['person'], 'diagnosis': 'OPR:KAAA001',
        'date_in': procedures['performed'], 'date_out': procedures['performed'],
        'register': 'somatic', 'c_types': procedures['procedure_kind'],
        'procedure_origin': procedures['procedure_origin'],
        'procedure_kind': procedures['procedure_kind'],
    })
    pd.concat([primary, fallback], ignore_index=True).to_csv(tmp_path / 'lpr.csv', index=False)
    (tmp_path / 'pheno.tsv').write_text('OP\tmain=OPR:KAAA001;rule_out=ICD10:G99\n')
    if not separate:
        options['file'] = ''
    if explicit_filters:
        options.update(include=['procedure_origin=LPR*'], exclude=['procedure_kind=M'])
    if indexed:
        gp.index_diag_file(str(tmp_path / 'lpr.csv'), ['iid', 'diagnosis'], separator=',')
    output = tmp_path / 'opr-policy.tsv'
    _run_getpheno(tmp_path, output, opr_options=options,
                  run_psyk=source_flag == 'PSYK', run_lpr=source_flag == 'LPR',
                  type_exclusions='M', type_inclusions='A', lowmem=lowmem, batchsize=2)
    result = _read_tsv(output).set_index('iid')
    assert set(result.index[result['OP'] == 'Case']) == ({'003'} if explicit_filters else {'002', '003'})
    assert result.loc['003', 'OP_earliest_date'].startswith('2017-07-01')
    assert not any(col.startswith('_opr') for col in result.columns)


def test_opr_filters_can_be_loaded_from_ini(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    ini = tmp_path / 'filters.ini'
    ini.write_text('[Default]\nopr_filters = person=003;procedure=*KAAA*\nopr_exclusions = person=002\n')
    output = tmp_path / 'filtered.tsv'
    _run_getpheno(tmp_path, output, opr_options=options, ini_file=ini,
                  default_args=['--oprFilter', '--oprExclude'])
    assert set(_read_tsv(output).query("OP == 'Case'")['iid']) == {'003'}


@pytest.mark.parametrize('lowmem', [False, True])
@pytest.mark.parametrize('date_settings', ['inherit', 'in_only', 'ini', 'legacy_ini'])
def test_opr_in_out_dates_reach_output(tmp_path, monkeypatch, lowmem, date_settings):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    procedures = pd.read_csv(options['file'], sep='\t', dtype=str)
    procedures['ended'] = ['2016-06-04', None, '2018-08-02', '2017-01-01']
    options['date_col'] = ''
    ini = None
    defaults = []
    if date_settings in ('inherit', 'in_only'):
        procedures.rename(columns={'performed': 'date_in', 'ended': 'date_out'}, inplace=True)
        if date_settings == 'in_only':
            procedures.drop(columns='date_out', inplace=True)
    else:
        ini = tmp_path / 'date-columns.ini'
        settings = ('opr_date_col = wrong_legacy\nopr_date_in_col = performed\n'
                    if date_settings == 'ini' else 'opr_date_col = performed\nopr_date_in_col =\n')
        ini.write_text('[Default]\n' + settings + 'opr_date_out_col = ended\n')
        defaults = ['--oprdatecol', '--oprdateoutcol']
    procedures.to_csv(options['file'], sep='\t', index=False)
    output = tmp_path / 'procedure-dates.tsv'
    _run_getpheno(tmp_path, output, opr_options=options, ini_file=ini,
                  default_args=defaults, lowmem=lowmem, batchsize=2)
    result = _read_tsv(output).set_index('iid')
    assert '2016-06-01' in result.loc['002', 'OP_In_Dates']
    assert ('2016-06-01' if date_settings == 'in_only' else '2016-06-04') in result.loc['002', 'OP_Out_Dates']
    assert '2017-07-01' in result.loc['003', 'OP_Out_Dates']


@pytest.mark.parametrize('in_flag', ['--oprdatecol', '--oprdateincol'])
@pytest.mark.parametrize('inherit_blank', [False, True])
def test_opr_date_cli_aliases_override_ini(tmp_path, in_flag, inherit_blank):
    options = _write_opr_inputs(tmp_path)
    procedures = pd.read_csv(options['file'], sep='\t', dtype=str)
    procedures['ended'] = '2019-01-01'
    if inherit_blank:
        procedures.rename(columns={'performed': 'date_in', 'ended': 'date_out'}, inplace=True)
    procedures.to_csv(options['file'], sep='\t', index=False)
    ini = tmp_path / 'date-cli.ini'
    ini.write_text('[Default]\nopr_date_col = wrong_old\nopr_date_in_col = wrong_in\nopr_date_out_col = wrong_out\n')
    output = tmp_path / 'cli-dates.tsv'
    completed = subprocess.run([
        sys.executable, str(REPO_ROOT / 'get_pheno.py'), '--ini', str(ini),
        '-i', str(tmp_path / 'stam.csv'), '-f', str(tmp_path / 'lpr.csv'),
        '-g', str(tmp_path / 'pheno.tsv'), '-o', str(output), '--iidcol', 'iid',
        '--fcol', 'diagnosis', '--din', 'date_in', '--don', 'date_out', '--DateFormat', '%Y-%m-%d',
        '--opr', options['file'], '--oprcol', 'procedure', '--opriidcol', 'person', '--oprsep', '\t',
        in_flag, '' if inherit_blank else 'performed',
        '--oprdateoutcol', '' if inherit_blank else 'ended', '--eM',
    ], cwd=tmp_path, text=True, capture_output=True)
    assert completed.returncode == 0, completed.stderr + completed.stdout[-2000:]
    assert '2019-01-01' in _read_tsv(output).set_index('iid').loc['002', 'OP_Out_Dates']


@pytest.mark.parametrize('separate', [False, True])
def test_explicit_opr_filter_can_leave_all_controls(tmp_path, monkeypatch, separate):
    monkeypatch.chdir(tmp_path)
    options = _write_opr_inputs(tmp_path)
    (tmp_path / 'pheno.tsv').write_text('OP\tmain=OPR:KAAA001;rule_out=ICD10:G99\n')
    options['include'] = ['person=nobody'] if separate else ['iid=nobody']
    if not separate:
        options['file'] = ''
    output = tmp_path / 'no-procedures.tsv'
    _run_getpheno(tmp_path, output, opr_options=options)
    assert _read_tsv(output)['OP'].tolist() == ['Control'] * 6


@pytest.mark.parametrize('namespace', ['ICD10', 'ICD10-CM'])
@pytest.mark.parametrize('indexed', [False, True])
def test_zero_padded_wildcard_ranges_do_not_include_other_categories(tmp_path, monkeypatch, namespace, indexed):
    monkeypatch.chdir(tmp_path)
    _write_opr_inputs(tmp_path)
    frame = pd.read_csv(tmp_path / 'lpr.csv', dtype=str)
    frame['diagnosis'] = [f'{namespace}:{code}' for code in ['F001', 'F019', 'F021', 'F201', 'F99']]
    frame.to_csv(tmp_path / 'lpr.csv', index=False)
    (tmp_path / 'pheno.tsv').write_text(f'LOW\t{namespace}:F00*-F02*\n')
    if indexed:
        gp.index_diag_file(str(tmp_path / 'lpr.csv'), ['iid', 'diagnosis'], separator=',')
    output = tmp_path / 'low-categories.tsv'
    _run_getpheno(tmp_path, output, lowmem=indexed, batchsize=2)
    assert set(_read_tsv(output).query("LOW == 'Case'")['iid']) == {'001', '002'}


def assert_generated_output_matches_fixture(actual_path, expected_path):
    actual = _read_tsv(actual_path)
    expected = _read_tsv(expected_path)

    pd.testing.assert_frame_equal(actual, expected, check_dtype=False)
    assert actual.columns.tolist() == expected.columns.tolist()
    assert "degen" not in actual.columns


def test_main_basic_multi_phenotype_output_matches_fixture(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    output_path = tmp_path / "actual.tsv"

    _run_getpheno(input_dir, output_path)

    assert_generated_output_matches_fixture(
        output_path,
        FIXTURE_ROOT / "basic_multi" / "expected" / "output.tsv",
    )


def test_main_lowmem_multi_phenotype_output_matches_fixture(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    output_path = tmp_path / "actual.lowmem.tsv"

    _run_getpheno(input_dir, output_path, lowmem=True, batchsize=2)

    assert_generated_output_matches_fixture(
        output_path,
        FIXTURE_ROOT / "basic_multi" / "expected" / "output.tsv",
    )


@pytest.mark.parametrize("variant", ["normalized_iids", "spaces", "multiline"])
def test_standard_and_batched_input_variants_match(tmp_path, monkeypatch, variant):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    if variant == "spaces":
        renamed = tmp_path / "input with spaces"
        input_dir.rename(renamed)
        input_dir = renamed
    else:
        source = input_dir / "lpr.csv"
        records = pd.read_csv(source, dtype=str, keep_default_na=False)
        if variant == "normalized_iids":
            records["iid"] = " " + records["iid"] + ".0 "
        else:
            records["note"] = "quoted\nmultiline value"
        records.to_csv(source, index=False)
    for lowmem in (False, True):
        output = tmp_path / f"actual.{lowmem}.tsv"
        _run_getpheno(input_dir, output, lowmem=lowmem, batchsize=2)
        assert_generated_output_matches_fixture(
            output, FIXTURE_ROOT / "basic_multi" / "expected" / "output.tsv"
        )


def test_main_iid_intersection_mode(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    processor = gp.process_pheno_and_exclusions

    def intersect(*args, **kwargs):
        kwargs["MatchFI"] = True
        return processor(*args, **kwargs)

    monkeypatch.setattr(gp, "process_pheno_and_exclusions", intersect)
    output = tmp_path / "intersection.tsv"
    _run_getpheno(input_dir, output)
    expected = set(pd.read_csv(input_dir / "lpr.csv", dtype=str)["iid"])
    expected &= set(pd.read_csv(input_dir / "stam.csv", dtype=str)["iid"])
    assert set(_read_tsv(output)["iid"]) == expected


@pytest.mark.parametrize("failure", ["duplicates", "provenance", "publication"])
def test_failed_omop_export_preserves_existing_outputs(tmp_path, monkeypatch, failure):
    inputs = _copy_fixture_inputs("basic_multi", tmp_path)
    result = tmp_path / "result.tsv"
    shutil.copyfile(FIXTURE_ROOT / "basic_multi" / "expected" / "output.tsv", result)
    output_dir = tmp_path / "omop"
    args = (result, "iid", inputs / "omop_person_map.tsv", inputs / "omop_cohorts.tsv", output_dir)
    paths = gp.export_omop_cohort(*args)
    before = {path: Path(path).read_bytes() for path in paths}
    if failure == "duplicates":
        records = _read_tsv(result)
        pd.concat([records, records]).to_csv(result, sep="\t", index=False)
    elif failure == "provenance":
        def fail_hash(path):
            raise OSError("provenance failure")
        monkeypatch.setattr(gp, "_sha256_file", fail_hash)
    else:
        manifest = _read_tsv(inputs / "omop_cohorts.tsv")
        manifest["cohort_definition_name"] += " updated"
        manifest.to_csv(inputs / "omop_cohorts.tsv", sep="\t", index=False)
        replace = gp.os.replace

        def fail_publication(source, destination):
            if str(destination) == paths[2]:
                raise OSError("publication failure")
            return replace(source, destination)

        monkeypatch.setattr(gp.os, "replace", fail_publication)
    with pytest.raises((ValueError, OSError), match="duplicate|failure"):
        gp.export_omop_cohort(*args)
    assert {path: Path(path).read_bytes() for path in paths} == before
    assert not list(output_dir.glob(".getpheno-omop-*"))


@pytest.mark.parametrize("lowmem", [False, True])
def test_main_omop_output_matches_fixture(tmp_path, monkeypatch, lowmem):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    output_path = tmp_path / f"actual.omop.{lowmem}.tsv"
    omop_dir = tmp_path / f"omop.{lowmem}"

    _run_getpheno(
        input_dir,
        output_path,
        lowmem=lowmem,
        batchsize=2,
        write_omop_cohort=True,
        omop_person_map=input_dir / "omop_person_map.tsv",
        omop_definitions=input_dir / "omop_cohorts.tsv",
        omop_output_dir=omop_dir,
    )

    assert_generated_output_matches_fixture(
        output_path,
        FIXTURE_ROOT / "basic_multi" / "expected" / "output.tsv",
    )
    assert_generated_output_matches_fixture(
        omop_dir / "COHORT.tsv",
        FIXTURE_ROOT / "basic_multi" / "expected" / "COHORT.tsv",
    )
    assert_generated_output_matches_fixture(
        omop_dir / "COHORT_DEFINITION.tsv",
        FIXTURE_ROOT / "basic_multi" / "expected" / "COHORT_DEFINITION.tsv",
    )
    assert (omop_dir / "GETPHENO_PROVENANCE.json").is_file()
    checksum_text = output_path.with_name(
        output_path.name.rpartition(".")[0] + ".checksums.sha256"
    ).read_text(encoding="utf-8")
    assert "COHORT.tsv" in checksum_text
    assert "COHORT_DEFINITION.tsv" in checksum_text
    assert "GETPHENO_PROVENANCE.json" in checksum_text


def test_omop_export_rejects_missing_person_mapping(tmp_path):
    result = tmp_path / "result.tsv"
    result.write_text(
        "iid\tMDD\tMDD_earliest_date\n1\tCase\t2000-01-01\n",
        encoding="utf-8",
    )
    person_map = tmp_path / "map.tsv"
    person_map.write_text("iid\tperson_id\n2\t102\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "phenotype_name\tcohort_definition_id\tcohort_definition_name\tdefinition_type_concept_id\tsubject_concept_id\n"
        "MDD\t1001\tMDD cases\t0\t0\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="no person_id for 1 exported cohort memberships"):
        gp.export_omop_cohort(
            result,
            "iid",
            person_map,
            manifest,
            tmp_path / "omop",
            person_map_iid_col="iid",
        )


def test_omop_export_can_write_explicit_control_cohort(tmp_path):
    result = tmp_path / "result.tsv"
    result.write_text(
        "iid\tMDD\tMDD_earliest_date\tbirthdate\n"
        "1\tCase\t2000-01-01\t1980-01-01\n"
        "2\tControl\t\t1981-02-03\n",
        encoding="utf-8",
    )
    person_map = tmp_path / "map.tsv"
    person_map.write_text("iid\tperson_id\n1\t101\n2\t102\n", encoding="utf-8")
    manifest = tmp_path / "manifest.tsv"
    manifest.write_text(
        "phenotype_name\tcohort_definition_id\tcohort_definition_name\tdefinition_type_concept_id\tsubject_concept_id\tcontrol_cohort_definition_id\tcontrol_cohort_definition_name\tcontrol_start_date_column\tcohort_initiation_date\n"
        "MDD\t1001\tMDD cases\t0\t0\t1002\tMDD controls\tbirthdate\t2026-01-01\n",
        encoding="utf-8",
    )

    paths = gp.export_omop_cohort(
        result,
        "iid",
        person_map,
        manifest,
        tmp_path / "omop",
        person_map_iid_col="iid",
    )

    cohort = _read_tsv(paths[0])
    assert cohort[["cohort_definition_id", "subject_id"]].values.tolist() == [
        ["1001", "101"], ["1002", "102"]
    ]


@pytest.mark.parametrize("lowmem", [False, True])
def test_main_multiple_diagnosis_files_matches_fixture(tmp_path, monkeypatch, lowmem):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    source = pd.read_csv(input_dir / "lpr.csv")
    source.iloc[[0, 3]].to_csv(
        input_dir / "diagnoses.csv", index=False
    )
    procedures = source.iloc[[1, 2]].rename(columns={"diagnosis": "code"})
    procedures.to_csv(input_dir / "procedures.csv", index=False)
    output_path = tmp_path / f"actual.multiple.{lowmem}.tsv"

    _run_getpheno(
        input_dir,
        output_path,
        lowmem=lowmem,
        batchsize=2,
        diagnosis_files="diagnoses.csv,procedures.csv",
        diagnostic_col="diagnosis,code",
        diagnosis_type_col="",
    )

    assert_generated_output_matches_fixture(
        output_path,
        FIXTURE_ROOT / "basic_multi" / "expected" / "output.tsv",
    )


@pytest.mark.parametrize("lowmem", [False, True])
def test_main_preserves_stam_qc_intersection_with_multiple_stam_files(
    tmp_path, monkeypatch, lowmem
):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    pd.DataFrame({
        "iid": ["1", "2"],
        "birthdate": ["1980-01-01", "1975-06-15"],
        "sex": ["F", "M"],
        "primary_value": ["p1", "p2"],
        "shared_value": ["first-1", "first-2"],
    }).to_csv(input_dir / "stam_primary.csv", index=False)
    pd.DataFrame({
        "iid": ["2", "3", "4"],
        "birthdate": ["", "1990-03-20", "1988-12-31"],
        "sex": ["", "F", "M"],
        "secondary_value": ["s2", "s3", "s4"],
        "shared_value": ["second-2", "second-3", "second-4"],
    }).to_csv(input_dir / "stam_secondary.csv", index=False)
    (input_dir / "qc.fam").write_text(
        "family1 1 0 0 0 -9\n"
        "family2 2 0 0 0 -9\n"
        "family3 3 0 0 0 -9\n"
        "family4 4 0 0 0 -9\n"
        "family5 5 0 0 0 -9\n",
        encoding="utf-8",
    )
    output_path = tmp_path / f"actual.stam-qc.{lowmem}.tsv"

    _run_getpheno(
        input_dir,
        output_path,
        lowmem=lowmem,
        batchsize=2,
        stam_files="stam_primary.csv,stam_secondary.csv",
        qced_file="qc.fam",
        build_entry_exit=True,
    )

    result = _read_tsv(output_path).set_index("iid")
    assert set(result.index) == {"1", "2", "3", "4"}
    assert len(result) == 4
    assert result.loc["2", "shared_value"] == "first-2"
    assert result.loc["3", "secondary_value"] == "s3"
    assert result.loc["4", "MDD"] == "Control"
    assert result.loc["4", "PSY"] == "Control"
    assert result.loc["4", "EntryDate"] == ""
    assert result.loc["4", "ExitDate"] == ""


def test_main_ini_false_source_flags_keep_all_diagnosis_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    source = pd.read_csv(input_dir / "lpr.csv")
    source.iloc[[0, 3]].to_csv(input_dir / "diagnoses.csv", index=False)
    source.iloc[[1, 2]].rename(columns={"diagnosis": "code"}).to_csv(
        input_dir / "procedures.csv", index=False
    )
    ini_file = tmp_path / "get_pheno.ini"
    ini_file.write_text(
        "\n".join([
            "[Default]",
            f"lpr_file = {input_dir / 'diagnoses.csv'},{input_dir / 'procedures.csv'}",
            "diagnostic_col = diagnosis,code",
            "runLPRonly = False",
            "runPSYKonly = False",
        ]),
        encoding="utf-8",
    )
    output_path = tmp_path / "actual.ini.tsv"

    _run_getpheno(
        input_dir,
        output_path,
        diagnosis_type_col="",
        ini_file=ini_file,
        default_args=["-f", "--fcol", "--LPR", "--PSYK"],
    )

    assert_generated_output_matches_fixture(
        output_path,
        FIXTURE_ROOT / "basic_multi" / "expected" / "output.tsv",
    )


@pytest.mark.parametrize("lowmem", [False, True])
def test_main_loads_death_filters_from_ini(tmp_path, monkeypatch, lowmem):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    pd.DataFrame({
        "iid": ["1", "2", "3", "4"],
        "D_DODSDATO": ["2000-01-01", "", "2020-01-01", ""],
    }).to_csv(input_dir / "death.csv", index=False)
    ini_file = tmp_path / "get_pheno.ini"
    ini_file.write_text(
        "\n".join([
            "[Default]",
            f"addition_information_file = {input_dir / 'death.csv'}",
            "deathdatecol = D_DODSDATO",
            "MinMaxAgeDeath = 20,40",
            "keepAlive = False",
            "onlyDead = True",
        ]),
        encoding="utf-8",
    )
    output_path = tmp_path / f"actual.death.ini.{lowmem}.tsv"

    _run_getpheno(
        input_dir,
        output_path,
        lowmem=lowmem,
        batchsize=1,
        ini_file=ini_file,
        default_args=[
            "-j", "--deathdatecol", "--MinMaxAgeDeath", "--keepAlive", "--onlyDead"
        ],
    )

    result = _read_tsv(output_path)
    assert result["iid"].tolist() == ["1", "3"]
    assert result["AgeAtDeath"].tolist() == ["20", "29"]


@pytest.mark.parametrize("lowmem", [False, True])
def test_main_multiple_additional_files_merge_by_iid(tmp_path, monkeypatch, lowmem):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("basic_multi", tmp_path)
    pd.DataFrame({
        "iid": ["1", "2", "3", "4"],
        "person_value": ["p1", "p2", "p3", "p4"],
    }).to_csv(input_dir / "t_person.csv", index=False)
    pd.DataFrame({
        "iid": ["1", "2"],
        "D_DODSDATO": ["2020-03-04", "2021-04-05"],
        "D_FINDEDATO": ["2020-03-05", "2021-04-06"],
    }).to_csv(input_dir / "t_dodsaarsag_2.csv", index=False)
    output_path = tmp_path / f"actual.additional.{lowmem}.tsv"

    _run_getpheno(
        input_dir,
        output_path,
        lowmem=lowmem,
        batchsize=2,
        additional_file="t_person.csv,t_dodsaarsag_2.csv",
    )

    result = _read_tsv(output_path).set_index("iid")
    assert result.loc["1", "person_value"] == "p1"
    assert result.loc["1", "D_DODSDATO"] == "2020-03-04"
    assert result.loc["2", "D_FINDEDATO"] == "2021-04-06"
    assert result.loc["3", "D_DODSDATO"] == ""


@pytest.mark.parametrize("lowmem", [False, True])
def test_main_ipsych_comparison_matches_fixture(tmp_path, monkeypatch, capsys, lowmem):
    monkeypatch.chdir(tmp_path)
    input_dir = _copy_fixture_inputs("ipsych_comparison", tmp_path)
    output_path = tmp_path / f"actual.ipsych.{lowmem}.tsv"

    _run_getpheno(
        input_dir,
        output_path,
        lowmem=lowmem,
        batchsize=2,
        additional_file="additional.csv",
        status_col="status_code",
        status_date_col="status_date",
        ipsych_options={
            "enabled": True,
            "dk_parents_col": "DKDK",
            "status_col": "status_code",
            "status_date_col": "status_date",
        },
    )

    assert "neither status column" not in capsys.readouterr().out

    assert_generated_output_matches_fixture(
        output_path,
        FIXTURE_ROOT / "ipsych_comparison" / "expected" / "output.tsv",
    )


@pytest.mark.parametrize("lowmem,verylowmem", [(False, False), (True, False), (False, True)])
@pytest.mark.parametrize("indexed", [False, True])
@pytest.mark.parametrize("date_format", ["%Y-%m-%d", "%d/%m/%Y"])
def test_main_ipsych_presets_without_definition_file(tmp_path, monkeypatch, lowmem, verylowmem, indexed, date_format):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gp.socket, "gethostname", lambda: "local-test")
    iids = [f"{number:03}" for number in range(1, 13)]
    stam = pd.DataFrame({
        "iid": iids, "birthdate": ["1985-01-01" if iid in ("007", "012") else "2000-01-01" for iid in iids],
        "sex": ["F"] * 12, "stat": [1] * 12, "fkode": [5100] * 12,
        "pnr_m": [f"mother-{iid}" for iid in iids], "DKDK": [True] * 12,
    })
    stam["birthdate"] = pd.to_datetime(stam["birthdate"]).dt.strftime(date_format)
    stam.to_csv(tmp_path / "stam.csv", index=False)
    events = [
        ("001", "ICD10:F32.1", "2015-12-31", "pSyK"),
        ("002", "ICD10:F20.0", "2015-01-01", "prefix_PsYc_suffix"),
        ("003", "ICD10:F31.1", "2015-01-01", "PSYK"),
        ("004", "ICD10:F84.5", "2001-01-01", "PSYC"),
        ("005", "ICD10:F90.0", "2001-01-01", "psyc_diag"),
        ("006", "ICD10:F50.1", "2016-12-31", "somatic_LpR3"),
        ("007", "ICD10:Z37.0", "2014-06-01", "LPR"),
        ("007", "ICD10:F32.0", "2015-06-01", "psyc_diag"),
        ("007", "ICD10:F33.0", "2015-06-02", "PSYK"),
        ("008", "ICD10:F32.0", "2015-01-01", "LPR"),
        ("009", "ICD10:F32.0", "2009-12-31", "PSYK"),
        ("010", "ICD10:F50.0", "2017-01-01", "LPR"),
        ("011", "ICD10:F90.0", "2000-12-31", "PSYK"),
        ("012", "ICD8:650", "2014-01-01", "LPR"),
        ("012", "ICD8:296", "2015-01-01", "PSYK"),
    ]
    diagnoses = pd.DataFrame(events, columns=["iid", "diagnosis", "date_in", "origin"])
    diagnoses["date_out"] = diagnoses["date_in"]
    for col in ("date_in", "date_out"):
        diagnoses[col] = pd.to_datetime(diagnoses[col]).dt.strftime(date_format)
    diagnoses["register"] = "unselected"
    diagnoses["source"] = "unselected"
    diagnoses.to_csv(tmp_path / "lpr.csv", index=False)
    if indexed:
        gp.index_diag_file(str(tmp_path / "lpr.csv"), index_columns=["iid", "diagnosis"], separator=",")
    output = tmp_path / "presets.tsv"
    _run_getpheno(
        tmp_path, output, phenotype_file="", register_col="origin",
        ipsych_options={"enabled": True}, lowmem=lowmem, verylowmem=verylowmem, batchsize=3,
        date_format=date_format,
    )
    result = _read_tsv(output).set_index("iid")
    assert result.index.tolist() == iids
    expected_cases = {
        "AFFEK": {"001", "003", "007"}, "SKIZO": {"002"},
        "SKIZOSPEK": {"002"}, "BIPOL": {"003"}, "ASD": {"004"},
        "ADHD": {"005"}, "ANO": {"006"}, "PPD": {"007", "012"},
    }
    for name, expected in expected_cases.items():
        assert set(result.index[result[name] == "Case"]) == expected
    assert float(result.loc["007", "PPD_n_diags"]) == 1


@pytest.mark.parametrize("register_filter", ["PSYK", "LPR"])
def test_main_register_flag_filters_rows_in_every_input_file(tmp_path, monkeypatch, register_filter):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(gp.socket, "gethostname", lambda: "local-test")
    pd.DataFrame({"iid": ["001", "002"], "birthdate": ["1985-01-01"] * 2,
                  "sex": ["F", "M"]}).to_csv(tmp_path / "stam.csv", index=False)
    for filename, col, sources in [("first.csv", "origin", ["psyc_diag", "lpr3"]),
                                   ("second.csv", "registry", ["LPR", "prefix_psyk"] )]:
        pd.DataFrame({"iid": ["001", "002"], "diagnosis": ["ICD10:F32"] * 2,
                      "date_in": ["2015-01-01"] * 2, "date_out": ["2015-01-02"] * 2,
                      col: sources}).to_csv(tmp_path / filename, index=False)
    (tmp_path / "pheno.tsv").write_text("CUSTOM\tICD10:F32\n")
    ini = tmp_path / "config.ini"
    ini.write_text("[Default]\nregistercol = origin,registry\n")
    output = tmp_path / "selected.tsv"
    _run_getpheno(
        tmp_path, output, diagnosis_files="first.csv,second.csv", ini_file=ini,
        default_args=["--registercol"], run_psyk=register_filter == "PSYK", run_lpr=register_filter == "LPR",
    )
    result = _read_tsv(output)
    assert result["CUSTOM"].tolist() == ["Case", "Case"]
    assert result["CUSTOM_n_diags"].astype(float).tolist() == [1, 1]


def test_cli_presets_index_dtype_and_psyk_preserves_livebirth_evidence(tmp_path):
    pd.DataFrame({"iid": ["001"], "birthdate": ["1985-01-01"], "sex": ["F"],
                  "stat": [1], "fkode": [5100], "pnr_m": ["mother"], "DKDK": [True]}).to_csv(
                      tmp_path / "stam.csv", index=False)
    pd.DataFrame({
        "iid": ["001"] * 2, "diagnosis": ["ICD10:Z37.0", "ICD10:F32.1"],
        "date_in": ["2014-06-01", "2015-06-01"], "date_out": ["2014-06-01", "2015-06-01"],
        "origin": ["LPR", "psyc_diag"],
    }).to_csv(tmp_path / "lpr.csv", index=False)
    command = [
        sys.executable, str(REPO_ROOT / "get_pheno.py"),
        "-f", str(tmp_path / "lpr.csv"), "-i", str(tmp_path / "stam.csv"),
        "-o", str(tmp_path / "output.tsv"), "--ini", str(tmp_path / "missing.ini"),
        "--iidcol", "iid", "--fcol", "diagnosis", "--din", "date_in", "--don", "date_out",
        "--registercol", "origin", "--IndexDtypes", '{"origin":"float"}',
        "--DateFormat", "%Y-%m-%d", "--iPSYCHComparison", "--PSYK",
    ]
    for extra_args in (["--BuildIndex"], []):
        completed = subprocess.run(command + extra_args, cwd=tmp_path, text=True, capture_output=True)
        assert completed.returncode == 0, completed.stderr + completed.stdout[-2000:]
    with pd.HDFStore(tmp_path / "lpr.h5", mode="r") as store:
        frame = store.select(store.keys()[0])
    assert frame["origin"].tolist() == ["LPR", "psyc_diag"]
    assert frame["iid"].tolist() == ["001", "001"]
    result = _read_tsv(tmp_path / "output.tsv")
    assert result["PPD"].tolist() == ["Case"]
