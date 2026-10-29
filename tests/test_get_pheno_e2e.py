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
        True,
        date_in_col,
        date_out_col,
        str(input_dir / qced_file) if qced_file else "",
        "",
        "",
        diagnosis_type_col,
        "",
        "",
        "",
        False,
        "",
        "",
        False,
        False,
        False,
        False,
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
    )


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
