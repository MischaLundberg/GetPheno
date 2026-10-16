import shutil
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
):
    gp.main(
        ",".join(str(input_dir / path) for path in diagnosis_files.split(",")),
        str(input_dir / "pheno.tsv"),
        str(input_dir / "stam.csv"),
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
        "",
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
        "%Y-%m-%d",
        status_col,
        status_date_col,
        "",
        False,
        1,
        "MainPheno",
        False,
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
        False,
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
        False,
        False,
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
