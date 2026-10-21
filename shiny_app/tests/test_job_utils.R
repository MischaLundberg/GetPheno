args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[startsWith(args, "--file=")][[1]])
app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(app_dir, "R", "job_utils.R"))

old_runtime_user <- Sys.getenv("GETPHENO_UNIX_USER", unset = NA_character_)
Sys.setenv(GETPHENO_UNIX_USER = "container_user")
stopifnot(current_unix_user() == "container_user")
if (is.na(old_runtime_user)) Sys.unsetenv("GETPHENO_UNIX_USER") else
  Sys.setenv(GETPHENO_UNIX_USER = old_runtime_user)

temporary <- tempfile("getpheno-shiny-test-")
dir.create(temporary)
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)

project_root <- file.path(temporary, "project")
output_root <- file.path(temporary, "output")
dir.create(project_root)
dir.create(output_root)
ini <- file.path(project_root, "get_pheno.ini")
python <- file.path(temporary, "python")
getpheno <- file.path(temporary, "get_pheno.py")
freeze_old <- file.path(project_root, "chb_oral_freeze_cohort_20240614_20260126.tsv")
freeze_new <- file.path(project_root, "chb_oral_freeze_cohort_20260123_20260818.tsv")
diagnosis <- file.path(project_root, "diagnoses.tsv")
procedures <- file.path(project_root, "procedures.tsv")
standard_info <- file.path(project_root, "t_person.tsv")
optional_info <- file.path(project_root, "laboratory.tsv")
genetics_root <- file.path(temporary, "genetics", "chb_oral_freeze_20260123")
general_exclusions <- file.path(genetics_root, "oral_exclusion_latest")
qced <- file.path(genetics_root, "01.Genotypes", "GSA_oral_qc.fam")
writeLines("[Default]", ini)
writeLines("#!/bin/sh", python)
writeLines("print('test')", getpheno)
dir.create(dirname(qced), recursive = TRUE)
for (path in c(freeze_old, freeze_new, diagnosis, procedures, standard_info, optional_info,
               general_exclusions, qced)) writeLines("test", path)
Sys.chmod(python, "0700")

registry_file <- file.path(temporary, "projects.conf")
writeLines(c(
  "Project: test",
  "ProjectLabel: Test project",
  "Protocol: default",
  "ProtocolLabel: Default protocol",
  paste("ProjectRoot:", project_root),
  paste("IniFile:", ini),
  paste("OutputRoot:", output_root),
  "AccessGroup:",
  paste("Python:", python),
  paste("GetPheno:", getpheno),
  paste("FreezePattern:", file.path(project_root, "chb_oral_freeze_cohort_*.tsv")),
  paste0("DiagnosisFiles: ", diagnosis, "|", procedures),
  "DiagnosticColumns: diagnosis,code",
  "IidColumns: cpr_enc",
  "DateInColumns: date_in",
  "DateOutColumns: date_out",
  "DiagnosisTypeColumns: diagnosis_type,",
  "DeathDateColumn: D_DODSDATO",
  paste("StandardAdditionalFiles:", standard_info),
  paste0("OptionalAdditionalFiles: Laboratory=", optional_info),
  paste("GeneralExclusionsTemplate:", file.path(temporary, "genetics", "chb_oral_freeze_{freeze_base}", "oral_exclusion_latest")),
  paste("QcedTemplate:", file.path(temporary, "genetics", "chb_oral_freeze_{freeze_base}", "01.Genotypes", "GSA_oral_qc.fam")),
  "Sbatch: sbatch",
  "Partition: normal",
  "Account: account",
  "Time: 01:00:00",
  "Memory: 4G",
  "Cpus: 2"
), registry_file)

registry <- read_project_registry(registry_file)
stopifnot(nrow(registry) == 1L)
stopifnot(nrow(accessible_registry(registry, groups = character())) == 1L)

freezes <- discover_freezes(registry[1, , drop = FALSE])
stopifnot(nrow(freezes) == 2L)
stopifnot(freezes$path[[1]] == freeze_new)
stopifnot(freezes$label[[1]] == "20260818 (cohort 20260123)")
additional_choices <- parse_optional_additional_files(registry$OptionalAdditionalFiles)
stopifnot(additional_choices$label[[1]] == "Laboratory")
protocol_files <- resolve_protocol_files(
  registry[1, , drop = FALSE], freeze_new, additional_choices$key[[1]]
)
stopifnot(protocol_files$stam == freeze_new)
stopifnot(identical(protocol_files$diagnosis, c(diagnosis, procedures)))
stopifnot(protocol_files$diagnostic_columns == "diagnosis,code")
stopifnot(protocol_files$death_date_column == "D_DODSDATO")
stopifnot(identical(protocol_files$additional, c(standard_info, optional_info)))
stopifnot(protocol_files$general_exclusions == general_exclusions)
stopifnot(protocol_files$qced == qced)

stopifnot(phenotype_count("F32\nF33") == 1L)
stopifnot(phenotype_count("MDD\tF32,F33\nSCZ\tF20") == 2L)
stopifnot(phenotype_count("Disorder\tDiagnosis\nMDD\tF32") == 1L)

run_id <- new_run_id(as.POSIXct("2026-10-02 12:00:00", tz = "UTC"))
stopifnot(is_valid_run_id(run_id))
run_dir <- file.path(output_root, "tester", run_id)
dir.create(run_dir, recursive = TRUE)
phenotype_file <- write_private_text(file.path(run_dir, paste0(run_id, ".phenotypes.tsv")), "MDD\tF32,F33")
iid_subset <- write_private_text(file.path(run_dir, paste0(run_id, ".iids.txt")), "1001\n1002")
output_file <- file.path(run_dir, paste0(run_id, ".tsv"))

command <- build_getpheno_command(
  registry[1, , drop = FALSE], run_id, phenotype_file, output_file,
  options = list(
    ipsych = TRUE,
    fastgwa = TRUE,
    plink2 = TRUE,
    exact_match = TRUE,
    gender = "F",
    fyob = "1981-05-01",
    min_max_age = "10,65",
    death_age_range = "20,80",
    keep_alive = TRUE,
    only_dead = FALSE,
    diag_exclusions = "H,M",
    diag_inclusions = "A,B",
    icd_mode = "ICD8",
    single_phenotype = TRUE,
    name = "MDD"
  ),
  protocol_files = protocol_files,
  iid_subset = iid_subset
)
stopifnot("--iPSYCHComparison" %in% command$args)
stopifnot("--write_fastGWA_format" %in% command$args)
stopifnot("--write_Plink2_format" %in% command$args)
stopifnot("--ICD8" %in% command$args)
stopifnot(command$args[match("-i", command$args) + 1L] == freeze_new)
stopifnot(command$args[match("-f", command$args) + 1L] == paste(diagnosis, procedures, sep = ","))
stopifnot(command$args[match("-j", command$args) + 1L] == paste(standard_info, optional_info, sep = ","))
stopifnot(command$args[match("--ge", command$args) + 1L] == general_exclusions)
stopifnot(command$args[match("--qced", command$args) + 1L] == qced)
stopifnot(command$args[match("--fcol", command$args) + 1L] == "diagnosis,code")
stopifnot(command$args[match("--iidcol", command$args) + 1L] == "cpr_enc")
stopifnot(command$args[match("--DiagTypecol", command$args) + 1L] == "diagnosis_type,")
stopifnot(command$args[match("--deathdatecol", command$args) + 1L] == "D_DODSDATO")
stopifnot(command$args[match("--MinMaxAgeDeath", command$args) + 1L] == "20,80")
stopifnot("--keepAlive" %in% command$args)
stopifnot(command$args[match("--selectIIDs", command$args) + 1L] == iid_subset)
stopifnot(command$args[match("--name", command$args) + 1L] == "MDD")

script <- write_slurm_script(registry[1, , drop = FALSE], run_id, run_dir, command)
script_text <- paste(readLines(script), collapse = "\n")
stopifnot(grepl("#SBATCH --partition=normal", script_text, fixed = TRUE))
stopifnot(grepl("--iPSYCHComparison", script_text, fixed = TRUE))
stopifnot(grepl("status.txt", script_text, fixed = TRUE))
stopifnot(system2("bash", script) == 0L)
stopifnot(readLines(file.path(run_dir, "status.txt"), warn = FALSE)[[1]] == "COMPLETED")

fake_sbatch <- file.path(temporary, "sbatch")
writeLines(c("#!/bin/sh", "printf '%s\\n' '12345;cluster'"), fake_sbatch)
Sys.chmod(fake_sbatch, "0700")
stopifnot(submit_slurm_job(script, fake_sbatch) == "12345")

metadata_path <- file.path(run_dir, "metadata.tsv")
write_run_metadata(metadata_path, list(run_id = run_id, job_id = "12345"))
metadata <- read_run_metadata(metadata_path)
stopifnot(metadata$run_id == run_id, metadata$job_id == "12345")

found <- find_run_directory(registry, "tester", run_id)
stopifnot(identical(found, normalizePath(run_dir)))

cat("All GetPheno Shiny utility tests passed.\n")
