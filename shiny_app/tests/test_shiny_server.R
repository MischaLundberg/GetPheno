if (!requireNamespace("shiny", quietly = TRUE)) {
  cat("Skipping Shiny server test because package 'shiny' is not installed.\n")
  quit(status = 0L)
}

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[startsWith(args, "--file=")][[1]])
app_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
repo_root <- normalizePath(file.path(app_dir, ".."), mustWork = TRUE)
temporary <- tempfile("getpheno-shiny-server-")
dir.create(file.path(temporary, "project"), recursive = TRUE)
dir.create(file.path(temporary, "output"))
on.exit(unlink(temporary, recursive = TRUE), add = TRUE)

ini <- file.path(temporary, "project", "get_pheno.ini")
death_file <- file.path(temporary, "project", "death.tsv")
fake_python <- file.path(temporary, "python")
fake_sbatch <- file.path(temporary, "sbatch")
writeLines("[Default]", ini)
writeLines("iid\tD_DODSDATO", death_file)
writeLines(c("#!/bin/sh", "exit 0"), fake_python)
writeLines(c("#!/bin/sh", "printf '%s\\n' '987654;test-cluster'"), fake_sbatch)
Sys.chmod(c(fake_python, fake_sbatch), "0700")

registry <- file.path(temporary, "projects.conf")
writeLines(c(
  "Project: test",
  "ProjectLabel: Test project",
  "Protocol: standard",
  "ProtocolLabel: Standard protocol",
  paste("ProjectRoot:", file.path(temporary, "project")),
  paste("IniFile:", ini),
  paste("OutputRoot:", file.path(temporary, "output")),
  "AccessGroup:",
  paste("Python:", fake_python),
  paste("GetPheno:", file.path(repo_root, "get_pheno.py")),
  paste("StandardAdditionalFiles:", death_file),
  "DeathDateColumn: D_DODSDATO",
  paste("Sbatch:", fake_sbatch),
  "Partition: normal",
  "Account: test",
  "Time: 01:00:00",
  "Memory: 4G",
  "Cpus: 2"
), registry)

old_config <- Sys.getenv("GETPHENO_SHINY_CONFIG", unset = NA_character_)
on.exit({
  if (is.na(old_config)) Sys.unsetenv("GETPHENO_SHINY_CONFIG") else
    Sys.setenv(GETPHENO_SHINY_CONFIG = old_config)
}, add = TRUE)
Sys.setenv(GETPHENO_SHINY_CONFIG = registry)

environment <- new.env(parent = globalenv())
sys.source(file.path(app_dir, "app.R"), envir = environment, chdir = TRUE)
ui_html <- as.character(environment$ui)
stopifnot(grepl("field-help", ui_html, fixed = TRUE))
stopifnot(grepl("main=ICD10:F32", ui_html, fixed = TRUE))
stopifnot(grepl("Lifetime exclusion", ui_html, fixed = TRUE))
stopifnot(grepl("IID subset", ui_html, fixed = TRUE))

shiny::testServer(environment$server, {
  session$setInputs(
    project = "test",
    protocol = "standard",
    phenotype_mode = "paste",
    phenotype_text = "MDD\tICD10:F32,ICD10:F33",
    phenotype_name_value = "MDD",
    iid_subset_mode = "paste",
    iid_subset_text = "1001\n1002",
    lifetime_mode = "none",
    post_mode = "none",
    prior_mode = "none",
    ipsych = TRUE,
    exact_match = TRUE,
    formats = c("fastgwa", "plink2"),
    icd_mode = "ICD8",
    gender = "F",
    fyob = "1981-05-01",
    min_max_age = "10,65",
    death_age_range = "20,80",
    keep_alive = TRUE,
    only_dead = FALSE,
    diag_exclusions = "H,M",
    diag_inclusions = "A,B"
  )
  session$setInputs(submit = 1)
  session$flushReact()
})

user_root <- file.path(temporary, "output", Sys.info()[["user"]])
runs <- list.dirs(user_root, recursive = FALSE, full.names = TRUE)
stopifnot(length(runs) == 1L)
run_id <- basename(runs[[1]])
metadata <- read.delim(
  file.path(runs[[1]], "metadata.tsv"), header = FALSE,
  col.names = c("key", "value"), sep = "\t", stringsAsFactors = FALSE
)
stopifnot(metadata[metadata$key == "job_id", "value"] == "987654")
stopifnot(grepl("--selectIIDs", metadata[metadata$key == "command", "value"], fixed = TRUE))
stopifnot(grepl("--deathdatecol", metadata[metadata$key == "command", "value"], fixed = TRUE))
stopifnot(grepl("--keepAlive", metadata[metadata$key == "command", "value"], fixed = TRUE))
stopifnot(file.exists(file.path(runs[[1]], paste0(run_id, ".iids.txt"))))

shiny::testServer(environment$server, {
  session$setInputs(
    project = "test", protocol = "standard",
    lookup_code = run_id, lookup = 1
  )
  session$flushReact()
  stopifnot(any(grepl(run_id, as.character(output$run_details), fixed = TRUE)))
  stopifnot(nrow(output$run_metadata) == 4L)
})

shiny::testServer(environment$server, {
  session$setInputs(
    project = "test", protocol = "standard",
    phenotype_mode = "none", iid_subset_mode = "none",
    lifetime_mode = "none", post_mode = "none", prior_mode = "none",
    ipsych = TRUE, formats = character(), gender = "", fyob = "",
    min_max_age = "0,0", death_age_range = "0,0"
  )
  session$setInputs(submit = 1)
  session$flushReact()
})
all_runs <- list.dirs(user_root, recursive = FALSE, full.names = TRUE)
preset_run <- setdiff(all_runs, runs)
stopifnot(length(preset_run) == 1L)
preset_metadata <- read.delim(
  file.path(preset_run, "metadata.tsv"), header = FALSE,
  col.names = c("key", "value"), sep = "\t", stringsAsFactors = FALSE
)
preset_command <- preset_metadata[preset_metadata$key == "command", "value"]
stopifnot(grepl("--iPSYCHComparison", preset_command, fixed = TRUE))
stopifnot(!grepl("'-g'", preset_command, fixed = TRUE))

cat("GetPheno Shiny submission and retrieval tests passed.\n")
