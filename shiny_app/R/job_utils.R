`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0L || is.na(value[[1]]) || !nzchar(value[[1]])) {
    fallback
  } else {
    value
  }
}

trim_value <- function(value) {
  trimws(as.character(value %||% "")[[1]])
}

current_unix_user <- function() {
  configured <- trim_value(Sys.getenv("GETPHENO_UNIX_USER"))
  if (nzchar(configured)) {
    if (!grepl("^[A-Za-z_][A-Za-z0-9_.-]*$", configured)) {
      stop("GETPHENO_UNIX_USER contains unsupported characters")
    }
    return(configured)
  }
  trim_value(Sys.info()[["user"]] %||% Sys.getenv("USER"))
}

unix_groups <- function() {
  output <- tryCatch(
    system2("id", "-Gn", stdout = TRUE, stderr = FALSE),
    error = function(error) character()
  )
  unique(strsplit(paste(output, collapse = " "), "[[:space:]]+")[[1]])
}

read_project_registry <- function(path) {
  if (!file.exists(path)) {
    stop("Project registry not found: ", path)
  }
  registry <- as.data.frame(read.dcf(path), stringsAsFactors = FALSE)
  required <- c(
    "Project", "ProjectLabel", "Protocol", "ProtocolLabel",
    "ProjectRoot", "IniFile", "OutputRoot", "Python", "GetPheno"
  )
  missing <- setdiff(required, names(registry))
  if (length(missing)) {
    stop("Project registry is missing fields: ", paste(missing, collapse = ", "))
  }
  optional <- c(
    "AccessGroup", "Sbatch", "Partition", "Account", "Time", "Memory", "Cpus",
    "FreezePattern", "DiagnosisFile", "DiagnosisFiles", "DiagnosticColumns",
    "IidColumns", "DateInColumns", "DateOutColumns", "DiagnosisTypeColumns",
    "RecordNumberColumns", "RegisterColumns", "DeathDateColumn", "StandardAdditionalFiles",
    "OptionalAdditionalFiles", "GeneralExclusionsTemplate", "QcedTemplate"
  )
  for (field in optional) {
    if (!field %in% names(registry)) registry[[field]] <- ""
  }
  registry[] <- lapply(registry, trimws)
  if (any(!grepl("^[A-Za-z0-9_.-]+$", registry$Project))) {
    stop("Project identifiers may only contain letters, numbers, dot, underscore, and hyphen")
  }
  if (any(!grepl("^[A-Za-z0-9_.-]+$", registry$Protocol))) {
    stop("Protocol identifiers may only contain letters, numbers, dot, underscore, and hyphen")
  }
  if (any(grepl("[[:space:]]", registry$OutputRoot))) {
    stop("OutputRoot paths must not contain whitespace")
  }
  duplicate_keys <- duplicated(paste(registry$Project, registry$Protocol, sep = "/"))
  if (any(duplicate_keys)) stop("Each project/protocol pair must be unique")
  registry
}

record_is_accessible <- function(record, groups = unix_groups()) {
  project_root <- trim_value(record$ProjectRoot)
  ini_file <- trim_value(record$IniFile)
  output_root <- trim_value(record$OutputRoot)
  access_group <- trim_value(record$AccessGroup)
  group_ok <- !nzchar(access_group) || access_group %in% groups
  group_ok &&
    dir.exists(project_root) && file.access(project_root, 4L) == 0L &&
    file.exists(ini_file) && file.access(ini_file, 4L) == 0L &&
    dir.exists(output_root) && file.access(output_root, 2L) == 0L
}

accessible_registry <- function(registry, groups = unix_groups()) {
  keep <- vapply(
    seq_len(nrow(registry)),
    function(index) record_is_accessible(registry[index, , drop = FALSE], groups),
    logical(1)
  )
  registry[keep, , drop = FALSE]
}

extract_freeze_dates <- function(path) {
  matches <- regmatches(basename(path), gregexpr("[0-9]{8}", basename(path)))[[1]]
  matches <- matches[grepl("^[0-9]{8}$", matches)]
  list(
    base = if (length(matches)) matches[[1]] else "",
    version = if (length(matches)) matches[[length(matches)]] else ""
  )
}

discover_freezes <- function(record) {
  pattern <- trim_value(record$FreezePattern)
  if (!nzchar(pattern)) {
    return(data.frame(path = character(), label = character(), stringsAsFactors = FALSE))
  }
  paths <- Sys.glob(pattern)
  paths <- paths[file.exists(paths) & file.access(paths, 4L) == 0L]
  if (!length(paths)) {
    return(data.frame(path = character(), label = character(), stringsAsFactors = FALSE))
  }
  dates <- lapply(paths, extract_freeze_dates)
  base_dates <- vapply(dates, `[[`, "", "base")
  versions <- vapply(dates, `[[`, "", "version")
  labels <- ifelse(
    nzchar(versions),
    paste0(versions, " (cohort ", base_dates, ")"),
    basename(paths)
  )
  order_index <- order(versions, base_dates, basename(paths), decreasing = TRUE)
  data.frame(path = paths[order_index], label = labels[order_index], stringsAsFactors = FALSE)
}

expand_freeze_template <- function(template, freeze_path = "") {
  template <- trim_value(template)
  if (!nzchar(template)) return("")
  dates <- extract_freeze_dates(freeze_path)
  replacements <- c(
    "{freeze_file}" = freeze_path,
    "{freeze_dir}" = if (nzchar(freeze_path)) dirname(freeze_path) else "",
    "{freeze_name}" = if (nzchar(freeze_path)) basename(freeze_path) else "",
    "{freeze_base}" = dates$base,
    "{freeze_version}" = dates$version
  )
  output <- template
  for (token in names(replacements)) {
    output <- gsub(token, replacements[[token]], output, fixed = TRUE)
  }
  if (grepl("\\{freeze_[^}]+\\}", output)) stop("Unsupported freeze template token")
  output
}

split_config_paths <- function(value) {
  value <- trim_value(value)
  if (!nzchar(value)) return(character())
  trimws(strsplit(value, "|", fixed = TRUE)[[1]])
}

parse_optional_additional_files <- function(value) {
  entries <- split_config_paths(value)
  if (!length(entries)) return(data.frame(key = character(), label = character(), path = character()))
  parsed <- lapply(seq_along(entries), function(index) {
    entry <- entries[[index]]
    separator <- regexpr("=", entry, fixed = TRUE)[[1]]
    if (separator < 1L) {
      label <- basename(entry)
      path <- entry
    } else {
      label <- trimws(substr(entry, 1L, separator - 1L))
      path <- trimws(substr(entry, separator + 1L, nchar(entry)))
    }
    data.frame(key = paste0("additional_", index), label = label, path = path)
  })
  do.call(rbind, parsed)
}

resolve_protocol_files <- function(record, freeze_path = "", selected_additional = character()) {
  freezes <- discover_freezes(record)
  if (nrow(freezes)) {
    if (!nzchar(trim_value(freeze_path)) || !freeze_path %in% freezes$path) {
      stop("Select an available cohort freeze")
    }
  } else if (nzchar(trim_value(record$FreezePattern))) {
    stop("No readable cohort freeze matches the configured pattern")
  }

  optional <- parse_optional_additional_files(record$OptionalAdditionalFiles)
  if (length(selected_additional) && any(!selected_additional %in% optional$key)) {
    stop("An unsupported additional-information file was selected")
  }
  optional_paths <- optional$path[optional$key %in% selected_additional]
  standard_paths <- split_config_paths(record$StandardAdditionalFiles)
  diagnosis_paths <- split_config_paths(record$DiagnosisFiles %||% record$DiagnosisFile)
  expand_all <- function(paths) unname(vapply(
    paths, expand_freeze_template, "", freeze_path = freeze_path
  ))

  files <- list(
    stam = freeze_path,
    diagnosis = expand_all(diagnosis_paths),
    additional = c(expand_all(standard_paths), expand_all(optional_paths)),
    general_exclusions = expand_freeze_template(record$GeneralExclusionsTemplate, freeze_path),
    qced = expand_freeze_template(record$QcedTemplate, freeze_path),
    diagnostic_columns = trim_value(record$DiagnosticColumns),
    iid_columns = trim_value(record$IidColumns),
    date_in_columns = trim_value(record$DateInColumns),
    date_out_columns = trim_value(record$DateOutColumns),
    diagnosis_type_columns = trim_value(record$DiagnosisTypeColumns),
    register_columns = trim_value(record$RegisterColumns),
    record_number_columns = trim_value(record$RecordNumberColumns),
    death_date_column = trim_value(record$DeathDateColumn)
  )
  configured <- c(files$stam, files$diagnosis, files$additional, files$general_exclusions, files$qced)
  configured <- configured[nzchar(configured)]
  if (any(grepl(",", configured, fixed = TRUE))) stop("Configured input paths must not contain commas")
  unavailable <- configured[!file.exists(configured) | file.access(configured, 4L) != 0L]
  if (length(unavailable)) {
    stop("Configured input file is unavailable or unreadable: ", unavailable[[1]])
  }
  files
}

random_hex <- function(bytes = 6L) {
  raw_value <- tryCatch(
    {
      connection <- file("/dev/urandom", open = "rb", raw = TRUE)
      on.exit(close(connection), add = TRUE)
      readBin(connection, what = "raw", n = bytes)
    },
    error = function(error) raw()
  )
  if (length(raw_value) != bytes) {
    return(paste(sample(c(0:9, letters[1:6]), bytes * 2L, replace = TRUE), collapse = ""))
  }
  paste(format(raw_value), collapse = "")
}

new_run_id <- function(now = Sys.time()) {
  timestamp <- format(now, "%Y%m%dT%H%M%SZ", tz = "UTC")
  paste0("GP-", timestamp, "-", random_hex())
}

is_valid_run_id <- function(value) {
  grepl("^GP-[0-9]{8}T[0-9]{6}Z-[0-9a-f]{12}$", trim_value(value))
}

normalize_definition_text <- function(text) {
  text <- paste(as.character(text %||% ""), collapse = "\n")
  gsub("\r\n?", "\n", text)
}

phenotype_count <- function(text) {
  lines <- trimws(strsplit(normalize_definition_text(text), "\n", fixed = TRUE)[[1]])
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(0L)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  two_column <- vapply(fields, length, integer(1)) >= 2L
  if (!any(two_column)) return(1L)
  data_lines <- lines[two_column]
  second <- tolower(trimws(vapply(strsplit(data_lines, "\t", fixed = TRUE), `[[`, "", 2L)))
  as.integer(length(data_lines) - as.integer(second[[1]] %in% c("diagnosis", "diagnoses")))
}

read_uploaded_text <- function(path, max_bytes = 1024L * 1024L) {
  if (!file.exists(path)) stop("Uploaded file is no longer available")
  size <- file.info(path)$size
  if (is.na(size) || size <= 0L) stop("Uploaded file is empty")
  if (size > max_bytes) stop("Definition files must not exceed 1 MiB")
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

write_private_text <- function(path, text) {
  writeLines(normalize_definition_text(text), path, useBytes = TRUE)
  Sys.chmod(path, mode = "0600")
  path
}

materialize_definition <- function(mode, uploaded_path, text, destination, label, required = FALSE) {
  mode <- trim_value(mode)
  if (identical(mode, "none")) {
    if (required) stop(label, " is required")
    return(NULL)
  }
  content <- if (identical(mode, "upload")) {
    if (!nzchar(trim_value(uploaded_path))) stop("Select a file for ", label)
    read_uploaded_text(uploaded_path)
  } else if (identical(mode, "paste")) {
    normalize_definition_text(text)
  } else {
    stop("Unsupported input mode for ", label)
  }
  if (!nzchar(trimws(content))) stop(label, " is empty")
  write_private_text(destination, content)
}

validate_date_or_blank <- function(value, label) {
  value <- trim_value(value)
  if (!nzchar(value)) return("")
  parsed <- as.Date(value, format = "%Y-%m-%d")
  if (is.na(parsed) || format(parsed, "%Y-%m-%d") != value) {
    stop(label, " must use YYYY-MM-DD")
  }
  value
}

validate_age_range <- function(value, label = "MinMaxAge") {
  value <- gsub("[[:space:]]+", "", trim_value(value))
  if (!grepl("^[0-9]+,[0-9]+$", value)) {
    stop(label, " must contain two non-negative integers separated by a comma")
  }
  limits <- as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
  if (limits[[2]] > 0L && limits[[1]] > limits[[2]]) stop(label, " minimum cannot exceed maximum")
  value
}

validate_diag_types <- function(value, label) {
  value <- toupper(gsub("[[:space:]]+", "", trim_value(value)))
  if (!nzchar(value)) return("")
  tokens <- strsplit(value, ",", fixed = TRUE)[[1]]
  if (any(!tokens %in% c("+", "A", "B", "C", "G", "H", "M"))) {
    stop(label, " contains unsupported diagnosis types")
  }
  paste(unique(tokens), collapse = ",")
}

validate_name <- function(value) {
  value <- trim_value(value)
  if (!nzchar(value)) return("")
  if (!grepl("^[A-Za-z][A-Za-z0-9_.-]{0,63}$", value)) {
    stop("Phenotype name must start with a letter and contain at most 64 safe characters")
  }
  value
}

build_getpheno_command <- function(record, run_id, phenotype_file, output_file,
                                   options, exclusions = list(), protocol_files = list(),
                                   iid_subset = "") {
  python <- trim_value(record$Python)
  script <- trim_value(record$GetPheno)
  ini_file <- trim_value(record$IniFile)
  phenotype_file <- trim_value(phenotype_file)
  if (!nzchar(phenotype_file) && !isTRUE(options$ipsych)) {
    stop("Phenotype definitions are required unless iPSYCH comparison is selected")
  }
  for (path in c(script, ini_file, if (nzchar(phenotype_file)) phenotype_file)) {
    if (!file.exists(path)) stop("Required run file is unavailable: ", path)
    if (file.access(path, 4L) != 0L) stop("Required run file is not readable: ", path)
  }
  iid_subset <- trim_value(iid_subset)
  if (nzchar(iid_subset) && (!file.exists(iid_subset) || file.access(iid_subset, 4L) != 0L)) {
    stop("IID subset file is unavailable or unreadable: ", iid_subset)
  }
  if (!file.exists(python) || file.access(python, 1L) != 0L) {
    stop("Configured Python executable is unavailable or not executable: ", python)
  }

  args <- c(script, "--ini", ini_file, if (nzchar(phenotype_file)) c("-g", phenotype_file), "-o", output_file)
  add_flag <- function(flag, enabled) if (isTRUE(enabled)) flag else character()
  add_value <- function(flag, value) if (nzchar(trim_value(value))) c(flag, trim_value(value)) else character()
  death_date_column <- trim_value(protocol_files$death_date_column %||% "")
  death_age_range <- validate_age_range(options$death_age_range %||% "0,0", "MinMaxAgeDeath")
  keep_alive <- isTRUE(options$keep_alive)
  only_dead <- isTRUE(options$only_dead)
  if (keep_alive && only_dead) stop("Keep living and only deceased cannot both be selected")
  if (!nzchar(death_date_column) && (!identical(death_age_range, "0,0") || keep_alive || only_dead)) {
    stop("The selected protocol does not configure a death-date column")
  }

  args <- c(
    args,
    add_value("-i", protocol_files$stam %||% ""),
    add_value("-f", paste(protocol_files$diagnosis %||% character(), collapse = ",")),
    add_value("-j", paste(protocol_files$additional %||% character(), collapse = ",")),
    add_value("--ge", protocol_files$general_exclusions %||% ""),
    add_value("--qced", protocol_files$qced %||% ""),
    add_value("--fcol", protocol_files$diagnostic_columns %||% ""),
    add_value("--iidcol", protocol_files$iid_columns %||% ""),
    add_value("--din", protocol_files$date_in_columns %||% ""),
    add_value("--don", protocol_files$date_out_columns %||% ""),
    add_value("--DiagTypecol", protocol_files$diagnosis_type_columns %||% ""),
    add_value("--registercol", protocol_files$register_columns %||% ""),
    add_value("--recnum", protocol_files$record_number_columns %||% ""),
    add_value("--deathdatecol", death_date_column),
    add_value("--selectIIDs", iid_subset),
    add_flag("--iPSYCHComparison", options$ipsych),
    add_flag("--write_fastGWA_format", options$fastgwa),
    add_flag("--write_Plink2_format", options$plink2),
    add_flag("--eM", options$exact_match),
    add_value("--Fgender", if (options$gender %in% c("F", "M")) options$gender else ""),
    add_value("--Fyob", validate_date_or_blank(options$fyob, "Fyob")),
    "--MinMaxAge", validate_age_range(options$min_max_age),
    if (nzchar(death_date_column)) c("--MinMaxAgeDeath", death_age_range) else character(),
    add_flag("--keepAlive", nzchar(death_date_column) && keep_alive),
    add_flag("--onlyDead", nzchar(death_date_column) && only_dead),
    add_value("--DiagTypeExclusions", validate_diag_types(options$diag_exclusions, "Diagnosis exclusions")),
    add_value("--DiagTypeInclusions", validate_diag_types(options$diag_inclusions, "Diagnosis inclusions"))
  )

  icd_mode <- trim_value(options$icd_mode)
  if (identical(icd_mode, "ICD8")) args <- c(args, "--ICD8")
  if (identical(icd_mode, "ICD10")) args <- c(args, "--ICD10")
  if (isTRUE(options$single_phenotype)) {
    args <- c(args, add_value("--name", validate_name(options$name)))
  }

  exclusion_flags <- c(
    lifetime = "--LifetimeExclusion",
    post = "--PostExclusion",
    prior = "--OneyPriorExclusion"
  )
  for (key in names(exclusion_flags)) {
    path <- exclusions[[key]] %||% ""
    if (nzchar(trim_value(path))) args <- c(args, exclusion_flags[[key]], path)
  }

  list(
    executable = python,
    args = args,
    display = paste(c(shQuote(python), shQuote(args)), collapse = " "),
    run_id = run_id,
    output_file = output_file
  )
}

safe_slurm_value <- function(value, label, pattern) {
  value <- trim_value(value)
  if (!nzchar(value)) return("")
  if (!grepl(pattern, value)) stop("Invalid ", label, " in project registry")
  value
}

write_slurm_script <- function(record, run_id, run_dir, command) {
  if (!is_valid_run_id(run_id)) stop("Invalid run identifier")
  if (grepl("[[:space:]]", run_dir)) {
    stop("Slurm run directories must not contain whitespace")
  }
  partition <- safe_slurm_value(record$Partition, "partition", "^[A-Za-z0-9_.-]+$")
  account <- safe_slurm_value(record$Account, "account", "^[A-Za-z0-9_.-]+$")
  time <- safe_slurm_value(record$Time, "time", "^[0-9]+(?:-[0-9]{1,2})?:[0-9]{2}:[0-9]{2}$|^[0-9]{1,2}:[0-9]{2}:[0-9]{2}$")
  memory <- safe_slurm_value(record$Memory, "memory", "^[0-9]+[KMGTP]?[Bb]?$|^[0-9]+$")
  cpus <- safe_slurm_value(record$Cpus %||% "1", "CPU count", "^[1-9][0-9]*$")

  directives <- c(
    paste0("#SBATCH --job-name=", substr(gsub("[^A-Za-z0-9_-]", "_", run_id), 1L, 64L)),
    paste0("#SBATCH --output=", file.path(run_dir, "slurm-%j.out")),
    paste0("#SBATCH --error=", file.path(run_dir, "slurm-%j.err")),
    if (nzchar(partition)) paste0("#SBATCH --partition=", partition),
    if (nzchar(account)) paste0("#SBATCH --account=", account),
    if (nzchar(time)) paste0("#SBATCH --time=", time),
    if (nzchar(memory)) paste0("#SBATCH --mem=", memory),
    paste0("#SBATCH --cpus-per-task=", cpus)
  )
  directives <- directives[!vapply(directives, is.null, logical(1))]

  command_line <- paste(c(shQuote(command$executable), shQuote(command$args)), collapse = " ")
  script <- c(
    "#!/bin/bash",
    directives,
    "",
    "set -uo pipefail",
    "umask 077",
    paste("cd", shQuote(run_dir)),
    paste("printf '%s\\n' RUNNING >", shQuote(file.path(run_dir, "status.txt"))),
    "set +e",
    command_line,
    "exit_code=$?",
    "set -e",
    paste0("if [ \"$exit_code\" -eq 0 ]; then printf '%s\\n' COMPLETED > ",
           shQuote(file.path(run_dir, "status.txt")),
           "; else printf '%s\\n' FAILED > ",
           shQuote(file.path(run_dir, "status.txt")), "; fi"),
    paste("printf '%s\\n' \"$exit_code\" >", shQuote(file.path(run_dir, "exit_code.txt"))),
    "exit \"$exit_code\""
  )
  path <- file.path(run_dir, paste0(run_id, ".sbatch"))
  writeLines(script, path, useBytes = TRUE)
  Sys.chmod(path, mode = "0700")
  path
}

submit_slurm_job <- function(script_path, sbatch = "sbatch") {
  sbatch <- trim_value(sbatch %||% "sbatch")
  executable <- Sys.which(sbatch)
  if (!nzchar(executable)) stop("sbatch executable was not found: ", sbatch)
  output <- suppressWarnings(system2(
    executable,
    args = c("--parsable", shQuote(script_path)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status") %||% 0L
  if (!identical(as.integer(status), 0L)) {
    stop("Slurm submission failed: ", paste(output, collapse = "\n"))
  }
  job_id <- sub(";.*$", "", trimws(output[[1]] %||% ""))
  if (!grepl("^[0-9]+$", job_id)) {
    stop("Could not parse the Slurm job ID from: ", paste(output, collapse = " "))
  }
  job_id
}

write_run_metadata <- function(path, values) {
  values <- vapply(values, function(value) gsub("[\t\r\n]", " ", as.character(value)), "")
  write.table(
    data.frame(key = names(values), value = unname(values)),
    path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  Sys.chmod(path, mode = "0600")
  invisible(path)
}

read_run_metadata <- function(path) {
  if (!file.exists(path)) return(list())
  table <- read.delim(path, header = FALSE, col.names = c("key", "value"),
                      quote = "", stringsAsFactors = FALSE)
  as.list(stats::setNames(table$value, table$key))
}

find_run_directory <- function(registry, user, run_id) {
  if (!is_valid_run_id(run_id)) return(NULL)
  roots <- unique(file.path(registry$OutputRoot, user))
  for (root in roots) {
    candidate <- file.path(root, run_id)
    if (!dir.exists(candidate)) next
    normalized_root <- normalizePath(root, mustWork = TRUE)
    normalized_candidate <- normalizePath(candidate, mustWork = TRUE)
    if (startsWith(normalized_candidate, paste0(normalized_root, .Platform$file.sep))) {
      return(normalized_candidate)
    }
  }
  NULL
}

query_job_state <- function(run_dir) {
  status_file <- file.path(run_dir, "status.txt")
  local_status <- ""
  if (file.exists(status_file)) {
    local_status <- trim_value(readLines(status_file, n = 1L, warn = FALSE))
    if (local_status %in% c("RUNNING", "COMPLETED", "FAILED", "SUBMISSION_FAILED")) {
      return(local_status)
    }
  }
  metadata <- read_run_metadata(file.path(run_dir, "metadata.tsv"))
  job_id <- trim_value(metadata$job_id)
  if (!grepl("^[0-9]+$", job_id)) return("SUBMISSION_UNKNOWN")
  squeue <- Sys.which("squeue")
  if (nzchar(squeue)) {
    state <- tryCatch(
      system2(squeue, c("-h", "-j", job_id, "-o", "%T"), stdout = TRUE, stderr = FALSE),
      error = function(error) character()
    )
    state <- trim_value(state)
    if (nzchar(state)) return(state)
  }
  sacct <- Sys.which("sacct")
  if (nzchar(sacct)) {
    state <- tryCatch(
      system2(
        sacct,
        c("-n", "-X", "-j", job_id, "--format=State", "--parsable2"),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(error) character()
    )
    state <- trim_value(state)
    state <- sub("[|+[:space:]].*$", "", state)
    if (nzchar(state)) return(state)
  }
  if (nzchar(local_status)) local_status else "NOT_IN_QUEUE"
}

list_run_files <- function(run_dir) {
  files <- list.files(run_dir, full.names = TRUE, recursive = FALSE, all.files = FALSE)
  files[file.info(files)$isdir %in% FALSE]
}
