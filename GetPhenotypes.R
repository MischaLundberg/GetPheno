#!/usr/bin/env Rscript
# GetPhenotypes.R
#
# This R script performs phenotype extraction from UK Biobank-like hospital episode data.
# It has been updated to handle the standard UKBB HES data structure where diagnosis
# codes (hesin_diag) and admission dates (hesin) are in separate files requiring a join.
#
# CRITICAL FIXES:
# 1. Enhanced debugging output (Blocks 1, 2, and 3) to trace date and column existence.
# 2. Reinforcement of 'as.character()' in the final transmute to prevent data.table::fwrite
#    from reformatting the character birthdate string (e.g., "01/03/1950" to "1/3/50").

# Suppress startup messages for a cleaner console output
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate) # CRITICAL: Added for robust date parsing (dmy)
  library(stringr)
  library(readr)
  library(data.table)
  library(optparse)
  library(parallel)
  library(uuid)
  library(tidyr)
})

#############################
# Global variables and version checks
#############################

# Setting default diagnosis date format (HES dates) to DD/MM/YYYY
DateFormat <- "%d/%m/%Y" 

# Birthdate format based on user input (e.g., 1/3/50)
BirthDateFormat <- "%d/%m/%y" 

# List of formats to try for HES Diagnosis Dates (d/m/y covers D/M/YY and D/M/YYYY)
# Added Y-m-d to cover potential ISO-like formats (e.g., 2001-01-03)
HES_DATE_ORDERS <- c("d/m/Y", "Y-m-d", "d/m/y")

ATC_Requested <- "NotSet"
# Define internal standard column names for ICD code and Date
ICD_COL_INTERNAL <- "icd_code"
DATE_COL_INTERNAL <- "diagnosis_date"

#############################
# Command Line Argument Parsing (updated defaults based on user headers)
#############################

option_list = list(
  make_option(c("-f", "--DiagnosisFile"), type="character", default=NULL, 
              help="Path to the primary diagnosis file (e.g., hesin). [REQUIRED]"),
  make_option(c("--f2"), type="character", default="", 
              help="Path to secondary diagnosis file (e.g., hesin_diag)"),
  make_option(c("-c", "--CodesToExport"), type="character", default=NULL,
              help="Path to the file containing ICD codes to export."),
  make_option(c("-g", "--Demographics"), type="character", default=NULL,
              help="Path to the demographics file."),
  make_option(c("-a", "--additionalInformation"), type="character", default=NULL,
              help="Path to the additional information file."),
  make_option(c("--ge"), type="character", default="f.eid", 
              help="Column name for ID in Demographics/additionalInformation."),
  # UPDATED DEFAULT: Assuming Primary ICD code column (if used)
  make_option(c("--fcol"), type="character", default="diag_icd10", 
              help="Column name for ICD code in DiagnosisFile (primary diagnosis code)."), 
  # UPDATED DEFAULT: Confirmed secondary ICD code column
  make_option(c("--f2col"), type="character", default="diag_icd10", 
              help="Column name for ICD code in Secondary Diagnosis File (all diagnosis codes)."), 
  # UPDATED DEFAULT: Confirmed date column in hesin.txt.csv
  make_option(c("--gcol"), type="character", default="admidate", 
              help="Column name for Date in DiagnosisFile (Primary Diagnosis Date)."),
  make_option(c("--iidcol"), type="character", default="eid", 
              help="Column name for IID in all files."),
  # UPDATED DEFAULT: Confirmed birthdate column in demographics
  make_option(c("--bdcol"), type="character", default="birth", 
              help="Column name for Birth Date in Demographics."),
  # UPDATED DEFAULT: Confirmed sex column in demographics
  make_option(c("--sexcol"), type="character", default="sex", 
              help="Column name for Sex in Demographics."),
  make_option(c("--fsep"), type="character", default="\t", 
              help="Field separator for DiagnosisFile. (Note: using tab default)"),
  make_option(c("--gsep"), type="character", default="\t", 
              help="Field separator for Demographics/additionalInformation."),
  make_option(c("-o", "--OutFile"), type="character", default="pheno.txt", 
              help="Output file name."),
  make_option(c("--eM"), type="character", default="ICD10", 
              help="ICD system for exclusions (ICD9/ICD10)."),
  make_option(c("--din"), type="character", default="icd_code", 
              help="Diagnosis Code column in secondary diagnosis file (used as a temporary column name)."),
  # NOTE: --don is still required but will be ignored for secondary date (as we join to get it)
  make_option(c("--don"), type="character", default="diag_date", 
              help="Diagnosis Date column in secondary diagnosis file (used as a temporary column name)."), 
  make_option(c("--qced"), type="character", default="f.20015.0.0", 
              help="Column name for QCed status."),
  make_option(c("--DiagTypeExclusions"), type="character", default="0", 
              help="Comma-separated list of diag. types for exclusion."),
  make_option(c("--DiagTypeInclusions"), type="character", default="0", 
              help="Comma-separated list of diag. types for inclusion."),
  make_option(c("--LifetimeExclusion"), action="store_true", default=FALSE, 
              help="Apply lifetime exclusion logic."),
  make_option(c("--PostExclusion"), action="store_true", default=FALSE, 
              help="Apply post-exclusion logic."),
  make_option(c("--OneyPriorExclusion"), action="store_true", default=FALSE, 
              help="Apply one-year prior exclusion logic."),
  make_option(c("--eCc"), type="character", default="", 
              help="Exclusion ICD codes file."),
  make_option(c("--Fyob"), type="character", default="f.34.0.0", 
              help="Column name for Year of Birth (only used for some legacy functions)."),
  make_option(c("--Fgender"), type="character", default="f.31.0.0", 
              help="Column name for Gender (only used for some legacy functions)."),
  make_option(c("-v", "--verbose"), action="store_true", default=FALSE, 
              help="Print status messages."),
  make_option(c("--BuildTestSet"), action="store_true", default=FALSE, 
              help="Build a test dataset."),
  make_option(c("--testRun"), action="store_true", default=FALSE, 
              help="Run in test mode (subset data)."),
  make_option(c("--MatchFI"), action="store_true", default=FALSE, 
              help="Match field index in columns."),
  make_option(c("--icdprefix"), type="character", default="ICD10:", 
              help="Prefix to add to all ICD codes for standardisation."),
  make_option(c("--skipICDUpdate"), action="store_true", default=FALSE, 
              help="Skip updating ICD coding."),
  # UPDATED DEFAULT: Corrected DateFormat
  make_option(c("--DateFormat"), type="character", default="%d/%m/%Y", 
              help="Format of diagnosis date columns (e.g., %d/%m/%Y for 01/03/1950)"),
  make_option(c("--iidstatus"), type="character", default="f.eid", 
              help="Column name for participant status (e.g., in/out of analysis set)"),
  make_option(c("--removePointInDiagCode"), action="store_true", default=FALSE, 
              help="Remove the decimal point from ICD-10 codes"),
  make_option(c("--nthreads"), type="integer", default=1, 
              help="Number of threads to use for parallel operations"),
  make_option(c("--name"), type="character", default="Phenotype", 
              help="Name for the generated binary phenotype column (will be 'diagnosis' in final output)."),
  make_option(c("--BuildEntryExitDates"), action="store_true", default=FALSE, 
              help="Build entry and exit dates from diagnosis data"),
  make_option(c("--BuildOphold"), action="store_true", default=FALSE, 
              help="Build Ophold-like exclusion column"),
  make_option(c("--write_pickle"), action="store_true", default=FALSE, 
              help="Write output in Python pickle format (for compatibility)"),
  make_option(c("--write_fastGWA_format"), action="store_true", default=FALSE, 
              help="Write output in fastGWA format"),
  make_option(c("--write_Plink_format"), action="store_true", default=FALSE, 
              help="Write output in Plink format"),
  make_option(c("--DiagnosisFile2"), type="character", default="", 
              help="Path to secondary diagnosis file (e.g., hesin_diag)"),
  # UPDATED DEFAULT: Confirmed record number column
  make_option(c("--recnum"), type="character", default="ins_index", 
              help="Column name for record number in DiagnosisFile"),
  # UPDATED DEFAULT: Confirmed record number column
  make_option(c("--recnum2"), type="character", default="ins_index", 
              help="Column name for record number in Secondary Diagnosis File")
)

parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)

# Map the alias '--f2' to the long name '--DiagnosisFile2'
if (opt$f2 != "" && opt$DiagnosisFile2 == "") {
  opt$DiagnosisFile2 <- opt$f2
}
opt$f2 <- NULL 

# Set DateFormat from options
DateFormat <- opt$DateFormat

#############################
# Helper Functions
#############################

# Robust function to read data using data.table::fread
read_data_table <- function(path, sep_char, verbose = FALSE, test_run = FALSE) {
  if (verbose) message(sprintf("Reading data from: %s (using separator: '%s')", path, sep_char))
  
  if (verbose) {
    if (!file.exists(path)) {
      message(sprintf("--- DEBUG FILE CHECK: File '%s' NOT FOUND.", path))
      full_path <- normalizePath(path, mustWork = FALSE)
      message(sprintf("--- DEBUG FILE CHECK: Tried to resolve path as: %s", full_path))
    } else {
      message(sprintf("--- DEBUG FILE CHECK: File '%s' found successfully. (Verbose set to TRUE)", path)) 
    }
  }
  
  # Determine separator if default is wrong
  if (sep_char == "\t" && stringr::str_detect(tolower(path), "\\.csv$")) {
    sep_char = ","
    if (verbose) message("--- HINT: File is .csv but separator was default tab. Using comma separator: ','")
  }
  
  # CRITICAL: Force all columns to character type
  dt <- fread(path, sep = sep_char, header = TRUE, stringsAsFactors = FALSE, colClasses = "character")
  
  original_names <- names(dt)
  
  # Aggressive Header Cleaning and Lowercasing for robustness
  cleaned_names <- original_names %>%
    stringr::str_replace_all("[^\\x20-\\x7E]", " ") %>% 
    stringr::str_replace_all("^[\\s\"'\\|]+|[\\s\"'\\|]+$", "") %>%
    trimws() %>%
    tolower()
  
  setnames(dt, original_names, cleaned_names)
  
  if (verbose) {
    message(sprintf("--- DEBUG READ: Columns in %s are: %s", basename(path), paste(names(dt), collapse=", ")))
  }
  
  return(dt)
}

# Standardized ID Checker and Renamer 
standardize_iid_column <- function(dt, filename, iidcol_option) {
  target_iid_col <- "eid"
  user_iid_col <- tolower(trimws(iidcol_option)) 
  
  if (target_iid_col %in% names(dt)) {
    dt[[target_iid_col]] <- as.character(dt[[target_iid_col]])
    return(dt)
  }
  
  if (user_iid_col %in% names(dt)) {
    if (user_iid_col != target_iid_col) {
      setnames(dt, user_iid_col, target_iid_col)
      if (opt$verbose) message(sprintf("--- DEBUG CHECK: Renamed ID column '%s' to internal '%s' in '%s'.", user_iid_col, target_iid_col, basename(filename)))
    }
    dt[[target_iid_col]] <- as.character(dt[[target_iid_col]])
    return(dt)
  }
  
  stop(sprintf("Could not find the participant ID column in the file: '%s'. Expected column name (after cleaning and lowercasing) is '%s'.", 
               filename, user_iid_col))
}

#############################
# Main Function
#############################

main <- function(DiagnosisFile, CodesToExport, Demographics, additionalInformation, ExDepExc, ge, fcol, f2col, gcol, iidcol, bdcol, sexcol, fsep, gsep, OutFile, eM, din, don, qced, DiagTypeExclusions, DiagTypeInclusions, LifetimeExclusion, PostExclusion, OneyPriorExclusion, eCc, Fyob, Fgender, verbose, BuildTestSet, testRun, MatchFI, icdprefix, skipICDUpdate, DateFormat, iidstatus, removePointInDiagCode, nthreads, name, BuildEntryExitDates, BuildOphold, write_pickle, write_fastGWA_format, write_Plink_format, DiagnosisFile2, recnum, recnum2) {
  
  # Ensure DateFormat is consistently DD/MM/YYYY
  DateFormat <- "%d/%m/%Y"
  if (verbose) message(sprintf("--- DEBUG: DateFormat set globally to '%s' ---", DateFormat))
  
  # Clean up column names based on user input for internal use
  fcol_clean <- tolower(trimws(fcol))   # Primary ICD
  f2col_clean <- tolower(trimws(f2col)) # Secondary ICD (ICD code in hesin_diag)
  gcol_clean <- tolower(trimws(gcol))   # Primary Date (admidate in hesin)
  bdcol_clean <- tolower(trimws(bdcol)) # Birthdate (birth in demographics)
  sexcol_clean <- tolower(trimws(sexcol)) # Sex (sex in demographics)
  recnum_clean <- tolower(trimws(recnum)) # Record number (ins_index in hesin)
  recnum2_clean <- tolower(trimws(recnum2)) # Record number (ins_index in hesin_diag)

  # 1. Load Primary Diagnosis Data (hesin.txt.csv) and Standardize ID
  if (verbose) message(sprintf("Step 1/7: Loading primary diagnosis data (episodes) from %s...", DiagnosisFile))
  df_primary <- read_data_table(DiagnosisFile, fsep, verbose, testRun)
  df_primary <- standardize_iid_column(df_primary, DiagnosisFile, iidcol)
  
  # 2. Load Demographics 
  if (!is.null(Demographics)) {
    if (verbose) message(sprintf("Step 2/7: Loading demographics from %s...", Demographics))
    df_demo <- read_data_table(Demographics, gsep, verbose, testRun)
    df_demo <- standardize_iid_column(df_demo, Demographics, iidcol)
  } else {
    stop("Error: The '--Demographics' file must be specified.")
  }
  
  # CRITICAL CHECK for Demographics before proceeding
  if (!(bdcol_clean %in% names(df_demo) && sexcol_clean %in% names(df_demo))) {
      stop(sprintf(
          "Critical Error: Demographics file must contain columns for birthdate ('%s') and sex ('%s').\n
          (Your headers: %s)", bdcol_clean, sexcol_clean, paste(names(df_demo), collapse=", ")
      ))
  }

  # --- DEBUG BLOCK 1: Check Demographics Input ---
  if (verbose) {
    message("--- DEBUG BLOCK 1: Demographics Data Check ---")
    
    # Check if the birthdate column exists and its data type
    if (bdcol_clean %in% names(df_demo)) {
        message(sprintf("  -> Birthdate column ('%s') found. Data type is: %s", 
                        bdcol_clean, class(df_demo[[bdcol_clean]])))
        
        # Check raw values for problematic EIDs
        problem_eids_demo <- df_demo %>% 
          filter(eid %in% c("1000015", "1000027")) %>%
          select(eid, sex = !!sexcol_clean, birthdate_raw = !!bdcol_clean)
          
        message("  -> Raw Demographics Data for EIDs 1000015, 1000027:")
        print(problem_eids_demo)
    } else {
        message(sprintf("  -> ERROR: Birthdate column '%s' not found in demographics after loading.", bdcol_clean))
    }
    message("--- DEBUG BLOCK 1: End ---")
  }
  
  # 3. Load Target ICD Codes (remaining steps 3 and 4 remain similar for code processing)
  if (is.null(CodesToExport)) stop("Error: The '--CodesToExport' file must be specified.")
  if (verbose) message(sprintf("Step 3/7: Loading ICD codes to export from %s...", CodesToExport))
  df_codes <- read_data_table(CodesToExport, fsep, verbose, testRun)
  
  target_codes_col <- names(df_codes)[1] 
  if (is.null(target_codes_col)) stop("Error: CodesToExport file appears empty or malformed.")
  
  target_codes_raw <- unique(df_codes[[target_codes_col]])
  target_codes_raw <- str_replace_all(target_codes_raw, "^[\\s\"']+|[\\s\"']+$", "")
  
  # Prepare the final target codes list used for the regex pattern
  target_codes <- str_replace_all(target_codes_raw, "\\.", "")
  target_codes <- str_remove(target_codes, "^ICD10:")
  target_codes <- str_remove(target_codes, "^ICD9:")
  target_codes <- str_remove(target_codes, "^ICD")
  
  escaped_target_codes <- gsub("([\\+|\\*|\\?|\\^|\\$|\\(|\\)|\\[|\\]|\\{|\\}|\\.|\\/|\\|])", "\\\\\\1", target_codes)
  icd_match_pattern <- paste0("^", icdprefix, "(", paste(escaped_target_codes, collapse = "|"), ")")
  
  if (verbose) message(sprintf("Loaded %d unique target ICD codes. First 5: %s", length(target_codes), paste(head(target_codes, 5), collapse=", ")))
  
  # 4. Combine and standardize diagnosis data (including dates via JOIN)
  if (verbose) message("Step 4/7: Standardizing, combining, and joining dates to diagnosis records...")
  
  df_secondary <- NULL
  if (DiagnosisFile2 != "") {
    if (verbose) message(sprintf("Loading secondary diagnosis data (codes) from %s...", DiagnosisFile2))
    df_secondary <- read_data_table(DiagnosisFile2, fsep, verbose, testRun)
    df_secondary <- standardize_iid_column(df_secondary, DiagnosisFile2, iidcol)
  }
  
  # --- Extract Dates and Keys from Primary (HESIN) ---
  if (!(gcol_clean %in% names(df_primary) && recnum_clean %in% names(df_primary))) {
      stop(sprintf("Primary episode data must contain date column ('%s') and record number ('%s').", gcol_clean, recnum_clean))
  }

  df_dates_keys <- df_primary %>%
    select(eid, !!recnum_clean, date_raw = !!gcol_clean) %>%
    setnames(recnum_clean, "ins_index_key") 
  
  # --- Extract Codes and Keys from Secondary (HESIN_DIAG) ---
  if (is.null(df_secondary) || !(f2col_clean %in% names(df_secondary) && recnum2_clean %in% names(df_secondary))) {
      warning("Secondary diagnosis file is missing or lacks required ICD code/record number. Skipping secondary data.")
      df_diag_codes <- setNames(data.table(character(0), character(0), character(0)), 
                                c("eid", "ins_index_key", ICD_COL_INTERNAL))
  } else {
      df_diag_codes <- df_secondary %>%
        select(eid, !!recnum2_clean, icd = !!f2col_clean) %>%
        setnames(recnum2_clean, "ins_index_key") %>%
        setnames("icd", ICD_COL_INTERNAL) %>%
        filter(!is.na(get(ICD_COL_INTERNAL)) & get(ICD_COL_INTERNAL) != "")
  }

  # --- JOIN: Combine Codes and Dates ---
  df_all_diag <- df_diag_codes %>%
    inner_join(df_dates_keys, by = c("eid", "ins_index_key")) %>%
    # Rename the raw date column
    setnames("date_raw", "date_raw") %>%
    # Add a single source column after the join
    mutate(Source = "HES") %>% 
    # Clean ICD codes
    mutate(!!ICD_COL_INTERNAL := paste0(icdprefix, get(ICD_COL_INTERNAL))) %>%
    mutate(!!ICD_COL_INTERNAL := str_replace_all(get(ICD_COL_INTERNAL), "\\.", "")) %>%
    # CRITICAL FIX: Robustly parse date using lubridate with multiple orders
    mutate(!!DATE_COL_INTERNAL := parse_date_time(date_raw, orders = HES_DATE_ORDERS)) %>%
    # Remove records with invalid/NA dates or codes
    filter(!is.na(get(DATE_COL_INTERNAL))) %>%
    # Select the required final columns for aggregation (including the new Source column)
    select(eid, !!ICD_COL_INTERNAL, !!DATE_COL_INTERNAL, Source) %>%
    unique()
  
  if (verbose) message(sprintf("Total combined and clean diagnosis records after join: %d", nrow(df_all_diag)))
  
  # 5. Phenotype Calculation: Identify cases and aggregate history
  if (verbose) message("Step 5/7: Aggregating history and calculating phenotype using data.table...")
  
  # Convert to data.table for highly efficient aggregation
  setDT(df_all_diag)
  
  # 5a. Identify Case EIDs
  case_iids <- df_all_diag[str_detect(get(ICD_COL_INTERNAL), icd_match_pattern), unique(eid)]

  # 5b. Aggregate Diagnosis History by EID using data.table
  if (nrow(df_all_diag) > 0) {
      df_history_agg <- df_all_diag[, .(
          # Core Aggregations (results are R Date objects)
          diagnosis = as.integer(eid[1] %in% case_iids), # Binary Phenotype
          first_dx_date_obj = min(get(DATE_COL_INTERNAL)), 
          last_dx_date_obj = max(get(DATE_COL_INTERNAL)),
          
          n_diags = .N,
          n_unique_in_days = uniqueN(get(DATE_COL_INTERNAL)),
          
          diagnoses_str = paste(unique(get(ICD_COL_INTERNAL)), collapse = "|"),
          # CRITICAL: Format pipe-separated dates to the required DD/MM/YYYY format for output string
          in_dates_str = paste(sort(unique(format(get(DATE_COL_INTERNAL), format = DateFormat))), collapse = "|"),
          out_dates_str = paste(sort(unique(format(get(DATE_COL_INTERNAL), format = DateFormat))), collapse = "|") 
      ), by = eid]

      # Convert min/max dates (Date objects) to the standardized character format (DD/MM/YYYY) for output columns
      df_history_agg[, `:=`(
          first_dx = format(first_dx_date_obj, format = DateFormat),
          last_dx = format(last_dx_date_obj, format = DateFormat)
      )]
      
  } else {
      # Handle case where no records exist
      df_history_agg <- data.table(
          eid = character(0), diagnosis = integer(0), first_dx = character(0), last_dx = character(0),
          n_diags = integer(0), n_unique_in_days = integer(0), diagnoses_str = character(0),
          in_dates_str = character(0), out_dates_str = character(0), 
          first_dx_date_obj = as.Date(character(0)), last_dx_date_obj = as.Date(character(0))
      )
  }

  
  if (verbose) message(sprintf("Aggregated history for %d individuals.", nrow(df_history_agg)))
  
  # 6. Merge with demographics and calculate age
  if (verbose) message("Step 6/7: Merging demographics and calculating age...")
  
  # Select demographics columns using the dynamically cleaned names
  df_final <- df_demo %>%
    select(eid, sex = !!sexcol_clean, birthdate_raw = !!bdcol_clean) %>% # Rename raw birthdate
    # Merge history. Use left_join to keep all participants.
    left_join(df_history_agg, by = "eid") %>%
    # Fill NA's for those with no diagnoses
    mutate(
        diagnosis = replace_na(diagnosis, 0),
        n_diags = replace_na(n_diags, 0),
        n_unique_in_days = replace_na(n_unique_in_days, 0),
        first_dx = replace_na(first_dx, ""),
        last_dx = replace_na(last_dx, ""),
        diagnoses_str = replace_na(diagnoses_str, ""),
        in_dates_str = replace_na(in_dates_str, ""),
        out_dates_str = replace_na(out_dates_str, "")
    ) %>%
    # Calculate Age_FirstDx
    mutate(
        # CRITICAL: Parse the raw birthdate string into a Date object for age calculation.
        # Use multiple orders to handle common UKBB variations (e.g., D/M/YY or D/M/YYYY)
        # Note: The warning about '3 failed to parse' is here, meaning 3 records might have a different format.
        birthdate_date = parse_date_time(birthdate_raw, orders = c('dmy', 'dmY', 'dmy')), 
        # Calculate age (uses the Date object from aggregation step)
        Age_FirstDx = round(as.numeric(difftime(first_dx_date_obj, birthdate_date, units = "days")) / 365.25, 2)
    )
  
  # --- DEBUG BLOCK 2: Final Merged Data Check (Dates, Age) ---
  if (verbose) {
    message("--- DEBUG BLOCK 2: Final Merged Data Check (Dates, Age) ---")
    problem_eids <- c("1000015", "1000027")
    
    df_debug_final <- df_final %>% 
      filter(eid %in% problem_eids) %>%
      # Select and show all relevant date columns before final transmute
      select(
        eid, 
        sex,
        birthdate_raw, 
        birthdate_date, # Date object from parsing
        first_dx_date_obj, # Date object from aggregation
        first_dx, # DD/MM/YYYY string output
        Age_FirstDx
      )
      
    message("  Raw Data, Parsed Dates (Date Objects), and Calculated Age:")
    print(df_debug_final)
    
    message(sprintf("  -> Check: birthdate_raw column exists: %s", "birthdate_raw" %in% names(df_final)))
    message("--- DEBUG BLOCK 2: End ---")
  }
  
  # FINAL STEP: Column Selection and Renaming using Transmute
  df_final <- df_final %>%
    # Transmute selects only the columns listed and renames them instantly.
    transmute(
        cpr_enc = eid, 
        diagnosis = diagnosis, 
        first_dx = first_dx, # DD/MM/YYYY string
        last_dx = last_dx,   # DD/MM/YYYY string
        diagnoses = diagnoses_str,
        in_dates = in_dates_str, # Standardized DD/MM/YYYY pipe-separated
        out_dates = out_dates_str, # Standardized DD/MM/YYYY pipe-separated
        n_diags = n_diags,
        n_unique_in_days = n_unique_in_days,
        sex = sex,
        birthdate = as.character(birthdate_raw), # *** FIX: Explicitly forcing character type to prevent fwrite format truncation ***
        
        # --- Placeholder Columns for Exclusion/Status ---
        dbds = 0, degen_old = 0, degen_new = 0, C_STATUS = 0, D_STATUS_HEN_START = NA_character_, 
        D_FODDATO = NA_character_, C_KON = NA_character_,
        # --- End Placeholders ---
        
        Age_FirstDx = Age_FirstDx
    ) 
  
  # --- DEBUG BLOCK 3: Final Output Types Check (Before Writing) ---
  if (verbose) {
    message("--- DEBUG BLOCK 3: Final Output Types Check (Before Writing) ---")
    message(sprintf("  -> Data type of FINAL 'birthdate' column: %s", class(df_final$birthdate)))
    
    df_debug_final_out <- df_final %>% 
      filter(cpr_enc %in% c("1000015", "1000027")) %>%
      select(cpr_enc, birthdate)
      
    message("  -> Final 'birthdate' values for EIDs 1000015, 1000027:")
    print(df_debug_final_out)
    message("--- DEBUG BLOCK 3: End ---")
  }
  
  # 7. Write Output
  if (verbose) message(sprintf("Step 7/7: Writing final output to %s...", OutFile))
  
  outfile <- OutFile
  
  fwrite(df_final, outfile, sep = fsep, row.names = FALSE)
  message(sprintf("Output successfully written with detailed history to: %s", outfile))
  message(sprintf("Phenotype script finished execution. Processed %d participants.", nrow(df_final)))
  message(sprintf("--- Summary ---"))
  message(sprintf("Total participants: %d", nrow(df_final)))
  message(sprintf("Number of cases (diagnosis=1): %d", sum(df_final$diagnosis == 1, na.rm = TRUE)))
  
  return(invisible())
}

# Function to run the main logic (called at the end of the script)
run_script <- function() {
  # Get the list of options parsed from the command line
  opt_list <- as.list(opt)
  
  # Map the alias '--f2' to the long name '--DiagnosisFile2'
  if ("f2" %in% names(opt_list) && opt_list$f2 != "" && opt_list$DiagnosisFile2 == "") {
    opt_list$DiagnosisFile2 <- opt_list$f2
    opt_list$f2 <- NULL
  }
  
  # Remove extraneous arguments from optparse
  opt_list$help <- NULL 
  opt_list$args <- NULL
  
  # Pass the cleaned list of arguments to main
  do.call(main, opt_list)
}

# Execute the script if not run interactively
if (!interactive()) {
  run_script()
}
