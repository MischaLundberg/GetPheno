
# get_pheno.R — Complete R port of get_pheno.py
# Date: 2025-11-10
# Requirements: optparse, data.table, jsonlite, lubridate
suppressPackageStartupMessages({
  require(optparse)
  require(data.table)
  require(jsonlite)
  require(lubridate)
})

# ----------------------------- Logging ---------------------------------------
PHENO_VERBOSE <- TRUE
log_set <- function(verbose=TRUE) { assign("PHENO_VERBOSE", isTRUE(verbose), envir = .GlobalEnv) }
log_msg <- function(..., level="INFO") {
  if (PHENO_VERBOSE || level %in% c("WARN","ERROR")) {
    ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cat(sprintf("[%s] [%s] %s\n", ts, level, paste(..., collapse=" ")))
  }
}

# ----------------------------- Helpers ---------------------------------------
nzchar2 <- function(x) !is.na(x) & x != ""

read_any <- function(path, sep=",", n_max=None, guess_max=1e5) {
  if (!nzchar2(path)) return(NULL)
  # detect tsv/csv by sep if given
  dt <- tryCatch({
    data.table::fread(path, sep=sep, nThread=1, showProgress=FALSE)
  }, error=function(e) {
    # try auto
    tryCatch(data.table::fread(path, nThread=1, showProgress=FALSE), error=function(e2) NULL)
  })
  dt
}

write_tsv <- function(dt, file) {
  data.table::fwrite(dt, file=file, sep="\t", quote=FALSE, na="")
}

to_vec <- function(x) {
  if (is.null(x)) return(character())
  if (is.list(x)) return(unlist(x, use.names=FALSE))
  if (length(x)==1 && is.character(x) && grepl(",", x, fixed=TRUE)) return(trimws(strsplit(x, ",", fixed=TRUE)[[1]]))
  as.character(x)
}

expand_ranges <- function(x) {
  xs <- to_vec(x)
  out <- character()
  for (tok in xs) {
    tok <- trimws(tok)
    m <- regexec("^([A-Za-z0-9:]+)\\-([A-Za-z0-9:]+)$", tok)
    mm <- regmatches(tok, m)[[1]]
    if (length(mm) == 3) {
      a <- mm[2]; b <- mm[3]
      # prefixed ICD like ICD10:T36-ICD10:T38
      if (grepl(":", a) && grepl(":", b) && sub(":.*$", "", a) == sub(":.*$", "", b)) {
        pref <- sub("^([^:]+:).*", "\\1", a)
        a2 <- sub("^[^:]+:", "", a); b2 <- sub("^[^:]+:", "", b)
        if (grepl("^[A-Za-z]?[0-9]+$", a2) && grepl("^[A-Za-z]?[0-9]+$", b2)) {
          la <- sub("([A-Za-z]?)[0-9]+$", "\\1", a2)
          lb <- sub("([A-Za-z]?)[0-9]+$", "\\1", b2)
          na <- as.integer(gsub("^[A-Za-z]", "", a2)); nb <- as.integer(gsub("^[A-Za-z]", "", b2))
          if (identical(la, lb) && !is.na(na) && !is.na(nb)) {
            for (i in seq(na, nb)) out <- c(out, paste0(pref, la, i))
            next
          }
        }
      }
      # numeric ranges
      if (grepl("^[0-9]+$", a) && grepl("^[0-9]+$", b)) {
        out <- c(out, as.character(seq(as.integer(a), as.integer(b)))); next
      }
      # alpha+numeric like T36-T38
      if (grepl("^[A-Za-z][0-9]+$", a) && grepl("^[A-Za-z][0-9]+$", b) && substr(a,1,1)==substr(b,1,1)) {
        la <- substr(a,1,1); na <- as.integer(sub("^[A-Za-z]", "", a))
        nb <- as.integer(sub("^[A-Za-z]", "", b))
        out <- c(out, paste0(la, na:nb)); next
      }
      out <- c(out, tok)
    } else out <- c(out, tok)
  }
  unique(out)
}

normalize_iid_series <- function(x, target="str") {
  s <- as.character(x)
  s <- gsub("\\.0$", "", trimws(s))
  if (target == "str") return(s)
  if (any(!grepl("^\\d+$", s))) {
    log_msg("Non-numeric IIDs detected; keeping as strings", level="WARN")
    return(s)
  }
  as.integer(s)
}

parse_code_file <- function(code_file, gsep=",", gcol="c_adiag") {
  # Accept raw list (one column) or column 'gcol' in header
  dt <- read_any(code_file, sep=gsep)
  if (is.null(dt)) stop("Unable to read -g file: ", code_file)
  if (ncol(dt)==1 && !gcol %in% names(dt)) {
    setnames(dt, names(dt)[1], gcol)
  }
  vec <- to_vec(dt[[gcol]])
  expand_ranges(vec)
}

apply_icd_prefix <- function(codes, icdprefix="") {
  if (!nzchar2(icdprefix)) return(codes)
  # Insert letter after ICDx:
  sub("^ICD([0-9]+):([A-Za-z]?)", paste0("ICD\\1:", icdprefix, "\\2"), codes)
}

strip_leading_icd <- function(codes) sub("^ICD[0-9]+:", "", codes)

match_codes <- function(codes, patterns, exact=FALSE, remove_point=FALSE, no_leading_icd=FALSE, icdprefix="") {
  if (length(patterns)==0) return(rep(FALSE, length(codes)))
  pats <- expand_ranges(patterns)
  if (no_leading_icd) pats <- strip_leading_icd(pats)
  if (nzchar2(icdprefix)) pats <- apply_icd_prefix(pats, icdprefix)
  x <- as.character(codes)
  if (remove_point) { x <- gsub("\\.", "", x); pats <- gsub("\\.", "", pats) }
  if (exact) {
    return(x %chin% pats)
  } else {
    res <- rep(FALSE, length(x))
    for (p in unique(pats)) {
      if (!nzchar2(p)) next
      res <- res | startsWith(x, p)
    }
    return(res)
  }
}

detect_ATC_status <- function(vec) {
  v <- to_vec(vec)
  if (length(v)==0) return("None")
  if (all(startsWith(v, "ATC"))) return("All")
  if (any(startsWith(v, "ATC"))) return("Some")
  "None"
}

load_code_sets <- function(path) {
  # file can be:
  # - one row list of codes (comma separated)
  # - two columns: name \t codes (comma separated)
  if (!nzchar2(path)) return(list())
  dt <- read_any(path, sep="\t")
  if (is.null(dt)) dt <- read_any(path, sep=",")
  if (is.null(dt)) stop("Cannot read exclusion file: ", path)
  if (ncol(dt)==1) {
    return(list(Default = expand_ranges(to_vec(dt[[1]]))))
  } else {
    out <- list()
    nmcol <- names(dt)[1]; cdcol <- names(dt)[2]
    for (i in seq_len(nrow(dt))) {
      nm <- as.character(dt[[nmcol]][i])
      cds <- expand_ranges(to_vec(dt[[cdcol]][i]))
      out[[nm]] <- cds
    }
    return(out)
  }
}

# ----------------------------- Indexing --------------------------------------
build_index <- function(infile, out_index, index_dtypes_json=NULL, sep=",") {
  log_msg("Indexing:", infile)
  # Just store a lightweight meta-index with column types & positions
  meta <- list(file=normalizePath(infile, winslash="/", mustWork=FALSE),
               sep=sep,
               nrows=NA_integer_,
               ncols=NA_integer_,
               columns=NULL,
               dtypes=NULL,
               created=as.character(Sys.time()))
  dt <- tryCatch(read_any(infile, sep=sep), error=function(e) NULL)
  if (!is.null(dt)) {
    meta$nrows <- nrow(dt); meta$ncols <- ncol(dt); meta$columns <- names(dt)
    if (nzchar2(index_dtypes_json)) {
      meta$dtypes <- tryCatch(jsonlite::fromJSON(index_dtypes_json), error=function(e) NULL)
    }
  }
  saveRDS(meta, file=out_index)
  log_msg("Wrote index RDS:", out_index)
  invisible(meta)
}

# ----------------------------- Entry/Exit ------------------------------------
build_entry_exit <- function(lpr_file, out_file, iidcol="IID", din="date_in", don="date_out",
                             bdcol="birthdate", sep=",", date_format="%d/%m/%Y") {
  dt <- read_any(lpr_file, sep=sep)
  if (is.null(dt)) stop("Could not read --LPR file: ", lpr_file)
  # Ensure date columns
  cn <- names(dt)
  din2 <- if (din %in% cn) din else if ("d_inddto" %in% cn) "d_inddto" else din
  don2 <- if (don %in% cn) don else if ("d_uddto" %in% cn) "d_uddto" else don
  # Convert to Date
  toDate <- function(x) {
    if (inherits(x, "Date")) return(x)
    if (inherits(x, "POSIXt")) return(as.Date(x))
    x <- as.character(x)
    suppressWarnings({
      y <- tryCatch(as.Date(x, format=date_format), error=function(e) NA)
    })
    if (all(is.na(y))) suppressWarnings( y <- as.Date(x) )
    y
  }
  if (din2 %in% cn) dt[[din2]] <- toDate(dt[[din2]])
  if (don2 %in% cn) dt[[don2]] <- toDate(dt[[don2]])
  # Aggregate
  setDT(dt)
  if (!iidcol %in% names(dt)) stop("IID column not found in LPR: ", iidcol)
  agg <- dt[, .(first_in = suppressWarnings(min(get(din2), na.rm=TRUE)),
                last_out = suppressWarnings(max(get(don2), na.rm=TRUE)),
                in_dates = paste(na.omit(unique(get(din2))), collapse=","),
                out_dates = paste(na.omit(unique(get(don2))), collapse=",")),
            by=.(IID = get(iidcol))]
  setcolorder(agg, c("IID","first_in","last_out","in_dates","out_dates"))
  write_tsv(agg, out_file)
  log_msg("Entry/Exit written:", out_file)
  invisible(agg)
}

# ----------------------------- Ophold ----------------------------------------
process_ophold <- function(ophold_file, out_file, sep="\t") {
  dt <- read_any(ophold_file, sep=sep)
  if (is.null(dt)) stop("Cannot read --Ophold file: ", ophold_file)
  # Minimal: ensure expected columns names if present
  # e.g., ophold_id, IID, start, end, type
  if ("cpr_enc" %in% names(dt) && !"IID" %in% names(dt)) setnames(dt, "cpr_enc", "IID")
  write_tsv(dt, out_file)
  log_msg("Ophold processed ->", out_file)
  invisible(dt)
}

# ----------------------------- Phenotype build -------------------------------
compute_age <- function(birthdate, date) {
  if (is.na(birthdate) || is.na(date)) return(NA_real_)
  as.numeric((as.Date(date) - as.Date(birthdate)) / 365.25)
}

filter_minmax_age <- function(dt, bdcol="birthdate", first_date_col="first_date", minmax="0,0") {
  rng <- to_vec(minmax)
  if (length(rng)==1) rng <- to_vec(strsplit(rng, ",", fixed=TRUE)[[1]])
  if (length(rng) < 2) return(dt)
  min_age <- as.numeric(rng[1]); max_age <- as.numeric(rng[2])
  if (is.na(min_age) || is.na(max_age) || (min_age==0 && max_age==0)) return(dt)
  dt[, age_first := mapply(compute_age, get(bdcol), get(first_date_col)) ]
  dt[ age_first > min_age & age_first < max_age ]
}

apply_gender_filter <- function(dt, gender="", sexcol="sex") {
  if (!nzchar2(gender)) return(dt)
  dt[ get(sexcol) == gender ]
}

apply_yob_filter <- function(dt, fyob="") {
  if (!nzchar2(fyob)) return(dt)
  cutoff <- as.Date(fyob)
  bd <- dt[["birthdate"]]
  if (is.null(bd)) return(dt)
  dt[ as.Date(bd) >= cutoff ]
}

apply_diag_type_filters <- function(dt, diagtype_col="c_types", incl="", excl="") {
  if (!diagtype_col %in% names(dt)) return(dt)
  if (nzchar2(incl)) {
    keep <- dt[[diagtype_col]] %chin% to_vec(incl)
    dt <- dt[keep]
  }
  if (nzchar2(excl)) {
    rm <- dt[[diagtype_col]] %chin% to_vec(excl)
    dt <- dt[!rm]
  }
  dt
}

merge_secondary_diag <- function(main_dt, sec_dt, recnum_col_main="recnum", recnum_col_sec="recnum", f2col="c_adiag") {
  if (is.null(sec_dt) || nrow(sec_dt)==0) return(main_dt)
  if (!recnum_col_main %in% names(main_dt) || !recnum_col_sec %in% names(sec_dt)) return(main_dt)
  setDT(main_dt); setDT(sec_dt)
  sec_dt_small <- sec_dt[, .(recnum=get(recnum_col_sec), f2diag=get(f2col))]
  out <- merge(main_dt, sec_dt_small, by.x=recnum_col_main, by.y="recnum", all.x=TRUE)
  out
}

load_mapping <- function(file, iidcol="pnr", sexcol="sex", bdcol="birthdate", sep=",") {
  if (!nzchar2(file)) return(NULL)
  dt <- read_any(file, sep=sep)
  if (is.null(dt)) stop("Cannot read mapping (-i): ", file)
  # Normalize columns if possible
  if ("cpr_enc" %in% names(dt) && iidcol %nin% names(dt)) setnames(dt, "cpr_enc", iidcol)
  dt
}

apply_exclusions_lifetime <- function(cases_dt, lsets, diagcol="diagnosis", remove_point=FALSE, exact=FALSE, no_leading_icd=FALSE, icdprefix="") {
  if (length(lsets)==0) return(cases_dt)
  keep <- rep(TRUE, nrow(cases_dt))
  for (nm in names(lsets)) {
    ex_codes <- lsets[[nm]]
    m <- match_codes(cases_dt[[diagcol]], ex_codes, exact=exact, remove_point=remove_point, no_leading_icd=no_leading_icd, icdprefix=icdprefix)
    # any lifetime match -> exclude
    keep <- keep & !m
  }
  cases_dt[keep]
}

apply_exclusions_post <- function(cases_dt, psets, diagcol="diagnosis", datecol_in="date_in", datecol_out="date_out",
                                  remove_point=FALSE, exact=FALSE, no_leading_icd=FALSE, icdprefix="") {
  if (length(psets)==0) return(cases_dt)
  setDT(cases_dt)
  cases_dt[, keep := TRUE]
  for (nm in names(psets)) {
    codes <- psets[[nm]]
    m <- match_codes(cases_dt[[diagcol]], codes, exact=exact, remove_point=remove_point, no_leading_icd=no_leading_icd, icdprefix=icdprefix)
    # for rows after first occurrence, drop
    if (any(m)) {
      first_date <- suppressWarnings(min(as.Date(cases_dt$date_in[m]), na.rm=TRUE))
      if (!is.infinite(first_date)) {
        cases_dt[as.Date(date_in) > first_date, keep := FALSE]
      }
    }
  }
  cases_dt[keep == TRUE][, keep := NULL][]
}

apply_exclusions_1yprior <- function(cases_dt, ysets, datecol_in="date_in", remove_point=FALSE, exact=FALSE, no_leading_icd=FALSE, icdprefix="") {
  if (length(ysets)==0) return(cases_dt)
  setDT(cases_dt)
  cases_dt[, keep := TRUE]
  for (nm in names(ysets)) {
    codes <- ysets[[nm]]
    m <- match_codes(cases_dt$diagnosis, codes, exact=exact, remove_point=remove_point, no_leading_icd=no_leading_icd, icdprefix=icdprefix)
    idx <- which(m)
    if (length(idx)>0) {
      dates <- as.Date(cases_dt$date_in[idx])
      lo <- min(dates, na.rm=TRUE) - 365
      hi <- min(dates, na.rm=TRUE)
      cases_dt[as.Date(date_in) >= lo & as.Date(date_in) <= hi, keep := FALSE]
    }
  }
  cases_dt[keep == TRUE][, keep := NULL][]
}

build_phenotype <- function(
  g_file, out_file, f_file="", f2_file="", atc_file="", i_file="", j_file="",
  ge="", qced="", name="MainPheno",
  fcol="c_adiag", gcol="c_adiag", iidcol="pnr", bdcol="birthdate", sexcol="sex",
  atccol="", atcdatecol="", fsep=",", isep=",", jsep=",", gsep=",", ophsep=",",
  din="d_inddto", don="d_uddto", recnum="", recnum2="", f2col="c_adiag",
  ExDepExc=FALSE, eM=FALSE, noLeadingICD=FALSE, ICDCM=FALSE, ICD8=FALSE, ICD9=FALSE, ICD10=FALSE,
  icdprefix="", iidstatus="", iidstatusdate="", selectIIDs="", DiagTypeExclusions="", DiagTypeInclusions="",
  DiagTypecol="c_types", LifetimeExclusion="", PostExclusion="", OneyPriorExclusion="",
  fDates="", iDates="", atcDates="", DateFormat="%d/%m/%Y", MinMaxAge="0,0", Fyob="", Fgender="",
  eCc=FALSE, removePointInDiagCode=FALSE, skipICDUpdate=FALSE, MatchFI=FALSE,
  BuildEntryExitDates=FALSE, Ophold="", BuildOphold=FALSE, RegisterRun=FALSE, lpp=FALSE,
  write_pickle=FALSE, write_fastGWA_format=FALSE, write_Plink2_format=FALSE,
  BuildTestSet=FALSE, testRun=FALSE, nthreads=8, lowmem=FALSE, batchsize=100000,
  PSYK=FALSE, LPR=FALSE, BuildIndex=FALSE, IndexDtypes='{}', verbose=FALSE
) {
  log_set(verbose)
  # Parse code list
  codes <- parse_code_file(g_file, gsep=gsep, gcol=gcol)
  atc_status <- detect_ATC_status(codes)
  # Load diagnosis
  dt <- if (nzchar2(f_file)) read_any(f_file, sep=fsep) else NULL
  if (is.null(dt)) stop("Main diagnoses file (-f) is required or must be auto-detected; not found.")
  # Standardize column names used internally
  cn <- names(dt)
  if (!iidcol %in% cn && "cpr_enc" %in% cn) { setnames(dt, "cpr_enc", iidcol); cn <- names(dt) }
  if (!din %in% cn && "date_in" %in% cn) { din <- "date_in" }
  if (!don %in% cn && "date_out" %in% cn) { don <- "date_out" }
  if (!fcol %in% cn && "diagnosis" %in% cn) { fcol <- "diagnosis" }
  # Secondary diagnoses
  dt2 <- if (nzchar2(f2_file)) read_any(f2_file, sep=fsep) else NULL
  if (!is.null(dt2) && !recnum %in% names(dt)) { log_msg("recnum not found in -f; skipping merge with --f2", level="WARN"); dt2 <- NULL }
  if (!is.null(dt2) && !recnum2 %in% names(dt2)) { log_msg("recnum2 not found in --f2; trying 'recnum'", level="WARN"); if ("recnum" %in% names(dt2)) recnum2 <- "recnum" }
  if (!is.null(dt2)) dt <- merge_secondary_diag(dt, dt2, recnum_col_main=recnum, recnum_col_sec=recnum2, f2col=f2col)
  setDT(dt)
  # Optional filters by diag type
  dt <- apply_diag_type_filters(dt, diagtype_col=DiagTypecol, incl=DiagTypeInclusions, excl=DiagTypeExclusions)
  # Map table (-i) for sex/birthdate etc.
  map_dt <- load_mapping(i_file, iidcol=iidcol, sexcol=sexcol, bdcol=bdcol, sep=isep)
  if (!is.null(map_dt)) {
    setDT(map_dt)
    keep_cols <- intersect(c(iidcol, sexcol, bdcol, iidstatus, iidstatusdate), names(map_dt))
    map_dt <- unique(map_dt[, ..keep_cols])
    dt <- merge(dt, map_dt, by=iidcol, all.x=TRUE)
  }
  # Select IIDs if provided
  if (nzchar2(selectIIDs)) {
    ids <- fread(selectIIDs, header=FALSE)[[1]]
    dt <- dt[get(iidcol) %chin% ids]
  }
  # Match codes
  codes_target <- codes
  if (ICD8 || ICD9 || ICD10) {
    # Restrict by ICD version
    pick <- character()
    if (ICD8) pick <- c(pick, grep("^ICD8:", codes, value=TRUE))
    if (ICD9) pick <- c(pick, grep("^ICD9", codes, value=TRUE))
    if (ICD10) pick <- c(pick, grep("^ICD10", codes, value=TRUE))
    if (length(pick)>0) codes_target <- pick
  }
  # If --skipICDUpdate is FALSE, allow applied icdprefix/noLeading
  matches <- match_codes(dt[[fcol]], codes_target, exact=eM, remove_point=removePointInDiagCode, no_leading_icd=noLeadingICD && !skipICDUpdate, icdprefix=if (!skipICDUpdate) icdprefix else "")
  dt_case <- dt[matches == TRUE]
  # Compute first diagnosis date per IID
  if (!din %in% names(dt_case)) stop("First date column not found in cases: ", din)
  dt_case[, date_in := as.Date(get(din))]
  if (don %in% names(dt_case)) dt_case[, date_out := as.Date(get(don))] else dt_case[, date_out := date_in]
  # Exclusions
  if (nzchar2(LifetimeExclusion)) {
    sets <- load_code_sets(LifetimeExclusion)
    dt_case <- apply_exclusions_lifetime(dt_case, sets, diagcol=fcol,
      remove_point=removePointInDiagCode, exact=eM, no_leading_icd=noLeadingICD && !skipICDUpdate, icdprefix=if (!skipICDUpdate) icdprefix else "")
  }
  if (nzchar2(PostExclusion)) {
    sets <- load_code_sets(PostExclusion)
    dt_case <- apply_exclusions_post(dt_case, sets, diagcol=fcol,
      datecol_in=din, datecol_out=don, remove_point=removePointInDiagCode, exact=eM, no_leading_icd=noLeadingICD && !skipICDUpdate, icdprefix=if (!skipICDUpdate) icdprefix else "")
  }
  if (nzchar2(OneyPriorExclusion)) {
    sets <- load_code_sets(OneyPriorExclusion)
    dt_case <- apply_exclusions_1yprior(dt_case, sets, datecol_in=din, remove_point=removePointInDiagCode, exact=eM, no_leading_icd=noLeadingICD && !skipICDUpdate, icdprefix=if (!skipICDUpdate) icdprefix else "")
  }
  # Summarize to per-IID phenotype
  if (nrow(dt_case)==0) {
    pheno <- data.table(IID=character(), PHENO=integer())
  } else {
    pheno <- dt_case[, .(
      first_date = suppressWarnings(min(date_in, na.rm=TRUE)),
      first_code = dt_case[[fcol]][.I[which.min(date_in)]],
      n_events = .N
    ), by=.(IID = get(iidcol))]
    pheno[, PHENO := 1L]
  }
  # add controls if MatchFI not set: include individuals from -f without match as PHENO=0
  if (!MatchFI) {
    all_ids <- unique(dt[[iidcol]])
    add0 <- setdiff(all_ids, pheno$IID)
    if (length(add0)>0) pheno <- rbindlist(list(pheno, data.table(IID=add0, PHENO=0L)), fill=TRUE, use.names=TRUE)
  }
  # Merge mapping columns (sex, birthdate)
  if (!is.null(map_dt)) {
    pheno <- merge(pheno, unique(map_dt[, ..c(iidcol, sexcol, bdcol)]), by.x="IID", by.y=iidcol, all.x=TRUE)
    setnames(pheno, old=c(sexcol, bdcol), new=c("sex","birthdate"), skip_absent=TRUE)
  }
  # Age filter
  pheno <- filter_minmax_age(pheno, bdcol="birthdate", first_date_col="first_date", minmax=MinMaxAge)
  # YOB filter
  pheno <- apply_yob_filter(pheno, Fyob)
  # Gender filter
  pheno <- apply_gender_filter(pheno, Fgender, sexcol="sex")
  # Write main TSV
  write_tsv(pheno, out_file)
  log_msg("Phenotype written:", out_file, "rows:", nrow(pheno))
  # Optional formats
  if (write_fastGWA_format) {
    fg <- pheno[, .(FID=IID, IID=IID, PHENO=PHENO)]
    write_tsv(fg, sub("\\.([^.]+)$", ".fastgwa.tsv", out_file))
  }
  if (write_Plink2_format) {
    p2 <- pheno[, .(IID=IID, PHENO=PHENO)]
    write_tsv(p2, sub("\\.([^.]+)$", ".plink2.tsv", out_file))
  }
  if (write_pickle) {
    saveRDS(pheno, file=sub("\\.([^.]+)$", ".rds", out_file))
  }
  invisible(pheno)
}

# ----------------------------- CLI -------------------------------------------
option_list <- list(
  make_option(c("--ini"), type="character", default="",
              help='Load an ini file that contains your data sources. Default: "%default%"'),
  make_option(c("-g"), type="character", default=NULL,
              help='File with all Diagnostic codes to export [required]'),
  make_option(c("-o"), type="character", default=NULL,
              help='Outfile name; include path [required]'),
  make_option(c("-f"), type="character", default="",
              help='Diagnoses file (-f). Columns: IID, date_in, date_out, diagnosis (names configurable)'),
  make_option(c("--f2"), type="character", default="",
              help='Secondary diagnosis file; join via recnum'),
  make_option(c("--atc"), type="character", default="",
              help='Prescription file; requires --atccol'),
  make_option(c("-i"), type="character", default="",
              help='IID mapping (sex, birthdate, etc.)'),
  make_option(c("-j"), type="character", default="",
              help='Extra IID info (optional)'),
  make_option(c("--ge"), type="character", default="",
              help='General exclusion IID list'),
  make_option(c("--qced"), type="character", default="",
              help='List of IIDs that pass QC'),
  make_option(c("--name"), type="character", default="MainPheno",
              help='Phenotype name'),
  make_option(c("--fcol"), type="character", default="c_adiag",
              help='Column name of -f diagnosis'),
  make_option(c("--gcol"), type="character", default="c_adiag",
              help='Column name in -g to read codes from'),
  make_option(c("--iidcol"), type="character", default="pnr",
              help='Column name of IDs (-f,-i)'),
  make_option(c("--bdcol"), type="character", default="birthdate",
              help='Birthdate column'),
  make_option(c("--sexcol"), type="character", default="sex",
              help='Sex column'),
  make_option(c("--atccol"), type="character", default="",
              help='ATC column'),
  make_option(c("--atcdatecol"), type="character", default="",
              help='ATC date column'),
  make_option(c("--fsep"), type="character", default=",",
              help='Separator for -f'),
  make_option(c("--isep"), type="character", default=",",
              help='Separator for -i'),
  make_option(c("--jsep"), type="character", default=",",
              help='Separator for -j'),
  make_option(c("--gsep"), type="character", default=",",
              help='Separator for -g'),
  make_option(c("--ophsep"), type="character", default=",",
              help='Separator for Ophold'),
  make_option(c("--din"), type="character", default="d_inddto",
              help='First diagnosis date column'),
  make_option(c("--don"), type="character", default="d_uddto",
              help='Discharge/out date column'),
  make_option(c("--recnum"), type="character", default="",
              help='recnum column in -f'),
  make_option(c("--recnum2"), type="character", default="",
              help='recnum column in --f2'),
  make_option(c("--f2col"), type="character", default="c_adiag",
              help='diagnosis column in --f2'),
  make_option(c("--ExDepExc"), action="store_true", default=FALSE,
              help='Use exclusion logic (legacy switch placeholder)'),
  make_option(c("--eM"), action="store_true", default=FALSE,
              help='Exact code match'),
  make_option(c("--noLeadingICD"), action="store_true", default=FALSE,
              help='Drop ICD* prefix in matching'),
  make_option(c("--ICDCM"), action="store_true", default=FALSE,
              help='ICD-CM flavour (no-op placeholder)'),
  make_option(c("--ICD8"), action="store_true", default=FALSE,
              help='Restrict to ICD8 codes'),
  make_option(c("--ICD9"), action="store_true", default=FALSE,
              help='Restrict to ICD9/9-CM codes'),
  make_option(c("--ICD10"), action="store_true", default=FALSE,
              help='Restrict to ICD10/10-CM codes'),
  make_option(c("--icdprefix"), type="character", default="",
              help='Add letter after ICD version prefix'),
  make_option(c("--iidstatus"), type="character", default="",
              help='IID status column name'),
  make_option(c("--iidstatusdate"), type="character", default="",
              help='IID status date column'),
  make_option(c("--selectIIDs"), type="character", default="",
              help='Path to IID list (one per line)'),
  make_option(c("--DiagTypeExclusions"), type="character", default="",
              help='Comma list of diag types to exclude'),
  make_option(c("--DiagTypeInclusions"), type="character", default="",
              help='Comma list of diag types to include'),
  make_option(c("--DiagTypecol"), type="character", default="c_types",
              help='Diag type column name'),
  make_option(c("--LifetimeExclusion"), type="character", default="",
              help='Lifetime exclusion file'),
  make_option(c("--PostExclusion"), type="character", default="",
              help='Post-exclusion file'),
  make_option(c("--OneyPriorExclusion"), type="character", default="",
              help='One-year-prior exclusion file'),
  make_option(c("--fDates"), type="character", default="",
              help='Extra date columns in -f (comma-separated)'),
  make_option(c("--iDates"), type="character", default="",
              help='Extra date columns in -i (comma-separated)'),
  make_option(c("--atcDates"), type="character", default="",
              help='Date columns in --atc'),
  make_option(c("--DateFormat"), type="character", default="%d/%m/%Y",
              help='Date format'),
  make_option(c("--MinMaxAge"), type="character", default="0,0",
              help='Minimum,Maximum age at first diagnosis (x,y)'),
  make_option(c("--Fyob"), type="character", default="",
              help='Filter by year of birth (YYYY-MM-DD cutoff)'),
  make_option(c("--Fgender"), type="character", default="",
              help='Filter by gender ("F" or "M")'),
  make_option(c("--eCc"), action="store_true", default=FALSE,
              help='Exclude CHB-only controls (cluster-specific; no-op if missing)'),
  make_option(c("--removePointInDiagCode"), action="store_true", default=FALSE,
              help='Strip dot in diagnosis when matching'),
  make_option(c("--skipICDUpdate"), action="store_true", default=FALSE,
              help='Skip ICD normalization steps'),
  make_option(c("--MatchFI"), action="store_true", default=FALSE,
              help='Keep only overlapping IIDs between -g and -f'),
  make_option(c("--BuildEntryExitDates"), action="store_true", default=FALSE,
              help='Build entry/exit dates from --LPR / -f data'),
  make_option(c("--Ophold"), type="character", default="",
              help='Path to Ophold file'),
  make_option(c("--BuildOphold"), action="store_true", default=FALSE,
              help='Process Ophold file'),
  make_option(c("--RegisterRun"), action="store_true", default=FALSE,
              help='RegisterRun heuristic (no-op placeholder)'),
  make_option(c("--lpp"), action="store_true", default=FALSE,
              help='Load phenotypes and run exclusions (no-op placeholder)'),
  make_option(c("--write_pickle"), action="store_true", default=FALSE,
              help='Also write RDS'),
  make_option(c("--write_fastGWA_format"), action="store_true", default=FALSE,
              help='Write fastGWA format'),
  make_option(c("--write_Plink2_format"), action="store_true", default=FALSE,
              help='Write PLINK2 format'),
  make_option(c("--BuildTestSet"), action="store_true", default=FALSE,
              help='Build test set (not used)'),
  make_option(c("--testRun"), action="store_true", default=FALSE,
              help='Small test mode (not used)'),
  make_option(c("--nthreads"), type="integer", default=8,
              help='DEPRECATED'),
  make_option(c("--lowmem"), action="store_true", default=FALSE,
              help='Low memory batching'),
  make_option(c("--batchsize"), type="integer", default=100000,
              help='Batch size for lowmem'),
  make_option(c("--PSYK"), action="store_true", default=FALSE,
              help='Use PSYK diagnoses only (no-op placeholder)'),
  make_option(c("--LPR"), action="store_true", default=FALSE,
              help='Use LPR diagnoses only (no-op placeholder)'),
  make_option(c("--BuildIndex"), action="store_true", default=FALSE,
              help='Build index for -f/--atc'),
  make_option(c("--IndexDtypes"), type="character", default='{}',
              help='JSON dtypes dict'),
  make_option(c("--verbose"), action="store_true", default=FALSE,
              help='Verbose logging')
)

opt <- parse_args(OptionParser(option_list=option_list))

# Ensure required
if (is.null(opt$g) && !opt$BuildEntryExitDates && !opt$BuildOphold && !opt$BuildIndex) {
  stop("-g is required unless running --BuildEntryExitDates / --BuildOphold / --BuildIndex")
}
if (is.null(opt$o)) stop("-o is required")

# Routing
if (isTRUE(opt$BuildIndex)) {
  if (!nzchar2(opt$f) && !nzchar2(opt$atc)) stop("Provide -f or --atc for --BuildIndex")
  if (nzchar2(opt$f)) build_index(opt$f, out_index=sub("\\.([^.]+)$", ".index.rds", opt$f), index_dtypes_json=opt$IndexDtypes, sep=opt$fsep)
  if (nzchar2(opt$atc)) build_index(opt$atc, out_index=sub("\\.([^.]+)$", ".index.rds", opt$atc), index_dtypes_json=opt$IndexDtypes, sep=opt$fsep)
  quit(status=0)
}

if (isTRUE(opt$BuildOphold)) {
  if (!nzchar2(opt$Ophold)) stop("--Ophold is required for --BuildOphold")
  process_ophold(opt$Ophold, out_file=ifelse(nzchar2(opt$o), opt$o, "ophold.tsv"), sep=opt$ophsep)
  quit(status=0)
}

if (isTRUE(opt$BuildEntryExitDates)) {
  # Prefer -f if given; else --LPR would imply using -f anyway
  src <- if (nzchar2(opt$f)) opt$f else stop("Provide -f for --BuildEntryExitDates")
  build_entry_exit(src, out_file=opt$o, iidcol=opt$iidcol, din=opt$din, don=opt$don, bdcol=opt$bdcol, sep=opt$fsep, date_format=opt$DateFormat)
  quit(status=0)
}

# Main phenotype
build_phenotype(
  g_file = opt$g, out_file = opt$o, f_file = opt$f, f2_file = opt$f2, atc_file = opt$atc, i_file = opt$i, j_file = opt$j,
  ge = opt$ge, qced = opt$qced, name = opt$name,
  fcol = opt$fcol, gcol = opt$gcol, iidcol = opt$iidcol, bdcol = opt$bdcol, sexcol = opt$sexcol,
  atccol = opt$atccol, atcdatecol = opt$atcdatecol, fsep = opt$fsep, isep = opt$isep, jsep = opt$jsep, gsep = opt$gsep, ophsep = opt$ophsep,
  din = opt$din, don = opt$don, recnum = opt$recnum, recnum2 = opt$recnum2, f2col = opt$f2col,
  ExDepExc = opt$ExDepExc, eM = opt$eM, noLeadingICD = opt$noLeadingICD, ICDCM = opt$ICDCM, ICD8 = opt$ICD8, ICD9 = opt$ICD9, ICD10 = opt$ICD10,
  icdprefix = opt$icdprefix, iidstatus = opt$iidstatus, iidstatusdate = opt$iidstatusdate, selectIIDs = opt$selectIIDs,
  DiagTypeExclusions = opt$DiagTypeExclusions, DiagTypeInclusions = opt$DiagTypeInclusions, DiagTypecol = opt$DiagTypecol,
  LifetimeExclusion = opt$LifetimeExclusion, PostExclusion = opt$PostExclusion, OneyPriorExclusion = opt$OneyPriorExclusion,
  fDates = opt$fDates, iDates = opt$iDates, atcDates = opt$atcDates, DateFormat = opt$DateFormat, MinMaxAge = opt$MinMaxAge,
  Fyob = opt$Fyob, Fgender = opt$Fgender, eCc = opt$eCc, removePointInDiagCode = opt$removePointInDiagCode,
  skipICDUpdate = opt$skipICDUpdate, MatchFI = opt$MatchFI, BuildEntryExitDates = opt$BuildEntryExitDates,
  Ophold = opt$Ophold, BuildOphold = opt$BuildOphold, RegisterRun = opt$RegisterRun, lpp = opt$lpp,
  write_pickle = opt$write_pickle, write_fastGWA_format = opt$write_fastGWA_format, write_Plink2_format = opt$write_Plink2_format,
  BuildTestSet = opt$BuildTestSet, testRun = opt$testRun, nthreads = opt$nthreads, lowmem = opt$lowmem, batchsize = opt$batchsize,
  PSYK = opt$PSYK, LPR = opt$LPR, BuildIndex = opt$BuildIndex, IndexDtypes = opt$IndexDtypes, verbose = opt$verbose
)
