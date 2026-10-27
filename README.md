[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21902133.svg)](https://doi.org/10.5281/zenodo.21902133)

# GetPheno
[Build Phenotypes](#1-build-phenotypes)

An optional low-dependency R Shiny interface can submit GetPheno runs to Slurm and retrieve results using a unique run code. It supports permission-aware project/protocol selection, newest-first cohort freezes, protocol-defined additional files, and clearly separated optional run settings. See [shiny_app/README.md](shiny_app/README.md).

A Docker image containing GetPheno, Python, R, the Shiny app, and Slurm client tools is provided in `Dockerfile`. Pushes to `main` publish `ghcr.io/mischalundberg/getpheno:edge`; release tags such as `V1.30` publish both the matching version tag and `latest`. Production deployments should pin a release tag. The per-user cluster deployment is documented under [Docker deployment](shiny_app/README.md#docker-deployment).

## Current capabilities

- Builds one or multiple phenotypes in one run from ICD8, ICD9, ICD10, ICD9-CM, ICD10-CM, and ATC code requests.
- Supports tab-separated multi-phenotype request files with one phenotype name per row and comma-separated codes per phenotype.
- Supports advanced phenotype definitions with `main=`, `sub=`, `rule_out=`, code ranges, and explicit wildcard codes.
- Diagnostic requests may always include trailing `*`. By default, all requested codes are matched as prefixes; with `--eM`, only codes ending in `*` remain prefix matches and all other codes are exact. You do not need to include the `--eM` flag for `*` to pick up prefix matches.
- Supports lifetime, one-year-prior, and post-onset exclusion files, including ATC-based exclusions.
- Supports IID subsets and optional death-date, age-at-death, living-person, and deceased-only filters.
- Supports `--lowmem` and `--verylowmem` batching with fixed output columns across batches, so appended TSV output keeps a stable schema.
- Supports `--iPSYCHComparison` with user-defined phenotypes or eight built-in presets and phenotype-specific eligibility rules.
- Restricts each batch to the current STAM/IID universe and applies QC (`--qced`) and general exclusions before case/control output is finalized.
- Writes per-phenotype status columns as `Case` or `Control`; phenotype detail columns remain blank when not applicable.
- Preserves STAM columns first in the output, followed by requested phenotypes, exclusions, covariates, and derived fields.
- Writes SHA256 checksums to `<output-stem>.checksums.sha256` and records that exact filename in the generated run README.
- Optionally writes OMOP CDM `COHORT` and `COHORT_DEFINITION` result tables with mapped `person_id` values and provenance.

# Getting Started (common)
In order to download `GetPheno`, you should clone this repository via the commands
```  
git clone https://github.com/MischaLundberg/GetPheno
cd GetPheno
```

In order to install the Python dependencies, you will need the [Anaconda](https://www.anaconda.com/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with the dependencies:

```
conda env create --file environment.yml
conda activate GetPheno

## to deactivate the environment, type
#conda deactivate
```

In case you are updating your current version of the pipeline, it is best practice to also update your environment to the updated prerequisites.
If you are using an Anaconda environment, you can do so by typing
```
conda env update --name GetPheno --file environment.yml
```


If you receive any errors while running scripts of the GetPheno pipeline, please ensure your versioning for the prerequisites is according to the tested versions.


## 1 Build Phenotypes:


### Input:
As input you need a file that contains a header, followed by your diagnostic codes of interest (one per row).
```
diagnosis
F33.0
F32.0
```
Alternatively, you can specify multiple phenotypes in one run. Put the phenotype name in the first column and comma-separated codes in the second column. The columns must be tab-separated. Codes can mix diagnostic classifications, for example ICD10 and ICD8.
```
MDD	ICD10:F33.0,ICD10:F32.0
SCZ	ICD10:F20,ICD8:295.3
BPD	ICD10:F30
```

Per-phenotype output columns are named after the requested phenotype, for example `MDD`, `SCZ`, and `BPD`, and contain `Case` or `Control`.

Advanced phenotype definitions can be used in the second column when a phenotype needs rule logic instead of a plain code list. Separate rule clauses with `;` and codes within a clause with `,`.
```
MD	main=ICD10:F32;rule_out=ICD10:F25
Recurrent_MD	main=ICD10:F33;sub=ICD10:F10,ICD10:F11;rule_out=ICD10:F25
Poisoning_Attempt	main=ICD10:F;sub=ICD10:T36-ICD10:T50,ICD10:T52-ICD10:T60
```

`main=` defines the main diagnosis codes for candidate cases. `sub=` defines supporting/secondary diagnosis codes. `rule_out=` removes IIDs with matching codes from that advanced phenotype case set. Code ranges such as `ICD10:T36-ICD10:T50` are expanded automatically before matching. A trailing `*`, for example `ICD10:F32*`, keeps that code prefix-based whether or not `--eM` is used.

Don't forget to specify your header `--gcol` if it differs from the default. It is otherwise set automatically to `diagnosis` or `c_adiag` depending on the cluster.

A Simple run would look like

```
./get_pheno.py -g diagnosis.request -o diagnosis.casecontrol
```
Keep in mind that there may be crucial flags to set, for example `--eM`.
The script will (if not told differently) update the ICD codes from i.e. F32.0 to ICD10:DF320 automatically on CHB/DBDS and IBP clusters.

Without `--eM`, all code requests are prefix matches. Operationally, this is equivalent to adding a trailing `*` to every requested code after ICD/ATC normalization, even if the user did not write `*`. `F25*,F33*,F32.2` is therefore interpreted as `F25*`, `F33*`, and `F32.2*`: it matches codes starting with `F25`, codes starting with `F33`, and codes starting with `F32.2`. For example, `ICD8:296.2` also matches subcodes such as `ICD8:296.29` after normalization.

With `--eM`, codes without trailing `*` are exact matches, while codes with trailing `*` stay prefix matches. `F25*,F33*,F32.2` is interpreted as `F25*`, `F33*`, and exact `F32.2`: it matches codes starting with `F25`, codes starting with `F33`, and only the exact code `F32.2`.

Ranges are expanded first and then matched using the same rules. For example, `ICD10:T36-ICD10:T38` becomes `ICD10:T36,ICD10:T37,ICD10:T38`. Without `--eM`, these expanded codes are prefix matches, so `ICD10:T360` would also match `ICD10:T36`. With `--eM`, the expanded codes are exact unless you explicitly write wildcard requests such as `ICD10:T36*,ICD10:T37*,ICD10:T38*`.

We also support a `get_pheno.ini` file. This allows running with the same base settings every time. `get_pheno.py` will use the `get_pheno.ini` in (1) the current working directory or (2) the directory where `get_pheno.py` is located. Manually set flags supersede settings in `get_pheno.ini`.

### Multiple diagnosis/event files

`-f` accepts comma-separated files. Column options accept either one name, which is expected in every file, or exactly one name per file. Per-file columns are renamed to the first name before the rows are combined.

```bash
./get_pheno.py \
  -f lpr3_diagnoses.tsv,lpr3_procedures.tsv \
  --fcol diagnosis,code \
  --iidcol cpr_enc \
  --din date_in --don date_out \
  --DiagTypecol diagnosis_type, \
  --fsep $'\t' \
  -i cohort.tsv -g phenotypes.tsv -o output.tsv
```

Here `diagnosis` is canonical and `code` from the second file is renamed to `diagnosis`. A single value such as `--din date_in` is checked in both files. If it is absent from one file, GetPheno warns once and fills canonical `date_in` with `NA` for rows from that file. Use per-file names such as `--din date_in,procedure_date` when headers differ. An explicit empty entry, as in `--DiagTypecol diagnosis_type,`, declares that the second file has no corresponding column and fills it with `NA` without a warning.

This rule applies to `--fcol`, `--iidcol`, `--din`, `--don`, `--DiagTypecol`, and `--recnum`. Supplying neither one name nor exactly one name per file is an error. The same settings can be placed in the INI file as comma-separated values.

### Multiple additional-information files

`-j` also accepts comma-separated files, for example:

```bash
-j t_person.tsv,t_dodsaarsag_2.tsv --iidcol cpr_enc --jsep $'\t'
```

These files are merged horizontally by IID rather than appended as rows. All IID columns are read as strings before pandas type inference, preserving encrypted or numeric-looking identifiers exactly. The first file determines row and column priority: new columns from later files are appended, while later values only fill missing values in overlapping columns. Duplicate IIDs within an individual metadata file retain their first row. GetPheno logs the number of unique and overlapping IIDs for each additional file and warns if a later file has no overlap with preceding `-j` files.

### iPSYCH-comparison cohort mode

`--iPSYCHComparison` creates phenotypes directly from the raw diagnosis file (`-f`). Supply a `-g` phenotype-request file to use your own definitions, or omit `-g` to use all eight built-in presets. For example:

```text
skiz	ICD10:F20
skizospek	ICD10:F20,ICD10:F21,ICD10:F22,ICD10:F23,ICD10:F24,ICD10:F25,ICD10:F28,ICD10:F29
bipol	ICD10:F30,ICD10:F31
affek	ICD10:F30,ICD10:F31,ICD10:F32,ICD10:F33,ICD10:F34,ICD10:F35,ICD10:F36,ICD10:F37,ICD10:F38,ICD10:F39
autism	ICD10:F84.0,ICD10:F84.1,ICD10:F84.5,ICD10:F84.8,ICD10:F84.9
adhd	ICD10:F90.0
anorek	ICD10:F50.0,ICD10:F50.1
```

Example:

```bash
./get_pheno.py \
  -f diagnoses.tsv \
  -i stam.tsv \
  -j additional_information.tsv \
  -g ipsych_phenotypes.tsv \
  -o ipsych_comparison.tsv \
  --fsep $'\t' --isep $'\t' --jsep $'\t' --gsep $'\t' \
  --iidcol pnr --bdcol birthdate \
  --deathdatecol D_DODSDATO \
  --fcol diagnosis --gcol diagnosis \
  --din date_in --don date_out \
  --iidstatus status_code --iidstatusdate status_date \
  --iPSYCHComparison
```

The default cohort rules mirror the supplied iPSYCH-comparison script:

- Birthdate is at least `1981-05-01` and earlier than `2009-01-01`. These bounds are applied to the STAM/IID universe first, so IIDs outside the birth window are never reported as cases or controls.
- Danish birthplace is determined from `fkode` using the code ranges from the comparison script.
- A nonmissing maternal CPR/identifier is required in `pnr_m`.
- The individual must be alive or resident in Denmark at age one: `statd >= birthdate + 365 days`, or `stat` is `1` or `3`. When `--deathdatecol` is supplied, a recorded death before the first birthday excludes the IID even if a status field would otherwise retain it.
- Repeated mother/birthdate combinations are treated as multiple births and all are excluded. Use `--iPSYCHKeepMultiples` to retain them.
- Parental birthplace is reported but is not an inclusion criterion, matching the supplied R code. The mode auto-detects `DKDK` or `both_parents_DK`, or derives the indicator from `fkode_m` and `fkode_f`.
- With a supplied phenotype file, diagnosis episodes must overlap the configured date window. By default there is no lower diagnosis boundary and the exclusive upper boundary is `2016-01-01`. Both `--din` and `--don` are parsed and used; an episode beginning before the upper boundary and ending after it is retained.

Without `-g`, the supplied `iPSYCH_*_Codes` definitions produce these columns:

| Preset | Register contains (case insensitive) | Minimum age at diagnosis | Last diagnosis date (inclusive) |
| --- | --- | --- | --- |
| AFFEK, BIPOL, SKIZOSPEK, SKIZO | PSYK | 10 years | 2015-12-31 |
| ADHD, ASD | PSYK | 1 year | 2015-12-31 |
| ANO | PSYK or LPR | 1 year | 2016-12-31 |
| PPD | PSYK | No additional age threshold | 2015-12-31 |

PPD additionally requires a supplied live-birth code on the mother's own IID, dated from the diagnosis date back through one calendar year. Birth day and first anniversary are included. All available birth events are considered, including LPR events when `--PSYK` selects psychiatric diagnoses. Presets match diagnosis-code prefixes, including subcodes, regardless of `--eM`. `--iPSYCHDxStart` optionally adds a lower diagnosis-date boundary; preset-specific upper dates take precedence over `--iPSYCHDxEnd`. Custom `-g` definitions retain their existing matching and date-window behavior.

`--registercol origin` selects the register/source column; `registercol = origin` configures it in the INI file. A comma-separated mapping supports different names across `-f` files. When filtering is needed and no name is supplied, each file is checked for `register`, then `source`. Missing required columns cause an error. `--PSYK` and `--LPR` filter rows by substring in this column across every file; file order no longer selects the source. `IndexDtypes` controls storage types only; configured register columns are stored as strings when indexing.

The cohort metadata used by these rules may be supplied in either the STAM file (`-i`) or the additional-information file (`-j`). GetPheno combines both tables by IID before applying the restrictions, then limits both tables to the retained IID universe. Birthdate remains required in `-i`. A death-date column is optional: datasets containing only `stat` continue to use the status-only rule.

`--deathdatecol D_DODSDATO` identifies the optional death date; blank values mean no recorded death. `--MinMaxAgeDeath 20,80` keeps deceased people whose completed age at death is within the inclusive range. Add `--keepAlive` to retain people with blank death dates alongside that range, or use `--onlyDead` to retain only recorded deaths. `--onlyDead` can be combined with the age range; `--keepAlive` and `--onlyDead` are mutually exclusive. `--selectIIDs` accepts a headerless file containing one IID per line.

If no DKDK/both-Danish-born-parents indicator and no `fkode_m`/`fkode_f` columns are available, GetPheno prints and logs a warning because the parental-birthplace indicator cannot be reported; no IID is excluded for this reason. Missing child birthplace, maternal identifier, or age-one status fields also produce explicit warnings when their corresponding criterion cannot be applied. Supplied phenotype files require both diagnosis date columns for episode overlap; presets use `--din` as the diagnosis date.

Cohort and custom-definition date defaults can be changed with `--iPSYCHBirthStart`, `--iPSYCHBirthEnd`, `--iPSYCHDxStart`, `--iPSYCHDxEnd`, `--iPSYCHDKParentsCol`, `--iPSYCHBirthplaceCol`, and `--iPSYCHMotherCol`. Use the existing `--iidstatus` and `--iidstatusdate` options when the age-one fields are not named `stat` and `statd`; these names may refer to columns in `-i` or `-j`. The same restricted IID universe and diagnosis window are used in standard, `--lowmem`, and `--verylowmem` execution.

### Options (Flags):
```
usage: get_pheno.py [-h] [--ini INI] [-g G] -o O [-f F] [--f2 F2] [--atc ATC] [-i I] [-j J] [--ge GE] [--qced QCED] [--name NAME] [--fcol FCOL] [--gcol GCOL] [--iidcol IIDCOL] [--bdcol BDCOL] [--sexcol SEXCOL] [--atccol ATCCOL] [--atcdatecol ATCDATECOL] [--fsep FSEP] [--gsep GSEP] [--ophsep OPHSEP] [--din DIN] [--don DON] [--recnum RECNUM] [--recnum2 RECNUM2] [--f2col F2COL] [--ExDepExc] [--eM]
                                [--noLeadingICD] [--ICDCM] [--ICD8] [--ICD9] [--ICD10] [--iidstatus IIDSTATUS] [--iidstatusdate IIDSTATUSDATE] [--DiagTypeExclusions DIAGTYPEEXCLUSIONS] [--DiagTypeInclusions DIAGTYPEINCLUSIONS] [--LifetimeExclusion LIFETIMEEXCLUSION] [--PostExclusion POSTEXCLUSION] [--OneyPriorExclusion ONEYPRIOREXCLUSION] [--fDates FDATES] [--iDates IDATES] [--atcDates ATCDATES] [--DateFormat DATEFORMAT] [--MinMaxAge MINMAXAGE] [--Fyob FYOB]
                                [--iPSYCHComparison] [--iPSYCHBirthStart IPSYCHBIRTHSTART] [--iPSYCHBirthEnd IPSYCHBIRTHEND] [--iPSYCHDxStart IPSYCHDXSTART] [--iPSYCHDxEnd IPSYCHDXEND] [--iPSYCHDKParentsCol IPSYCHDKPARENTSCOL] [--iPSYCHBirthplaceCol IPSYCHBIRTHPLACECOL] [--iPSYCHMotherCol IPSYCHMOTHERCOL] [--iPSYCHKeepMultiples]
                                [--Fgender FGENDER] [--eCc] [--removePointInDiagCode] [--skipICDUpdate] [--MatchFI] [--BuildEntryExitDates] [--Ophold OPHOLD] [--BuildOphold] [--RegisterRun] [--lpp] [--write_pickle] [--write_fastGWA_format] [--write_Plink2_format] [--BuildTestSet] [--testRun] [--nthreads NTHREADS] [--lowmem] [--verylowmem] [--batchsize BATCHSIZE] [--registercol REGISTERCOL] [--PSYK] [--LPR] [--BuildIndex] [--IndexDtypes INDEXDTYPES] [--verbose]

Extracts a Phenotype from input files based on IIDs and diagnostic codes. The best way to start, is to generate a test dataset: 'python get_phenotype.py -g "" -o ./ --BuildTestSet' and then run 'python get_phenotype.py -g "" -o testrun.tsv --eM --ExDepExc --testRun'

optional arguments:
  -h, --help            show this help message and exit
  --ini INI             Load an ini file that contains your data sources. Default: ""
  -g G                  File with all Diagnostic codes to export
  -o O                  Outfile name; dont forget to add the location, otherwise it will be saved in the local dir.
  -f F                  Diagnosis file. Should have at least "IID", "date_in", "date_out", and "diagnosis" as columns. Names can be defined using --iidcol, --din, --don, and --fcol. No default. Will be automatically determined if not entered for GenomeDK and CHB/DBDS (DEGEN protocol)
  --f2 F2               Secondary diagnosis file. This is meant to be used if you have files with secondary diagnosis information as i.e. in Denmark. Here we would usually use the 'recnum' to merge these entries with '-f'. These files also contain less information as the '-f' files. Should have at least "recnum", and "diagnosis" as columns. Names can be defined using --recnum2 and --f2col. The recnum column in the
                        '-f' files can be specified using --recnum. Default: ""
  --atc ATC             PRescription information file. Will be handeld simillar to -f but needs additional information to be specified: --atccol. Default: ""
  -i I                  This file should adds based on a mapping using the supplied "IID" column information about Gender/Sex (needed), Brithdate (needed) ... e.g., information from stamdata file. No default. Will be automatically determined if not entered for GenomeDK and CHB/DBDS (DEGEN protocol)
  -j J                  This file should have additional "IID" information e.g. date of entry to cohort. This is not needed and can also be skipped. No default. Will be automatically determined if not entered for GenomeDK and CHB/DBDS (DEGEN protocol)
  --ge GE               General exclusion list. This referrs to a list of IIDs that should be excluded from the study. Within CHB/DBDS it will default to /data/preprocessed/genetics/chb_degen_freeze_20210503/degen_exclusion_latest
  --qced QCED           List with all IIDs that pass initial QC. Default: "" (or '/data/preprocessed/genetics/chb_degen_freeze_20210503/DEGEN_GSA_FINAL.fam' in CHB/DBDS)
  --name NAME           Define your main Phenotype name. Defaults to "MainPheno"
  --fcol FCOL           Columname of -f and -i file to be mapped. Defaulting to "c_adiag" (or diagnosis in CHB/DBDS)
  --gcol GCOL           Columname of -g file to be mapped against. Defaulting to "c_adiag" (or diagnosis in CHB/DBDS)
  --iidcol IIDCOL       Columname of IDs in -f and -i file. Defaulting to "pnr" (or cpr_enc in CHB/DBDS)
  --bdcol BDCOL         Columname of Birthdate in files. Defaults to "birthdate"
  --sexcol SEXCOL       Columname of Sex/Gender in files. Defaults to "sex"
  --atccol ATCCOL       Columname of ATC codes in files. Defaults to ""
  --atcdatecol ATCDATECOL
                        Columname of ATC prescription date column in the --atc file(s). Defaults to ""
  --fsep FSEP           Separator of -f i.e. tab; default "," (or "\t" in CHB/DBDS)
  --gsep GSEP           Separator of -g i.e. tab; default "," (or "\t" in CHB/DBDS)
  --ophsep OPHSEP       Separator of Ophold file - currently only available on CHB/DBDS,DST i.e. "\t"; default "," (or "\t" in CHB/DBDS)
  --din DIN             Columname of first diagnosis date. e.g. 'd_inddto' or 'date_in'. Default: "d_inddto" (or 'date_in' in CHB/DBDS)
  --don DON             Columname of first diagnosis date. e.g. 'd_uddto' or 'date_out'. Default: "d_uddto" (or 'date_out' in CHB/DBDS)
  --recnum RECNUM       Columname of the recnum field in -f files. Default: "" (or 'recnum' on NCRR)
  --recnum2 RECNUM2     Columname of the recnum field in -f files. Default: "" (or 'recnum' on NCRR)
  --f2col F2COL         Columname of first diagnosis date. e.g. 'd_uddto' or 'date_out'. Default: "c_adiag" (or 'date_out' in CHB/DBDS)
  --ExDepExc            List of diagnostic codes to exclude (e.g. Cases holding also these codes will be excluded). This is currently precoded and will be changed to allowing to specify a file with codes.
  --eM                  Exact Match. Codes without trailing * must match exactly. Codes with trailing * remain prefix matches.
  --noLeadingICD        Set this, if your diagnostic codes in "-g" do have a leading ICD*:, but your "-f" does not. This will only take affect together with --ExDepExc.
  --ICDCM               Set this, if your Cohort is based on ICD10/9-CM codes and not base WHO. This will only take affect together with --ExDepExc.
  --ICD8                Restrict phenotype export to ICD8 codes.
  --ICD9                Restrict phenotype export to ICD9/ICD9-CM codes.
  --ICD10               Restrict phenotype export to ICD10/ICD10-CM codes.
  --iidstatus IIDSTATUS
                        Status column, including the age-one status used by --iPSYCHComparison. The column may be in -i or -j. In iPSYCH mode, blank defaults to "stat".
  --iidstatusdate IIDSTATUSDATE
                        Date column for --iidstatus. The column may be in -i or -j. In iPSYCH mode, blank defaults to "statd".
  --DiagTypeExclusions DIAGTYPEEXCLUSIONS
                        List of diagnostic types to exclude. e.g. "H,M". Potential c_types may be +,A,B,C,G,M,H. A (main diagnosis), B (secondary diagnosis), G (grundmorbus), H (referral), + (secondary diagnosis), C (complication), and M (temporary) are all possible diagnosis types. From 1995 to 1998, registering a referral diagnosis (H) was optional; thereafter, it became mandatory for certain referral methods.
                        Grundmorbus (G) was optional from 1995 to 2001, used exclusively for psychiatric patients in 2002–2003, and then discontinued entirely. Complication (C) and temporary diagnosis (M) were both discontinued at the end of 2013. Default: ""
  --DiagTypeInclusions DIAGTYPEINCLUSIONS
                        List of diagnostic types to include. e.g. "A,B". Potential c_types may be +,A,B,C,G,M,H. A (main diagnosis), B (secondary diagnosis), G (grundmorbus), H (referral), + (secondary diagnosis), C (complication), and M (temporary) are all possible diagnosis types. From 1995 to 1998, registering a referral diagnosis (H) was optional; thereafter, it became mandatory for certain referral methods.
                        Grundmorbus (G) was optional from 1995 to 2001, used exclusively for psychiatric patients in 2002–2003, and then discontinued entirely. Complication (C) and temporary diagnosis (M) were both discontinued at the end of 2013. Default: ""
  --LifetimeExclusion LIFETIMEEXCLUSION
                        Define Lifetime exclusions (if a case has any of the listed codes, it will be excluded). This should be a file containing either one row with all codes listed (comma separated) or per row a exclusion name, followed by a list of diagnostic codes (similar to the input). e.g. "BPD
                        ICD10:F30,ICD10:F30.0,ICD10:F30.1,ICD10:F30.2,ICD10:F30.8,ICD10:F30.9,ICD10:F31,ICD10:F31.0,ICD10:F31.1,ICD10:F31.2,ICD10:F31.3,ICD10:F31.4,ICD10:F31.5,ICD10:F31.6,ICD10:F31.7,ICD10:F31.8,ICD10:F31.9" If you are using CHB/DBDS data, please remember, that the ICD10 codes have to start with ICD10:D***, e.g. ICD10:DF33.0. Default: ""
  --PostExclusion POSTEXCLUSION
                        Define Post exclusions (if a case has any of the listed codes, all main diagnoses after the first occuring date of the listed codes will be excluded). This should be a file containing either one row with all codes listed (comma separated) or per row a exclusion name, followed by a list of diagnostic codes (similar to the input). Default: ""
  --OneyPriorExclusion ONEYPRIOREXCLUSION
                        Define One Year Prio exclusions (all entries for a case that happen within one year prior to any date of the listed codes, these entries will be excluded). This should be a file containing either one row with all codes listed (comma separated) or per row a exclusion name, followed by a list of diagnostic codes (similar to the input). Default: ""
  --fDates FDATES       Which columns in your -f file are Dates? This automatically uses the columns you supplied under --bdcol, --don and --din. You only need to specify this, if you have additonal Dates that should be reformatted. Specify as comma separated list without spaces. Default: "[]"
  --iDates IDATES       Which columns in your -i file are Dates? This automatically uses the columns you supplied under --bdcol, --don and --din. You only need to specify this, if you have additonal Dates that should be reformatted. Specify as comma separated list without spaces. Default: "[]"
  --atcDates ATCDATES   Which columns in your --atc file are Dates? This automatically uses the columns you supplied under --bdcol, --don and --din. You only need to specify this, if you have additonal Dates that should be reformatted. Specify as comma separated list without spaces. Default: "[]"
  --DateFormat DATEFORMAT
                        Format of Dates in your data. Default "%d/%m/%Y" --> 31/01/1980
  --MinMaxAge MINMAXAGE
                        Filter by Min and Max Age at first diagnosis. This should be given as comma separated numerics in the form of x,y; i.e. 18,50 or 17.99,50.01, 0,0 (no exclusion) and is interpreted as inclusion as case if first diagnosis age >x and <y. Default: "0,0"
  --Fyob FYOB           Filter by Year of Birth. Everyone before this date will be excluded. Date needs to be given in the following format: "YYYY-MM-DD". Default: ""
  --iPSYCHComparison, --iPsychComparison, --ipsych-comparison
                        Build the user-supplied -g phenotypes within the iPSYCH-comparison cohort.
  --iPSYCHBirthStart IPSYCHBIRTHSTART
                        Inclusive birthdate boundary. Default: "1981-05-01"
  --iPSYCHBirthEnd IPSYCHBIRTHEND
                        Exclusive birthdate boundary. Default: "2009-01-01"
  --iPSYCHDxStart IPSYCHDXSTART
                        Inclusive diagnosis-episode boundary; blank means no lower bound.
  --iPSYCHDxEnd IPSYCHDXEND
                        Exclusive diagnosis-episode boundary for supplied definitions; presets use their own upper dates. Default: "2016-01-01"
  --iPSYCHDKParentsCol IPSYCHDKPARENTSCOL
                        Optional DKDK reporting column; otherwise auto-detected or derived from fkode_m/fkode_f. This does not filter IIDs.
  --iPSYCHBirthplaceCol IPSYCHBIRTHPLACECOL
                        Child birthplace-code column. Default: "fkode"
  --iPSYCHMotherCol IPSYCHMOTHERCOL
                        Maternal identifier column. Default: "pnr_m"
  --iPSYCHKeepMultiples
                        Keep repeated mother/birthdate combinations; the default excludes all such IIDs.
  --Fgender FGENDER     Filter by gender. Only Individuals with the selected Gender ("F" or "M") will be included in the output. Default: ""
  --eCc                 Exclude Individuals that are identified as Controls and are only part of CHB and not DBDS. This flag only works when ran on CHB/DBDS cluster.
  --removePointInDiagCode
                        If your supplied diagnostic codes contain a point but your -f file does not
  --skipICDUpdate       If your supplied diagnostic codes are already in line with your datastructure, you can use this toggle to skip the updating (highly recommended if your input is in a correct format!)
  --MatchFI             If you want to keep only the IIDs overlapping between -g and -f use this flag.
  --BuildEntryExitDates
                        If you want to extract additional information about the controls. This basically adds first and last date of observation (any diagnostic code) and writes it to first in date, last out date, in_dates, and out_dates. NB: This can be used if you don't have an Entry and Exit date known. Be aware, this will increase time needed for calculation.
  --Ophold OPHOLD       If you want to load and process the ophold file (on IBP Cluster only). Default ""
  --BuildOphold         If you want to update the ophold file used (on IBP Cluster only).
  --RegisterRun         We will try to determine this automatically based on known Servers. If you are using this method on an unknown Server and your diagnostic codes are in the follwing format DYXXxx; where Y stands for the letter of the group and XX for the main and xx for the subcode.
  --lpp                 Set this if you want to load phenotypes (one phenotype per file) and run our exclusions on them. Keep in mind, that -g will only allow one file to be supplied. All others will need to be supplied through the appropriate exclusion flags.: "False" (or 'date_out' in CHB/DBDS)
  --write_pickle        If you want to write the results to a pickle format in addition to the standard output.
  --write_fastGWA_format
                        If you want to write the phenotype into fastGWA format.
  --write_Plink2_format
                        If you want to write the phenotype into PLINK2 format.
  --write_omop_cohort   Write OMOP CDM COHORT and COHORT_DEFINITION result tables.
  --omopPersonMap       Protected IID-to-person_id mapping file.
  --omopPersonMapIidCol Source-IID column in the mapping; blank uses --iidcol.
  --omopPersonIdCol     OMOP person-ID column; default person_id.
  --omopCohortDefinitions
                        OMOP cohort-definition manifest.
  --omopOutputDir       Result directory; default <output-stem>.omop.
  --omopsep             Input separator for both OMOP files; default tab.
  --BuildTestSet        Build the test set that can be used to see an example of the input data or to test if your setup runs smoothly
  --testRun             Run only on smaller test data input (only available for the test dataset)
  --nthreads NTHREADS   DEPRECATED! - How many threads should be used. Default: "8"
  --lowmem              Experimental! - Processes plaintext input in IID batches and appends each batch to one output file using a fixed precomputed output schema. For HDF5 input, use --verylowmem for batched extraction.
  --batchsize BATCHSIZE
                        Experimental! - This will set the batches (when --lowmem is set) to the desired value. Default: "100000"
  --verylowmem          Experimental! - Applies low-memory batching to HDF5-indexed input.
  --BuildIndex          Build an indexed HDF5 file for faster repeated low-memory runs.
  --IndexDtypes         JSON dtype mapping used while building HDF5 indexes.
  --registercol         Register/source column; shared name or one name per -f file (also INI: registercol).
  --PSYK                Keep rows whose register/source contains PSYK (case insensitive).
  --LPR                 Keep rows whose register/source contains LPR (case insensitive).
  --verbose             Verbose output
```

### OMOP cohort result tables

Use `--write_omop_cohort` with a protected IID mapping and a definition manifest:

```bash
./get_pheno.py ... --write_omop_cohort \
  --omopPersonMap iid_to_person.tsv \
  --omopCohortDefinitions omop_cohorts.tsv
```

The person map requires the source IID column (defaults to `--iidcol`) and `person_id`. Source IIDs are not written to the OMOP files, and the mapping is not copied. Missing mappings, duplicate IDs, invalid dates, and undefined cohort IDs stop the export.

The definition manifest requires:

```text
phenotype_name  cohort_definition_id  cohort_definition_name  definition_type_concept_id  subject_concept_id
```

Optional columns are `status_column`, `case_value`, `cohort_start_policy` (`earliest`, `latest`, or `column`), `cohort_start_date_column`, `cohort_end_policy` (`start`, `latest`, or `column`), `cohort_end_date_column`, description, syntax, and initiation date. Defaults are the phenotype status column, `Case`, earliest diagnosis, and an end date equal to the start date.

Controls are absent from a case cohort by default. To export a distinct control cohort, supply `control_cohort_definition_id`, `control_cohort_definition_name`, and `control_start_date_column`; optional control fields include `control_value`, `control_end_date_column`, description, syntax, and concept-ID overrides.

The output directory contains `COHORT.tsv`, `COHORT_DEFINITION.tsv`, and `GETPHENO_PROVENANCE.json`. These files are also included in `<output-stem>.checksums.sha256`. The exporter validates the required fields, identifier references, dates, and duplicate periods against the [OMOP CDM 5.4 cohort-results structure](https://ohdsi.github.io/CommonDataModel/cdm54.html#cohort). Vocabulary validation and the OHDSI Data Quality Dashboard still require loading the files into the target OMOP environment.

### Minimum output in file based on standard settings:
|Column|Explanation|
| ------------- | ------------- |
|cpr_enc|IID of the individual|
|diagnosis|Diagnosis of phenotype: "Case" or "Control"; within DBDS for CHB individuals "Non-Case" instead of "Control"|
|first_dx|Date of first Diagnosis|
|last_dx|Date of last Diagnosis|
|diagnoses|All diagnostic codes for this IID on the phenotype|
|in_dates|All in_dates for the corresponding "diagnses" (Start of the current diagnosis; NOT INPATIENT DATE)|
|out_dates|All out_dates for the corresponding "diagnses" (End of the current diagnosis; NOT OUTPATIENT DATE)|
|n_diags|Number of diagnoses identified|
|n_unique_in_days|Number of unique days of diagnoses|
|sex|Gender|
|birthdate|Birthdate|
|Age_FirstDx|Age at first diagnosis|

When multiple phenotypes are requested, each phenotype gets its own status and detail columns:

|Column|Explanation|
| ------------- | ------------- |
|PhenoName|Case/Control status for that phenotype|
|PhenoName_Codes|Codes identified for that phenotype|
|PhenoName_In_Dates|Input dates for the corresponding phenotype codes|
|PhenoName_Out_Dates|Output dates for the corresponding phenotype codes|
|PhenoName_earliest_date|Earliest date identified for that phenotype|
|PhenoName_latest_date|Latest date identified for that phenotype|
|PhenoName_n_diags|Number of diagnoses identified for that phenotype|
|PhenoName_n_unique_in_days|Number of unique diagnosis days for that phenotype|

### Additional output in file based on Covariate or ExDEP settings:
|Column|Explanation|
| ------------- | ------------- |
|Covariate [i.e. ADHD]|Name of the additional covariate phenotype requested as columnname (this can include many covariates). The values represent the diagnostic codes identified for this IID|
|Covariate_Dates|Dates of diagnoses of the corresponding Covariate|
|Covariate_earliest_date|Earliest date of diagnosis of the corresponding Covariate|
|diagnoses_Level2_modifier|Diagnostic code modifying the main phenotype diagnosis|
|disorder_Level2_modifier|The disorder name for the corresponding diagnostic code|
|date_Level2_modifier|Date for the corresponding diagnostic code|
|Level2_diagnoses|Diagnostic codes left after the current exclusions|
|Level2_dates|Dates for all corresponding diagnostic codes|
|Level2_AgeExclusion|If this IID would be excluded due to Age restrictions|
|Level2_ExclusionReason|Which ecxclusion reasons used|
|Level2_FirstDx|Age at first diagnosis after these exclusions|
|Pheno_diagnoses_in_percent_lost_due_to_Covariate|Percentage of phenotype diagnoses lost due to current covariate|
|Covariate_Inflicted_changes|Per Covariate information about which main phenotype diagnoses were inflicted through the current covariate|
|Covariate_Inflicted_changes_MDD_dates|Correspondig dates|
|Covariate_Inflicted_changes_MDD_codes|Corresponding codes|
|Covariate_Inflicted_changes_MDD_dXnumber|Number of corresponding changes|
|Level3_CaseControl|Case (1) Control (0) status after all but Age exclusions|
|Level3_Age_FirstDx|Age at corresponding first diagnosis|
|Level3_CaseControl_AgeExclusions|Case (1) Control (0) status after all (incl. Age) exclusions|
|Level3_Sankey|Sankey information, otherwise not needed (for plots)|
|Level3_Sankey_data|Sankey information, otherwise not needed (for plots)|
|Level3_Sankey_source|Sankey information, otherwise not needed (for plots)|
|Level3_Sankey_target|Sankey information, otherwise not needed (for plots)|
|Level3_Sankey_value|Sankey information, otherwise not needed (for plots)|

### Output ordering and integrity files

The output columns start with the columns from the STAM/input IID file (`-i`), followed by optional additional IID information (`-j`), phenotype columns, exclusion/covariate columns, and derived output fields. Cohort membership columns such as `chb`, `dbds`, `degen`, `degen_new`, `degen_old`, `dcb`, `drb`, and `sle` are only written when they are present in the `-i` or `-j` input headers; when multiple are present, they are kept together in the input-column block. Cohort columns are metadata and do not implicitly remove IIDs with no active cohort assignment; exclusions require an explicit filter such as `--ge` or `--exclCHBcontrols`.

Each run writes a checksum file named `<output-stem>.checksums.sha256`. The generated run README contains the matching command:

```
sha256sum -c ./<output-stem>.checksums.sha256
```
