# ExDep resources
[Build Phenotypes](#1-build-phenotypes)


# Getting Started (Delphi - IBP Cluster)
As Delphi is a pain to get anything started up, you need to load many modules to even get the easiest things to run... Here an example (this may be different as soon as something gets updated on the cluster).
```
module load tools
module load shared
module load shared #needs to be loaded twice
module load anaconda3/2021.11
module load gcc/9.4.0 intel/perflibs/64 apache-arrow/9.0.0_LIBS
module load R/4.3.0
module load rstudio/1.1.423
```
The script is located in the ```aSchorkLab\scripts``` directory


# Getting Started (common)
In order to download `ExDEP_Plotting`, you should clone this repository via the commands
```  
git clone https://github.com/MischaLundberg/ExDEP_Plotting
cd ExDEP_Plotting
```

In order to install the Python dependencies, you will need the [Anaconda](https://store.continuum.io/cshop/anaconda/) Python distribution and package manager. After installing Anaconda, run the following commands to create an environment with the dependencies:

```
conda env create --file environment.yml
conda activate ExDEP

## to deactivate the environment, type
#conda deactivate
```

In case you are updating your current version of the pipeline, it would be best practice to also update your environment to the updated prerequisetes.
If you are using a Anaconda environment, you can do so by typing
```
conda env update --name ExDEP --file environment.yml
```


If you receive any errors while running scripts of the ExDEP_Plotting pipeline, please ensure your versioning for the prerequisites is according to the tested versions.


## 1 Build Phenotypes:


### Input:
As Input you need a file that contains a header, followed by your diagnostic codes of interest (one per row).
```
diagnosis
F33.0
F32.0
```
Don't forget to specify your header --gcol (will be otherwise set automatically to diagnosis or c_adiag depending on the cluster; see below) 

A Simple run would look like

```
./get_phenotype.py -g diagnosis.request -o diagnosis.casecontrol
```
Keep in mind, that there may be some crucial flags to set: i.e. --eM
The script will (if not told differently) update the ICD codes from i.e. F32.0 to ICD10:DF320 automatically on CHB/DBDS and IBP clusters.

### Options (Flags):
```
usage: get_phenotype.py  [-h] [-f F] -g G [-i I] [-j J] [--name NAME] [--ExDepExc] [--ge GE] [--fcol FCOL] [--gcol GCOL]
                            [--iidcol IIDCOL] [--bdcol BDCOL] [--sexcol SEXCOL] -o O [--fsep FSEP] [--gsep GSEP] [--eM] [--din DIN]
                            [--don DON] [--DateFormat DATEFORMAT] [--qced QCED] [--testRun] [--eCc] [--iidstatus IIDSTATUS]
                            [--DiagTypeExclusions DIAGTYPEEXCLUSIONS] [--DiagTypeInclusions DIAGTYPEINCLUSIONS]
                            [--LifetimeExclusion LIFETIMEEXCLUSION] [--PostExclusion POSTEXCLUSION]
                            [--OneyPriorExclusion ONEYPRIOREXCLUSION] [--Fyob FYOB] [--Fgender FGENDER] [--nthreads NTHREADS]
                            [--BuildTestSet] [--removePointInDiagCode] [--skipICDUpdate] [--MatchFI] [--CasesOnly] [--verbose]

Extracts a Phenotype from input files based on IIDs and diagnostic codes. The best way to start, is to generate a test dataset: 'python
get_phenotype.py -g "" -o ./ --BuildTestSet' and then run 'python get_phenotype.py -g "" -o testrun.tsv --eM --ExDepExc --testRun'

optional arguments:
  -h, --help            show this help message and exit
  -f F                  Diagnosis file. Should have at least "IID", "date_in", "date_out", and "diagnosis" as columns. Names can be
                        defined using --iidcol, --din, --don, and --fcol. No default. Will be automatically determined if not entered for
                        GenomeDK and CHB/DBDS (DEGEN protocol)
  -g G                  File with all Diagnostic codes to export
  -i I                  This file should adds based on a mapping using the supplied "IID" column information about Gender/Sex (needed),
                        Brithdate (needed) ... e.g., information from stamdata file. No default. Will be automatically determined if not
                        entered for GenomeDK and CHB/DBDS (DEGEN protocol)
  -j J                  This file should have additional "IID" information e.g. date of entry to cohort. This is not needed and can also
                        be skipped. No default. Will be automatically determined if not entered for GenomeDK and CHB/DBDS (DEGEN
                        protocol)
  --name NAME           Define your main Phenotype name. Defaults to "MainPheno"
  --ExDepExc            List of diagnostic codes to exclude (e.g. Cases holding also these codes will be excluded). This is currently
                        precoded and will be changed to allowing to specify a file with codes.
  --ge GE               General exclusion list. This referrs to a list of IIDs that should be excluded from the study. Within CHB/DBDS it
                        will default to /data/preprocessed/genetics/chb_degen_freeze_20210503/degen_exclusion_latest
  --fcol FCOL           Columname of -f and -i file to be mapped. Defaulting to "c_adiag" (or diagnosis in CHB/DBDS)
  --gcol GCOL           Columname of -g file to be mapped against. Defaulting to "c_adiag" (or diagnosis in CHB/DBDS)
  --iidcol IIDCOL       Columname of IDs in -f and -i file. Defaulting to "pnr" (or cpr_enc in CHB/DBDS)
  --bdcol BDCOL         Columname of Birthdate in files. Defaults to "birthdate"
  --sexcol SEXCOL       Columname of Sex/Gender in files. Defaults to "sex"
  -o O                  Outfile name; dont forget to add the location, otherwise it will be saved in the local dir.
  --fsep FSEP           Separator of -f i.e. tab; default "," (or "\t" in DBDS/CHB)
  --gsep GSEP           Separator of -g i.e. tab; default "," (or "\t" in DBDS/CHB)
  --eM                  Exact Match. If the ICD or diagnosis code should max exactly and not e.g. searching for ICD:F33 and receiving
                        ICD:F331, ICD:F33, ICD:F339 and so on.
  --din DIN             Columname of first diagnosis date. e.g. 'd_inddto' or 'date_in'. Default: "d_inddto" (or 'date_in' in DBDS/CHB)
  --don DON             Columname of first diagnosis date. e.g. 'd_uddto' or 'date_out'. Default: "d_uddto" (or 'date_out' in DBDS/CHB)
  --DateFormat DATEFORMAT
                        Format of Dates in your data. Default "%d/%m/%Y" --> 31/01/1980
  --qced QCED           List with all IIDs that pass initial QC. Default: "" (or
                        '/data/preprocessed/genetics/chb_degen_freeze_20210503/DEGEN_GSA_FINAL.fam' in DBDS/CHB)
  --testRun             Run only on smaller test data input (only available for the test dataset)
  --eCc                 Exclude Individuals that are identified as Controls and are only part of CHB and not DBDS. This flag only works
                        when ran on CHB/DBDS cluster.
  --iidstatus IIDSTATUS
                        Information about the status of the IID, i.e. if they withdrew their consent (), moved outside the country (), or
                        died (). Default: ""
  --DiagTypeExclusions DIAGTYPEEXCLUSIONS
                        List of diagnostic types to exclude. e.g. "H,M". Potential c_types may be +,A,B,C,G,M,H. Default: ""
  --DiagTypeInclusions DIAGTYPEINCLUSIONS
                        List of diagnostic types to include. e.g. "A,B". Potential c_types may be +,A,B,C,G,M,H. Default: ""
  --LifetimeExclusion LIFETIMEEXCLUSION
                        Define Lifetime exclusions (if a case has any of the listed codes, it will be excluded). This should be a file
                        containing either one row with all codes listed (comma separated) or per row a exclusion name, followed by a list
                        of diagnostic codes (similar to the input). e.g. "BPD ICD10:F30,ICD10:F30.0,ICD10:F30.1,ICD10:F30.2,ICD10:F30.8,I
                        CD10:F30.9,ICD10:F31,ICD10:F31.0,ICD10:F31.1,ICD10:F31.2,ICD10:F31.3,ICD10:F31.4,ICD10:F31.5,ICD10:F31.6,ICD10:F3
                        1.7,ICD10:F31.8,ICD10:F31.9" If you are using CHB/DBDS data, please remember, that the ICD10 codes have to start
                        with ICD10:D***, e.g. ICD10:DF33.0. Default: ""
  --PostExclusion POSTEXCLUSION
                        Define Post exclusions (if a case has any of the listed codes, all main diagnoses after the first occuring date
                        of the listed codes will be excluded). This should be a file containing either one row with all codes listed
                        (comma separated) or per row a exclusion name, followed by a list of diagnostic codes (similar to the input).
                        Default: ""
  --OneyPriorExclusion ONEYPRIOREXCLUSION
                        Define One Year Prio exclusions (all entries for a case that happen within one year prior to any date of the
                        listed codes, these entries will be excluded). This should be a file containing either one row with all codes
                        listed (comma separated) or per row a exclusion name, followed by a list of diagnostic codes (similar to the
                        input). Default: ""
  --Fyob FYOB           Filter by Year of Birth. Everyone before this date will be excluded. Date needs to be given in the following
                        format: "YYYY-MM-DD". Default: ""
  --Fgender FGENDER     Filter by gender. Only Individuals with the selected Gender ("F" or "M") will be included in the output. Default:
                        ""
  --nthreads NTHREADS   How many threads should be used. Default: "8"
  --BuildTestSet        Build the test set that can be used to see an example of the input data or to test if your setup runs smoothly
  --removePointInDiagCode
                        If your supplied diagnostic codes contain a point but your -f file does not
  --skipICDUpdate       If your supplied diagnostic codes are already in line with your datastructure, you can use this toggle to skip
                        the updating (highly recommended if your input is in a correct format!)
  --MatchFI             If you want to keep only the IIDs overlapping between -g and -f use this flag.
  --CasesOnly           If you want to only extract cases use this flag.
  --verbose             Verbose output
```

### Minimum output in file based on standard settings:
|Column|Explanation|
| ------------- | ------------- |
|cpr_enc|IID of the individual|
|diagnosis|Diagnosis of phenotype: "Case" or "Control"; within DBDS for CHB individuals "Non-Case" instead of "Control"|
|first_dx|Date of first Diagnosis|
|last_dx|Date of last Diagnosis|
|diagnoses|All diagnostic codes for this IID on the phenotype|
|in_dates|All in_dates for the corresponding "diagnses"|
|out_dates|All out_dates for the corresponding "diagnses"|
|n_diags|Number of diagnoses identified|
|n_unique_in_days|Number of unique days of diagnoses|
|sex|Gender|
|birthdate|Birthdate|
|Age_FirstDx|Age at first diagnosis|

#### Further examples based on CHB/DBDS cluster (might be different on other sources; this only adds information from the input files without processing it):
|Column|Explanation|
| ------------- | ------------- |
|dbds|If on CHB/DBDS cluster, TRUE/FALSE if this IID is from DBDS|
|degen_old|If on CHB/DBDS cluster, TRUE/FALSE if this IID is from old CHB|
|degen_new|If on CHB/DBDS cluster, TRUE/FALSE if this IID is from curren CHB|
|C_STATUS||
|D_STATUS_HEN_START||
|D_FODDATO||
|C_KON||
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
|Level3_CaseControl|Case (1) Control (0) status after all exclusions|
|Level3_Age_FirstDx|Age at corresponding first diagnosis|
|Level3_CaseControl_AgeExclusions||
|Level3_Sankey|Sankey information, otherwise not needed (for plots)|
|Level3_Sankey_data|Sankey information, otherwise not needed (for plots)|
|Level3_Sankey_source|Sankey information, otherwise not needed (for plots)|
|Level3_Sankey_target|Sankey information, otherwise not needed (for plots)|
|Level3_Sankey_value|Sankey information, otherwise not needed (for plots)|
