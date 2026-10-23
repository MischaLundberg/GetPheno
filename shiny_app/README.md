# GetPheno Shiny launcher

This application submits GetPheno jobs to Slurm and retrieves their files by a generated run code. It has one runtime dependency beyond base R: `shiny`.

## Security model

Run one Shiny process per authenticated Unix user. The app then calls `sbatch` as that user, and project access is enforced by the cluster's Unix groups and filesystem ACLs. The configured project is displayed only when all of the following are true:

- The optional `AccessGroup` matches one of the user's Unix groups.
- `ProjectRoot` and `IniFile` are readable.
- `OutputRoot` is writable.

Every run is stored below `OutputRoot/<unix-user>/<run-code>/`, with user and run directories set to mode `0700`. The run code is a locator, not an authentication credential.

Do not expose a single Shiny process running as a shared service account and assume that a login form changes its filesystem identity. For a shared HTTPS service, use an existing HPC portal that starts the app under the authenticated user's UID, such as Open OnDemand. A custom password database is intentionally not included.

## Configure projects

Copy the example registry and edit it as an administrator:

```bash
cd shiny_app
cp projects.conf.example projects.conf
chmod 0644 projects.conf
```

The registry uses base R's DCF format. Add one record per project/protocol combination, separated by a blank line:

```text
Project: oral
ProjectLabel: IBP Oral project
Protocol: chb_dbds
ProtocolLabel: CHB/DBDS
ProjectRoot: /data/preprocessed/project_data_oral
IniFile: /dprhdbds/shared/data/projects_chb/ibp_oral/config/get_pheno.ini
OutputRoot: /dprhdbds/shared/data/projects_chb/ibp_oral/getpheno_runs
AccessGroup: ibp_oral
Python: /services/tools/conda/envs/GetPheno/bin/python
GetPheno: /dprhdbds/shared/apps/GetPheno/get_pheno.py
FreezePattern: /data/preprocessed/project_data_oral/chb_oral_freeze_cohort_*.tsv
DiagnosisFile: /data/preprocessed/project_data_oral/oral_registry_lpr_lpr3_diagnoses_pat_type.tsv
StandardAdditionalFiles: /data/preprocessed/project_data_oral/t_person.tsv|/data/preprocessed/project_data_oral/t_dodsaarsag_2.tsv
OptionalAdditionalFiles: Epikur=/data/preprocessed/project_data_oral/epikur.tsv|Laboratory=/data/preprocessed/project_data_oral/lab_forsker.tsv|SSIK=/data/preprocessed/project_data_oral/t_ssik.tsv|SSR=/data/preprocessed/project_data_oral/t_ssr.tsv
DeathDateColumn: D_DODSDATO
GeneralExclusionsTemplate: /data/preprocessed/genetics/chb_oral_freeze_{freeze_base}/oral_exclusion_latest
QcedTemplate: /data/preprocessed/genetics/chb_oral_freeze_{freeze_base}/01.Genotypes/GSA_oral_qc.fam
Sbatch: /usr/bin/sbatch
Partition: normal
Account: ibp_oral
Time: 1-00:00:00
Memory: 32G
Cpus: 4
```

This oral record maps the newest selected freeze to `-i`, the diagnosis file to `-f`, `t_person.tsv` plus selected optional files to `-j`, and the freeze-specific genetics paths to `--ge` and `--qced`. `{freeze_base}` is the first eight-digit date in the selected cohort filename; `{freeze_version}` is the last. The app lists readable files matching `FreezePattern`, sorts their dates newest first, and selects the newest by default.

Protocols with multiple diagnosis/event files use `DiagnosisFiles` with pipe-separated paths and optional per-file column mappings:

```text
DiagnosisFiles: /data/project/lpr3_diagnoses.tsv|/data/project/lpr3_procedures.tsv
DiagnosticColumns: diagnosis,code
IidColumns: cpr_enc
DateInColumns: date_in
DateOutColumns: date_out
DiagnosisTypeColumns: diagnosis_type,
```

The launcher converts the path separator to the comma-separated `-f` argument and passes the mappings as `--fcol`, `--iidcol`, `--din`, `--don`, `--DiagTypecol`, and `--recnum` (`RecordNumberColumns`). `DiagnosisFile` remains supported for a single file. Each mapping uses one shared name or exactly one name per file; a trailing empty entry explicitly represents an unavailable column.

`OptionalAdditionalFiles` is an administrator-controlled allowlist in `Label=/path|Label=/path` form. Users can select any available entries in addition to `StandardAdditionalFiles`; they cannot enter arbitrary server paths. All `-j` files must use the protocol's configured separator and IID column. GetPheno merges their columns by IID in listed order, with the standard file taking precedence for overlapping nonmissing values.

Set `DeathDateColumn` only when the configured standard inputs provide that column. The app then shows age-at-death, retain-living, and deceased-only controls and passes the column through `--deathdatecol`. The example keeps `t_dodsaarsag_2.tsv` in `StandardAdditionalFiles` so the configured `D_DODSDATO` column is always available.

Use a separate record and INI file for each selectable protocol. Add records for the available project roots as required: `project_data_brain`, `project_data_cid`, `project_data_cvd`, `project_data_degen`, `project_data_repro`, `project_data_cancer`, `project_data_covid`, `project_data_dbds`, and `project_data_oral`, all below `/data/preprocessed/`. Their exact protocol files and access groups must be administrator-configured; they are not inferred from directory names.

Create each output root with permissions appropriate for its project group. For example:

```bash
mkdir -p /cluster/projects/oral/getpheno_runs
chgrp oral_project /cluster/projects/oral/getpheno_runs
chmod 2770 /cluster/projects/oral/getpheno_runs
```

Keep the registry and INI files administrator-owned. Paths containing whitespace are rejected because Slurm output directives do not handle them consistently.

## Install and run

```bash
Rscript shiny_app/install.R
Rscript shiny_app/run_app.R \
  --config=/absolute/path/to/projects.conf \
  --host=127.0.0.1 \
  --port=3838
```

For a per-user session, connect through SSH instead of exposing the port:

```bash
ssh -L 3838:127.0.0.1:3838 USER@CLUSTER
```

Then open `http://127.0.0.1:3838`. Run the application only on a node where cluster policy permits persistent interactive services and Slurm client commands.

## Docker deployment

The repository-level `Dockerfile` contains the `get_pheno.py` revision present in the Docker build context, Python 3, the Python runtime dependencies, R, Shiny, and Slurm client commands. It does not download source code when the container starts. The GitHub workflow publishes every `main` revision as `ghcr.io/mischalundberg/getpheno:edge`. Pushing a release tag such as `V1.30` publishes `:V1.30` and advances `:latest`; every build also receives a commit-SHA tag. It uses the fixed `rocker/r-ver:4.3.3` base. `compose.yaml` pulls the configured image, starts one non-root container for one authenticated Unix user, and publishes Shiny only on `127.0.0.1`.

For production, set an explicit release in `.env.docker` rather than relying on the moving `latest` tag:

```text
GETPHENO_IMAGE=ghcr.io/mischalundberg/getpheno:V1.30
```

No Dockerfile edit is needed for the next release: create and push the next Git tag after the release commit is finalized. The workflow builds directly from that checked-out tag and records its full commit SHA in the image metadata.

Docker does not replace the app's Unix authorization model. Start one container per Unix user with that user's numeric UID, primary GID, login name, and project GIDs. The kernel then applies the mounted filesystem's normal ownership, group, and ACL checks.

Create the deployment environment:

```bash
cp .env.docker.example .env.docker
sed -i "s/^GETPHENO_UID=.*/GETPHENO_UID=$(id -u)/" .env.docker
sed -i "s/^GETPHENO_GID=.*/GETPHENO_GID=$(id -g)/" .env.docker
sed -i "s/^GETPHENO_UNIX_USER=.*/GETPHENO_UNIX_USER=$(id -un)/" .env.docker
```

Set `GETPHENO_ACCESS_GID` to one required project group's numeric GID. For access to multiple project groups, copy `compose.groups.yaml.example` to `compose.groups.yaml`, list every required numeric GID, and include both Compose files below.

Set these absolute paths in `.env.docker`:

- `GETPHENO_PROJECTS_CONFIG`: the administrator-owned project registry.
- `GETPHENO_PREPROCESSED_ROOT`: read-only source-data root, mounted at the same absolute path.
- `GETPHENO_SHARED_ROOT`: shared project/output root, mounted read-write and protected by Unix permissions.
- `GETPHENO_SERVICES_ROOT`: read-only software root containing the configured Python interpreter.
- `GETPHENO_SLURM_CONF_DIR`: cluster Slurm configuration directory.
- `GETPHENO_MUNGE_DIR`: directory containing the cluster Munge socket.

Validate and start:

```bash
docker compose --env-file .env.docker \
  -f compose.yaml -f compose.groups.yaml config
docker compose --env-file .env.docker \
  -f compose.yaml -f compose.groups.yaml pull
docker compose --env-file .env.docker \
  -f compose.yaml -f compose.groups.yaml up -d
docker compose --env-file .env.docker \
  -f compose.yaml -f compose.groups.yaml logs -f getpheno
```

When only one project GID is needed, omit `-f compose.groups.yaml`. Open `http://127.0.0.1:3838`, or forward that address over SSH as described above.

If the GHCR package is private, authenticate first with a GitHub token that has `read:packages`:

```bash
echo "$GITHUB_TOKEN" | docker login ghcr.io -u YOUR_GITHUB_USER --password-stdin
```

To build the current local checkout instead of pulling GHCR, include `compose.build.yaml`:

```bash
docker compose --env-file .env.docker \
  -f compose.yaml -f compose.build.yaml -f compose.groups.yaml build
docker compose --env-file .env.docker \
  -f compose.yaml -f compose.build.yaml -f compose.groups.yaml up -d
```

The image records its source commit in `/opt/getpheno/IMAGE_REVISION` and the OCI revision label:

```bash
docker image inspect "$GETPHENO_IMAGE" \
  --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}'
```

### Slurm requirements

The registry's `Python`, `GetPheno`, input, and output paths are written into the submitted batch script. They must therefore exist at the same absolute paths inside the Shiny container and on every relevant Slurm compute node. Consequently, pulling a newer Shiny image does not update a separately configured shared `GetPheno:` script. Do not set `GetPheno: /opt/getpheno/get_pheno.py` unless `/opt/getpheno` is also a shared compute-node path; use or update the shared installation for submitted jobs. Running the exact containerized revision on compute nodes requires a cluster-supported container runtime such as Apptainer and corresponding Slurm integration.

The container mounts Slurm configuration and the Munge socket but does not run Slurm or Munge daemons. Confirm before accepting submissions:

```bash
docker compose --env-file .env.docker exec getpheno id
docker compose --env-file .env.docker exec getpheno sbatch --version
docker compose --env-file .env.docker exec getpheno sinfo
```

The displayed UID and groups must match the intended user. The Slurm client version must be compatible with the cluster. Some clusters also require host networking or additional authentication/configuration mounts; these are cluster-specific and must be approved by the Slurm administrator.

The image cannot be built in an environment without Docker/Podman access. On HPC systems where Docker is prohibited, use Apptainer or run the documented native R process under Open OnDemand instead.

## Submitted options

The phenotype definition, project, protocol, and any protocol-configured cohort freeze are the run inputs. All other form groups are explicitly marked **Optional** and start collapsed. The form supports:

- Uploaded or pasted phenotype definitions.
- An optional pasted or uploaded IID subset containing one IID per line without a header.
- Uploaded or pasted lifetime, post-onset, and one-year-prior exclusions.
- iPSYCH comparison, exact matching, ICD8-only, and ICD10-only modes.
- fastGWA and PLINK2 output.
- Gender, earliest birth date, age range, and diagnosis-type filters.
- Protocol-dependent age-at-death, retain-living, and deceased-only filters.
- `--name` when exactly one phenotype is detected.
- A newest-first cohort-freeze selector when `FreezePattern` is configured.
- Protocol-standard and allowlisted optional additional-information files.

Hover over or focus an information icon beside any option to see its accepted format and effect. The phenotype help includes simple, multiple-phenotype, advanced `main=`/`sub=`/`rule_out=`, wildcard, and exclusion examples.

Each submission creates a unique identifier such as `GP-20261002T120000Z-a1b2c3d4e5f6`. The identifier is used in the run directory, phenotype request, Slurm script, and result filename. The Retrieve run tab shows Slurm state and allows files to be downloaded after reconnecting.

## Test

The utility test does not require Shiny or Slurm. The server test uses Shiny when it is installed and still uses a fake scheduler:

```bash
Rscript shiny_app/tests/test_job_utils.R
Rscript shiny_app/tests/test_shiny_server.R
```

They validate registry parsing, project access, phenotype counting, safe command construction, Slurm-script execution, job-ID parsing, metadata, submission, and run retrieval. The server test skips cleanly when Shiny is not installed.
