#!/usr/bin/env bash
set -euo pipefail

config="${GETPHENO_SHINY_CONFIG:-/config/projects.conf}"
host="${GETPHENO_HOST:-0.0.0.0}"
port="${GETPHENO_PORT:-3838}"
runtime_user="${GETPHENO_UNIX_USER:-}"

if [[ ! -r "$config" ]]; then
    printf 'GetPheno project registry is not readable: %s\n' "$config" >&2
    exit 1
fi
if [[ ! "$port" =~ ^[0-9]+$ ]] || (( port < 1 || port > 65535 )); then
    printf 'GETPHENO_PORT must be an integer from 1 to 65535.\n' >&2
    exit 1
fi
if [[ -n "$runtime_user" && ! "$runtime_user" =~ ^[A-Za-z_][A-Za-z0-9_.-]*$ ]]; then
    printf 'GETPHENO_UNIX_USER contains unsupported characters.\n' >&2
    exit 1
fi

mkdir -p "$HOME"
chmod 0700 "$HOME"

exec Rscript /opt/getpheno/shiny_app/run_app.R \
    "--config=$config" \
    "--host=$host" \
    "--port=$port"
