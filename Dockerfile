FROM rocker/r-ver:4.3.3

ARG DEBIAN_FRONTEND=noninteractive
ARG GETPHENO_REVISION=unknown

LABEL org.opencontainers.image.source="https://github.com/MischaLundberg/GetPheno" \
      org.opencontainers.image.revision="${GETPHENO_REVISION}"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        python3 \
        python3-pip \
        slurm-client \
    && rm -rf /var/lib/apt/lists/*

COPY docker/requirements.txt /tmp/getpheno-requirements.txt
RUN python3 -m pip install --no-cache-dir -r /tmp/getpheno-requirements.txt \
    && rm /tmp/getpheno-requirements.txt \
    && Rscript -e 'install.packages("shiny", repos = getOption("repos"))'

WORKDIR /opt/getpheno
COPY . /opt/getpheno
RUN chmod 0755 /opt/getpheno/get_pheno.py /opt/getpheno/docker/entrypoint.sh \
    && python3 -m py_compile /opt/getpheno/get_pheno.py \
    && printf '%s\n' "${GETPHENO_REVISION}" > /opt/getpheno/IMAGE_REVISION

ENV GETPHENO_SHINY_CONFIG=/config/projects.conf \
    GETPHENO_HOST=0.0.0.0 \
    GETPHENO_PORT=3838 \
    HOME=/tmp/getpheno-home \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

EXPOSE 3838

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD curl --fail --silent --show-error "http://127.0.0.1:${GETPHENO_PORT}/" >/dev/null || exit 1

ENTRYPOINT ["/opt/getpheno/docker/entrypoint.sh"]
