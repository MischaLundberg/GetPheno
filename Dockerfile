FROM rocker/r-ver:4.3.3

ARG DEBIAN_FRONTEND=noninteractive
ARG GETPHENO_VERSION=development
ARG GETPHENO_REVISION=unknown

LABEL org.opencontainers.image.title="GetPheno" \
      org.opencontainers.image.description="Extract phenotypes from biobank-scale files" \
      org.opencontainers.image.source="https://github.com/MischaLundberg/GetPheno" \
      org.opencontainers.image.version="${GETPHENO_VERSION}" \
      org.opencontainers.image.revision="${GETPHENO_REVISION}"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        python3 \
        python3-pip \
        python3-dev \
        build-essential \
        pkg-config \
        libhdf5-dev \
        liblzo2-dev \
        libbz2-dev \
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
    && printf '%s\n' "${GETPHENO_REVISION}" > /opt/getpheno/IMAGE_REVISION \
    && printf 'GetPheno %s (%s)\n' "${GETPHENO_VERSION}" "${GETPHENO_REVISION}" > /etc/getpheno-version \
    && printf '#!/bin/sh\nexec python3 /opt/getpheno/get_pheno.py "$@"\n' > /usr/local/bin/getpheno \
    && printf '#!/bin/sh\nexec /opt/getpheno/docker/entrypoint.sh shiny "$@"\n' > /usr/local/bin/getpheno-shiny \
    && chmod 0755 /usr/local/bin/getpheno /usr/local/bin/getpheno-shiny

ENV GETPHENO_SHINY_CONFIG=/config/projects.conf \
    GETPHENO_HOST=0.0.0.0 \
    GETPHENO_PORT=3838 \
    HOME=/tmp/getpheno-home \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

WORKDIR /data

RUN groupadd --system getpheno \
    && useradd --system --gid getpheno --home-dir /tmp/getpheno-home getpheno \
    && mkdir -p /tmp/getpheno-home /data \
    && chown -R getpheno:getpheno /tmp/getpheno-home /data /opt/getpheno

USER getpheno

EXPOSE 3838

HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD curl --fail --silent --show-error "http://127.0.0.1:${GETPHENO_PORT}/" >/dev/null || exit 1

ENTRYPOINT ["/opt/getpheno/docker/entrypoint.sh"]
CMD ["shiny"]
