# Container with cdtest and the ATAC-seq pipeline tools (environment-atac.yml).
#   docker build -t cdtest .
#   docker run --rm -v "$PWD":/work -w /work cdtest pytest -q
#   docker run --rm -v "$PWD":/work -w /work cdtest python scripts/reproduce_from_intermediates.py
FROM mambaorg/micromamba:1.5.10
COPY --chown=$MAMBA_USER:$MAMBA_USER environment-atac.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1
WORKDIR /opt/cdtest
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/cdtest
RUN pip install --no-cache-dir -e ".[full,test]"
WORKDIR /work
