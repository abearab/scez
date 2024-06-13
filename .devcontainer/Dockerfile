FROM continuumio/miniconda3

# Install system dependencies
USER root
RUN apt-get update \
  && apt-get install -y build-essential git tree curl \
  && apt-get clean \
  && apt-get purge 

# Install conda env
COPY environment.yml /tmp/env.yaml
RUN conda env update -n base -f /tmp/env.yaml \
    && conda clean --all --yes \
    && rm -f /tmp/env.yaml