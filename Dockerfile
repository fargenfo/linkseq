FROM nfcore/base:latest

LABEL \
    authors="olavur@fargen.fo" \
    description="Image with tools used in exolink" \
    maintainer="Ólavur Mortensen <olavur@fargen.fo>"

RUN apt update -yqq && \
    apt install -yqq \
    unzip \
    ttf-dejavu

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/exolink-alpha/bin:$PATH
