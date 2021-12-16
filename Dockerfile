FROM nfcore/base
LABEL description="Docker image containing all requirements for pfitmap-nextflow workflow"
COPY environment.yml /

RUN conda install -c conda-forge mamba
RUN mamba clean -a
RUN mamba env create -f /environment.yml && mamba clean -a

ENV PATH /opt/conda/envs/pfitmap-nextflow-env/bin:$PATH
