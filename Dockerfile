FROM nfcore/base
LABEL description="Docker image containing all requirements for pfitmap-nextflow workflow"
COPY environment.yml /

RUN conda update conda
RUN conda env create -f /environment.yml && conda clean -a

ENV PATH /opt/conda/envs/pfitmap-nextflow-env/bin:$PATH
