FROM continuumio/miniconda3
LABEL authors="Michael J. Bale" \
      description="Docker image containing all software requirements for the michaeljbale/CnRFlow pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/CnRFlow-0.9.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name CnRFlow-0.9.0 > CnRFlow-0.9.0.yml