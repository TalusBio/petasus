FROM mambaorg/micromamba:latest
MAINTAINER William E Fondrie <wfondrie@talus.bio>
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER dist/*.whl /tmp/

WORKDIR /app

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba install -c conda-forge -n base gcc gxx_impl_linux-64 gxx_linux-64 binutils_impl_linux-64 && \
    micromamba install -c conda-forge -n base binutils && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN WHEEL=$(ls /tmp/*.whl) && pip install ${WHEEL}

RUN pip install petasus

# FROM mambaorg/micromamba:latest
# MAINTAINER William E Fondrie <wfondrie@talus.bio>
# COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
# COPY --chown=$MAMBA_USER:$MAMBA_USER dist/*.whl /tmp/

# WORKDIR /app

# RUN micromamba install -y -n base -f /tmp/environment.yml && \
#     micromamba install -c conda-forge -n base gcc gxx_impl_linux-64 gxx_linux-64 binutils_impl_linux-64 && \
#     micromamba install -c conda-forge -n base binutils && \
#     micromamba clean --all --yes

# ARG MAMBA_DOCKERFILE_ACTIVATE=1
# RUN WHEEL=$(ls /tmp/*.whl) && pip install ${WHEEL}

# RUN echo "conda activate petasus" >> ~/.bashrc
# SHELL ["/bin/bash", "--login", "-c"]

# RUN pip install git+https://github.com/TalusBio/petasus
