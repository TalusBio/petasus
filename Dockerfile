# Base image for building
FROM python:3.11 as build

# Install essentials:
RUN apt-get update && apt-get install -y build-essential libhdf5-dev

# Install uv:
ADD https://astral.sh/uv/install.sh /install.sh
RUN chmod -R 655 /install.sh  && /install.sh && rm /install.sh

# Install our package:
ENV VIRTUAL_ENV=/opt/venv \
    PATH="/opt/venv/bin:$PATH"

ADD . /src
RUN /root/.cargo/bin/uv venv --no-cache /opt/venv && \
    /root/.cargo/bin/uv pip install --no-cache /src

# The app stage (slim debian build):
FROM python:3.11-slim-bookworm
COPY --from=build /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
