FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# --- Dependencies ---
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    # Python dependencies
    python3 \
    python3-pip \
    python3-dev \
    # R dependencies
    r-base \
    r-base-dev \
    gfortran \
    # General build tools
    build-essential \
    cargo \
    wget \
    pkg-config \
    # C-library dependencies for Python/R packages
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libtirpc-dev \
    # C-library dependencies for R graphics
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev     \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev 


# ---  Packages ---
COPY install_r_packages.r .
RUN Rscript install_r_packages.r
COPY python_requirements.txt .
RUN pip3 install --break-system-packages --no-cache-dir -r python_requirements.txt 


# Security in a world where no one cares tho
COPY --chown=1000:1000 app /app

WORKDIR /app
USER 1000

RUN chmod +x start.sh

EXPOSE 8501


