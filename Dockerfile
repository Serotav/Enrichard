FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PYTHONUNBUFFERED=1

# --- Dependencies ---
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    r-base \
    r-base-dev \
    python3-dev \
    gfortran \
    build-essential \
    cargo \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libtirpc-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# ---  Packages ---
COPY python_requirements.txt .
RUN pip3 install --break-system-packages --no-cache-dir -r python_requirements.txt 
COPY install_r_packages.r .
RUN Rscript install_r_packages.r


WORKDIR /app
COPY ./app .
COPY ./start.sh /start.sh

RUN chmod +x /start.sh

EXPOSE 8501
