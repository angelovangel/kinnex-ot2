FROM --platform=linux/amd64 rocker/shiny:4.6.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libuv1-dev \
    git \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Install Opentrons CLI in a virtual environment and symlink opentrons_simulate into PATH
RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir opentrons==8.8.1 && \
    ln -s /opt/venv/bin/opentrons_simulate /usr/local/bin/opentrons_simulate

ENV R_REPOS="https://packagemanager.posit.co/cran/__linux__/jammy/latest"

WORKDIR /app

COPY renv.lock .Rprofile ./
COPY renv/activate.R renv/

RUN R -e "options(repos = c(CRAN = Sys.getenv('R_REPOS'))); install.packages('renv'); renv::restore()"

COPY . .

EXPOSE 3838

CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
