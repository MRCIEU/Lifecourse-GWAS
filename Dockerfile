FROM rocker/r-ver:latest AS renv_base

WORKDIR /usr/local/src/myscripts
RUN apt-get update && \
    apt-get install -yyy \
        build-essential \
        libcurl4-gnutls-dev \
        libxml2-dev \
        libssl-dev \
        libgmp3-dev \
        cmake \
        libcairo2-dev \
        libxt-dev \
        libharfbuzz-dev \
        libtiff-dev \
        libzstd-dev \
        git \
        libgit2-dev \
        libfribidi-dev \
        pandoc \
        wget && \
        wget https://launchpad.net/ubuntu/+source/icu/70.1-2/+build/23145450/+files/libicu70_70.1-2_amd64.deb && \
        dpkg -i libicu70_70.1-2_amd64.deb && \
        ln -s /usr/lib/x86_64-linux-gnu/libgit2.so.1.7 /usr/lib/x86_64-linux-gnu/libgit2.so.1.1

WORKDIR /project
COPY renv.lock renv.lock

RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

RUN R -e 'options( \
    repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"), \
    HTTPUserAgent = sprintf( \
        "R/%s R (%s)", \
        getRversion(), \
        paste(getRversion(), \
          R.version["platform"], \
          R.version["arch"], \
          R.version["os"]))); \
    renv::install(c("stringi", "git2r")); \
    renv::restore(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"))'

COPY . .