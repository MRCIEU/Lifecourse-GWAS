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
        libtiff5

# Wait until MRCIEU R-Universe has built latest version of the TwoSampleMR binary
# Should be 1 hour or maybe overnight
# Check https://mrcieu.r-universe.dev/TwoSampleMR
# I set the HTTPUserAgent option in order to obtain binary packages from the
# Public Posit Package Manager (otherwise source packages which will have to 
# be built will be obtained).

WORKDIR /project
COPY renv.lock renv.lock

RUN mkdir -p renv
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

RUN R -e 'options( \
    repos = c(universe = "https://mrcieu.r-universe.dev/bin/linux/jammy/4.3/", \
        CRAN = "https://packagemanager.posit.co/cran/__linux__/jammy/latest"), \
    HTTPUserAgent = sprintf( \
        "R/%s R (%s)", \
        getRversion(), \
        paste(getRversion(), \
          R.version["platform"], \
          R.version["arch"], \
          R.version["os"]))); \
    install.packages("renv", dependencies = TRUE); \
    renv::restore()'

COPY . .