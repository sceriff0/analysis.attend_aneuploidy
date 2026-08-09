# syntax=docker/dockerfile:1
# =============================================================================
# ATTEND — reproducible R environment
# =============================================================================
# Base: official Bioconductor image for the RELEASE_3_18 line. This pins
# R 4.3.2 + Bioconductor 3.18 with BiocManager already configured to that
# release — the exact contract in CLAUDE.md. (Bioconductor 3.18 ships R 4.3.x.)
#
# Build:   docker build -t attend:3.18 .
# Run an analysis shell:
#   docker run --rm -it -v "$PWD":/home/rstudio/ATTEND -w /home/rstudio/ATTEND attend:3.18 R
# Knit the whole site:
#   docker run --rm -v "$PWD":/home/rstudio/ATTEND -w /home/rstudio/ATTEND attend:3.18 \
#     Rscript -e 'workflowr::wflow_build()'
# RStudio Server (browse http://localhost:8787, user rstudio / pass attend):
#   docker run --rm -p 8787:8787 -e PASSWORD=attend \
#     -v "$PWD":/home/rstudio/ATTEND attend:3.18
# =============================================================================
FROM bioconductor/bioconductor_docker:RELEASE_3_18

LABEL org.opencontainers.image.title="ATTEND" \
      org.opencontainers.image.description="Multi-omic cancer-cohort QC & analysis pipeline (workflowr + renv, R 4.3.2 / Bioc 3.18)" \
      org.opencontainers.image.source="https://github.com/bolt3x/ATTEND"

# -----------------------------------------------------------------------------
# System libraries.
# The Bioconductor base already carries most of these, but the spatial stack
# (sf, spatstat.geom) needs the GDAL/GEOS/PROJ/udunits headers explicitly,
# and arrow benefits from the bundled-build prerequisites.
# -----------------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgdal-dev \
        libgeos-dev \
        libproj-dev \
        libudunits2-dev \
        libgsl-dev \
        libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /home/rstudio/ATTEND

# -----------------------------------------------------------------------------
# Package installation.
#
# renv is NOT yet initialized in this project (no renv.lock). Until it is, we
# install the dependency set explicitly from the loader below. The Dockerfile
# is written to flip to the reproducible path automatically the moment a
# renv.lock appears:
#
#   * renv.lock present  -> renv::restore() (exact pinned versions)
#   * renv.lock absent   -> install the curated CRAN + Bioconductor list
#
# To switch permanently once you run `renv::init(bioconductor = "3.18")`:
# just rebuild — renv.lock will be copied in and restore() takes over.
# -----------------------------------------------------------------------------

# Copy only the dependency manifests first so this layer caches across code edits.
# install_deps.R is always present, which keeps the COPY valid even though the
# optional renv.lock / .Rprofile globs currently match nothing (a glob matching
# zero files would otherwise fail the build). Once `renv::init()` creates them,
# they are picked up automatically and install_deps.R switches to renv::restore().
COPY install_deps.R renv.lock* .Rprofile* ./

RUN Rscript install_deps.R

# -----------------------------------------------------------------------------
# Project source. Mount-over at runtime with -v for live editing; the COPY
# below makes the image self-contained for CI / one-shot knits.
# -----------------------------------------------------------------------------
COPY . .

# Default to an R session. Override with the run commands documented in the header.
CMD ["R"]
