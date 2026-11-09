FROM rocker/r-ver:4.4.3

ARG CRAN_REPOSITORY=https://packagemanager.posit.co/cran/2025-04-04
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
    && apt-get clean
RUN Rscript -e 'packages <- c("shiny", "validate", "httr2", "processx", "zip"); install.packages(packages, repos = commandArgs(TRUE)[1], Ncpus = 1); missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))' "${CRAN_REPOSITORY}"
RUN Rscript -e 'if (!requireNamespace("openssl", quietly = TRUE)) install.packages("openssl", repos = commandArgs(TRUE)[1]); stopifnot(requireNamespace("openssl", quietly = TRUE))' "${CRAN_REPOSITORY}"

COPY DESCRIPTION NAMESPACE /opt/omicsqc/
COPY R/ /opt/omicsqc/R/
COPY inst/ /opt/omicsqc/inst/
RUN R CMD INSTALL /opt/omicsqc && \
    useradd --system --create-home omicsqc && mkdir -p /tmp/omicsqc
COPY shiny/ /opt/omicsqc-shiny/

ENV TMPDIR=/tmp/omicsqc
USER omicsqc
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/opt/omicsqc-shiny', host='0.0.0.0', port=3838)"]
