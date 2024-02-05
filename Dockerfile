####################################################
# 
### BUILD IMAGE
# docker build -t icrsc/icr-8q24 .
#
### RUN IMAGE
# docker run -p 3838:3838 --rm --name icr-8q24 icrsc/icr-8q24
#
# Should then be visible on http://localhost:3838/
#
# PUSH IMAGE
# docker push icrsc/icr-8q24
####################################################

FROM rocker/r-ver:4.3.2

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \ 
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    apt-file \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy relevant app files
RUN mkdir www
COPY www ./www
RUN mkdir Data
COPY Data ./Data
COPY global.R LocusExplorerFiles.R server.R ui.R README.md .

RUN Rscript -e "install.packages('rspm') ; rspm::enable() ; install.packages(c('shiny', 'tidyverse', 'plotly', 'shiny', 'shinydashboard', 'shinythemes', 'shinyWidgets','remotes','visNetwork','ggplot2','ggrepel','DT','markdown'))" 
RUN Rscript -e "install.packages('rspm') ; rspm::enable() ; install.packages('igraph')" 
RUN installGithub.r oncogenetics/oncofunco

# Run app (change port accordingly if required)
EXPOSE 3838
CMD Rscript -e "shiny::runApp('/app', port = 3838, host = '0.0.0.0')"

