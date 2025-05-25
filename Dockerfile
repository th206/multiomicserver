# Base image with R and Shiny Server
FROM rocker/shiny:4.2.3

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && apt-get clean

# Install R packages
RUN R -e "install.packages(c('shiny', 'shinydashboard', 'readxl', 'dplyr', 'tidyr', 'visNetwork', 'igraph', 'DT', 'shinyWidgets', 'markdown', 'rlang'), repos='https://cloud.r-project.org/')"

# Copy app to Shiny server
COPY . /srv/shiny-server/
RUN chown -R shiny:shiny /srv/shiny-server

# Expose Shiny port
EXPOSE 3838

# Run the app
CMD ["/usr/bin/shiny-server"]
