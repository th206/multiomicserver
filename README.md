# Multi-omics Shiny App
A Shiny app for visualizing multi-omics networks.

## Running the server locally in RStudio:
All files required to run the app locally using RStudio are in this GitHub repository. Start-up rstudio and then launch app.R as a shiny app.

Note: you need to have the following libraries installed in RStudio:

```bash
install.packages(c("shiny", "shinydashboard", "readxl", "dplyr", "tidyr", "visNetwork", "DT", "shinyWidgets", "igraph","rlang"))
```
## Usage
Run:
```bash
shiny:runApp()
```
## Running via a web server
The app is available following this link: http://128.84.40.245 

## Running the Dockerized Shiny App
1. Pull the image from GitHub Container Registry:
   ```bash
   docker pull ghcr.io/th206/multiomicserver:latest
   ```
2. Run the container:
   ```bash
   docker run -p 3838:3838 ghcr.io/th206/multiomicserver:latest
   ```
   Navigate to your browser and open the following page: http://localhost:3838

Here is a screenshot of the server:
![My Image](www/about.png)
