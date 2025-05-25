## Running the Dockerized Shiny App
1. Pull the image from GitHub Container Registry:
   ```bash
   docker pull ghcr.io/th206/multiomicserver:latest
   ```
2. Run the container:
   ```bash
   docker run -p 3838:3838 ghcr.io/th206/multiomicserver:latest
