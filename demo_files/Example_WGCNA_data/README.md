
---
  
  # 🚀 WGCNA Pipeline (Docker)
  
  This repository provides a reproducible WGCNA analysis workflow using Docker.  
The pipeline can be executed either via an R script (`run_wgcna.R`) .

---
  
  ## 🧱 Build the Docker Image
  
  ```bash
nohup docker build -t wgcna-r433 . > docker.log &
  ```

This builds the image in the background and writes logs to `docker.log`.

---
  
  ## 🐳 Start the Container
  
  Mount your working directory into the container:
  
  ```bash
docker run -it -v /path_to_mount:/mnt/work wgcna-r433 bash
```

  ## ▶️ Run the WGCNA Pipeline (R Script Version)
  
  Show help and copy the default config:
  
  ```bash
Rscript run_wgcna.R
```

Run the full pipeline with a config file:
  
  ```bash
Rscript run_wgcna.R config.yml
```

---
  

