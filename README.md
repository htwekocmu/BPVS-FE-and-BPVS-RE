# BPVS-FE-and-BPVS-RE — quick install & use

This repository now includes a minimal R package-like layout and helper functions so others can load and use the BPVS-FE and BPVS-RE models easily.

Quick install (as a package) from GitHub (branch add-r-package-and-shiny):

```r
# install remotes if you don't have it
install.packages("remotes")
remotes::install_github("htwekocmu/BPVS-FE-and-BPVS-RE", ref = "add-r-package-and-shiny")
library(BPVS)

# load FE model and predict
model_fe <- BPVS::load_bpvs_model("FE")
res_fe <- BPVS::predict_bpvs(newdata = your_data_frame, type = "FE", model = model_fe)
# res_fe$predictions and res_fe$selected

# load RE model and predict
model_re <- BPVS::load_bpvs_model("RE")
res_re <- BPVS::predict_bpvs(newdata = your_data_frame, type = "RE", model = model_re)

# Backwards-compatible wrappers:
# model <- BPVS::load_model() # defaults to FE
# res <- BPVS::predict_model(newdata = your_data_frame)
```

Adding your models to the repo
- Train your models locally and run data-raw/save_bpvs_models.R (edit it) to write the .rds files to inst/extdata/.
- Commit inst/extdata/bpvs_fe_model.rds and inst/extdata/bpvs_re_model.rds to the repository. If files are large, consider Git LFS.

Shiny app
- A simple Shiny app has been added at shiny/app.R. It allows users to pick FE vs RE, upload a CSV of features, and download predictions + selected covariates.

Notes
- If your models require preprocessing (scales, factor mappings), save the preprocessing object together with the model or implement preprocessing inside predict_bpvs().
