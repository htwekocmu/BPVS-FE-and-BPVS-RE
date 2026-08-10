# Helper: save both BPVS FE and RE fitted models to inst/extdata
#
# Edit and run this script locally after you've trained your models in R.
# Replace `fitted_fe` and `fitted_re` with your model object names.

if (!dir.exists("inst/extdata")) dir.create("inst/extdata", recursive = TRUE)

# Example placeholders (replace with your trained models):
# fitted_fe <- <your trained BPVS-FE model object>
# fitted_re <- <your trained BPVS-RE model object>

# Save to the package extdata location
# saveRDS(fitted_fe, file = "inst/extdata/bpvs_fe_model.rds")
# saveRDS(fitted_re, file = "inst/extdata/bpvs_re_model.rds")

message("Edit this file to save your real BPVS-FE and BPVS-RE fitted models to inst/extdata/")
