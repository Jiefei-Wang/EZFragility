# Script executed on Mac Server
library(glmnet)

pipeline_xls <- read.csv("/Volumes/bigbrain/Multipatient/patient_data_all.csv")[1, ]
subject_code <- pipeline_xls$subject_code
project <- pipeline_xls$project_name
electrodes <- dipsaus::parse_svec(pipeline_xls$load_electrodes)
display <- electrodes
epoch_name <- pipeline_xls$epoch_file_name
reference_name <- pipeline_xls$reference_name
condition <- pipeline_xls$condition
soz <- dipsaus::parse_svec(pipeline_xls$SOZ)
sozc <- electrodes[!(electrodes %in% soz)]

fl <- \(x) x[x > 30 & x < 40]
electrodes <- fl(electrodes)
display <- fl(display)
soz <- fl(soz)
sozc <- fl(sozc)

subject <- raveio::as_rave_subject(paste0(project, "/", subject_code))
trial_num <- match(condition, subject$get_epoch(epoch_name)$table$Condition)
fragility_pipeline <- raveio::pipeline("karaslab_fragility", paths = "../rave_pipeline/modules/")
source("../rave_pipeline/modules/karaslab_fragility/R/shared-functions.R")

fragility_pipeline$set_settings(
  project_name = project,
  subject_code = subject_code,
  epoch_name = epoch_name,
  epoch_time_window = c(-10, 20),
  reference_name = reference_name,
  load_electrodes = electrodes,
  display_electrodes = display,
  condition = condition,
  trial_num = trial_num,
  t_window = 100,
  t_step = 50,
  sz_onset = 0,
  lambda = FALSE,
  threshold_start = 0,
  threshold_end = 20,
  threshold = 0.5,
  soz = soz,
  sozc = sozc
)

results <- c(fragility_pipeline$run(c("repository", "adj_frag_info")))

repository <- results$repository
t_window <- fragility_pipeline$get_settings("t_window")
t_step <- fragility_pipeline$get_settings("t_step")
soz <- fragility_pipeline$get_settings("soz")
sozc <- fragility_pipeline$get_settings("sozc")
f <- results$adj_frag_info$frag
adj_frag_info <- results$adj_frag_info

save(repository, t_window, t_step, soz, sozc, f, adj_frag_info, file = "pt1.rda")

# pt1.rda copied to local machine package "data-raw" directory
#
# Then...
#
# load("data-raw/pt1.rda")
# use_data(repository, t_window, t_step, soz, sozc, f, adj_frag_info,
#          internal = TRUE,
#          overwrite = TRUE)

