# Load manually
load("~/github/sablefishData/data/obs_released_3d_mon_10a_pool_sml.rda")
load("~/github/sablefishData/data/obs_recovered_5d_mon_10a_pool_sml.rda")
load("~/github/sablefishData/data/obs_capture_rate_2d_mon_10a.rda")
load("~/github/sablefishData/data/obs_report_ratio_2d_mon_10a.rda")
load("~/github/sablefishData/data/obs_template_2d_mon_10a.rda")

# Load fit
load("~/github/sablefishFit/data/fit_obs_mon_10a_size_sml_static_size_bias_nb1.rda")
fit_list <- fit_obs_mon_10a_size_sml_static_size_bias_nb1

# Define arguments for testing
released_3d <- obs_released_3d_mon_10a_pool_sml
recovered_5d <- obs_recovered_5d_mon_10a_pool_sml
capture_rate_2d <- obs_capture_rate_2d_mon_10a
report_ratio_2d <- obs_report_ratio_2d_mon_10a
template_2d <- obs_template_2d_mon_10a
tag_loss_rate <- 0.02 / 12 # Divide by 12
imm_loss_ratio <- 0.1
recapture_delay <- 1
error_family <- 1
result_units <- 12 # Check
time_process <- 0
time_pattern <- 0
pattern_size <- 0
newton_steps <- 0
nlminb_loops <- 5
openmp_cores <- floor(parallel::detectCores() / 2)
capture_map_2d <- array(rep(1, 10), dim = c(10, 1))
structure_list <- NULL
parameter_list <- NULL
optimizer_list <- NULL
nlminb_control <- mmmTMBcontrol()

# Load Meghan obs data for testing
#------------------------------------------------------------------------------#
# Load manually
load("~/github/meghanData/data/obs_released_3d_mon_12a_pool_sml.rda")
load("~/github/meghanData/data/obs_recovered_5d_mon_12a_pool_sml.rda")
load("~/github/meghanData/data/obs_released_3d_mon_12a_size_sml.rda")
load("~/github/meghanData/data/obs_recovered_5d_mon_12a_size_sml.rda")
load("~/github/meghanData/data/obs_capture_rate_2d_mon_12a.rda")
load("~/github/meghanData/data/obs_report_ratio_2d_mon_12a.rda")
load("~/github/meghanData/data/obs_template_2d_mon_12a.rda")

# Pool sml
released_3d <- obs_released_3d_mon_12a_pool_sml
recovered_5d <- obs_recovered_5d_mon_12a_pool_sml
capture_map_2d <- array(rep(1, 12), dim = c(12, 1))

# Size sml
released_3d <- obs_released_3d_mon_12a_size_sml
recovered_5d <- obs_recovered_5d_mon_12a_size_sml
capture_map_2d <- array(rep(1:3, each = 12), dim = c(12, 3))

# General
capture_rate_2d <- obs_capture_rate_2d_mon_12a
report_ratio_2d <- obs_report_ratio_2d_mon_12a
template_2d <- obs_template_2d_mon_12a
tag_loss_rate <- 0.02 / 12 # Divide by 12
imm_loss_ratio <- 0.1
recapture_delay <- 1
error_family <- 1
result_units <- 12 # Check
time_process <- 0
time_pattern <- 0
pattern_size <- 0
newton_steps <- 0
nlminb_loops <- 5
openmp_cores <- floor(parallel::detectCores() / 2)
structure_list <- NULL
parameter_list <- NULL
optimizer_list <- NULL
nlminb_control <- mmmTMBcontrol()

