# ############################################################################ #
# ############################################################################ #
#                                                                              #
# ---------------------- WP1 Data wrangling ---------------------------------- #
# Script_01_XLSX_to_RData.R                                                    #
# Author: Patrizia Maier                                                       #
#                                                                              #
# ############################################################################ #
# ############################################################################ #


# ------------------------------------------------------------------------------
# ::: LOAD PACKAGES ::: #
# ------------------------------------------------------------------------------

library(readxl)
library(foreign)
library(tidyverse)


# ------------------------------------------------------------------------------
# ::: STARMAZE NAVIGATION DATA ::: #
# ------------------------------------------------------------------------------

# read-in data 
date <- readline(prompt = "Please enter the date string of the result file ")
file_path <- paste("../WP1_data/WP1_results/wp1_navigation_data_", date, ".xlsx", sep="")
data <- read_xlsx(file_path, col_names = T, na = "999")
rm(date, file_path)

# tidy data 
data <- data %>% 
  mutate_at(c("session", "goal_vis", "arrow_vis", "goal_identity"), factor) %>% 
  mutate(group=factor(group, levels=c("Delay1H", "Delay1D", "Delay2W")),
         trial_type=factor(trial_type, levels=c("learn", "testN", "testE", "testA"))) %>% 
  mutate(dtw_to_testE=dtw_to_testE*10,
         initial_rotation_velocity=initial_rotation_velocity*10)

# mark low-performing subjects according to fixed criterion
temp_performance <- data %>%
  filter(session==1, trial_type=="testN", block==4) %>% 
  group_by(id) %>%
  summarise_at(c("correct_alley"), mean, na.rm=T)

crit <- 0.625
low_performer <- temp_performance$id[which(temp_performance$correct_alley<crit)]

data <- data %>% 
  mutate(exclude_low_performer = id %in% low_performer)


# save as RData
file_name <- "../WP1_data/WP1_results/wp1_navigation_data.RData"
save(data, file=file_name)
rm(temp_performance, low_performer, crit, file_name)


# ------------------------------------------------------------------------------
# ::: NON-NAVIGATION DATA ::: #
# ------------------------------------------------------------------------------

# read-in data 
in_file <- "../WP1_data/WP1_results/wp1_questionnaire_data.sav"
questionnaire_data <- read.spss(in_file, use.value.labels=T, to.data.frame=T)
rm(in_file)

file_path <- "../WP1_data/WP1_results/wp1_results_map_scoring.txt"
geometry_data <- read.delim(file_path, header = T, sep = ",")
rm(file_path)

file_path <- "../WP1_data/WP1_results/wp1_results_lm_scoring.txt"
landmark_data <- read.delim(file_path, header = T, sep = ",") %>% 
  filter(gmda_measure %in% c("SCanOrg", "CanAcc", "NLandmarks")) %>% 
  pivot_wider(names_from=gmda_measure, values_from=score, names_prefix="gmda_")
rm(file_path)

# combine data
demo_data <- questionnaire_data %>%
  left_join(geometry_data,) %>% 
  left_join(landmark_data)
rm(questionnaire_data, geometry_data, landmark_data)

# update exclusion marker
participants <- data %>% filter(exclude_low_performer==FALSE) %>% pull(id) %>% unique()
demo_data <- demo_data %>%
  mutate(include = ifelse(id %in% participants, "include", "exclude"))
rm(participants)

# save as RData
file_name <- "../WP1_data/WP1_results/wp1_questionnaire_data.Rdata"
save(demo_data, file=file_name)


# ------------------------------------------------------------------------------

# clear workspace
rm(list = ls())
