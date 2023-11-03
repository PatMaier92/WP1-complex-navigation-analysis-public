# ############################################################################ #
# ############################################################################ #
#                                                                              #
# ---------------------- WP1 Landmark scoring data --------------------------- #
# Script_01_GMDA_to_txt.R                                                      #
# Author: Patrizia Maier                                                       #
#                                                                              #
# ############################################################################ #
# ############################################################################ #


# ------------------------------------------------------------------------------
# ::: LOAD PACKAGES ::: #
# ------------------------------------------------------------------------------

library(tidyverse)


# ------------------------------------------------------------------------------
# ::: LOAD GMDA RESULT FILES ::: #
# ------------------------------------------------------------------------------

#  get file list
file_list <- paste("Data", list.files(path="Data", pattern=".*template_9.*Summary\\.csv$"), sep="/")

# preprocess gmda data 
gmda_data <- file_list %>%
  purrr::map_df(read_csv,
                col_names=c("Measure Type", "Filename", "Measure", "Score", "Score_2"), 
                col_type=cols(`Measure Type`=col_character(),
                              Filename=col_character(),
                              Measure=col_character(),
                              Score=col_double(),
                              Score_2=col_double()),
                skip=9, 
                n_max=8) %>% 
  # correct for delimiter error in raw data 
  unite(Score, Score_2, col="Score", sep=".", na.rm=T) %>%
  rename("id"="Filename")

# preprocess bdr data 
bdr_data <- file_list %>%
  purrr::map_df(read_csv, 
                col_names=c("Measure Type", "Filename", "Measure", "Score"), 
                col_type=cols(`Measure Type`=col_character(),
                              Filename=col_character(),
                              Measure=col_character(),
                              Score=col_double()),
                skip=21, 
                n_max=10) %>% 
  select(!c(`Measure Type`)) %>% 
  filter(Measure %in% c("theta")) %>% 
  rename("id"="Filename")

# # check for "bad" theta values
# threshold <- as.numeric(readline("Enter theta threshold: "))
# bad_theta <- bdr_data %>% filter(Measure=="theta") %>% filter(abs(Score) > threshold) 
# bad_theta
# rm(bad_theta)


# ------------------------------------------------------------------------------
# ::: SELECT AND COMBINE DATA ::: #
# ------------------------------------------------------------------------------

# select and combine 
data <- gmda_data %>% 
  select(!c(`Measure Type`)) %>% 
  filter(!Measure %in% c("Canonical Organization", "Rotational Bias", "Scaling Bias")) %>% 
  mutate(Score=as.numeric(Score),
         Measure=case_when(Measure=="Num Landmarks Missing" ~ "NLandmarks",
                           Measure=="SQRT(Canonical Organization)" ~ "SCanOrg", 
                           Measure=="Canonical Accuracy" ~ "CanAcc",
                           Measure=="Distance Accuracy" ~ "DistAcc",
                           Measure=="Angle Accuracy" ~ "AngleAcc",
                           TRUE ~ "NA")) %>% 
  full_join(bdr_data, by=c("id", "Measure", "Score")) %>% 
  mutate(id=as.numeric(id),
         Score=case_when(Measure=="NLandmarks" ~ (1-(Score/9)), T ~ Score)) %>%
  arrange(id) %>% 
  rename(gmda_measure=Measure, score=Score)
rm(bdr_data, gmda_data)

# save as .txt  
out_file <-  "../../WP1_data/WP1_results/wp1_results_lm_scoring.txt"
write.table(data, file=out_file, sep=",", row.names=F, quote=F)


# ------------------------------------------------------------------------------

# clear workspace
rm(list = ls())