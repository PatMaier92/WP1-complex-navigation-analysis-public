# ############################################################################ #
# ############################################################################ #
#                                                                              #
# ---------------------- WP1 Statistical analysis ---------------------------- #
# Script_02_Analyzer.R                                                         #
# Author: Patrizia Maier                                                       #
#                                                                              #
# ############################################################################ #
# ############################################################################ #


# ------------------------------------------------------------------------------
# ::: LOAD PACKAGES ::: #
# ------------------------------------------------------------------------------

## ---- load_analysis_packages
library(tidyverse)
library(performance)
library(car)
library(janitor)
library(gtsummary)
library(lme4)
library(afex)
library(emmeans)
library(ggsignif)
library(effectsize)
library(corx)
library(ggpol)
library(colorspace)
library(patchwork)
library(papaja)
## ----


# ------------------------------------------------------------------------------
# ::: DATA SETUP ::: #
# ------------------------------------------------------------------------------

# read-in NAVIGATION data 
in_file <- "../WP1_data/WP1_results/wp1_navigation_data.RData"
load(in_file)
rm(in_file)

original_data <- data
data <- data %>% 
  filter(exclude_trial_matlab==0, exclude_low_performer==FALSE) %>% 
  mutate(version=factor(if_else(id %% 2 == 0, "A", "B")))

# check sample size 
length(unique(data$id[data$group=="Delay1H"]))
length(unique(data$id[data$group=="Delay1D"]))
length(unique(data$id[data$group=="Delay2W"]))


# read-in NON-NAVIGATION data 
in_file <- "../WP1_data/WP1_results/WP1_questionnaire_data.RData"
load(in_file)
rm(in_file)

## ---- data_prep
demo_data <- demo_data %>% 
  filter(include=="include") %>% 
  select(id, group, dfb_q01_sex, dfb_q02_age, dfb_q03_eduyrs_total_clean, 
         dfb_q05_language_german, dfb_q06_handiness, dfb_q16_comp_expertise, dfb_q17_comp_freq, sbsds_score, 
         prokrustes_distance, gmda_NLandmarks, gmda_CanAcc)

pre_data <- data %>%
  select(id, group, version, session, block, trial_type, goal_i, memory_score) %>% 
  filter(session==1, trial_type %in% c("testN")) %>%
  mutate(trial_type1=case_when(trial_type=="testN" & block==2 ~ "testN1", T ~ "testN2"),
         trial_type1=factor(trial_type1, levels=c("testN1", "testN2"))) %>% 
  droplevels()

pre_post_data <- data %>%
  select(id, group, version, session, block, trial_type, goal_i, memory_score, correct_object) %>% 
  filter(trial_type %in% c("testA", "testE")) %>% 
  droplevels()

covariates <- data %>% 
  filter(session==1, (trial_type=="testN" & block==4) | trial_type=="testE" | trial_type=="testA") %>% 
  select(id, group, session, trial_type, trial, memory_score, correct_object) %>% 
  group_by(id, trial_type) %>% 
  summarise_at(c("memory_score", "correct_object"), mean, na.rm=T) %>% 
  pivot_wider(names_from=trial_type, values_from=c(memory_score, correct_object)) %>% 
  janitor::remove_empty(which="cols") %>% 
  rename(base_t1_testN=memory_score_testN, base_t1_testE=memory_score_testE, base_t1_testA=memory_score_testA)

main_data <- data %>%
  filter(session==2) %>%
  select(id, group, session, trial, block, goal_i, memory_score, strategy_score_allo, strategy_score_ego,
         time, excess_path_length, dtw_to_testE, initial_rotation_velocity, version) %>%
  left_join(covariates, by=c("id")) %>%
  mutate(block_f=factor(block),
         goal_f=factor(goal_i)) %>% 
  mutate(base_t1_testN=base_t1_testN-mean(base_t1_testN, na.rm=T),
         base_t1_testE=base_t1_testE-mean(base_t1_testE, na.rm=T),
         base_t1_testA=base_t1_testA-mean(base_t1_testA, na.rm=T))

main_data_long <- main_data %>% 
  pivot_longer(cols=c("strategy_score_allo", "strategy_score_ego")) %>% 
  separate(name, into=c("variable", "reference"), sep="_(?=[^_]+$)") %>% 
  pivot_wider(names_from="variable", values_from=value)
## ----

## ---- plot_settings 
# labels 
group_labels <- c("Delay1H"="1h", "Delay1D"="1d", "Delay2W"="2wk")
type_labels <- c("allo"="place", "ego"="response")
session_labels <- c("1"="Learning", "2"="Retrieval", "3"="Retrieval")

# colors
# scales::show_col()
group_colors_c <- c("#fdbf02", "#003399", "#D64457") 
group_colors_f <- lighten(group_colors_c, 0.3)
type_colors_c <- c("#784421", "#216778")
type_colors_f <- lighten(type_colors_c, 0.3) 

# variable labels
l_memory_score <- "memory score [%]"
l_strategy_score <- "strategy score [%]"
l_time <- "latency [s]"
l_excess_path_length <- "excess path length [vu]"
l_dynamic_time_warp_ego <- "norm. trajectory distance [vu]"
l_initial_rotation_velocity <- "initial rotation velocity [rad/s]"

# boxplot wrapper 
afex_boxplot_wrapper <- function(model, xv, tv, pv, ylabel, xlabel=l_session, ymin=0, ymax=1, ybreaks=waiver(), 
                                 tracevis=1, colcol=group_colors_c, fillcol=group_colors_f, boxwidth=0.5, add_points=FALSE) {
  
  if(add_points) {
    p <- afex_plot(model, x=xv, trace=tv, panel=pv, id="id", 
                   error="model", dodge=0.8,
                   mapping=c("shape", "fill", "color"),
                   factor_levels=list(group=group_labels, reference=type_labels),
                   legend_title=NULL, 
                   data_geom=geom_boxjitter,
                   data_arg=list(width=boxwidth, jitter.params=list(seed=100), jitter.size=1, 
                                 outlier.size=1, outlier.intersect=T, show.legend=FALSE),
                   point_arg=list(size=3), 
                   line_arg=list(size=1.25, linetype=tracevis),
                   error_arg=list(size=1.25, width=0)) 
  } 
  else {
    p <- afex_plot(model, x=xv, trace=tv, panel=pv, id="id", 
                   error="model", dodge=0.8,
                   mapping=c("shape", "fill", "color"),
                   factor_levels=list(group=group_labels, session=session_labels),
                   legend_title=NULL, 
                   data_geom=geom_boxplot, 
                   data_arg=list(width=boxwidth, outlier.colour="lightgrey", show.legend=FALSE),
                   point_arg=list(size=3), 
                   line_arg=list(size=1.25, linetype=tracevis),
                   error_arg=list(size=1.25, width=0))
  }
  
  p <- p + 
    scale_fill_manual(values=fillcol) + 
    scale_color_manual(values=colcol) +
    scale_y_continuous(breaks=ybreaks, expand=expansion(mult=c(0, 0.3))) + 
    coord_cartesian(ylim=c(ymin, ymax)) + 
    theme_classic(base_size=14) + 
    theme(legend.position="top", legend.justification=c(0,0),
          strip.background=element_rect(color=NA, fill=NA)) +
    labs(x=xlabel, y=ylabel)
  
  return(p)
}
## ---- 

## ---- analysis_settings
# options("contrasts")
options(contrasts=c(unordered="contr.sum", ordered="contr.poly"))
options(emmeans=list(lmerTest.limit=15000))
options(afex=list(lmerTest.limit=15000))
## ----

## ---- papaja_output_helper
# fix for papaja bug in apa-style emmeans output when using Bonferroni correction
apa_bonferroni_fix <- function(list) {
  list <- list %>% modify_depth(2, str_replace, pattern="\\\\scriptsize ", replacement="")
  
  return(list)
}

# apa-style table for random effects 
apa_random_table <- function(varcor, LRT=NULL) {
  
  # base table
  table <- varcor %>% 
    as.data.frame() %>% 
    mutate(SD=if_else(is.na(var2), sdcor, NaN), 
           r=if_else(!is.na(var2), sdcor, NaN)) %>% 
    mutate_at(vars(SD, r), round, 3) %>% 
    select(-vcov, -sdcor) %>% 
    unite('Random effect', var1:var2, sep=" x ", remove=T, na.rm=T) %>% 
    mutate_at(vars(`Random effect`), str_replace_all, pattern="re1.", replacement="") %>% 
    mutate_at(vars(`Random effect`), str_replace_all, pattern="_", replacement=" ") %>% 
    mutate_at(vars(`grp`), str_replace_all, pattern=".1", replacement="") %>% 
    mutate_at(vars(`grp`), str_replace_all, pattern=".2", replacement="") %>% 
    mutate_at(vars(-SD, -r), str_to_title) %>% 
    rename(`Grouping`=grp) %>% 
    label_variable(SD="$SD$", r="$r$")
  
  # optional: add LRT results 
  if (!is.null(LRT)) {
    table <- table %>%    
      full_join(LRT, by=c("Grouping", "Random effect")) %>% 
      label_variable(p.value="$p$") %>%
      arrange(Grouping)
  }
  
  return(table)
}
## ----


# ############################################################################ #
# ############################################################################ #

# ------------------------------------------------------------------------------
# ::: DEMOGRAPHICS TABLE ::: #
# ------------------------------------------------------------------------------
## ---- demo_table
demo_data_temp <- demo_data %>% 
  select(id, group, dfb_q01_sex, dfb_q02_age, dfb_q03_eduyrs_total_clean,
         dfb_q05_language_german, dfb_q06_handiness, dfb_q16_comp_expertise, dfb_q17_comp_freq, sbsds_score) %>% 
  mutate(dfb_q01_sex=fct_recode(dfb_q01_sex, female="weiblich", male="männlich"),
         dfb_q05_language_german=fct_recode(dfb_q05_language_german, yes="Deutsch ist Muttersprache", no="Deutsch ist nicht Muttersprache"),
         dfb_q06_handiness=fct_recode(dfb_q06_handiness, right="rechtshändig", left="linkshändig", both="beidhändig"),
         dfb_q16_comp_expertise=as.numeric(dfb_q16_comp_expertise),
         dfb_q17_comp_freq=as.numeric(dfb_q17_comp_freq)) %>% 
  droplevels()

demo_table <- demo_data_temp %>% 
  select(-c(id)) %>% 
  tbl_summary(by=group, 
              label=list(dfb_q01_sex ~ "Gender", dfb_q02_age ~ "Age", 
                         dfb_q03_eduyrs_total_clean ~ "Years of education", 
                         dfb_q05_language_german ~ "German native speaker", 
                         dfb_q06_handiness ~ "Handedness", 
                         dfb_q16_comp_expertise ~ "Self-rated computer expertise", 
                         dfb_q17_comp_freq ~ "Self-rated computer use frequency", 
                         sbsds_score ~ "Self-rated spatial abilities (SBSDS)"),
              type=list(dfb_q16_comp_expertise ~ 'continuous', 
                        dfb_q17_comp_freq  ~ 'continuous'),
              statistic=list(all_continuous() ~ "{mean} ({sd})", all_categorical() ~ "{n} ({p}%)"),
              digits=list(all_continuous() ~ c(1, 1)),
              missing="no") %>% 
  add_p(test=list(all_continuous() ~ "aov", all_categorical() ~  "fisher.test")) %>% 
  modify_header(label="")
## ----
rm(demo_data_temp, demo_table)


# ------------------------------------------------------------------------------
# ::: ANALYSIS: BASELINE LEARNING SESSION ::: #
# ------------------------------------------------------------------------------

## ---- base_memory_score
# model
model.base_ms <- mixed(memory_score ~ group + (1 | id) + (1| goal_i),
                       data=pre_data %>% filter(trial_type1=="testN2") , expand_re=T)

omega.base_ms <- omega_squared(model.base_ms, partial=T)
## ----

# # check model
# check_model(model.base_ms$full_model)

# random effects
VarCorr(model.base_ms$full_model)

# fixed effects
model.base_ms
rm(model.base_ms, omega.base_ms)

# ------------------------------------------------------------------------------

## ---- base_memory_score_ego
model.base_ms_ego <- mixed(memory_score ~ group + (1 | id) + (1 | goal_i),
                           data=pre_post_data %>% filter(trial_type=="testE", session==1), expand_re=T)

omega.base_ms_ego <- omega_squared(model.base_ms_ego, partial=T)
## ----
# random effects
VarCorr(model.base_ms_ego$full_model)

# fixed effects
model.base_ms_ego
rm(model.base_ms_ego, omega.base_ms_ego)

# ------------------------------------------------------------------------------

## ---- base_memory_score_allo
model.base_ms_allo <- mixed(memory_score ~ group + (1 | id) + (1 | goal_i),
                            data=pre_post_data %>% filter(trial_type=="testA", session==1), expand_re=T)

omega.base_ms_allo <- omega_squared(model.base_ms_allo, partial=T)
## ----
# random effects
VarCorr(model.base_ms_allo$full_model)

# fixed effects
model.base_ms_allo

rm(model.base_ms_allo, omega.base_ms_allo)


# ## ---- base_correct_object_allo
# model.base_co_allo <- mixed(correct_object ~ group + (1 | id) + (1 | goal_i),
#                             family=binomial(), method="LRT",
#                             data=pre_post_data %>% filter(trial_type=="testA", session==1), expand_re=T)
# ## ----
# # random effects
# VarCorr(model.base_co_allo$full_model)
# 
# # fixed effects
# model.base_co_allo
# 
# rm(model.base_co_allo)


# ------------------------------------------------------------------------------
# ::: ANALYSIS: PRE-POST CHANGE ::: #
# ------------------------------------------------------------------------------

## ---- change_memory_score_ego
model.change_ms_ego <- mixed(memory_score ~ group*session + (session | id) + (1 | goal_i),
                             data=pre_post_data %>% filter(trial_type=="testE"), expand_re=T)

omega.change_ms_ego <- omega_squared(model.change_ms_ego, partial=T)

# plot
p.values <- model.change_ms_ego$anova_table %>% slice(2) %>% pull("Pr(>F)") %>% apa_p(add_equals=T)
plot.change_ms_ego <- afex_boxplot_wrapper(model.change_ms_ego, "session", "group", NULL, l_memory_score, NULL) +
  geom_hline(yintercept=0.5, color="red", linetype="dotted") + 
  geom_signif(textsize=4, xmin=1, xmax=2, y_position=1.12, annotation=c(p.values[1]), color="black", tip_length=0.03) +
  geom_signif(xmin=c(0.75, 1.75), xmax=c(1.25, 2.25), y_position=c(1.1, 1.1), annotation="", color="black", tip_length=0.05)
rm(p.values)
## ----
# # check model
# check_model(model.change_ms_ego$full_model)

# random effects
VarCorr(model.change_ms_ego$full_model)

# fixed effects
model.change_ms_ego

rm(model.change_ms_ego, omega.change_ms_ego, plot.change_ms_ego)


## ---- change_memory_score_allo
model.change_ms_allo <- mixed(memory_score ~ group*session + (session | id) + (1 | goal_i),
                              data=pre_post_data %>% filter(trial_type=="testA"), expand_re=T)

omega.change_ms_allo <- omega_squared(model.change_ms_allo, partial=T)

# plot 
plot.change_ms_allo <- afex_boxplot_wrapper(model.change_ms_allo, "session", "group", NULL, l_memory_score, NULL)
## ---- 
# # check model
# check_model(model.change_ms_allo$full_model)

# random effects
VarCorr(model.change_ms_allo$full_model)

# fixed effects
model.change_ms_allo

## ---- post_change_memory_score_allo
emm <- emmeans(model.change_ms_allo, ~ group*session)
post.change_ms_allo_group_session <- summary(rbind(pairs(emm, simple="group"), pairs(emm, simple="session")), by=NULL, adjust="bonferroni")
rm(emm)

p.values <- post.change_ms_allo_group_session %>% pull(p.value) %>% apa_p(add_equals=T) %>% str_replace("= ", "")
plot.change_ms_allo <- plot.change_ms_allo + 
  geom_hline(yintercept=0.5, color="red", linetype="dotted") + 
  geom_signif(textsize=4, xmin=c(1.25, 1.75, 2), xmax=c(2.25, 2.25, 2.25), y_position=c(1.22, 1.15, 1.05), 
              annotation=c(p.values[9], p.values[5], p.values[6]), color="black", tip_length=0.05)
rm(p.values)
## ----
rm(model.change_ms_allo, omega.change_ms_allo, plot.change_ms_allo, post.change_ms_allo_group_session)


# ## ---- change_correct_object_allo
# model.change_co_allo <- mixed(correct_object ~ group*session + (session | id) + (1 | goal_i),
#                               family=binomial(), method="LRT",
#                               data=pre_post_data %>% filter(trial_type=="testA"), expand_re=T)
# ## ----
# # random effects
# VarCorr(model.change_co_allo$full_model)
# 
# # fixed effects
# model.change_co_allo
# emm <- emmeans(model.change_co_allo, ~ group*session)
# post.change_co_allo_group_session <- summary(rbind(pairs(emm, simple="group"), pairs(emm, simple="session")), by=NULL, adjust="bonferroni")
# rm(emm)
# 
# rm(model.change_co_allo, post.change_co_allo_group_session)


# ------------------------------------------------------------------------------
# ::: ANALYSIS: STRATEGY SCORE DATA RETRIEVAL SESSION ::: #
# ------------------------------------------------------------------------------

## ---- joint_strategy_score
# model
model.joint_ps <- mixed(strategy_score ~ group*reference + base_t1_testN + base_t1_testA + base_t1_testE +
                          (reference | id) + (reference | goal_i) + (reference | block), data=main_data_long, expand_re=T) 

omega.joint_ps <- omega_squared(model.joint_ps, partial=T)

# plot 
plot.joint_ps_avg <- afex_boxplot_wrapper(model.joint_ps, "group", "reference", NULL, l_strategy_score, NULL, colcol=type_colors_c, fillcol=type_colors_f, tracevis=0, boxwidth=0.35)
# plot.joint_ps_avg <- afex_boxplot_wrapper(model.joint_ps, "group", "reference", NULL, l_strategy_score, NULL, colcol=type_colors_c, fillcol=type_colors_f, tracevis=0, add_points=T)
## ---- 
# # check model
# check_model(model.joint_ps$full_model)

# random effects
VarCorr(model.joint_ps$full_model)
# lattice::dotplot(ranef(model.joint_ps$full_model))

# fixed effects
model.joint_ps

## ---- post_joint_strategy_score
emm1 <- emmeans(model.joint_ps, ~ group*reference, lmer.df="satterthwaite")
post.joint_ps_group_reference <- summary(rbind(pairs(emm1, simple="group"), pairs(emm1, simple="reference")), by=NULL, adjust="bonferroni")
post.joint_ps_group_reference_chance <- summary(emm1, null=0.5, adjust="bonferroni", infer=c(T,T)) # test against chance level
rm(emm1)

p.values <- post.joint_ps_group_reference %>% pull(p.value) %>% apa_p(add_equals=T) %>% str_replace("= ", "")
p.valuesc <- post.joint_ps_group_reference_chance %>% pull(p.value) %>% apa_p(add_equals=T) %>% str_replace("= ", "")

plot.joint_ps_avg <- plot.joint_ps_avg + 
  geom_hline(yintercept=0.5, color="red", linetype="dotted") + 
  geom_signif(textsize=3.5, xmin=c(0.8, 1.8), xmax=c(2.8, 2.8), y_position=c(1.23, 1.17), 
              annotation=c(p.values[2], p.values[3]), color="black", tip_length=0.05) + 
  geom_signif(textsize=3.5, xmin=c(0.8, 1.8), xmax=c(1.2, 2.2), y_position=c(1.08, 1.08), 
              annotation=c(p.values[7], p.values[8]), color="black", tip_length=0.05) +
  annotate("segment", x=c(0.625, 1.625, 2.625, 3.425), xend=c(0.625, 1.625, 2.625, 3.425), 
           y=c(0.5, 0.5, 0.5, 0.5), yend=c(0.8, 0.8, 0.8, 0.8), colour="black") + 
  annotate("segment", x=c(0.625, 1.625, 2.625, 3.395), xend=c(0.625, 1.625, 2.625, 3.395)+0.03, 
           y=c(0.5, 0.5, 0.5, 0.5), yend=c(0.5, 0.5, 0.5, 0.5), colour="black") + 
  annotate("segment", x=c(0.625, 1.625, 2.625, 3.395), xend=c(0.625, 1.625, 2.625, 3.395)+0.03, 
           y=c(0.8, 0.8, 0.8, 0.8), yend=c(0.8, 0.8, 0.8, 0.8), colour="black") + 
  annotate("text", angle=90, size=3.5, x=c(0.625, 1.625, 2.625, 3.575)-0.1, y=c(0.6, 0.6, 0.6, 0.6), 
           label=c(p.valuesc[1], p.valuesc[2], p.valuesc[3], p.valuesc[6]), colour="black")

rm(p.values, p.valuesc)
## ----
rm(model.joint_ps, omega.joint_ps, plot.joint_ps_avg, post.joint_ps_group_reference, post.joint_ps_group_reference_chance)


# ## ---- joint_strategy_score_block
# # model
# model.joint_ps <- mixed(strategy_score ~ group*block_f*reference + base_t1_testN + base_t1_testA + base_t1_testE +
#                           (reference | id) + (reference | goal_i), data=main_data_long, expand_re=T)
# 
# # plot
# plot.joint_ps <- afex_boxplot_wrapper(model.joint_ps, "block_f", "reference", "group", l_strategy_score, "block", colcol=type_colors_c, fillcol=type_colors_f)
# plot.joint_ps_avg <- afex_boxplot_wrapper(model.joint_ps, "reference", NULL, "group", l_strategy_score, "block", colcol=type_colors_c, fillcol=type_colors_f)
# ## ----
# # # check model
# # check_model(model.joint_ps$full_model)
# 
# # random effects
# VarCorr(model.joint_ps$full_model)
# # lattice::dotplot(ranef(model.joint_ps$full_model))
# 
# # fixed effects
# model.joint_ps
# 
# ## ---- post_joint_strategy_score_block
# # 3-way interaction
# # option 1
# emm1 <- emmeans(model.joint_ps, ~ group*reference*block_f, lmer.df="satterthwaite")
# post.joint_ps_group_reference_block_change <- summary(rbind(pairs(emm1, interaction=c("consec"), by=c("reference", "group"), exclude=c("2", "3", "4"))), adjust="bonferroni")
# e1 <- pairs(emm1, interaction=c("pairwise", "consec"), by="reference", exclude=c("2", "3", "4"))
# e2 <- pairs(emm1, interaction=c("pairwise", "consec"), by="group", exclude=c("2", "3", "4"))
# post.joint_ps_group_reference_block <- summary(rbind(e1, e2), by=NULL, adjust="bonferroni")
# rm(e1, e2)
# 
# # option 2
# post.joint_ps_group_reference_block_change2 <- summary(rbind(pairs(emm1, interaction=c("poly"), max.degree=1, by=c("reference", "group"))), adjust="bonferroni")
# e1 <- pairs(emm1, interaction=c("pairwise", "poly"), max.degree=1, by="reference")
# e2 <- pairs(emm1, interaction=c("pairwise", "poly"), max.degree=1, by="group")
# post.joint_ps_group_reference_block2 <- summary(rbind(e1, e2), by=NULL, adjust="bonferroni")
# rm(e1, e2, emm1)
# 
# # 2-way interaction
# emm2 <- emmeans(model.joint_ps, ~ reference*block_f, lmer.df="satterthwaite")
# post.joint_ps_reference_block_change <- summary(rbind(pairs(emm2, interaction=c("poly"), max.degree=1, by="reference")), adjust="bonferroni")
# post.joint_ps_reference_block <- pairs(emm2, interaction=c("pairwise", "poly"), max.degree=1)
# rm(emm2)
# 
# # 2-way interaction
# emm3 <- emmeans(model.joint_ps, ~ group*reference, lmer.df="satterthwaite")
# post.joint_ps_group_reference <- summary(rbind(pairs(emm3, simple="group"), pairs(emm3, simple="reference")), by=NULL, adjust="bonferroni")
# post.joint_ps_group_reference_chance <- summary(emm3, null=0.5, adjust="bonferroni", infer=c(T,T)) # test against chance level
# rm(emm3)
# ## ----
# rm(model.joint_ps, plot.joint_ps, plot.joint_ps_avg, post.joint_ps_group_reference_block_change, post.joint_ps_group_reference_block_change2,
#    post.joint_ps_group_reference_block, post.joint_ps_group_reference_block2, post.joint_ps_reference_block_change, post.joint_ps_reference_block,
#    post.joint_ps_group_reference, post.joint_ps_group_reference_chance)


# ------------------------------------------------------------------------------
# ::: ANALYSIS: NAVIGATION DATA IN RETRIEVAL SESSION ::: #
# ------------------------------------------------------------------------------

## ---- dynamic_time_warping
# model
model.dtw_to_testE <- mixed(dtw_to_testE ~ group + (1 | id) + (1 | goal_i) + (1 | block), 
                       data=main_data, expand_re=T) 

omega.dtw_to_testE <- omega_squared(model.dtw_to_testE, partial=T)

# plot 
plot.dtw_to_testE <- afex_boxplot_wrapper(model.dtw_to_testE, "group", NULL, NULL, l_dynamic_time_warp_ego, NULL, ymin=0, ymax=3.5, boxwidth=0.35)
## ---- 
# # check model
# check_model(model.dtw_to_testE$full_model)

# random effects
VarCorr(model.dtw_to_testE$full_model)
# lattice::dotplot(ranef(model.dtw_to_testE$full_model))

# fixed effects
model.dtw_to_testE

## ---- post_dynamic_time_warping
post.dtw_to_testE_group <- emmeans(model.dtw_to_testE, pairwise ~ group, lmer.df="satterthwaite", adjust="bonferroni")$contrasts

p.values <- post.dtw_to_testE_group %>% as.data.frame %>% pull(p.value) %>% apa_p(add_equals=T) %>% str_replace("= ", "")
plot.dtw_to_testE <- plot.dtw_to_testE + 
  geom_signif(textsize=4, xmin=c(1, 2), xmax=c(3, 3), y_position=c(4, 3.5), 
              annotation=c(p.values[2], p.values[3]), color="black", tip_length=0.05)
rm(p.values)
## ----
rm(model.dtw_to_testE, omega.dtw_to_testE, plot.dtw_to_testE, post.dtw_to_testE_group)

# cor.test(main_data$strategy_score_allo, main_data$dtw_to_testE)
# cor.test(main_data$strategy_score_ego, main_data$dtw_to_testE)

# ------------------------------------------------------------------------------

## ---- initial_rotation_velocity 
# model
model.initial_rotation_velocity <- mixed(initial_rotation_velocity ~ group + (1 | id) + (1 | goal_i) + (1 | block), 
                                data=main_data, expand_re=T) 

omega.initial_rotation_velocity <- omega_squared(model.initial_rotation_velocity, partial=T)

# plot 
plot.initial_rotation_velocity <- afex_boxplot_wrapper(model.initial_rotation_velocity, "group", NULL, NULL, l_initial_rotation_velocity, NULL, ymin=0, ymax=4, boxwidth=0.35)
## ---- 
# # check model
# check_model(model.initial_rotation_velocity$full_model)

# random effects
VarCorr(model.initial_rotation_velocity$full_model)
# lattice::dotplot(ranef(model.initial_rotation_velocity$full_model))

# fixed effects
model.initial_rotation_velocity

## ---- post_initial_rotation_velocity
post.initial_rotation_velocity_group <- emmeans(model.initial_rotation_velocity, pairwise ~ group, lmer.df="satterthwaite", adjust="bonferroni")$contrasts

p.values <- post.initial_rotation_velocity_group %>% as.data.frame %>% pull(p.value) %>% apa_p(add_equals=T) %>% str_replace("= ", "")
plot.initial_rotation_velocity <- plot.initial_rotation_velocity + 
  geom_signif(textsize=4, xmin=c(1), xmax=c(3), y_position=c(4.5), 
              annotation=c(p.values[2]), color="black", tip_length=0.05)
rm(p.values)
## ----
rm(model.initial_rotation_velocity, omega.initial_rotation_velocity, plot.initial_rotation_velocity, post.initial_rotation_velocity_group)

# ------------------------------------------------------------------------------

## ---- time
# model
model.time <- mixed(time ~ group + (1 | id) + (1 | goal_i) + (1 | block), 
                    data=main_data, expand_re=T) 

omega.time <- omega_squared(model.time, partial=T)

# plot 
plot.time <- afex_boxplot_wrapper(model.time, "group", NULL, NULL, l_time, NULL, ymin=0, ymax=30)
## ---- 
rm(model.time, omega.time, plot.time)

# ------------------------------------------------------------------------------

## ---- excess_path_length
# model
model.path <- mixed(excess_path_length ~ group + (1 | id) + (1 | goal_i) + (1 | block), 
                    data=main_data, expand_re=T) 

omega.path <- omega_squared(model.path, partial=T)

# plot 
plot.path <- afex_boxplot_wrapper(model.path, "group", NULL, NULL, l_excess_path_length, NULL, ymin=0, ymax=0.5)
## ---- 
rm(model.path, omega.path, plot.path)


# ############################################################################ #
# ############################################################################ #

# ------------------------------------------------------------------------------
# ::: ANALYSIS: CORRELATIONS ::: #
# ------------------------------------------------------------------------------

## ---- correlations
demo_data_temp <- demo_data %>% 
  select(id, group, dfb_q01_sex, dfb_q02_age, dfb_q03_eduyrs_total_clean, sbsds_score) %>% 
  rename(sex=dfb_q01_sex, age=dfb_q02_age, education=dfb_q03_eduyrs_total_clean, SBSOD=sbsds_score) %>% 
  mutate(sex=as.numeric(sex))

post_data_temp <- pre_post_data %>% 
  filter(session==3) %>% 
  select(-correct_object, -version, -session, -block, -group) %>% 
  group_by(id, trial_type) %>% 
  summarise_at(c("memory_score"), mean, na.rm=T) %>% 
  pivot_wider(names_from=trial_type, values_from=memory_score) %>% 
  relocate("testE", .after="testA")

corr_data <- main_data %>%
  select(id, strategy_score_allo, strategy_score_ego, initial_rotation_velocity, dtw_to_testE, time, excess_path_length) %>% 
  group_by(id) %>% 
  summarise_all(mean, na.rm=T) %>% 
  ungroup() %>% 
  left_join(post_data_temp) %>%
  left_join(covariates) %>% 
  left_join(demo_data_temp) %>%
  select(-correct_object_testA, -id, -group) %>% 
  drop_na() %>% 
  relocate(any_of(c("testA", "testE")), .after="dtw_to_testE") %>% 
  relocate(any_of(c("base_t1_testE")), .after="base_t1_testA") %>% 
  relocate(any_of(c("base_t1_testN", "base_t1_testA", "base_t1_testE")), .after="testE") %>% 
  rename("Place Strategy Score"="strategy_score_allo", "Response Strategy Score"="strategy_score_ego",
         "Place Memory Score"="testA", "Response Memory Score"="testE", 
         "Sex"="sex", "Age"="age", "Years of Education"="education",
         "Baseline Memory Score"="base_t1_testN", "Baseline Place Memory Score"="base_t1_testA", "Baseline Response Memory Score"="base_t1_testE",
         "Initial Rotation Velocity"="initial_rotation_velocity", "Norm. Trajectory Distance"="dtw_to_testE", 
         "Latency"="time", "Excess Path Length"="excess_path_length") %>% 
  as.data.frame()

n_bonferroni=69
table.correlation <- corx(corr_data, method="spearman", stars=c(0.05/n_bonferroni, 0.01/n_bonferroni, 0.001/n_bonferroni), p_adjust="none",
                   triangle="lower", describe=c(`$M$`=mean, `$SD$`=sd))
table.correlation_adj <- cbind(table.correlation$apa[,1:6], table.correlation$apa[,15:16])
## ----
rm(demo_data_temp, corr_data, table.correlation, table.correlation_adj)


# ############################################################################ #
# ############################################################################ #

# ------------------------------------------------------------------------------
# ::: SUPPLEMENT ANALYSIS: COGNITIVE MAP DRAWINGS ::: #
# ------------------------------------------------------------------------------

## ---- boundary
model.boundary <- aov_ez("id", "prokrustes_distance", demo_data, between=c("group")) 

# plot.boundary <-- afex_boxplot_wrapper(model.boundary, "group", NULL, NULL, "prokrustes distance [vu]", NULL, ymin=0, ymax=0.25)
## ---- 
rm(model.boundary, plot.boundary)


## ---- number_landmarks
model.landmarks <- aov_ez("id", "gmda_NLandmarks", demo_data, between=c("group")) 

# plot.landmarks <- afex_boxplot_wrapper(model.landmarks, "group", NULL, NULL, "landmark identity accuracy [%]", NULL)
## ---- 
rm(model.landmarks, plot.landmarks)


## ---- position_landmarks
model.positioning <- aov_ez("id", "gmda_CanAcc", demo_data, between=c("group")) 

# plot.positioning <- afex_boxplot_wrapper(model.positioning, "group", NULL, NULL, "landmark canonical accuracy [%]", NULL)
## ----
rm(model.positioning, plot.positioning)


# ############################################################################ #
# ############################################################################ #

rm(list = ls())

# ############################################################################ #
# ############################################################################ #