# This script was a modification of Gina Lewin's 2023 and 2024 PNAS Papers
#script calculates AS2 to compare a model to the human transcriptome and then graphs output using sunburst plots
# The script also outputs the AS2 and Z scores for each run that can be used to make any other type of plot
#for this paper, we used heatmaps to represent the AS2 scores 

#   1. Calculates AS2 or zscores for all samples in a model and plots sunburst graph
#   2. Calculates AS2 and z-scores after subsampling to keep input numbers consistent across different conditions


#set working directory
#install.packages("tidyverse")
library(dplyr)
library(tidyverse)

#install.packages("cowplot")
library(cowplot)

#install.packages("zeallot")
library(zeallot)

#install.packages("devtools")
library(devtools)

#install.packages("reshape2")
library(reshape2)

#install_github("didacs/ggsunburst")
library(ggsunburst)
library(readxl)

#install.packages("scico")
library(scico)

library(stringr)
library(purrr)
library(tibble)

score_target_vs_model_SA <- function(counts_DF, Target_samplenames, Model_samplenames) # given a counts DF, lists of target samples, and lists of model samples, will return a score for each gene
{
  
  all_tested_samples_list <- names(counts_DF)[-1]
  names(counts_DF)
  
  conditions_list <- all_tested_samples_list
  dim(counts_DF %>% select(COGID, all_tested_samples_list) %>% tibble::column_to_rownames("COGID"))
  
  df1 <- counts_DF
  setdiff(Target_samplenames, names(df1))
  DF_target <- df1  %>% select(COGID, Target_samplenames) %>% tidyr::gather(key=sample_name, value=expression, -COGID) %>% group_by(COGID) %>% dplyr::summarize(target_mean = mean(expression), target_SD = sd(expression))
  
  df_allzscores <- DF_target %>% inner_join(df1 %>% select(COGID, Model_samplenames)) %>% mutate_at(.vars=vars(-COGID, -target_mean, -target_SD), .funs=funs((.-target_mean)/target_SD))
  mean_modelZscoreDF <- df_allzscores %>% select(COGID, Model_samplenames)  %>% transmute(COGID, penalty_temp = purrr::pmap_dbl(.[c(-1)], function(...)  mean(c(...)))) %>% mutate(penalty = round(penalty_temp, digits = 4)) %>% select(COGID, penalty) 
  
  # below should get rid of NAs
  
  if(length(which(is.na(mean_modelZscoreDF$penalty)))>0)
  {
    mean_modelZscoreDF[which(is.na(mean_modelZscoreDF$penalty)),]$penalty <- 0
  }  
  final_DF <- mean_modelZscoreDF  
  
  return(final_DF)
}


#for dummy model
score_target_vs_model <- function(Target_samplenames, Model_samplenames)
{
  DF_target <- counts_normalized  %>% select(COGID, Target_samplenames) %>% gather(key=sample_name, value=expression, -COGID) %>% group_by(COGID) %>% dplyr::summarize(target_mean = mean(expression), target_SD = sd(expression))
  
  df_allzscores <- DF_target %>% inner_join(counts_normalized %>% select(COGID, Model_samplenames)) %>% mutate_at(.vars=vars(-COGID, -target_mean, -target_SD), .funs=funs((.-target_mean)/target_SD))
  
  mean_modelZscoreDF <- df_allzscores %>% select(COGID, Model_samplenames)  %>% transmute(COGID, penalty_temp = pmap_dbl(.[c(-1)], function(...)  mean(c(...)))) %>% mutate(penalty = round(penalty_temp, digits = 4)) %>% select(COGID, penalty) 
  
  # below should get rid of NAs
  
  if(length(which(is.na(mean_modelZscoreDF$penalty)))>0)
  {
    mean_modelZscoreDF[which(is.na(mean_modelZscoreDF$penalty)),]$penalty <- 0
  }  
  final_DF <- mean_modelZscoreDF  
  
  return(final_DF)
}

model_selfvalidation_leaveout <- function(sample_DF1, left_out_number)
{
  test_wound <- sample(1:length(human_list), left_out_number, replace = FALSE)
  train_wound <- setdiff(1:length(human_list), test_wound)
  
  dummy1 <- score_target_vs_model(human_list[train_wound], human_list[test_wound])
  return(dummy1)
}



#draw sunburst with color brewer pallette-----------------
library(RColorBrewer)
draw_sunburst_STANDARD_SA <- function(SA_df_counts, COG_SA, sb, target_list, model_list, maxlim, node_labs)
{
  
  COG_SA <- COG_SA  %>% mutate(Main = if_else(Main == "Other", paste0(Main, "@", Meta), Main))
  SA_DE <- score_target_vs_model_SA(SA_df_counts, target_list, model_list)
  
  Main_penalty <- COG_SA %>% left_join(SA_DE) %>% group_by(Main) %>% select(COGID, Main, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  Meta_penalty <- COG_SA %>% left_join(SA_DE) %>% group_by(Meta) %>% select(COGID, Meta, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  Main_penalty %>% arrange(penalty)
  sum(abs(SA_DE$penalty) < 2)
  COG_SA %>% left_join(SA_DE)  %>% arrange(abs(penalty)) #%>% filter(Sub == "Sulfur metabolism")
  
  Main_penalty %>% filter()
  
  names(Meta_penalty)[1] <- "role"
  names(Main_penalty)[1] <- "role"
  
  all_penalty_levels <- rbind(Meta_penalty, Main_penalty)
  all_penalty_levels$role <- all_penalty_levels$role %>% as.character
  
  sb$rects$name <- sb$rects$name %>% as.character
  all_penalty_levels <- na.omit(all_penalty_levels)
  temp1 <- sb$rects %>% left_join(all_penalty_levels, by = c("name"="role"))
  
  sb$rects$color <- temp1$penalty
  #the above was not working, so I made the below one up
  #sb$rects$color <- temp1$penalty[1:length(sb$rects$color)]
  sb$rects$leaf <- TRUE
  
  
  tempinnerplot <- sunburst(sb, node_labels = node_labs, node_labels.min = .2, node_labels.size = .8, leaf_labels = FALSE, rects.fill.aes = "color", rects.size  = 0.1) + scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = c(0.5, 1))
  d=data.frame(x1=.5, x2=1, y1=.5, y2=1)
  backplot_1 <- ggplot() + geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill= sum(SA_DE$penalty %>% abs < maxlim, na.rm = TRUE)/length(SA_DE$penalty))) + scale_fill_distiller(palette = "YlGnBu",direction = 1, limits = c(0.5, 1)) + theme(legend.position="none")
  prow <- plot_grid(tempinnerplot,
                    backplot_1,
                    align = 'v',
                    labels = c("inner", "middle"),
                    ncol = 1
  )
  
  plot_grid(prow)
  p <- plot_grid(prow)
  
  return(list(p, SA_DE))
}



##draw sunburst without text labels---------
library(RColorBrewer)
draw_sunburst_NO_TEXT <- function(SA_df_counts, COG_SA, sb, target_list, model_list, maxlim, node_labs)
{
  COG_SA <- COG_SA  %>% mutate(Main = if_else(Main == "Other", paste0(Main, "@", Meta), Main))
  SA_DE <- score_target_vs_model_SA(SA_df_counts, target_list, model_list)
  
  Main_penalty <- COG_SA %>% left_join(SA_DE) %>% group_by(Main) %>% select(COGID, Main, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  Meta_penalty <- COG_SA %>% left_join(SA_DE) %>% group_by(Meta) %>% select(COGID, Meta, penalty) %>% distinct %>% dplyr::summarize(penalty = mean(abs(penalty)<maxlim, na.rm=TRUE))
  Main_penalty %>% arrange(penalty)
  sum(abs(SA_DE$penalty) < 2)
  COG_SA %>% left_join(SA_DE)  %>% arrange(abs(penalty)) #%>% filter(Sub == "Sulfur metabolism")
  
  Main_penalty %>% filter()
  
  names(Meta_penalty)[1] <- "role"
  names(Main_penalty)[1] <- "role"
  
  all_penalty_levels <- rbind(Meta_penalty, Main_penalty)
  all_penalty_levels$role <- all_penalty_levels$role %>% as.character
  
  sb$rects$name <- sb$rects$name %>% as.character
  all_penalty_levels <- na.omit(all_penalty_levels)
  temp1 <- sb$rects %>% left_join(all_penalty_levels, by = c("name"="role"))
  
  sb$rects$color <- temp1$penalty
  #the above was not working, so I made the below one up
  #sb$rects$color <- temp1$penalty[1:length(sb$rects$color)]
  sb$rects$leaf <- TRUE
  
  
  tempinnerplot <- sunburst(sb, node_labels = FALSE, node_labels.min = .2, node_labels.size = .8, leaf_labels = FALSE, rects.fill.aes = "color", rects.size  = 0.1) + scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = c(0.5, 1))
  d=data.frame(x1=.5, x2=1, y1=.5, y2=1)
  backplot_1 <- ggplot() + geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill= sum(SA_DE$penalty %>% abs < maxlim, na.rm = TRUE)/length(SA_DE$penalty))) + scale_fill_distiller(palette = "YlGnBu", direction = 1, limits = c(0.5, 1)) + theme(legend.position = "none")
  prow <- plot_grid(tempinnerplot,
                    backplot_1,
                    align = 'v',
                    labels = c("inner", "middle"),
                    ncol = 1
  )
  
  plot_grid(prow)
  p <- plot_grid(prow)
  
  return(list(p, SA_DE))
}

