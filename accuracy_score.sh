###########

excel_sheets("humann_jointtable_eggnog_regrouped_unstratified.xlsx")
counts_normalized_with_IVW <- read_excel("humann_jointtable_eggnog_regrouped_unstratified.xlsx", sheet = 8)
dim(counts_normalized_with_IVW) 
names(counts_normalized_with_IVW) #"COGID" should be the name of column1

#read in the sunburstmetadata file 
excel_sheets("humann_jointtable_eggnog_regrouped_unstratified.xlsx")
metadata_file <- read_excel("humann_jointtable_eggnog_regrouped_unstratified.xlsx", sheet = 9)
dim(metadata_file) #
#set metadata: change list names and filter condition to match dataset for comparison
#Here I want the human samples to be clearly defined.
#define sample lists

human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")
length(human_list) #67
#set count normalized in function to be counts thati have here 
counts_normalized <- counts_normalized_with_IVW

#set conditions_list without COG ID
conditions_list <- names(counts_normalized_with_IVW)[-1]
conditions_list #samples only


##########################
#Now do the human resampling. 
#run the human samples subsamplig just for clarity
#setsession
iterations <- 1200 #set number of times to resample
leaveout <- 2 #set number of samples to use as test model
output_version <- "run1"

#run multiple iterations of AS2 using X samples as the gold-standard and Y samples as the "model." For example, for my 66 human samples, I run 64 as the human "gold-standard" and 2 as the test across 100 iterations.
#define functions below
set.seed(123)

counts_normalized_double <- c(1:iterations) %>% map(function(x) {dummy1 <- model_selfvalidation_leaveout(self_validation, leaveout)
dummy1$penalty <- abs(dummy1$penalty)
names(dummy1)[2] <- paste0("penalty_", x)
return(dummy1)})  %>% purrr::reduce(left_join)


write.table(counts_normalized_double, file = paste(output_version, "all_human_samples_resampled_2025", iterations, "_leaveout", leaveout, ".txt", sep = ""), sep ="\t", row.names = FALSE) 
#output is penalties (z-scores). In R or excel, calculate AS2 for each iteration then average iterations





########do this for the 4member community
metadata_file <- read_excel("humann_jointtable_eggnog_regrouped_unstratified.xlsx", sheet = 9)
#first make sample ID the rownames of the metadata
metadata_rows <- as.data.frame(metadata_file)
rownames(metadata_rows) <- metadata_rows$filename

#extract for fourmember_48Hwlm
# get the unique grousp in the metadata
unique(metadata_file$Condition) 

#get samples other than Human and four_mem_wlm
fourmember_48Hwlm <- rownames(metadata_rows)[metadata_rows$Condition %in% c("7member_48Hwlm_onetoone", "4mem_4dpi", "7mem_4dpi", "7mem_gluc_48H", "7member_48Hwlm")]

# select the columns that are not in in fourmember_48Hwlm
fourmember_48Hwlm_columns <- colnames(counts_normalized)[!colnames(counts_normalized) %in% fourmember_48Hwlm]

#get the counts
fourmember_48Hwlm_counts <- counts_normalized[, fourmember_48Hwlm_columns]

head(fourmember_48Hwlm_counts) #this has the samples and the COGIDs


#get the metadata
fourmember_48Hwlm_metadata_column <- rownames(metadata_rows)[metadata_rows$Condition %in% c("Human", "4member_48Hwlm")]

fourmember_48Hwlm_metadata <- metadata_rows[rownames(metadata_rows) %in% fourmember_48Hwlm_metadata_column, ]
head(fourmember_48Hwlm_metadata)
#remove the rownames
rownames(fourmember_48Hwlm_metadata) <- NULL
head(fourmember_48Hwlm_metadata)

dim(fourmember_48Hwlm_counts) 
dim(fourmember_48Hwlm_metadata)


#Now begin calculating for fourmember_48Hwlm
metadata_file <- fourmember_48Hwlm_metadata
human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")
invitro_list <- metadata_file %>% filter(type1 == "IVW") %>% .$filename %>% str_replace_all("-", "_") #adjust filter (currently "type1") to match correct column of metadata_file

#Calculate mean AS scores
#set session
reps <- 2 #number of replicates to run at once if subsampling
outputversion <- paste("run1.sub", reps, sep = "")
invitro_dataset <- "IVW"

#JUST OUTPUT AS2: outputs average AS2 value across iterations
c(1:choose(length(invitro_list), reps)) %>% purrr::map(function(.x) {dummy1 <- score_target_vs_model_SA(fourmember_48Hwlm_counts, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
dummy2 <- abs(dummy1) < 2
dummy3 <- mean(dummy2)
return(dummy3)}) %>% purrr::reduce(c) %>% mean




#########
########do for 4member in the murine data
metadata_file <- read_excel("humann_jointtable_eggnog_regrouped_unstratified.xlsx", sheet = 9)
#first make sample ID the rownames of the metadata
metadata_rows <- as.data.frame(metadata_file)
rownames(metadata_rows) <- metadata_rows$filename

#extract for fourmember_4dpi
#get samples other than Human and 4mem 4dpi
fourmember_4dpi <- rownames(metadata_rows)[metadata_rows$Condition %in% c("4member_48Hwlm", "7member_48Hwlm_onetoone", "7mem_4dpi", "7mem_gluc_48H", "7member_48Hwlm"  )]

# select the columns that are not in in fourmember_4dpi
fourmember_4dpi_columns <- colnames(counts_normalized)[!colnames(counts_normalized) %in% fourmember_4dpi]

#get the counts
fourmember_4dpi_counts <- counts_normalized[, fourmember_4dpi_columns]

head(fourmember_4dpi_counts) #this has the samples and the COGIDs


#get the metadata
fourmember_4dpi_metadata_column <- rownames(metadata_rows)[metadata_rows$Condition %in% c("Human", "4mem_4dpi")]

fourmember_4dpi_metadata <- metadata_rows[rownames(metadata_rows) %in% fourmember_4dpi_metadata_column, ]
head(fourmember_4dpi_metadata)
#remove the rownames
rownames(fourmember_4dpi_metadata) <- NULL
head(fourmember_4dpi_metadata)

dim(fourmember_4dpi_counts) 
dim(fourmember_4dpi_metadata) 


#Now begin calculating for fourmember_4dpi
metadata_file <- fourmember_4dpi_metadata
human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")
invitro_list <- metadata_file %>% filter(type1 == "IVW") %>% .$filename %>% str_replace_all("-", "_") #adjust filter (currently "type1") to match correct column of metadata_file

#Calculate mean AS scores
#set session
reps <- 2 #number of replicates to run at once if subsampling
outputversion <- paste("run1.sub", reps, sep = "")
invitro_dataset <- "IVW"

#JUST OUTPUT AS2: outputs average AS2 value across iterations
c(1:choose(length(invitro_list), reps)) %>% purrr::map(function(.x) {dummy1 <- score_target_vs_model_SA(fourmember_4dpi_counts, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
dummy2 <- abs(dummy1) < 2
dummy3 <- mean(dummy2)
return(dummy3)}) %>% purrr::reduce(c) %>% mean



###############
###do for the seven member community with reduced PA ratio
metadata_file <- read_excel("humann_jointtable_eggnog_regrouped_unstratified.xlsx", sheet = 9)
#first make sample ID the rownames of the metadata
metadata_rows <- as.data.frame(metadata_file)
rownames(metadata_rows) <- metadata_rows$filename

#extract for sevenmember_48Hwlm
#get samples other than Human and 7mem 48H WLM
sevenmember_48Hwlm <- rownames(metadata_rows)[metadata_rows$Condition %in% c("7member_48Hwlm_onetoone",  "4member_48Hwlm", "4mem_4dpi", "7mem_4dpi", "7mem_gluc_48H" )]

# select the columns that are not in in sevenmember_48Hwlm
sevenmember_48Hwlm_columns <- colnames(counts_normalized)[!colnames(counts_normalized) %in% sevenmember_48Hwlm]

#get the counts
sevenmember_48Hwlm_counts <- counts_normalized[, sevenmember_48Hwlm_columns]

head(sevenmember_48Hwlm_counts) #this has the samples and the COGIDs


#get the metadata
sevenmember_48Hwlm_metadata_column <- rownames(metadata_rows)[metadata_rows$Condition %in% c("Human", "7member_48Hwlm")]

sevenmember_48Hwlm_metadata <- metadata_rows[rownames(metadata_rows) %in% sevenmember_48Hwlm_metadata_column, ]
head(sevenmember_48Hwlm_metadata)
#remove the rownames
rownames(sevenmember_48Hwlm_metadata) <- NULL
head(sevenmember_48Hwlm_metadata)

dim(sevenmember_48Hwlm_counts) 
dim(sevenmember_48Hwlm_metadata) 


#Now begin calculating for sevenmember_48Hwlm
metadata_file <- sevenmember_48Hwlm_metadata
human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")
invitro_list <- metadata_file %>% filter(type1 == "IVW") %>% .$filename %>% str_replace_all("-", "_") #adjust filter (currently "type1") to match correct column of metadata_file

#Calculate mean AS scores
#set session
reps <- 2 #number of replicates to run at once if subsampling
outputversion <- paste("run1.sub", reps, sep = "")
invitro_dataset <- "IVW"

#JUST OUTPUT AS2: outputs average AS2 value across iterations
c(1:choose(length(invitro_list), reps)) %>% purrr::map(function(.x) {dummy1 <- score_target_vs_model_SA(sevenmember_48Hwlm_counts, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
dummy2 <- abs(dummy1) < 2
dummy3 <- mean(dummy2)
return(dummy3)}) %>% purrr::reduce(c) %>% mean



##########
####do for 7member one to one with the base ratio of species
metadata_file <- read_excel("humann_jointtable_eggnog_regrouped_unstratified.xlsx", sheet = 9)
#first make sample ID the rownames of the metadata
metadata_rows <- as.data.frame(metadata_file)
rownames(metadata_rows) <- metadata_rows$filename

#extract for sevenmember_48Hwlm_onetoone
#get samples other than Human and 7mem 1:1
sevenmember_48Hwlm_onetoone <- rownames(metadata_rows)[metadata_rows$Condition %in% c("7member_48Hwlm",  "4member_48Hwlm", "4mem_4dpi", "7mem_4dpi", "7mem_gluc_48H" )]

# select the columns that are not in in sevenmember_48Hwlm_onetoone
sevenmember_48Hwlm_onetoone_columns <- colnames(counts_normalized)[!colnames(counts_normalized) %in% sevenmember_48Hwlm_onetoone]

#get the counts
sevenmember_48Hwlm_onetoone_counts <- counts_normalized[, sevenmember_48Hwlm_onetoone_columns]

head(sevenmember_48Hwlm_onetoone_counts) #this has the samples and the COGIDs


#get the metadata
sevenmember_48Hwlm_onetoone_metadata_column <- rownames(metadata_rows)[metadata_rows$Condition %in% c("Human", "7member_48Hwlm_onetoone")]

sevenmember_48Hwlm_onetoone_metadata <- metadata_rows[rownames(metadata_rows) %in% sevenmember_48Hwlm_onetoone_metadata_column, ]
head(sevenmember_48Hwlm_onetoone_metadata)
#remove the rownames
rownames(sevenmember_48Hwlm_onetoone_metadata) <- NULL
head(sevenmember_48Hwlm_onetoone_metadata)

dim(sevenmember_48Hwlm_onetoone_counts) 
dim(sevenmember_48Hwlm_onetoone_metadata) 


#Now begin calculating for sevenmember_48Hwlm_onetoone
metadata_file <- sevenmember_48Hwlm_onetoone_metadata
human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")
invitro_list <- metadata_file %>% filter(type1 == "IVW") %>% .$filename %>% str_replace_all("-", "_") #adjust filter (currently "type1") to match correct column of metadata_file

#Calculate mean AS scores
#set session
reps <- 2 #number of replicates to run at once if subsampling
outputversion <- paste("run1.sub", reps, sep = "")
invitro_dataset <- "IVW"

#JUST OUTPUT AS2: outputs average AS2 value across iterations
c(1:choose(length(invitro_list), reps)) %>% purrr::map(function(.x) {dummy1 <- score_target_vs_model_SA(sevenmember_48Hwlm_onetoone_counts, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
dummy2 <- abs(dummy1) < 2
dummy3 <- mean(dummy2)
return(dummy3)}) %>% purrr::reduce(c) %>% mean





############
######do it for the 7member in hyperglycemia
metadata_file <- read_excel("humann_jointtable_eggnog_regrouped_unstratified.xlsx", sheet = 9)
#first make sample ID the rownames of the metadata
metadata_rows <- as.data.frame(metadata_file)
rownames(metadata_rows) <- metadata_rows$filename

#extract for sevenmember_gluc_48H


#get samples other than Human and 7mem gluc
sevenmember_gluc_48H <- rownames(metadata_rows)[metadata_rows$Condition %in% c("7member_48Hwlm_onetoone",  "4member_48Hwlm", "7mem_4dpi", "7member_48Hwlm", "4mem_4dpi" )]

# select the columns that are not in in sevenmember_gluc_48H
sevenmember_gluc_48H_columns <- colnames(counts_normalized)[!colnames(counts_normalized) %in% sevenmember_gluc_48H]

#get the counts
sevenmember_gluc_48H_counts <- counts_normalized[, sevenmember_gluc_48H_columns]

head(sevenmember_gluc_48H_counts) #this has the samples and the COGIDs


#get the metadata
sevenmember_gluc_48H_metadata_column <- rownames(metadata_rows)[metadata_rows$Condition %in% c("Human", "7mem_gluc_48H")]

sevenmember_gluc_48H_metadata <- metadata_rows[rownames(metadata_rows) %in% sevenmember_gluc_48H_metadata_column, ]
head(sevenmember_gluc_48H_metadata)
#remove the rownames
rownames(sevenmember_gluc_48H_metadata) <- NULL
head(sevenmember_gluc_48H_metadata)

dim(sevenmember_gluc_48H_counts) 
dim(sevenmember_gluc_48H_metadata) 


#Now begin calculating for sevenmember_gluc_48H
metadata_file <- sevenmember_gluc_48H_metadata
human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")
invitro_list <- metadata_file %>% filter(type1 == "IVW") %>% .$filename %>% str_replace_all("-", "_") #adjust filter (currently "type1") to match correct column of metadata_file

#Calculate mean AS scores
#set session
reps <- 2 #number of replicates to run at once if subsampling
outputversion <- paste("run1.sub", reps, sep = "")
invitro_dataset <- "IVW"

#JUST OUTPUT AS2: outputs average AS2 value across iterations
c(1:choose(length(invitro_list), reps)) %>% purrr::map(function(.x) {dummy1 <- score_target_vs_model_SA(sevenmember_gluc_48H_counts, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
dummy2 <- abs(dummy1) < 2
dummy3 <- mean(dummy2)
return(dummy3)}) %>% purrr::reduce(c) %>% mean


############
########do 7 member in the murine wounds 
metadata_file <- read_excel("humann_jointtable_eggnog_regrouped_unstratified.xlsx", sheet = 9)
#first make sample ID the rownames of the metadata
metadata_rows <- as.data.frame(metadata_file)
rownames(metadata_rows) <- metadata_rows$filename

#extract for sevenmember_4dpi

#get samples other than Human and 7mem murine
sevenmember_4dpi <- rownames(metadata_rows)[metadata_rows$Condition %in% c("7member_48Hwlm_onetoone",  "4member_48Hwlm", "7mem_gluc_48H", "7member_48Hwlm", "4mem_4dpi"  )]

# select the columns that are not in in sevenmember_4dpi
sevenmember_4dpi_columns <- colnames(counts_normalized)[!colnames(counts_normalized) %in% sevenmember_4dpi]

#get the counts
sevenmember_4dpi_counts <- counts_normalized[, sevenmember_4dpi_columns]

head(sevenmember_4dpi_counts) #this has the samples and the COGIDs


#get the metadata
sevenmember_4dpi_metadata_column <- rownames(metadata_rows)[metadata_rows$Condition %in% c("Human", "7mem_4dpi")]

sevenmember_4dpi_metadata <- metadata_rows[rownames(metadata_rows) %in% sevenmember_4dpi_metadata_column, ]
head(sevenmember_4dpi_metadata)
#remove the rownames
rownames(sevenmember_4dpi_metadata) <- NULL
head(sevenmember_4dpi_metadata)

dim(sevenmember_4dpi_counts) 
dim(sevenmember_4dpi_metadata) 


#Now begin calculating for sevenmember_4dpi
metadata_file <- sevenmember_4dpi_metadata
human_list <- metadata_file %>% filter(type1 == "human") %>% .$filename %>% str_replace_all("-", "_")
invitro_list <- metadata_file %>% filter(type1 == "IVW") %>% .$filename %>% str_replace_all("-", "_") #adjust filter (currently "type1") to match correct column of metadata_file

#Calculate mean AS scores
#set session
reps <- 2 #number of replicates to run at once if subsampling
outputversion <- paste("run1.sub", reps, sep = "")
invitro_dataset <- "IVW"

#JUST OUTPUT AS2: outputs average AS2 value across iterations
c(1:choose(length(invitro_list), reps)) %>% purrr::map(function(.x) {dummy1 <- score_target_vs_model_SA(sevenmember_4dpi_counts, human_list, combn(invitro_list, reps)[,.x]) %>% .$penalty
dummy2 <- abs(dummy1) < 2
dummy3 <- mean(dummy2)
return(dummy3)}) %>% purrr::reduce(c) %>% mean



