## 'Time variable' in metastatic dissemination
## Along with the Tropism project
## 
## The longer the disease course, the more metastatic sites detected?
## 
## 01/11/2022
## chris-kreitzer

clean()
setwd('~/Documents/GitHub/Tropism-Timing/')


## Data
data = read_excel('Data/tableS2_data.xlsx', skip = 2)
data = data[, c('cancer_type', 'sample_type', 'met_site_mapped', 'met_event_age')]

#' split the data
ncols_max = max(stringr::str_count(data$met_event_age, '\\/'), na.rm = T) + 1
column_names = paste0('Rank', 1:ncols_max)
data_split = tidyr::separate(data = data,
                             col = 'met_site_mapped',
                             sep = '\\/',
                             into = column_names,
                             remove = F,
                             convert = T)

ncols_time = max(stringr::str_count(data_split$met_event_age, '\\/'), na.rm = T) + 1
column_times = paste0('Time', 1:ncols_time)
data_split = tidyr::separate(data = data_split,
                             col = 'met_event_age',
                             sep = '\\/',
                             into = column_times,
                             remove = F,
                             convert = T)

#' modies
data_split = data_split[!is.na(data_split$met_site_mapped), ]
data_split = as.data.frame(data_split)
for(i in 40:ncol(data_split)){
  data_split[, i] = as.numeric(as.integer(data_split[, i]))
}


#' just keep samples < 12 metastatic sites
data_split = data_split[which(data_split$metCount <= 12), ]


###############################################################################
#' Prostate
Prostate = data_split[which(data_split$cancer_type == 'prostate cancer'), ]
Prostate = Prostate[which(Prostate$Rank1 != 'na'), ]

#' delete columns with no entry
for(i in 1:13){
  if(all(is.na(Prostate[, i]))){
    Prostate[, i] = NULL
  }
}

#' delete time for now
Prostate = Prostate[,c(4:15, 40:51)]






















