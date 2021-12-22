## What happens if we detect two metastatic sites at the same time point?
## look into an example;

setup(working.path = '~/Documents/GitHub/Tropism-Timing/')
clean()

# install.packages('readxl')
library(readxl)
# install.packages("fitdistrplus")
library(fitdistrplus)
library(tidyr)
# install.packages('BradleyTerry2')
library(BradleyTerry2)
library(ggplot2)
library(qvcalc)

## Data
data = read_excel('Data/tableS2_data.xlsx', skip = 2)


# starting with Prostate first
data_prostate = data[which(data$cancer_type == 'prostate cancer'), ]
data_PP = data_prostate[which(data_prostate$sample_type == 'Primary'), ]
data_PP_death = data_PP[which(data_PP$os_status == 'dead'), ]

# sites
ncols_max = max(stringr::str_count(data_PP_death$met_event_age, '\\/'), na.rm = T) + 1
sites_PP = data_PP_death[, c('met_site_mapped', 'met_event_age')]
column_names = paste0("Time", 1:ncols_max)

#' split the dataframe
data_pp_split = tidyr::separate(sites_PP, 
             col = met_event_age, 
             sep = '\\/', 
             into = column_names,
             remove = F, convert = T)


#' how many sites are detected at the first inspection (date)
data_pp_split = as.data.frame(data_pp_split)
for(i in 3:length(data_pp_split)){
  data_pp_split[, i] = as.integer(data_pp_split[, i])
}
data_pp_split[is.na(data_pp_split)] = 0
data_pp_split = data_pp_split[which(data_pp_split$Time1 != 0), ]
data_pp_split$First = rowSums(data_pp_split[, c(3:ncol(data_pp_split))] == data_pp_split[, 3], na.rm = T)

#' expand the metastatic sites
ncols_max_sites = max(stringr::str_count(data_pp_split$met_site_mapped, '\\/'), na.rm = T) + 1
column_sites = paste0('Site', 1:ncols_max_sites)
data_pp_split = tidyr::separate(data_pp_split,
                                col = met_site_mapped,
                                sep = '\\/',
                                into = column_sites,
                                remove = F,
                                convert = T)

#' extract sites involved
Diagnosis1 = data.frame()
for(i in 1:nrow(data_pp_split)){
  print(i)
  nu = data_pp_split$First[i]
  sites = paste(as.character(data_pp_split[i, c(2: (nu+1))]), collapse = '/', sep = '\\/')
  out = data.frame(Site = sites,
                   Time = data_pp_split$Time1[i],
                   number = nu)
  Diagnosis1 = rbind(Diagnosis1, out)
}

#' make a plot
plot(NA, xlim = c(1, 7), ylim = c(0, 110), xaxt = 'n')
barplot(sort(table(Diagnosis1$number), decreasing = T), 
        las = 1, main = 'Number of metastatic sites / date')

Diagnosis1



