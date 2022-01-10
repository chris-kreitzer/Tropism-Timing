## Tour through Prostate Cancer
## 
## 01/09/2022
## chris-kreitzer


clean()
setwd('~/Documents/GitHub/Tropism-Timing/')


#install.packages('PlackettLuce')
library(BradleyTerry2)


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

#' how many sites are we seeing; how is the frequency of these sites
unique_sites = as.character(unique(unlist(data_split[, 4:38])))
data_split$metCount = rowSums(!is.na(data_split[, 4:38]))

barplot(sort(data_split$metCount))
abline(h = quantile(data_split$metCount, probs = 0.9),
       col = 'red',
       lwd = 2,
       lty = 'dashed')


#' just keep samples < 12 metastatic sites
data_split = data_split[which(data_split$metCount <= 12), ]


###############################################################################
###############################################################################
#' Prostate
Prostate = data_split[which(data_split$cancer_type == 'prostate cancer'), ]
Prostate = Prostate[which(Prostate$Rank1 != 'na'), ]

barplot(sort(table(Prostate$metCount), decreasing = T),
        xlab = 'Metastases', ylab = 'Count', main = 'Metastatic counts for prostate cancer')
par(mar = c(8, 4,4,3))
barplot(sort(table(Prostate$Rank1), decreasing = T), las = 2, cex.names = 0.65)
barplot(sort(table(Prostate$Rank2), decreasing = T), las = 2, cex.names = 0.65)
barplot(sort(table(Prostate$Rank3), decreasing = T), las = 2, cex.names = 0.65)





