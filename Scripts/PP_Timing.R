## Tropism. Pattern of metastatic dissemination;
## We start with Prostate Cancer (primaries)
## We will start in creating a proper table which is 
## suitable for the BT model and then move on to the MLE-
## algorithm, I just established
## 
## chris-kreitzer
## 12/17/2021


setup(working.path = '~/Documents/GitHub/Tropism-Timing/')
clean()

# install.packages('readxl')
library(readxl)
# install.packages("fitdistrplus")
library(fitdistrplus)
library(tidyr)
install.packages('BradleyTerry2')
library(BradleyTerry2)


## Data
data = read_excel('Data/tableS2_data.xlsx', skip = 2)


# starting with Prostate first
data_prostate = data[which(data$cancer_type == 'prostate cancer'), ]
data_PP = data_prostate[which(data_prostate$sample_type == 'Primary'), ]
data_PP_death = data_PP[which(data_PP$os_status == 'dead'), ]

#' start with simple example
plot.new()
hist(data_PP_death$met_count, nclass = nclass.scott(data_PP_death$met_count), xaxt = 'n')
axis(1, at = seq(0, 22, 2))
summary(data_PP_death$met_count)

#' OS for primary prostate cancer patients
plot.new()
hist(data_PP_death$os_days, 
     nclass = nclass.scott(data_PP_death$os_days), 
     xaxt = 'n',
     main = paste0('Median OS = ', summary(data_PP_death$os_days)[['Median']], ' days'))
axis(1, at = seq(0, 2200, by = 500))

#' is there a association with metastatic sites and OS_days?
plot(data_PP_death$met_count, data_PP_death$os_days, main = 'Association with met. counts and OS', xlab = 'met count', ylab = 'OS')
abline(lm(data_PP_death$os_days ~ data_PP_death$met_count), lwd = 2.5, col = 'red')
model.simple = lm(data_PP_death$os_days ~ data_PP_death$met_count)
summary(lm(data_PP_death$os_days ~ data_PP_death$met_count))
plot(model.simple)

#' site distribution
sites = as.data.frame(data_PP_death[, c('met_site_mapped')])
sites_split = strsplit(as.character(sites$met_site_mapped), '/', fixed = T)
max.length = max(sapply(sites_split, length))
sites_new = lapply(sites_split, function(x){c(x, rep(NA, max.length - length(x)))})
sites_new = as.data.frame(do.call(rbind, sites_new))
View(sites_new)

barplot(sort(table(sites_new$V1), decreasing = T), las = 2, main = 'Counts: metastatic site 1')
barplot(sort(table(sites_new$V2), decreasing = T), las = 2)
barplot(sort(table(sites_new$V3), decreasing = T), las = 2)


#' modify the sites data frame into a suitable format for BT modeling
#' exclude entries where every site is NA in table
sites_new = sites_new[!is.na(sites_new$V1), ]
unique.sites = length(unique(sites_new$V1))
sites_names = unique(sites_new$V1)
occurrence = table(sites_new$V1)

#' modify the table
m = matrix(nrow = unique.sites, 
           ncol = unique.sites, 
           dimnames = list(sites_names, sites_names))
m = as.data.frame(m)

for(i in 1:nrow(m)){
  m[row.names(m)[i], ] = occurrence[which(names(occurrence) == row.names(m)[i])]
}














