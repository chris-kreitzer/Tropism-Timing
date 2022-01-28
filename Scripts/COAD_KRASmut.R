## Metastatic seeding pattern of Colorectal Cancers
## Is there any difference whether we have a KRAS mut vs KRAS wt
## in the metastatic seeding pattern
## 
## 01/18/22
## chris-kreitzer
## see: https://docs.google.com/document/d/1rBbO7hGpNjw3okbAkLUsmFwhwYj9pyIXx1kh7DQ25BA/edit
## 

setwd('~/Documents/GitHub/Tropism-Timing/')
clean()

## Libraries
library(readxl)
library(stringr)

## Input:
#' COAD
#' MSS (stable)
#' all: queried cBIO 01/18/22
COAD_all = read.csv('Data/COAD_MSS_all.tsv', sep = '\t')
COAD_all = COAD_all[, c('Patient.ID', 'Cancer.Type.Detailed', 'MSI.Type')]

#' COAD 
#' MSS (stable)
#' KRAS mutation
COAD = as.data.frame(t(read.csv('Data/COAD_MSS_KRASmut.tsv', sep = '\t', header = T)))
colnames(COAD) = COAD[1, ]
COAD = COAD[c(-1,-2), ]
COAD$ID = row.names(COAD)
COAD$ID = gsub(pattern = '\\.', replacement = '-', x = COAD$ID)
row.names(COAD) = NULL
COAD_MSS = COAD[which(COAD$`MSI Type` == 'Stable'), ]

#' annotation for KRAS mutants vs wild-type cases
COAD_all$KRAS = NA
for(i in 1:nrow(COAD_all)){
  COAD_all$KRAS[i] = ifelse(COAD_all$Patient.ID[i] %in% COAD_MSS$ID, 'mutant', 'wt')
}

## Input:
#' Mets (supplementary table 2)
#' met_site_mapped
#' met_event_age
mets = read_excel('Data/tableS2_data.xlsx', skip = 2)
mets = mets[, c('patient_id', 'cancer_type', 'sample_type', 'met_site_mapped', 'met_event_age')]

## Input:
#' Diagnosis age (age_dx)
#' Birth date
#' Anisha forwarded me the data
mets_time = read.csv('Data/tr_diagnosis.csv')


## Merge the data
#' mets (chronological) with KRAS mutation status
mets_COAD = merge(COAD_all, mets, by.x = 'Patient.ID', by.y = 'patient_id', all.x = T)
mets_COAD = mets_COAD[!is.na(mets_COAD$met_site_mapped), ]
#' time (diagnosis)
mets_COAD = merge(mets_COAD, mets_time[, c('DMP_ID', 'age_dx')], 
                  by.x = 'Patient.ID', 
                  by.y = 'DMP_ID', 
                  all.x = T)

mets_COAD$cancer_type = NULL


#' basic exploration:
boxplot(as.numeric(mets_COAD$age_dx) ~ mets_COAD$KRAS, 
        xlab = '', ylab = 'Days [Diagnosis]', main = 'Colorectal Cancer Patients')
title(sub = 'p-value = 0.01272 [t.test]')


## data modifications:
#' split the data
ncols_max = max(stringr::str_count(mets_COAD$met_event_age, '\\/'), na.rm = T) + 1
column_names = paste0('Rank', 1:ncols_max)
data_split = tidyr::separate(data = mets_COAD,
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

#' how many metastatic sites do we see:
data_split$metCount = rowSums(!is.na(data_split[, 7:37]))
barplot(sort(data_split$metCount, decreasing = T))
abline(h = quantile(data_split$metCount, probs = 0.9),
       lty = 'dashed', col = 'red')
text(x = 3900, y = 9.5, label = '<= 12 mets\n90% quantile')


#' exclude patients with more than 12 metastatic sites
data_split = data_split[, -(19:37)]
data_split = data_split[, -(32:50)]

cols.num = 20:32
data_split[cols.num] = sapply(data_split[cols.num], as.numeric)


## Which sites are over represented among the two groups?
mut = data_split[which(data_split$KRAS == 'mutant'), 7:18]
mut_table = table(t(mut))
par(mfrow = c(2,1))
pdf("MetFrequency_COAD_KRAS.pdf", width = 10, height = 18) 
barplot(sort(mut_table, decreasing = T), 
        las = 2, 
        cex.names = 0.65,
        ylab = '#metastatic lesions',
        yaxt = 'n',
        ylim = c(0, 2100))
axis(2, at = c(0, 500, 1000, 1500, 2001), labels = c(0, 500, 1000, 1500, 2000), line = -1.5, 
     lwd = 2, lwd.ticks = 2, hadj = 1, padj = 0.5, las = 1)
title(main = 'KRAS mutant Colorectal Cancers\nn=1,686')
box(which = 'figure', lwd = 2)

#' wild type data
wt = data_split[which(data_split$KRAS == 'wt'), 7:18]
wt_table = table(t(wt))
barplot(sort(wt_table, decreasing = T), 
        las = 2, 
        cex.names = 0.65,
        ylab = '#metastatic lesions',
        yaxt = 'n',
        ylim = c(0, 2100))
axis(2, at = c(0, 500, 1000, 1500, 2001), labels = c(0, 500, 1000, 1500, 2000), line = -1.5, 
     lwd = 2, lwd.ticks = 2, hadj = 1, padj = 0.5, las = 1)
title(main = 'KRAS wild-type Colorectal Cancers\nn=2,063')
box(which = 'figure', lwd = 2)

dev.off(which = dev.cur())

#' Are lung metastasis recurring in patients with KRASmut?
lung_mut = as.data.frame(apply(mut, 1, function(x) sum(str_count(string = x, pattern = 'lung'), na.rm = T)))
colnames(lung_mut) = 'count'

par(mfrow = c(2,1))
v = lung_mut$count
a = min(v)
b = max(v)
lung_hist = hist(v+0.001, 
           breaks = b-a, 
           xaxt = "n",
           las = 1,
           lwd = 2,
           hadj = 1,
           col = "orange",
           panel.first = grid(), 
           main = "Lung metastasis/patient\nKRAS-mutants", 
           xlab = "")
axis(side = 1, at = lung_hist$mids, labels = seq(a,b), lwd = 2)

#' Lung metastasis at KRAS wt patient
lung_wt = as.data.frame(apply(wt, 1, function(x) sum(str_count(string = x, pattern = 'lung'), na.rm = T)))
colnames(lung_wt) = 'count'

v = lung_wt$count
a = min(v)
b = max(v)
lung_hist = hist(v + 0.001, 
                 breaks = b - a, 
                 xaxt = "n",
                 las = 1,
                 lwd = 2,
                 hadj = 1,
                 col = "orange",
                 panel.first = grid(), 
                 main = "Lung metastasis/patient\nKRAS-wt", 
                 xlab = "")
axis(side = 1, at = lung_hist$mids, labels = seq(a,b), lwd = 2)

#' Lung metastatic recurrence attributable to KRASmut status?
data_split$lungMets = apply(data_split[,7:18], 1, 
                            function(x) sum(str_count(string = x, pattern = 'lung'), na.rm = T))
data_split$lungMetsrel = ifelse(data_split$met_event_age != 1 & !data_split$lungMets %in% c(0,1), 
                                data_split$lungMets / data_split$metCount, NA)

par(oma = c(5,10,5,10))
boxplot(data_split$lungMetsrel ~ data_split$KRAS,
        las = 2,
        xaxt = 'n',
        width = c(1,1), 
        col = 'orange', 
        xlab = '',
        ylab = 'avg.lung.met.fraction')
axis(side = 1, at = c(1,2), labels = c('KRAS\nmutant', 'KRAS\nwt'),
     font = 2, tick = F, outer = NA, line = 0.5)
text(x = 1.5, y = 0.9, label = 'p-value = 0.598')
box(lwd = 2)


## Finding repetitive tropism patterns in the data:
mut_all = data_split[which(data_split$KRAS == 'mutant'), 1:18]
mut_rle = rle(mut_all$met_site_mapped)

#' which pattern occurs most frequent?
mut_rle$values[which.max(mut_rle$lengths)]

met_seq = list()
for(i in 3:max(mut_rle$lengths)){
  seq_organgs = mut_rle$values[which(mut_rle$lengths == i)]
  met_seq[[i]] = seq_organgs
  names(met_seq)[[i]] = i
}

met_seq









