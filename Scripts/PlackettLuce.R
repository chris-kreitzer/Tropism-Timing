## Introduction of the PlackettLuce Model
## 
## 01/04/2022
## chris-kreitzer

clean()
setwd('~/Documents/GitHub/Tropism-Timing/')


#install.packages('PlackettLuce')
library(PlackettLuce)



## Data
data = read_excel('Data/tableS2_data.xlsx', skip = 2)

# starting with Prostate first
data_prostate = data[which(data$cancer_type == 'prostate cancer'), ]
data_PP = data_prostate[which(data_prostate$sample_type == 'Primary'), ]
data_PP_death = data_PP[which(data_PP$os_status == 'dead'), ]

sites = as.data.frame(data_PP_death[, c('met_site_mapped')])
sites_split = strsplit(as.character(sites$met_site_mapped), '/', fixed = T)
max.length = max(sapply(sites_split, length))
sites_new = lapply(sites_split, function(x){c(x, rep(NA, max.length - length(x)))})
sites_new = as.data.frame(do.call(rbind, sites_new))
sites_new = sites_new[!is.na(sites_new$V1), ]
colnames(sites_new) = paste0('rank', seq(1, length(sites_new)))
sites_new = sites_new[!is.na(sites_new$rank2), ]
sites_new = sites_new[which(sites_new$rank1 != sites_new$rank2), ]

#' make a ranking object
x = as.matrix(sites_new)
y = as.rankings(x, input = 'orderings')
mod = PlackettLuce(y)
avRank <- apply(y, 2, function(x) mean(x[x > 0]))
sort(avRank)


#' dummy example with 10 rankings and 11 items to choose;
#' see google doc for example

# dummy = matrix(c(1,2,3,4,5, NA, NA, NA, NA, NA, NA,
#                  4, NA, 5, 2,1,3, NA, NA, NA, NA, NA,
#                  5, 1, NA, 3, 2, 4, NA, NA, NA, NA, NA,
#                  3, 4, 5, 1, NA, NA, 2, NA, NA, NA, NA,
#                  NA, NA, NA, NA, 1, NA, NA, NA, NA, NA, NA,
#                  NA, NA, NA, NA, 2, 1, NA, NA, NA, NA, NA,
#                  2, NA, NA, 4, 3, 1, NA, 5, NA, NA, 3,
#                  2, NA, NA, NA, 1, 3, NA, NA, 4, NA, NA,
#                  NA, 1, 2, 3, 5, NA, NA, NA, NA, 4, NA,
#                  2, NA, NA, 3, NA, 1, NA, NA, NA, NA, NA),
#                nrow = 10, byrow = T)

# colnames(dummy) = c('regio', 'bladder', 'gential', 'other', 'bone',
#                     'dist', 'adrenal', 'cns', 'peripheral', 'breast', 'liver')

dummy = matrix(c(1,2,3,4,5,
                 5,4,6,1,3,
                 2,5,4,6,1,
                 4,7,1,2,3,
                 5, NA, NA, NA, NA,
                 6,5,11, NA, NA,
                 6,1,5,4,8,
                 5,1,6,9,NA,
                 2,3,4,10,5,
                 6,1,4, NA, NA),
               nrow = 10, byrow = T)
colnames(dummy) = paste0('rank', seq(1, 5, 1))

sites = c('regio', 'bladder', 'gential', 'other', 'bone',
          'dist', 'adrenal', 'cns', 'peripheral', 'breast', 'liver')
names(sites) = seq(1,11, 1)
attr(dummy, which = 'sites') = sites

#' make a ranking object
dummy.ranking = as.rankings(x = dummy, 
                            input = 'orderings', 
                            items = attr(dummy, 'sites'))

#' exploratory function
avRank = apply(dummy.ranking, 2, function(x) mean(x[x > 0]))
barplot(sort(avRank), las = 2, ylab = 'average rank')

#' statistically modelling the outcome with PlackettLuce
mod = PlackettLuce(rankings = dummy.ranking)
coef(summary(mod))
plot(sort(coef(mod)))
abline(h = 0)

#' quasi-variance errors; if the error intervals overlap with another
#' we can say that those are significantly different
qv = qvcalc(mod)
qv$qvframe = qv$qvframe[order(coef(mod)),]
plot(qv, xlab = NULL, ylab = "Ability (log)", main = NULL,
     xaxt = "n", xlim = c(1, 11))
axis(1, at = seq_len(11), labels = rownames(qv$qvframe), las = 2, cex.axis = 0.6)





###############################################################################
###############################################################################
## PanCancer: Starting with Francisco's Table
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
#' taking prostate cancer as an example; where we compare the BradleyTerry model with 
#' the PlackettLuce model (moreover, we can reference this with the PanCancer cohort)
Prostate = data_split[which(data_split$cancer_type == 'prostate cancer'), ]

#' delete columns with no entry
for(i in 4:26){
  print(i)
  if(all(is.na(Prostate[, i]))){
    Prostate[, i] = NULL
  }
}
Prostate[,16:20] = NULL
Prostate[,29:39] = NULL

#' how many unique sites are there in the data?
unique_sites_Prostate = as.character(unique(unlist(Prostate[, 4:15])))
names(unique_sites_Prostate) = seq(1, length(unique(unique_sites_Prostate)), 1)

Prostate_rankings = Prostate[, 4:15]

#' recode the named characters into numeric values
Prostate_rankings[Prostate_rankings == "regional_lymph"] = 1L
Prostate_rankings[Prostate_rankings == "genital_male"] = 4L
Prostate_rankings[Prostate_rankings == "lung"] = 7L
Prostate_rankings[Prostate_rankings == "peritoneum"] = 10L
Prostate_rankings[Prostate_rankings == "pleura"] = 13L
Prostate_rankings[Prostate_rankings == "adrenal_gland"] = 16L
Prostate_rankings[Prostate_rankings == "breast"] = 19L
Prostate_rankings[Prostate_rankings == "mediastinum"] = 22L
Prostate_rankings[Prostate_rankings == "regional_lymph"] = 1L
Prostate_rankings[Prostate_rankings == "other"] = 2L
Prostate_rankings[Prostate_rankings == "dist_lymph"] = 5L
Prostate_rankings[Prostate_rankings == "liver"] = 8L
Prostate_rankings[Prostate_rankings == "bowel"] = 11L
Prostate_rankings[Prostate_rankings == "kidney"] = 14L
Prostate_rankings[Prostate_rankings == "skin"] = 17L
Prostate_rankings[Prostate_rankings == "na"] = 20L
Prostate_rankings[Prostate_rankings == "head_and_neck"] = 23L
Prostate_rankings[Prostate_rankings == "bone"] = 3L
Prostate_rankings[Prostate_rankings == "bladder_or_urinary_tract"] = 6L
Prostate_rankings[Prostate_rankings == "lymph"] = 9L
Prostate_rankings[Prostate_rankings == "cns_brain"] = 12L
Prostate_rankings[Prostate_rankings == "peripheral_nervous_system"] = 15L
Prostate_rankings[Prostate_rankings == "biliary_tract"] = 18L

#' conversion to matrix
x = as.matrix(Prostate_rankings, rownames.force = F)
mode(x) = 'numeric'
attr(x, which = 'sites') = unique_sites_Prostate

#' make ranking of the whole data-frame
x_ranked = as.rankings(x = x, input = 'orderings', items = attr(x, 'sites'))
