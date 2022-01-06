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


#' dummy example with 5 rankings and 10 items to choose
dummy = data.frame(rank1 = c('regio', 'bone', 'bladder', 'other', 'bone', 'dist', 'dist', 'bone', 'bladder', 'dist'),
                   rank2 = c('bladder', 'other', 'bone', 'adrenal', 'bone', 'bone', 'regio', 'regio', 'gential', 'regio'),
                   rank3 = c('gential', 'dist', 'other', 'regio', NA, 'liver', 'bone', 'bone', 'other', 'other'),
                   rank4 = c('other', 'regio', 'dist', 'bladder', NA, NA, 'other', 'dist', 'breast', 'dist'),
                   rank5 = c('bone', 'genital', 'regio', 'genital', NA, NA, 'cns_brain', 'periphalNS', 'bone', 'regio'))

dummy = as.matrix(dummy)

dummy.ranking = as.rankings(x = dummy, input = 'orderings')
dummy.ranking[1]


mod = PlackettLuce(rankings = dummy.ranking)
coefs2 = round(coef(mod), 2)
coefs2[names(coefs2)[1:3]]





library(stringr)
sum(str_count(string = 'regio', dummy), na.rm = T)



