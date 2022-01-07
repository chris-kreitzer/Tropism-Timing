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


adjacency(dummy.ranking)

summary(mod, ref = NULL)














