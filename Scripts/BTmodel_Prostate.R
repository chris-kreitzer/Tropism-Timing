## Tour through Prostate Cancer;
## We start with the BradleyTerry Model and move on to the PlackettLuce model
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


########
#' modify the data for the BT model:
unique.sites = length(unique(Prostate$Rank1))
sites_names = unique(Prostate$Rank1)
occurrence = table(Prostate$Rank1)

#' modify the table
m = matrix(nrow = unique.sites, 
           ncol = unique.sites, 
           dimnames = list(sites_names, sites_names))
m = as.data.frame(m)

for(i in 1:nrow(m)){
  m[row.names(m)[i], ] = occurrence[which(names(occurrence) == row.names(m)[i])]
}

m = as.table(as.matrix(m))

#' counts to binomial
c.binom = BradleyTerry2::countsToBinomial(m)

pp.model = BradleyTerry2::BTm(cbind(win1, win2),
                              player1, player2,
                              ~player, id = 'player',
                              data = c.binom)
#' update the model
pp.model.update = update(pp.model, refcat = 'breast')
summary(pp.model.update)

#' setting the 'weakest' as reference. In this case we are using liver;
#' occurred just once

#' now lets extract the model estimates and make a plot
plot.pp = as.data.frame(coefficients(pp.model.update))
plot.pp$sig = coef(summary(pp.model.update))[,4]
plot.pp$variable = row.names(plot.pp)
plot.pp$variable = substr(plot.pp$variable, start = 7, nchar(plot.pp$variable))
colnames(plot.pp)[1] = 'estimate'

ggplot(plot.pp, aes(x = reorder(variable, estimate), y = estimate)) +
  geom_point(col = ifelse(plot.pp$sig < 0.05, 'red', 'blue')) +
  coord_flip() +
  scale_y_continuous(
    sec.axis = dup_axis()) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  labs(x = '', y = 'BT model estimates')


#' Intervals based on quasi standard errors
pp.model.qv = qvcalc(BTabilities(pp.model.update))
plot(pp.model.qv, las = 2, cex.axis = 0.65)
abline(h = 1)


#' ML estimation
ML_Prostate = BT_ML(m)
#' this convergence actually to the same as the simply frequency of occurrence


###############################################################################
###############################################################################
#' PlackettLuce model
Prostate = Prostate[, 4:ncol(Prostate)]

#' delete columns with no entry
for(i in 4:26){
  print(i)
  if(all(is.na(Prostate[, i]))){
    Prostate[, i] = NULL
  }
}
Prostate[,25:ncol(Prostate)] = NULL
Prostate[,13:ncol(Prostate)] = NULL

#' unique sites
unique_sites_Prostate = as.character(unique(unlist(Prostate[, 1:ncol(Prostate)])))
unique_sites_Prostate = unique_sites_Prostate[!is.na(unique_sites_Prostate)]
names(unique_sites_Prostate) = seq(1, length(unique(unique_sites_Prostate)), 1)

#' recode the named characters into numeric values
Prostate[Prostate == "regional_lymph"] = 1L
Prostate[Prostate == "genital_male"] = 4L
Prostate[Prostate == "lung"] = 7L
Prostate[Prostate == "peritoneum"] = 10L
Prostate[Prostate == "pleura"] = 13L
Prostate[Prostate == "adrenal_gland"] = 16L
Prostate[Prostate == "breast"] = 19L
Prostate[Prostate == "mediastinum"] = 20L
Prostate[Prostate == "other"] = 2L
Prostate[Prostate == "dist_lymph"] = 5L
Prostate[Prostate == "liver"] = 8L
Prostate[Prostate == "bowel"] = 11L
Prostate[Prostate == "kidney"] = 14L
Prostate[Prostate == "skin"] = 17L
Prostate[Prostate == "head_and_neck"] = 21L
Prostate[Prostate == "bone"] = 3L
Prostate[Prostate == "bladder_or_urinary_tract"] = 6L
Prostate[Prostate == "lymph"] = 9L
Prostate[Prostate == "cns_brain"] = 12L
Prostate[Prostate == "peripheral_nervous_system"] = 15L
Prostate[Prostate == "biliary_tract"] = 18L

#' convert
x_matrix = as.matrix(Prostate, rownames.force = F)
mode(x_matrix) = 'numeric'
attr(x_matrix, which = 'sites') = unique_sites_Prostate

#' make the modeling:
x_ranked = as.rankings(x_matrix, 
                       input = 'orderings', 
                       items = attr(x_matrix, 'sites'))


#' average rank of metastatic site
avRank = apply(x_ranked, 2, function(x) print(x))
barplot(sort(apply(avRank, 2, function(x) mean(x[x>0])), decreasing = F),  
        las = 2, cex.axis = 0.85, cex.names = 0.65, ylab = 'average Rank')


#' modeling with PlacketLuce
model_rank = PlackettLuce(rankings = x_ranked, npseudo = 1)

sort(coef(model_rank))
summary(model_rank)



## hyper2 package
test = data_split[1:10, c('Rank1', 'Rank2', 'Rank3', 'Rank4')]
r1 = as.character(test[1, ])
r2 = as.character(test[2, ])
r3 = as.character(test[3, ])
r4 = as.character(test[4, ])
r5 = as.character(test[5, ])

H = rankvec_likelihood(r1)
H = H + rankvec_likelihood(r2)
H = H + rankvec_likelihood(r4)
H = H + rankvec_likelihood(r5)

maxp(H)

#####################################
for(i in 1:nrow(test)){
  row_i = as.character(test[i, ])
  row_i = row_i[!duplicated(row_i)]
  hy_single = rankvec_likelihood(row_i)
  #hy = hy_single
  hy = hy + hy_single
}

o = maxp(hy)
dotchart(o)

a = as.character(test[1,])
a = a[!duplicated(a)]



## test

test = Prostate[1:20, 1:10]
test_ranked = as.rankings(test, input = 'orderings')

b = apply(test_ranked, 2, function(x) print(x))
apply(b, 2, function(x) mean(x[x>0]))
mean(b[,3])
u = b[,3] 
u = as.numeric(u)
weighted.mean(u, w = c(1,1,1,1,1,1,1,0,1,1,1,0,0,1,1,1,1,1,1,1))
weighted.mean(u, w = c(rep(0, 16), 10, 0, 0, 0))

b[which(b[,2] == 0), 2] = NA

mean(b[,2], na.rm = T)
weighted.mean(b[, 5], w = rep(length(b[which(b[, 5] != 0), 5]) / nrow(b), nrow(b)))

mean(u)

