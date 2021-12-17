## Tropism. Pattern of metastatic dissemination;
## We start with Prostate Cancer (primaries)
## 


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

barplot(sort(table(sites_new$V1), decreasing = T), las = 2)
barplot(sort(table(sites_new$V2), decreasing = T), las = 2)
barplot(sort(table(sites_new$V3), decreasing = T), las = 2)

devtools::install_github("hturner/BradleyTerry2")


###############################################################################
###############################################################################
#' We start with the BT model 2
#' We try to reproduce the example provided on the WikiPage (https://en.wikipedia.org/wiki/Bradley%E2%80%93Terry_model)

wiki = matrix(c(0, 2, 0, 1, 3, 0, 5 , 0, 0, 3, 0, 1, 4, 0, 3, 0), ncol = 4, byrow = T)
colnames(wiki) = c('A', 'B', 'C', 'D')
row.names(wiki) = colnames(wiki)
wiki = as.table(wiki)

wiki_counts = BradleyTerry2::countsToBinomial(wiki)
wiki_model = BradleyTerry2::BTm(cbind(win1, win2), player1, player2, ~player, id = 'player', data = wiki_counts)


## Model formulation for BradleyTerry - ML estimation;
## note that this function requires a fully-connected network model otherwise the ML cannot be estimated

BT_ML = function(data){
  if(!inherits(x = data, what = 'table')) stop('Input data must be in table format', call. = F)
  if(!any(grepl(pattern = '*BradleyTerry2', x = search()))){
    require(BradleyTerry2)
  } else cat('BradleyTerry2 was loaded previously\n')
  
  #' modify the input data
  wins_all = rowSums(data, na.rm = T)
  data_counts = BradleyTerry2::countsToBinomial(data)
  priors = as.numeric(rep(1, ncol(data)))
  names(priors) = colnames(data)

  #' for-loop to start ML approach
  iteration = function(priors, wins_all, data_counts){
    j = c()
    for(i in names(wins_all)){
      sub = data_counts[which(data_counts$player1 == i | data_counts$player2 == i), ]
      sub$prior1 = NA
      sub$prior2 = NA
      for(row in 1:nrow(sub)){
        sub$prior1[row] = priors[which(sub$player1[row] == names(priors))]
        sub$prior2[row] = priors[which(sub$player2[row] == names(priors))]
      }
      
      #' calculate the ratios
      sub$ratio = (sub$win1 + sub$win2) / (sub$prior1 + sub$prior2)
      sub.ratio = sum(sub$ratio)
      p.prior = wins_all[[i]] / sub.ratio
      names(p.prior) = i
      j = c(j, p.prior)
    }
    
    #' post MLE normalization
    prior_sum = sum(j)
    post_p = j / prior_sum
    return(post_p)
  }
  
  #' iterate ML until convergence
  posterior_ML = data.frame()
  p0 = iteration(priors = priors, wins_all = wins_all, data_counts = data_counts)
  
  p1 = iteration(priors = p0, wins_all = wins_all, data_counts = data_counts)
  p2 = iteration(priors = p1, wins_all = wins_all, data_counts = data_counts)
  p3 = iteration(priors = p2, wins_all = wins_all, data_counts = data_counts)
  p4 = iteration(priors = p3, wins_all = wins_all, data_counts = data_counts)
  p5 = iteration(priors = p4, wins_all = wins_all, data_counts = data_counts)
  p6 = iteration(priors = p5, wins_all = wins_all, data_counts = data_counts)
  p7 = iteration(priors = p6, wins_all = wins_all, data_counts = data_counts)
  p8 = iteration(priors = p7, wins_all = wins_all, data_counts = data_counts)
  p9 = iteration(priors = p8, wins_all = wins_all, data_counts = data_counts)
  p10 = iteration(priors = p9, wins_all = wins_all, data_counts = data_counts)
  p11 = iteration(priors = p10, wins_all = wins_all, data_counts = data_counts)
  p12 = iteration(priors = p11, wins_all = wins_all, data_counts = data_counts)
  p13 = iteration(priors = p12, wins_all = wins_all, data_counts = data_counts)
  p14 = iteration(priors = p13, wins_all = wins_all, data_counts = data_counts)
  p15 = iteration(priors = p14, wins_all = wins_all, data_counts = data_counts)
  p16 = iteration(priors = p15, wins_all = wins_all, data_counts = data_counts)
  p17 = iteration(priors = p16, wins_all = wins_all, data_counts = data_counts)
  p18 = iteration(priors = p17, wins_all = wins_all, data_counts = data_counts)
  p19 = iteration(priors = p18, wins_all = wins_all, data_counts = data_counts)
  p20 = iteration(priors = p19, wins_all = wins_all, data_counts = data_counts)
  
  posterior_ML = rbind(p0, p1, p2, p3, p4, p5, p6,
                       p7, p8, p9, p10, p11, p12,
                       p13, p14, p15, p16, p17,
                       p18, p19, p20)
  
  return(posterior_ML)
  
}

x = BT_ML(data = test)
x = data.frame(x, row.names = dimnames(x) [[1]])


#' Visualization; where the ML converges to 1
#' looking into the behavior of regional lymph nodes
plot(x$regio,  col = 'red', lwd = 1.5, pch = 1, ylim = c(0.08, 0.150), 
     yaxt = 'n', ylab = 'estimate', xlab = 'iteration', main = 'MLE for regional lymph nodes')
lines(x$regio)
points(x$regio,  col = 'red', lwd = 1.5, pch = 1, ylim = c(0.08, 0.150), 
     yaxt = 'n', ylab = 'estimate', xlab = 'iteration', main = 'MLE for regional lymph nodes')
box(lwd = 2)

#' looking into bone
plot(x$bone,  col = 'blue', lwd = 1.5, pch = 1, ylim = c(0.2, 0.350), 
     yaxt = 'n', ylab = 'estimate', xlab = 'iteration', main = 'MLE for bone metastasis')
lines(x$bone)
points(x$bone,  col = 'blue', lwd = 1.5, pch = 1, ylim = c(0.2, 0.35), 
       yaxt = 'n', ylab = 'estimate', xlab = 'iteration', main = 'MLE for bone metastasis')
box(lwd = 2)


points(x$bone)
lines(x$bladder)
points(x$bladder, col = 'red')
lines(x$other)
points(x$other, col = 'green')
lines(x$dist)
points(x$dist, col = 'blue')


plot.new()












#' test example
# test = matrix(c(0, 1, 1,1,1,3,0,3,3,3,2,2,0,2,2,1,1,1,0,1,3,3,3,3,0), ncol = 5, byrow = T)
# colnames(test) = c('regio', 'bone', 'bladder', 'other', 'dist')
# row.names(test) = colnames(test)
# test = as.table(test)
# priors = c(1, 1, 1, 1, 1)
# names(priors) = c('bone', 'regio', 'bladder', 'other', 'dist')


