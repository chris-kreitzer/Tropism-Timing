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
  } else cat('BradleyTerry2 was loaded previously')
  
  #' modify the input data
  wins_all = rowSums(data, na.rm = T)
  data_counts = BradleyTerry2::countsToBinomial(data)
  
}

BT_ML(data = wiki)


a = rowSums(wiki)
iteration = function(priors){
  j = vector(mode = 'numeric')
  for(i in names(a)){
    sub = wiki_counts[which(wiki_counts$player1 == i | wiki_counts$player2 == i), ]
    sub$prior1[which(sub$player1 %in% names(priors))] = priors[which(names(priors) %in% sub$player1)]
    sub$prior2[which(sub$player2 %in% names(priors))] = priors[which(names(priors) %in% sub$player2)]

    sub$ratio = (sub$win1 + sub$win2) / (sub$prior1 + sub$prior2)
    sub.ratio = sum(sub$ratio)
    pi = a[[i]] / sub.ratio
    names(pi) = i
    print(sub)
    j[i] = pi
  }
  prior = sum(j)
  post = j / prior
  return(post)
}

p1 = iteration(priors = b)

p1$prior1[which(p1$player1 %in% names(b))] = b[which(names(b) %in% p1$player1)]
p1$prior2[which(p1$player2 %in% names(b))] = b[which(names(b) %in% p1$player2)]
b
str(p1)
str(p1)

b[which(names(b) == 'D')][[1]]
str(b)
names(b)
names(a)
names(b)




b = c(1,1,1,1)
names(b) = c('A', 'B', 'C', 'D')
b

b[which(names(b) == 'A')][[1]]
a
b




iter = 0
while (iter < 20) {

  p1 = iteration(priors = 1))
  p2 = iteration(priors = p1)
  iter = iter + 1
  
}

p = iteration()
sum(p)

names(p)




