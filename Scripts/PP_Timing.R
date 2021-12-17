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
library(ggplot2)
library(qvcalc)

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

m = as.table(m)
#' counts to binomial
c.binom = BradleyTerry2::countsToBinomial(m)

pp.model = BradleyTerry2::BTm(cbind(win1, win2),
                              player1, player2,
                              ~player, id = 'player',
                              data = c.binom)
#' update the model
pp.model.update = update(pp.model, refcat = 'liver')
summary(pp.model.update)

#' setting the 'weakest' as reference. In this case we are using liver;
#' occurred just once

#' now lets extract the model estimates and make a plot
plot.pp = as.data.frame(coefficients(pp.model.update))
plot.pp$sig = coef(summary(pp.model.update))[,4]
plot.pp$variable = row.names(plot.pp)
colnames(plot.pp)[1] = 'estimate'

ggplot(plot.pp, aes(x = reorder(variable, estimate), y = estimate)) +
  geom_point() +
  coord_flip() +
  scale_y_continuous(
                     sec.axis = dup_axis()) +
  geom_hline(yintercept = 1, linetype = 'dashed') +
  labs(x = 'Metastatic sites', y = 'BT model estimates')
  

#' Intervals based on quasi standard errors
pp.model.qv = qvcalc(BTabilities(pp.model.update))
plot(pp.model.qv, las = 2)
abline(h = 1)


#' ML estimation of strength parameters; among each other 
x = BT_ML(data = m) #' defined in Example.BTmodel.R script
x$iteration = factor(x$iteration, 
                     levels = c('p0', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 
                                'p7', 'p8', 'p9', 'p10', 'p11', 'p12', 
                                'p13', 'p14', 'p15', 'p16', 'p17', 'p18', 
                                'p19', 'p20'))

#' Visualization
ggplot(x, aes(x = iteration)) +
  geom_point(aes(y = regional_lymph), col = 'red') +
  geom_line(aes(y = regional_lymph, group = 1), col = 'red') +
  geom_point(aes(y = bone), col = 'blue') +
  geom_line(aes(y = bone, group = 1), col = 'blue') +
  geom_point(aes(y = bladder_or_urinary_tract), col = 'green') +
  geom_line(aes(y = bladder_or_urinary_tract, group = 1), col = 'green') +
  geom_point(aes(y = other), col = 'brown') +
  geom_line(aes(y = other, group = 1), col = 'brown') +
  geom_point(aes(y = dist_lymph), col = 'grey') +
  geom_line(aes(y = dist_lymph, group = 1), col = 'grey') +
  geom_point(aes(y = genital_male), col = 'yellow') +
  geom_line(aes(y = genital_male, group = 1), col = 'yellow') +
  theme_bw() +
  labs(x = '', y = 'ML estimate')





