## What happens if we detect two metastatic sites at the same time point?
## look into an example;

setup(working.path = '~/Documents/GitHub/Tropism-Timing/')
clean()

# install.packages('readxl')
library(readxl)
# install.packages("fitdistrplus")
library(fitdistrplus)
library(tidyr)
# install.packages('BradleyTerry2')
library(BradleyTerry2)
library(ggplot2)
library(qvcalc)

## Data
data = read_excel('Data/tableS2_data.xlsx', skip = 2)


# starting with Prostate first
data_prostate = data[which(data$cancer_type == 'prostate cancer'), ]
data_PP = data_prostate[which(data_prostate$sample_type == 'Primary'), ]
data_PP_death = data_PP[which(data_PP$os_status == 'dead'), ]

# sites
ncols_max = max(stringr::str_count(sites_PP$met_event_age, '\\/'), na.rm = T) + 1
sites_PP = data_PP_death[, c('met_site_mapped', 'met_event_age')]
column_names = paste0("Time", 1:ncols_max)

#' split the dataframe
data_pp_split = tidyr::separate(sites_PP, 
             col = met_event_age, 
             sep = '\\/', 
             into = column_names,
             remove = F)

View(data_pp_split)


df %>%
  as_tibble() %>%
  mutate(dup = pmap_dbl(list(V1, V2, V3), ~ n_distinct(c(...)))) %>%
  filter(dup == 3) %>%
  select(-dup)

x = data_pp_split %>% 
  as_tibble() %>%
  mutate(dup = pmap_dbl(list(Time1, Time2, Time3), ~ n_distinct(c(...))))
View(x)



head(data_pp_split)
dup_rows = apply(data_pp_split[, 3:ncol(data_pp_split)], 1, FUN = function(x) ifelse(max(table(x, na.rm = TRUE)) > 1 ,TRUE, FALSE))
u = data_pp_split[1, 3:ncol(data_pp_split)]
apply(u, 1, FUN = function(x) ifelse(max(table(x)) > 1, TRUE, FALSE))
View(u)





x = data.frame(one = c('a', 'a', 'x'),
               two = c('a', 'a', 'c'),
               third = c('a', 'b', 'y'))

x = apply(x, 1, FUN = function(x) table(x))
u = data.table::rbindlist(x)

j = do.call(rbind.data.frame, x)














