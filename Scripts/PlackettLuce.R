## Introduction of the PlackettLuce Model
## 
## 01/04/2022
## chris-kreitzer



install.packages('PlackettLuce')
library(PlackettLuce)

pudding = PlackettLuce::pudding

i_wins <- data.frame(Winner = pudding$i, Loser = pudding$j)
j_wins <- data.frame(Winner = pudding$j, Loser = pudding$i)

ties <- data.frame(Winner = asplit(pudding[c("i", "j")], 1),
                   Loser = rep(NA, 15))
R <- as.rankings(rbind(i_wins, j_wins, ties),
                 input = "orderings")

str(R)
head(unclass(R), 2)
