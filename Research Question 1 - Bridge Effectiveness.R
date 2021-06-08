# Research Question: Implication of "Bridge"
library(ggplot2)
library(tidyverse)
library(knitr)
require(tidyverse)
require(knitr)
install.packages("broom")
install.packages("ggpubr")
library(ggpubr)

###### SET UP ######

# run 100 simulations at each bridge size with four different scales of communitities: 10, 50, 100, 500
# the corresponding generation time should be 10 turnover times
# so for n = 10, gen = 100, turn over time is 100/10 = 10 times
# n   turnover   gen needed
# 10  10         100
# 50  10         500
# 100 10         1000
# 500 10         5000

###### Simulations - Species Count ######

# create a dataframe to store results
Q1 <- NULL

# N = 10
n=10
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(SPECMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q1 <- rbind(Q1,result)
}

# N = 50
n=50
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(SPECMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q1 <- rbind(Q1,result)
}

# N = 100
n=100
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(SPECMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q1 <- rbind(Q1,result)
}

# N = 500
n=500
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(SPECMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q1 <- rbind(Q1,result)
}

# clean dataframe
colnames(Q1) <- c("n","bridge","species_count")
rownames(Q1) <- NULL

Q1

write.csv(Q1,"H:/ENGE817/Report/Q1newversion.csv",row.names=FALSE)


###### Analysis - Species Count ######

spec_bridge <- read.csv("/cloud/project/Q1_Spec_Count.csv")

spec_bridgep <- pivot_wider(spec_bridge,names_from=bridge,values_from=species_count)
spec_bridgep %>% 
  kable(caption="Expected Species Count with Different Bridge Sizes")



###### Simulation - BC Index #####

# create a dataframe to store results
Q1 <- NULL

# N = 10
n=10
gen=n*10
for (i in seq(from=0,to=0.1,by=0.1)){
  df <- t(BCMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q1 <- rbind(Q1,result)
}

# N = 50
n=50
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(BCMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q1 <- rbind(Q1,result)
}

# N = 100
n=100
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(BCMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q1 <- rbind(Q1,result)
}

# N = 500
n=500
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(BCMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q1 <- rbind(Q1,result)
}

# clean dataframe
colnames(Q1) <- c("n","bridge","bc_index")
rownames(Q1) <- NULL

Q1

write.csv(Q1,"cloud/project/Q1_BC_Index.csv",row.names=FALSE)


###### Analysis - BC Index ######

BC_bridge <- read.csv("/cloud/project/Q1_BC_Index.csv")
BC_bridgep <- pivot_wider(BC_bridge,names_from=bridge,values_from=bc_index)
BC_bridgep %>% 
  kable(caption="Expected Bray Curtis Dissimilarity After Ten Turnovers", digit=2)


##### ANCOVA test on bridge / n / BC Index #####
BC_ANCOVA <- read.csv("/cloud/project/Q1_BC_Ancova.csv")

result <- summary(aov(bc_index ~ n*bridge, data = BC_ANCOVA))
result %>% 
  tidy() %>%
  kable(caption="Linear Regression Analysis")


###### PLOT Species Count vs BC Index ####

plot1 <- spec_bridge %>% 
  ggplot(aes(x=n,y=species_count,group=bridge,color=bridge)) + geom_line() +
  ggtitle("Prediction of Species Count across Bridge Sizes") +
  theme(plot.title = element_text(hjust = 0.5))

plot2 <- BC_bridge %>% 
  ggplot(aes(x=n,y=bc_index,group=bridge,color=bridge)) + geom_line() +
  ggtitle("Prediction of BC Index across Bridge Sizes") +
  theme(plot.title = element_text(hjust = 0.5))

ggarrange(plot1, plot2, ncol = 2, nrow = 1)

# END OF ANALYSIS