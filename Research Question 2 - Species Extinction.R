install.packages("ggplot2")
library(ggplot2)
library(tidyverse)
library(dplyr)
library(knitr)
library(ggpubr)

#### Research Question: HOW LONG DOES IT TAKE TO LOSE 50% BIODIVERSITY? ####

##### FIRST PLOT EXAMINATION #####

par(mfrow=c(2,2))
drawplot(10,1,500)
drawplot(50,1,500)
drawplot(100,1, 500)
drawplot(200,1,500)


#### STEP 1: estimate the timestep to lose 50% species ###  
#### 1000 simulations each community for n 10 to 500
#### 100 simulations each community 

# N 10 to 50 by 10
for (i in seq(from=10,to=50,by=10)){
  max = GENMC(bio_min,i,1,100,"species_count",i*0.5,1000) %>% 
  arrange(desc(Freq)) %>% slice(1)
  print(i)
  print(max)
}

# N 60 to 100 by 10
for (i in seq(from=60,to=100,by=10)){
  max = GENMC(bio_min,i,1,200,"species_count",i*0.5,1000) %>% 
    arrange(desc(Freq)) %>% slice(1)
  print(i)
  print(max)
}

# N 100 to 500 by 100
for (i in seq(from=100,to=500,by=100)){
  max = GENMC(bio_min,i,1,1000,"species_count",i*0.5,1000) %>% 
    arrange(desc(Freq)) %>% slice(1)
  print(i)
  print(max)
}

# due to computational power limit, community above 500 individuals are examined by the mean of 10 trails

# N = 1000
N1000 <- GENMC(bio_min,1000,1,2000,"species_count",1000*0.5,10)
N1000
(1004+1011+1074+1081+948+962+970+971+990+996)/10

# N = 2000
N2000 <- GENMC(bio_min,2000,1,4000,"species_count",2000*0.5,10)
N2000
(1911+1973+1995+2001+2004+2013+2030+2062+2070+2144)/10

# N = 5000
N5000 <- GENMC(bio_min,5000,1,6000,"species_count",5000*0.5,10)
N5000
(4907+4924+4955+4958+4976+4981+5001+5008+5031+5090)/10

# N = 10000
N10000 <- GENMC(bio_min,10000,1,20000,"species_count",10000*0.5,10)
N10000 %>% arrange(desc(Freq)) %>% slice(1)
N10000
(10011+10039+10043+10150+10203+9793+9817+9904+9933+9978)/10


######## PLOT turnover time vs time step #######

spec_loss <-  read.csv("/cloud/project/specloss.csv")
colnames(spec_loss) <- c("n","bridge","timestep","species_left","turnover_time")
#table
spec_loss

# plot
plotq2_1 <- ggplot(spec_loss, aes(x=n,y=turnover_time)) + 
  geom_line(color="blueviolet", size=1)
plotq2_2 <- ggplot(spec_loss, aes(x=n,y=timestep)) +
  geom_line(color="coral1", size=1.5)

figure <- ggarrange(plotq2_1, plotq2_2, ncol = 2, nrow = 1.5) 

annotate_figure(figure,
                top = text_grob("Turnover Time vs Timestep of Reaching 50% Species Loss", 
                                color = "black", face = "bold", size = 20))


######### DOES SPECIES EXTINCT SLOWER WITHIN THE SMALL COMMUNITIES??? ##########

xn <- c(10,50,100,500)
spec <- c(1.402,4.987,9.516,45.688)
gent <- c(100,500,1000,5000)
turnover <- c(10,10,10,10)
df <- data.frame(xn,spec,gent,turnover)

pct <- df$spec/df$xn
df <- cbind(df,pct)
df

# after 10 turnover times, it's expected to have a 9% left of the species count for community larger than 50
# community n = 10 since to retain more species than the larger ones
# how about 20?

Q2 <- NULL

# N = 20
n=20
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(SPECMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q2 <- rbind(Q2,result)
}

mean(Q2[,3])/200
# 0.01165909

Q2a <- NULL

# N = 6
n=6
gen=n*10
for (i in seq(from=0,to=1,by=0.1)){
  df <- t(SPECMC(n,i,gen,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  Q2a <- rbind(Q2a,result)
}

mean(Q2a[,3])/60
#0.02066667

# how about i repeat the experiment with different turnover times?
# start from a midium size community:

# N 100 with different gens
# turnover gen
# 1        100
# 2        200
# 3        300
# ...      ...
# 10       1000

# count the speices left with turnover 1 to 10 times

turn_test <- NULL

# N = 10
n=10
for (i in seq(from=10,to=100,by=10)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 20
n=20
for (i in seq(from=20,to=200,by=20)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 50
n=50
for (i in seq(from=50,to=500,by=50)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 100
n=100
for (i in seq(from=100,to=1000,by=100)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 500
n=500
for (i in seq(from=500,to=5000,by=500)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 1000
n=1000
for (i in seq(from=1000,to=10000,by=1000)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

##### Plot turnover against species count ####

turnover <-  read.csv("/cloud/project/q2turnover_examine.csv")
colnames(turnover) <- c("n","gen","turnover","species_left","pct")
#table
turnover$pctloss = 1- turnover$pct


plot2 <- turnover %>% 
  ggplot(aes(x=turnover,y=pctloss,group=n)) + geom_line(size=1.5,color="blueviolet") + 
  geom_hline(aes(yintercept=0.9),linetype = 'dashed', color='coral1',size=1) +
  facet_wrap(~n) +
  ggtitle("Prediction in Species Loss by Turnover Time") +
  theme(plot.title = element_text(hjust = 0.5,size=20), legend.position = "none")
plot2


#### STEP 2: estimate the timestep to lose 50% dissimilarity ###  
#### 100 simulations for each bridge ##
### N = 20

# testing a n=20 community at 1 turnover
turn_test <- NULL
n=20
g=20
for (i in seq(from=0.1,to=1,by=0.1)){
  df <- t(BCMC(n,i,g,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# testing a n=20 community at 10 turnovers
n=10
g=200
for (i in seq(from=0.1,to=1,by=0.1)){
  df <- t(BCMC(n,i,g,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# testing a n=20 community at 10 turnovers
n=20
g=200
for (i in seq(from=0.1,to=1,by=0.1)){
  df <- t(BCMC(n,i,g,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# examine results
BCturn <-  read.csv("/cloud/project/q2_bc_turnover.csv")

plot(BCturn$turnover,BCturn$BC_index)

plotBC <- BCturn %>% 
  ggplot(aes(x=turnover,y=BC_index,group=b,colour=b)) + geom_line(size=0.5) + 
  geom_hline(aes(yintercept=0.5),linetype = 'dashed', color='coral1',size=1) +
  ggtitle("Prediction in BC Dissimilarity Change by Bridge Size and Turnover Time") +
  theme(plot.title = element_text(hjust = 0.5))
plotBC

# END OF ANALYSIS