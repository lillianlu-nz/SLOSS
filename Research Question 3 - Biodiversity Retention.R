library(ggplot2)
library(tidyverse)
library(dplyr)
library(knitr)
library(ggpubr)


#### STEP 1: estimate the timestep to hold 75% species ###  
#### 1000 simulations each community for n 10 to 500
#### 100 simulations each community 

# N 10 to 50 by 10
for (i in seq(from=10,to=50,by=10)){
  max = GENMC(bio_max,i,1,100,"species_count",i*0.75,1000) %>% 
  arrange(desc(Freq)) %>% slice(1)
  print(i)
  print(max)
}

# N 60 to 100 by 10
for (i in seq(from=60,to=100,by=10)){
  max = GENMC(bio_max,i,1,200,"species_count",i*0.75,1000) %>% 
    arrange(desc(Freq)) %>% slice(1)
  print(i)
  print(max)
}

# N 100 to 500 by 100
for (i in seq(from=100,to=500,by=100)){
  max = GENMC(bio_max,i,1,200,"species_count",i*0.75,1000) %>% 
    arrange(desc(Freq)) %>% slice(1)
  print(i)
  print(max)
}

# due to computational power limit, community above 500 individuals are examined by the mean of 10 trails

# N = 1000
N1000 <- GENMC(bio_max,1000,1,2000,"species_count",1000*0.75,10)
N1000
(315+324+326+315+328+336+337+337+339+340)/10

# N = 2000
N2000 <- GENMC(bio_max,2000,1,4000,"species_count",2000*0.75,10)
N2000
(650+662+672+676+678+679+694+696+679+700)/10

# N = 5000
N5000 <- GENMC(bio_max,5000,1,6000,"species_count",5000*0.75,10)
N5000
(1629+1648+1651+1666+1672+1681+1681+1682+1683+1684)/10

# N = 10000
N10000 <- GENMC(bio_max,10000,1,20000,"species_count",10000*0.75,10)
N10000 %>% arrange(desc(Freq)) %>% slice(1)
N10000
(3277+3300+3322+3325+3330+3331+3344+3351+3361+3395)/10


######## PLOT turnover time vs time step #######

spec_retain <-  read.csv("/cloud/project/specretain.csv")
colnames(spec_retain) <- c("n","bridge","timestep","species_main","turnover_time")
#table
spec_retain

# plot
plotq3_1 <- ggplot(spec_retain, aes(x=n,y=turnover_time)) + 
  geom_line(color="blueviolet", size=1)
plotq3_2 <- ggplot(spec_retain, aes(x=n,y=timestep)) +
  geom_line(color="coral1", size=1.5)

figure <- ggarrange(plotq3_1, plotq3_2, ncol = 2, nrow = 1.5) 

annotate_figure(figure,
                top = text_grob("Turnover Time vs Timestep of Retaining 75% Biodiversity", 
                                color = "black", face = "bold", size = 20))


##### we looked at species loss, let's check species retention #####
# count the speices left with turnover 1 to 10 times
# for ONE turnover time, we expect around 50% loss #
# let's test the speed between 0 to 1

turn_test <- NULL

# N = 10
n=10
for (i in seq(from=1,to=10,by=1)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 20
n=20
for (i in seq(from=2,to=20,by=2)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 50
n=50
for (i in seq(from=5,to=50,by=5)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 100
n=100
for (i in seq(from=10,to=100,by=10)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 500
n=500
for (i in seq(from=50,to=500,by=50)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

# N = 1000
n=1000
for (i in seq(from=100,to=1000,by=100)){
  df <- t(SPECMC(n,1,i,100))
  result <- c(n,i,mean(df[,ncol(df)]))
  turn_test <- rbind(turn_test,result)
}

colnames(turn_test) <- c("n","gen","species_left")
write.csv(turn_test, "H:/ENGE817/Report/q3turnover_examine.csv", row.names = TRUE)


##### Plot turnover against species count ####

turnover <-  read.csv("/cloud/project/q3turnover_examine.csv")
colnames(turnover) <- c("n","gen","turnover","species_left","pct","change")
#table
turnover

plot3 <- turnover %>% 
  ggplot(aes(x=turnover,y=pct,group=n)) + geom_line(size=1.5,color="dodgerblue1") + 
  geom_hline(aes(yintercept=0.75),linetype = 'dashed', color='coral1',size=1) +
  facet_wrap(~n) +
  ggtitle("Prediction in Species Retention by Turnover Time") +
  theme(plot.title = element_text(hjust = 0.5,size=20), legend.position = "none")

plot3

#### variation of change by community size ####
q3table <- turnover %>%
  group_by(n) %>%
  summarise(mean = mean(change), variation = var(change))

q3table %>%
  kable(caption="Change Statistics of Different Communities")

# END OF ANALYSIS
