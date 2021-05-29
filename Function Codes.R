#---------------------------------------#
#           FUNCTION: island()          #
#---------------------------------------#

island <- function(n=10, b=0.5, g=100){
  
  # package needed
  require(tidyverse)
  
  # create a metacommunity with N singletons
  prefix <- "J"
  suffix <- seq(1:n)
  gen <- paste(prefix, suffix, sep = "") 
  mc <- data.frame(gen) 
  
  # randomly & equally assign the singletons to two islands
  rangi <- mc %>% slice_sample(prop = .5, replace = FALSE)
  motu <- mc %>%  filter(!gen %in% rangi$gen)
  
  # duplicate dataframes for calculations
  r <- rangi
  m <- motu
  
  # create another dataframe to record the biodiversity status at each step
  rm_step <- rbind(rangi,motu)
  
  # set community size
  j = nrow(mc)
  
  # set probability of individual crossing b 
  prob = ((j-1)^-1)*(j/2)*b
  
  # set empty vectors for for loop
  island_die <- NULL
  island_fill <- NULL
  
  # loop begins
  for (i in 1:g) {
    island_die[i] <- sample(0:1,1)                   # choose an island, a site on this island dies
    
    if (island_die[i] == 1){                         # 1 = island r
      island_fill[i] <- rbinom(1,size=1,prob=prob) # choose an island, an individual on this island fills the site 
      if (island_fill[i] == 1){                      # 1 = counter island - m
        r[sample(nrow(r),1),] <- sample(m$gen,1)      
      }
      else {                                         # 0 = same island - r
        r[sample(nrow(r),1),] <- sample(r$gen,1)    
      }
      singlestep <- rbind(r,m)                       # abundance change at one single step
      rm_step <- cbind(rm_step, singlestep)        
    }
    
    else {                                           # 0 = island m
      island_fill[i] <- rbinom(1,size=1,prob=prob) # choose an island, an individual on this island fills the site
      if (island_fill[i] == 1){                      # 1 = counter island - r
        m[sample(nrow(r),1),] <- sample(r$gen,1)     
      }
      else {                                         # 0 = same island - m
        m[sample(nrow(r),1),] <- sample(m$gen,1)    
      }
      singlestep <- rbind(r,m)                      # abundance change at one single step
      rm_step <- cbind(rm_step, singlestep)    
    }
  }
  
  # rename columns
  gen <- length(rm_step)-1                        
  colnames(rm_step) <- paste(colnames(rm_step),0:gen,sep="")  
  
  rm_step <- t(rm_step) # transpose data
  
  colnames(rm_step) <- paste('island',1:n,sep="") # assign
  colnames(rm_step)[1:(n/2)] <- "rangi"
  colnames(rm_step)[(n/2+1):n] <- "motu"
  
  # calculate biodiversity
  species_count <- NULL
  rangi_count <- NULL
  motu_count <-  NULL
  BrayCurtis <- NULL
  
  for (i in 0:g+1){
    species_count[i] <- (length(unique(rm_step[i,])))
    rangi_count[i] <- (length(unique(rm_step[i,1:(n/2)])))
    motu_count[i] <- (length(unique(rm_step[i,((n/2)+1):n])))
  }
  
  rm_step <- as.data.frame(cbind(rm_step,species_count,rangi_count,motu_count))
  rm_step$species_count <- as.numeric(rm_step$species_count)
  rm_step$rangi_count <- as.numeric(rm_step$rangi_count)
  rm_step$motu_count <- as.numeric(rm_step$motu_count)
  
  df <- rm_step
  result_set <- NULL
  
  for (i in 1:nrow(df)) {
    
    t1 <- data.frame(table(t(df[i,1:(n/2)])))
    t2 <- data.frame(table(t(df[i,((n/2)+1):n])))
    
    temp_set <- NULL
    
    for (j in 1:nrow(t1)) {
      
      if(sum(as.character(t1$Var1[j]) == as.character(t2$Var1[1:nrow(t2)]))){
        
        frequency_in_t1 <- t1$Freq[j] # count freq on island A
        
        frequency_in_t2 <- t2 %>%
          filter(Var1 == as.character(t1$Var1[j])) %>%
          pull(Freq)                # count freq on island B
        
        # comparison
        if(frequency_in_t1 >= frequency_in_t2){
          comparison_result <- frequency_in_t2
        }
        else{
          comparison_result <- frequency_in_t1
        }
        temp_set <- c(temp_set, comparison_result) 
      }
    }
    result_set <- c(result_set, sum(temp_set)) # result_set stays here, length is correct
  }
  
  BC <- 1 - (2 * result_set) / n
  
  rm_step$BrayCurtis <- round(BC,2)
  
  return(rm_step)
}

if(FALSE){
  print("Error: please examine your codes")
}



#----------------------------------------#
#    FUNCTION: bio_max(), bio_min()      #
#----------------------------------------#

# create a function to identify the maximum generation to reach any % of diversity
bio_max <- function(df,column,number){
  df[df[,column]>=number,][nrow(df[df[,column]>=number,]),]
}

# create a function to identify the minimum generation to reach any % of diversity
bio_min <- function(df,column,number){
  df[df[,column]<=number,][1,]
}


#-------------------------------------------------#
#               FUNCTION: GENMC()                 #
#-------------------------------------------------#


GENMC <- function(func,n,b,g,column,d,t=100){
  
  diversegen <- NULL
  
  for (i in 1:t){
    diversegen[i] <- rownames(func(island(n,b,g),column,d)) # check generation
  }                       #bio_max(island(20,1,100),"BrayCurtis",0.9)
  as.data.frame(table(diversegen))
}



#-------------------------------------------------#
#                FUNCTION: BCMC()                 #
#-------------------------------------------------#

BCMC <- function(n=n,b=b,g=g,t=t){
  simulation <- NULL
  
  for (i in 1:t){
    step <- as.data.frame(island(n,b,g)[,"BrayCurtis"])
    simulation[i] <- step
  }
  simulation <- as.data.frame(simulation)
  
  # rename columns
  path <- ncol(simulation)                        
  colnames(simulation) <- paste("path",1:path,sep="")
  
  return(simulation)
}



#-------------------------------------------------#
#             FUNCTION: SPECMC()                 #
#-------------------------------------------------#

SPECMC <- function(n=n,b=b,g=g,t=t){
  simulation <- NULL
  
  for (i in 1:t){
    step <- as.data.frame(island(n,b,g)[,"species_count"])
    simulation[i] <- step
  }
  simulation <- as.data.frame(simulation)
  
  # rename columns
  path <- ncol(simulation)                        
  colnames(simulation) <- paste("path",1:path,sep="")
  
  return(simulation)
}


#################################################################################################
#################################################################################################
#################################################################################################
#################################################################################################



#------------------------------------------#
#           FUNCTION: BCplot()             #
#------------------------------------------#

BCplot <- function(n=n,b=b,g=g){
  for (i in 1){
    parameter <- island(n,b,g)
    plot(parameter$BrayCurtis,type='l',lwd=2, xlab="generation", ylab="Bray Curtis Dissimilarity",ylim=c(0,1))
    parameter <- island(n,b,g)
    points(parameter$BrayCurtis,type='l',lwd=2,col='darkgoldenrod1')
    parameter <- island(n,b,g)
    points(parameter$BrayCurtis,type='l',lwd=2,col='aquamarine3')
    parameter <- island(n,b,g)
    points(parameter$BrayCurtis,type='l',lwd=2,col='dodgerblue')
    parameter <- island(n,b,g)
    points(parameter$BrayCurtis,type='l',lwd=2,col='coral1')
    abline(h=c(0.1,0.5,0.9),lty='dashed', col="gray")
    abline(v=c(g*0.1,g*0.5,g*0.9),lty='dashed', col="grey")
  }
}

#--------------------------------------------#
#           FUNCTION: specplot()             #
#--------------------------------------------#

specplot <- function(n=n,b=b,g=g){
  for (i in 1){
    parameter <- island(n,b,g)
    plot(parameter$species_count,type='l',lwd=2, xlab="generation", ylab="species count",ylim=c(0,n))
    parameter <- island(n,b,g)
    points(parameter$species_count,type='l',lwd=2,col='lightskyblue4')
    parameter <- island(n,b,g)
    points(parameter$species_count,type='l',lwd=2,col='midnightblue')
    parameter <- island(n,b,g)
    points(parameter$species_count,type='l',lwd=2,col='mediumvioletred')
    parameter <- island(n,b,g)
    points(parameter$species_count,type='l',lwd=2,col='orange2')
    abline(h=c(n*0.1,n*0.5,n*0.9),lty='dashed', col="gray")
    abline(v=c(g*0.1,g*0.5,g*0.9),lty='dashed', col="grey")
  }
}

