# define a function to do replication
# this function initializes the population, pop, an array to follow
# the presence of an allele. In this case we have N elements
# one of then has an allele active.

# The presence of the allele is marked by a 1 and 0 means no allele.
# In this way it's very easy to know how many people has the allele,
# We only have to sum the pop variable

# Parameters: numreplications
#             N: Num. of individuals, of the population
#             days  Days or generations during which to evolve
#             botlleneck  
#                         1 -> no botellneck
#                         2 -> botellneck 1/2  -> 0.5
#                         4 -> botellneck 1/4  -> 0.25

# Returns the result of the replications
replication <- function(numreplications, N, days, bottleneck=1, delta_days=1) {
  # We initialise the list in which we will return the result 
  # of the replications.
  mresult = list() 
  
  # We initialise the vector with which we will track the population.
  # we will mark 1 when the allele is present and 0 if it is not.
  pop<-rep(0,N)                         # Population Vector

  #The initial condition of the exercise indicates that one and only one 
  #has the active allele and that there is no mutation of any kind.
  #It is only multiplied or extinguished by normal birth/death.

  # We activate an allele
  mutate <- sample(1:N, 1)
  pop[mutate] <- 1
  
  # For each of the replications we wish to carry out, we are going to simulate 
  # neutral evolution.To do this we have a function that will "evolve"
  # this population with a single individual with the allele 
  # The result of each replication shall be stored in the mresult variable
  for (i in 1:numreplications){
    resultado <- neutral_evolution(pop, N, days, bottleneck, delta_days)
    mresult <- rbind(mresult, list(fixed = resultado$fixed, day=resultado$day, S1=resultado$S1))
  }
  
  # Once the replications have been completed, we return the results.
  return (mresult)
}

# define a function to do neutral_evolution
neutral_evolution <- function(pop, N, days, bottleneck=1, delta_days=1) {
  obs <- 0
  result = list(
    fixed = 0,
    day = 0,
    S1 = list()
  )
  


  # the algorithm (Evolution)
  # We choose the algorithm used in class, it has been modified 
  #to add the bottleneck
  # different versions have been tested
  #   1 -> Perform a single bottleneck in half of the generations in our case 
  #        5000 (No good, very rarely the allele is active at that time, or has 
  #        become extinct or fixed much earlier).
  #   2 -> The bottleneck occurs in generation 50. As we do not know the status 
  #        of the allele at that time, nothing significant 
  #        has been detected either.
  #   3 -> The bottleneck occurs when the allele is present in half of the 
  #        population, in this case slight deviations from the normal model are 
  #        found due to the randomness of selection. 
  #   4 -> As the alleleo usually dies out very soon, the bottleneck was chosen 
  #        for all generations. The sensation is that as the bottleneck 
  #        increases (greater population decline) the system behaves strangely, 
  #        maintaining Pfix and decreasing Tfix.
  
  for (day in 1:days){
    #  reproduction 

    # In this approach we choose to apply the bottleneck in all generations. 
    # we choose to reproduce the 50% of the population that will keep 
    # their alleles
    # choose N/(2*bottleneck) individuals among the population
    offspring<-sample(1:N, ceiling(N/(2*bottleneck)), replace=TRUE)     
    # Here we will apply the bottleneck
    
    # We make sure to keep the population stable
    # N/(rateoffspring*bn) can give decimals with small populations, so we check 
    # for missing or excess individuals and add or remove accordingly.
    pop<-rep(pop[offspring], 2*bottleneck)       # update pop matrix
    if (length(pop) > N) {
      pop <- pop[1:N]
    } else if (length(pop) < N) {
      pop <- c(pop, rep(0, N - length(pop)))
    }

    # We already have a new generation, we count the number of individuals
    # with the active allele. 
    # We add the number of alleles in result$S1 which is a vector
    # that will contain the count of active alleles day by day, generation by generation.
    #  data retrieving
    if (day %in% seq(1,days,by=delta_days)) {     # to save memory and time, don't record all generations 
      obs<-obs+1
      #result$S1[obs]<-sum(result$pop)
      result$S1 <- append(result$S1, sum(pop))
      #print(day)
    }
    
    # Check if we are done
    # check for extinction, if there is no allele we are done.
    result$day <- day
    if (sum(pop) == 0) {
      result$fixed <- FALSE
      return(result)  # break
    }
    
    # Check for fixation
    # If all individuals have it, fixation has occurred, we are done.
    if (sum(pop) >= N) {
      result$fixed <- TRUE
      return(result)  # break
    }
    
  }
  # We have finished all iterations without fixation nor extinction
  result$fixed <- FALSE
  return(result)  
}



######  UTILITY FUNCTIONS
# Utility functions for saving, plotting, etc


# get data in an array for saving, ploting, etc.
# only fixed iterations
get_fixed_data <- function(resultado, N){
  # We convert the list from list to vector
  # We need to know whether a given iteration
  # has become fixed or extinct on the day it happened.
 
  fixed <- unlist(resultado[,1])
  fixday <- unlist(resultado[,2])
  
  # If there are no fixings, we go out
  if (sum(fixed)==0) {return(NA)}
  # We create a matrix with the fixation data and the day of the event.
  
  replication_result <- cbind(fixed, fixday)

  # We assign names, not necessary, but just in case we save it to a a file
  colnames(replication_result) <- c("fixed", "fixday")
  
  # We obtain a vector with the iterations that have been set
  # which will be the ones of interest to us.
  fixed <- which(replication_result[,"fixed"]==TRUE)
  
  # We get the maxfixday, to generate an array of mafixdaycol
  maxfixday <- max(fixday[fixed])
  
  # We generate the array to store the data of the iterations
  # that have been set.
  fixedallelefreq <- matrix(N, length(fixed), maxfixday)
  
  # To speed up the loading and facilitate the handling of the data for its 
  # representation we generate a matrix by filling the gaps
  # with the maximum value, with N
  # We lose space and gain agility.
  # Ideally, we would like to save only the data, as we do with all data.
  for (i in 1:length(fixed)){
    fixedallelefreq[i, 1:fixday[fixed[i]]] <- unlist(resultado[fixed[i],]$S1) #resultado[fixed[i],]$S1[1,1:fixday[fixed[i]]] 
  }
  
  data<-cbind(fixed, replication_result[fixed,2], fixedallelefreq)
  #return(list(fixed=fixed, fixday=fixday, maxfixday=maxfixday, 
  #            fixedallelefreq=fixedallelefreq, data=data))
  return(data)

}


# save all data to file. in csv format
# simulation takes so much time we need to save data for later analysis
save_all_data <- function(lst,  filename, ncols=10000){
  if (dim(lst)[1] < 1) {return()}
  i=1
  # each line -> i, fixed?, dayofevent(fix or extinc), sum(pop) day1, sum(pop) day2, sum(pop) day3, ... (until fix or extinc)
  write(c(i, as.integer(lst[i,1]), as.integer(lst[i,2]), as.integer(unlist(lst[i,3]))), filename, ncolumns=ncols, sep = ",", append = FALSE)
  for (i in 2:dim(lst)[1]){
    write(c(i, as.integer(lst[i,1]), as.integer(lst[i,2]), as.integer(unlist(lst[i,3]))), filename, ncolumns=ncols, sep = ",", append = TRUE)
  }  
}

# plot data
plot_data<-function(data, N, b, numreplications) {
  x <- length(data[1,])-2
  y <- 3:length(data[1,])
  plot(1:x, 
       #data[1,3:length(data[1,])], 
       data[1,y], 
       log="y",
       type="l", 
       col=1,
       xlim=c(1,max(data[,2])), 
       ylim=c(2,N),
       xlab="days",
       ylab="allele frequency"
  )
  if (length(data[,1])>1) {
    for (i in 2:length(data[,1])){
      lines(1:x, 
            data[i,y],
            col=i
      )  
    }
  }
  
  title(paste("Allele frequency of a population of", N, "individuals", 1/b, "bottleneck", 
              "\nshowing", length(data[,1]), "fixed replications from", numreplications, "replications"),
        paste("Days to fixation, Min.", min(data[,2]), "  Mean", round(mean(data[,2]),1), "  Max.", max(data[,2]), "  Pfix", length(data[,1])/N)
  )

}




###################  PROGRAM  ################

# We run the programme, calling the different functions

N=c(100,300,1000,3000) # number of individuals
#replicates = 100000
bottleneck=c(1,2,4)

days=10000             # length of experiment

numreplications = 5000               # length of experiment
delta_days=1

dirresult <- format(Sys.time(), "./results/%Y%m%d%H%M_bn_N_2/")
dir.create(dirresult, recursive = TRUE)
sink(paste(dirresult,"run.log", sep = ""), append = TRUE, split = TRUE )
print(paste("Resultados en dir", getwd(), dirresult))

iniTime <- Sys.time()
print(paste("Ini time:", iniTime))
for (n in N){
  for (b in bottleneck) {
    print(paste(numreplications, "replications of", n, "individuals with bottleneck of", 1/b))
    exectime = system.time({
      resultado <- replication(numreplications, n, days, b, delta_days)
    
    
      # We check that we have at least one fixation
      # This happens in tests with few iterations
      numfixed = sum(unlist(resultado[,1]))
      if ( numfixed > 0) {
        data<-get_fixed_data(resultado, n)
      
        write.csv(data, paste(dirresult, "fix_dat_sim_rep_", numreplications, "_pop_", n,"_bn_", b, ".csv", sep=""))
        save_all_data(resultado, paste(dirresult, "all_dat_sim_rep_", numreplications, "_pop_", n,"_bn_", b, ".csv", sep=""))
        plot_data(data, n, b, numreplications)
        
        #next
      } else {
        print(paste("No fixed replications from", numreplications, "total replications"))
        
      }
    })
    print(paste("Replications: ", numreplications, "   N: ", n, "    days: ", days,  " Fixed: ", sum(unlist(resultado[,1]))))
    print(exectime)
  }
  
}

endTime <- Sys.time()
print("Execution terminated")
print(paste("Ini time:", iniTime))
print(paste("End time:", endTime, "Elapsed Time:", endTime-iniTime))
sink()


