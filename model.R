library(plyr)
library(doMC)
library(ggplot2)
registerDoMC(cores=detectCores())

# command line argument, if present, indicates the results folder
args <- commandArgs(trailingOnly = T)
if (length(args) != 0) {
  res.folder <- args[1]
} else {
  res.folder <- './'
}
res.folder = '/Users/sestari/Google\ Drive/SimulationPE/NEW/P2/NewSimulator'



# Calcolate MU
mu <- function(packet.size=(1460+32)/2) {
  #mean packtes tx = (rate bit/s /packet.size bits) 
  (8*1024^2)/(packet.size*8)
  
}

# total offered load in bits per second
offered.load <- function(lambda, n.nodes, packet.size=(1460+32)/2) {
  lambda*n.nodes*packet.size*8/1024/1024
}



model_1 = function(lamda,mu,N,tr = NULL,mean_tr =NULL){
  
  return = c();
  result = c();
  total = c();
  
  v = matrix(ncol = 5, nrow = length(lamda))
  for (i in 1:length(lamda)){
    
    # Arrival rate 
    la = lamda[i]
    #Infnitesimal generator matrix
    m = matrix(c(
      1,10*la,0,0,0,0,0,0,0,0,0,0,
      1,-(9*la+mu),9*la,0,0,0,0,0,0,0,0,0,
      1,0,-(8*la+2*mu),8*la,0,0,0,0,0,0,0,2*mu,
      1,0,3*mu,-(3*mu+7*la),7*la,0,0,0,0,0,0,0,
      1,0,0,4*mu,-(4*mu+6*la),6*la,0,0,0,0,0,0,
      1,0,0,0,5*mu,-(5*mu+5*la),5*la,0,0,0,0,0,
      1,0,0,0,0,6*mu,-(6*mu+4*la),4*la,0,0,0,0,
      1,0,0,0,0,0,7*mu,-(7*mu+3*la),3*la,0,0,0,
      1,0,0,0,0,0,0,8*mu,-(8*mu+2*la),2*la,0,0,
      1,0,0,0,0,0,0,0,8*mu,-(8*mu+la),la,0,
      1,0,0,0,0,0,0,0,0,10*mu,-(10*mu),0,
      1,0,9*la,0,0,0,0,0,0,0,0,-(9*la+mu)
    ), byrow=T, ncol=12)
    
    B <- c(1, 0, 0,0,0,0,0,0,0,0,0,0)
    pi = solve(t(m), B)
    #average theory
    arrivals = la * pi[2]
    
    v[i,1] = 11
    v[i,2] = la
    v[i,3] = offered.load(arrivals,10) 
    v[i,4] = offered.load(la,10) 
    v[i,5] = 2
    result[i] =offered.load(arrivals,10) 
    total[i] =offered.load(la,10) 
    
    if ((!is.null(mean_tr)  && !is.null(mean_tr)) == TRUE){
      tr<- as.vector(rbind(tr,c(v[i,1],v[i,2],v[i,3],v[i,4])))  
      mean_tr = as.vector(rbind(mean_tr,c(v[i,4],v[i,3],v[i,5])))  
    }
  }
  plot(total,result)
  if ((!is.null(mean_tr)  && !is.null(mean_tr)) == TRUE){
    c(tr,mean_tr)
  }
}


model_2 = function(lamda,mu,N,tr = NULL,mean_tr =NULL){
  
  return = c();
  result = c();
  total = c();
  
  v = matrix(ncol = 5, nrow = length(lamda))
  for (i in 1:length(lamda)){
    
    # Arrival rate 
    la = lamda[i]
    #Infnitesimal generator matrix
    m = matrix(c(
      1,10*la,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      1,-(mu + 9*la),9*la,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      1,0,-(2*mu + 8*la),8*la,0,0,0,0,0,0,0,2*mu,0,0,0,0,0,0,0,0,
      1,0,0,-(3*mu + 7*la),7*la,0,0,0,0,0,0,0,3*mu,0,0,0,0,0,0,0,
      1,0,0,0,-(4*mu + 6*la),6*la,0,0,0,0,0,0,0,4*mu,0,0,0,0,0,0,
      1,0,0,0,0,-(5*mu + 5*la),5*la,0,0,0,0,0,0,0,5*mu,0,0,0,0,0,
      1,0,0,0,0,0,-(6*mu + 4*la),4*la,0,0,0,0,0,0,0,6*mu,0,0,0,0,
      1,0,0,0,0,0,0,-(7*mu + 3*la),3*la,0,0,0,0,0,0,0,7*mu,0,0,0,
      1,0,0,0,0,0,0,0,-(8*mu + 2*la),2*la,0,0,0,0,0,0,0,8*mu,0,0,
      1,0,0,0,0,0,0,0,0,-(9*mu + 1*la),1*la,0,0,0,0,0,0,0,9*mu,0,
      1,0,0,0,0,0,0,0,0,0,-10*mu,0,0,0,0,0,0,0,0,10*mu,
      1,0,9*la,0,0,0,0,0,0,0,0,-(mu + 9*la),0,0,0,0,0,0,0,0,
      1,0,0,8*la,0,0,0,0,0,0,0,1*mu,-(1*mu + 8*la),0,0,0,0,0,0,0,
      1,0,0,0,7*la,0,0,0,0,0,0,0,1*mu,-(1*mu + 7*la),0,0,0,0,0,0,
      1,0,0,0,0,6*la,0,0,0,0,0,0,0,1*mu,-(1*mu + 6*la),0,0,0,0,0,
      1,0,0,0,0,0,5*la,0,0,0,0,0,0,0,1*mu,-(1*mu + 5*la),0,0,0,0,
      1,0,0,0,0,0,0,4*la,0,0,0,0,0,0,0,1*mu,-(1*mu + 4*la),0,0,0,
      1,0,0,0,0,0,0,0,3*la,0,0,0,0,0,0,0,1*mu,-(1*mu + 3*la),0,0,
      1,0,0,0,0,0,0,0,0,2*la,0,0,0,0,0,0,0,1*mu,-(1*mu + 2*la),0,
      1,0,0,0,0,0,0,0,0,0,1*la,0,0,0,0,0,0,0,1*mu,-(1*mu + 1*la)
      
    ), byrow=T, ncol=20)
    B <- c(1, 0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
    pi = solve(t(m), B)
    #average theory
    arrivals = la * pi[2]
    
    v[i,1] = 12
    v[i,2] = la
    v[i,3] = offered.load(arrivals,10) 
    v[i,4] = offered.load(la,10) 
    v[i,5] = 3
    if ((!is.null(mean_tr)  && !is.null(mean_tr)) == TRUE){
      tr<- as.vector(rbind(tr,c(v[i,1],v[i,2],v[i,3],v[i,4])))  
      mean_tr = as.vector(rbind(mean_tr,c(v[i,4],v[i,3],v[i,5])))  
    }
    
    result[i] =offered.load(arrivals,10) 
    total[i] =offered.load(la,10) 
  }
  plot(total,result)
  if ((!is.null(mean_tr)  && !is.null(mean_tr)) == TRUE){
    c(tr,mean_tr)
  }
}


#print models
lamda= c(10,60,110,160,210,260,310,360,410,460,510,560,610,660,710,760,810,860,910,960,1010,1060,1110,1160,1210,1260,1310,1360,1410,1460,1510);
model_1(lamda=lamda ,mu = mu(), N=10)
model_2(lamda=lamda ,mu = mu(), N=10)



# possible packet states
PKT_RECEIVING = 0
PKT_RECEIVED = 1
PKT_CORRUPTED = 2
PKT_GENERATED = 3
PKT_QUEUE_DROPPED = 4

# determine whether a string contains a parsable number"
is.number <- function(string) {
  if (length(grep("^[[:digit:]]*$", string)) == 1)
    return (T)
  else
    return (F)
}

# gets the list of files with a certain prefix and suffix in a folder
get.data.files <- function(folder, suffix=".csv") {
  if (strsplit(suffix, '')[[1]][1] == '.')
    suffix <- paste('\\', suffix, sep='')
  return(list.files(folder, pattern=paste('.*', suffix, sep='')))
}

# splits the name of an output file by _ and extracts the values of simulation parameters
get.params <- function(filename, fields) {
  p <- strsplit(gsub(".csv", "", basename(filename)), "_")[[1]]
  #to add a column, we need to have something in the dataframe, so we add a
  #fake column which we remove at the end
  d <- data.frame(todelete=1)
  for (f in 1:length(fields)) {
    v <- p[f]
    if (is.number(v))
      d[[fields[[f]]]] <- as.numeric(v)
    else
      d[[fields[[f]]]] <- v
  }
  d$todelete <- NULL
  return (d)
}



# compute throughput: total bits received / simulation time
compute.throughput <- function(d, data.rate, sim.time, group=F) {
  fields <- c('lambda')
  if (!group)
    fields <- c('dst', fields)
  throughput <- ddply(d, fields, function(x) {
    received.packets <- subset(x, event == PKT_RECEIVED)
    return(data.frame(tr=sum(received.packets$size*8)/sim.time/(1024**2)))
  }, .parallel=T)
  return(throughput)
}


# if there is no aggregated file, load all csv files into a single one
aggregated.file <- paste(res.folder, 'alld.Rdata', sep='/')
if (!file.exists(aggregated.file)) {
  alld <- data.frame()
  # find all csv in current folder
  data.files <- get.data.files(res.folder, '.csv')
  for (f in data.files) {
    full.path <- paste(res.folder, f, sep='/')
    print(full.path)
    pars <- get.params(full.path, c('prefix', 'lambda', 'seed'))
    d <- read.csv(full.path)
    d <- cbind(d, pars)
    alld <- rbind(d, alld)
  }
  save(alld, file=aggregated.file)
} else {
  # otherwise simply load the aggregated file
  load(aggregated.file)
}

if (is.null(alld$time) == FALSE){
  
  # get simulation time and number of nodes from the simulation data
  sim.time <- max(alld$time)
  n.nodes <- length(unique(alld$src))
  
  # compute the statistics
  tr <- compute.throughput(alld, 8e6, sim.time)
  tr$ol <- offered.load(tr$lambda, n.nodes=n.nodes)
  #Create the throughput mean dataset
  mean_tr = aggregate(tr$tr, by=list(tr$ol), FUN=mean)
  mean_tr[3] = c(1)
  
  dataset =model_1(tr$lambda,mu(),n.nodes,tr,mean_tr)
  tr = data.frame(dataset[1],dataset[2],dataset[3],dataset[4])
  mean_tr =data.frame(dataset[5],dataset[6],dataset[7])
  
  dataset =model_2(tr$lambda,mu(),n.nodes,tr,mean_tr)
  tr = data.frame(dataset[1],dataset[2],dataset[3],dataset[4])
  mean_tr =data.frame(dataset[5],dataset[6],dataset[7])
  # and plot the results
  div <- 3
  
  
  ggplot(mean_tr, aes(x=Group.1, y=x, color = factor(V3, labels = c("Simulator", "Base model","Model improvement")))) +
    geom_line() +
    geom_point() +
    scale_linetype_manual(name="Type", values = c(1,5,3)) +
    scale_color_brewer(name="Type", palette = "Set1") +
    xlab('total offered load (Mbps)') +
    ylab('throughput at receiver (Mbps)') +
    ylim(c(0, 8))  +
    xlim(c(0, 90))
  ggsave(paste(res.folder, 'modelMEANthr.pdf', sep='/'), width=16/div, height=9/div)
  print(mean_tr)
}



