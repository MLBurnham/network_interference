suppressPackageStartupMessages({
  library(igraph)
  library(dplyr)
  library(usmap)
  library(ggplot2)
})

###################
## Data Generation
###################
set.seed(123)

edges <- read.csv('state_edges.csv', stringsAsFactors = FALSE)
states <- state.abb[!state.abb %in% list('AK','HI')]
z <- rbinom(length(states), 1, .50)

df <- data.frame(state = states, treatment = z)

g <- graph.data.frame(edges, directed = FALSE, vertices = df)

# get the degree for each node
degree <- degree(g)

# get the number of treated and control neighbors for each node
cn <- sapply(V(g), function(x) {allneighbors = neighbors(g, x) ; length(allneighbors[allneighbors$treatment == 0])})
tn <- sapply(V(g), function(x) {allneighbors = neighbors(g, x) ; length(allneighbors[allneighbors$treatment == 1])})

# proportion of neighbors treated
ptn <- tn/degree
#ptn[is.na(ptn) == TRUE] <- 1 # observations with no neighbors return na, so replace with 1

# function to generate the observed values. We'll pretend we don't know this function.
# Beta is the treatment effect and it slows down the doubling rate
# tau controls the spillover effect
# z is the binary treatment variable
# x is the number of neighbors not treated
# treated turns off spillover effects
natural_growth <- function(beta, z){
  2^(10)/2^(beta * z)
}

spill_over <- function(n, x, tau){
  n * (x * tau)
}

# create the starting value
init <- floor(runif(48, min = 1, max = 5))

# create the observed value. Assumes a doubling time of 3 days with no treatment.
natural_outcome <- init * natural_growth(beta = 2, z = z)
spill_outcome <- spill_over(natural_outcome, x = cn, tau = .05)
y <- natural_outcome + spill_outcome


df <- data.frame(init_cases = init, 
                 cases_1month = round(y), 
                 treatment = z, 
                 neighbors = degree, 
                 treat_neighbors = tn, 
                 control_neighbors = cn,
                 perc_treat_neighbor = ptn)

###################
## Uniformity Trial
###################
# estimate parameters
B1 <- 2.5
B2 <- .1

# if an observation had no spill over, return the original value. Else, remove the spillover effect.
no_spill <- ifelse(df$control_neighbors == 0, df$cases_1month, df$cases_1month/((df$control_neighbors * B2) + 1))
# remove the treatment direct effect
no_treat <- no_spill * 2^(B1 * df$treatment)
# round because there are no partially infected people. Add spill over from non treated states
df$y0 <- round(no_treat + (no_treat * df$neighbors * B2))

#calculate the test statistic
mu_treat <- mean(df[df$treatment == 1, 'y0'])
mu_control <- mean(df[df$treatment == 0, 'y0'])
y0_tstat <- mu_control - mu_treat
y0_tstat

#####################
## Random assignment
#####################
# get 100 unique permutations of the treatment vector
# ideally you want more trials than this, but I kept the number relatively low for speed
perms <- data.frame(sample(z, 48))
nperms <- 1
while(nperms < 100){
  perms <- cbind(perms, sample(z, 48))
  perms[!(duplicated(perms)), ]
  nperms <- length(perms)
}

random_trials <- c()
for(i in 1:length(perms)){ 
  
  # now generate a randomized treatment vector and assign it to the random samples
  fdf <- df
  fdf$fake_treat <- perms[,i]
  fdf <- tibble::rownames_to_column(fdf, "state") 
  
  # create the new graph object to find the number of adjacent treated nodes
  fakeg <- graph.data.frame(edges, directed = FALSE, vertices = fdf)
  
  # for each node, find the proportion of neighbors that were not treated
  fdf$fntn <- sapply(V(fakeg), function(x) {allneighbors = neighbors(fakeg, x) ; length(allneighbors[allneighbors$fake_treat == 0])})
   
  # Now that we have the randomized treatment and the randomized 
  # proportion of treated neighbors, we can calculate the randomized uniformity trial
  # remove spillover
  fno_spill <- ifelse(fdf$fntn == 0, fdf$cases_1month, fdf$cases_1month/((fdf$fntn * B2) + 1))
  # remove treatment
  fno_treat <- fno_spill * 2^(B1 * fdf$fake_treat)
  # add spillover back in
  fdf$fy0 <- round(fno_treat + (fno_treat * fdf$neighbors * B2))
  
  # And now we calculate our test statistic for the randomized uniformity trial and return it
  fmu_treat <- mean(fdf[fdf$fake_treat == 1, 'fy0'])
  fmu_control <- mean(fdf[fdf$fake_treat == 0, 'fy0'])
  fy0_tstat <- fmu_control - fmu_treat
  
  #rtstat <- sum(residuals(lm(fy0 ~ fake_treat + fptn, data = fdf))^2) # use RSS as test stat instead of diff of means
  random_trials <- append(random_trials, fy0_tstat)
}

# get the number of random trials less than than what we observed
num <- sum(abs(y0_tstat) < abs(random_trials))
p <- num/100
p

############################
## Testing for Interference
############################
# Set the estimated value of the parameter
B2 <- .1
# if an observation had no spill over, return the original value. Else, remove the spillover effect.
df$y0 <- ifelse(df$control_neighbors == 0, df$cases_1month, df$cases_1month/((df$control_neighbors * B2) + 1))

y0_rss <- sum(residuals(lm(y0 ~ control_neighbors, data = df[df$treatment == 0,]))^2)
y0_rss

random_trials <- c()
for(i in 1:length(perms)){
  
  # now generate a randomized treatment vector and assign it to the random samples
  fdf <- df
  fdf$fake_treat <- perms[,i]
  fdf <- tibble::rownames_to_column(fdf, "state") 
  
  # create the new graph object to find the number of adjacent treated nodes
  fakeg <- graph.data.frame(edges, directed = FALSE, vertices = fdf)
  
  # for each node, find the proportion of neighbors that were not treated
  fdf$fntn <- sapply(V(fakeg), function(x) {allneighbors = neighbors(fakeg, x) ; length(allneighbors[allneighbors$fake_treat == 0])})

  # Now that we have the randomized treatment and the randomized 
  # proportion of treated neighbors, we can calculate the randomized uniformity trial
  # remove spillover
  fdf$y0 <- ifelse(fdf$fntn == 0, fdf$cases_1month, fdf$cases_1month/((fdf$fntn * B2) + 1))
  
  # And now we calculate our test statistic for the randomized uniformity trial and return it
  fy0_tstat <- sum(residuals(lm(fdf$y0 ~ fdf$control_neighbors, data = fdf[fdf$treatment == 0,]))^2)
  
  #rtstat <- sum(residuals(lm(fy0 ~ fake_treat + fptn, data = fdf))^2) # use RSS as test stat instead of diff of means
  random_trials <- append(random_trials, fy0_tstat)
}

num <- sum(y0_rss < random_trials)
p <- num/100
p