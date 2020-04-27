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

# function to generate the observed values. We'll pretend we don't know this function.

# create the starting value
init <- floor(runif(48, min = 1, max = 5))

# add a random error
u <- abs(round(rnorm(48, mean = 20, sd = 10)))
# generate the data
y <- init + init * 1000 * abs(z-1) + 200 * cn + u

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

summary(lm(cases_1month ~ treatment + treat_neighbors, data = df))

# estimate parameters
B1 <- -2513
B2 <- -8

# remove the estimated spillover effect
df$y0 <- round(df$cases_1month - B1*df$treatment - df$treat_neighbors*B2)

#calculate the test statistic
y0_rss <- lm(y0 ~ treatment + treat_neighbors, data = df) %>%
            residuals %>%
            `^`(2) %>%
            sum

y0_rss

#####################
## Random assignment
#####################
# get 100 unique permutations of the treatment vector
# ideally you want more trials than this, but I kept the number relatively low for speed
set.seed(123)
perms <- data.frame(sample(z, 48)) # start with a single permuted vector of the original treatments
nperms <- 1 # initiate a counter for the number of unique permutations
while(nperms < 100){
    perms <- cbind(perms, sample(z, 48)) # permute the treatment
    perms[!(duplicated(perms)), ] # if the permutation is not unique, drop it
    nperms <- length(perms) # increase the counter
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
    fdf$ftn <- sapply(V(fakeg), function(x) {allneighbors = neighbors(fakeg, x) ; length(allneighbors[allneighbors$fake_treat == 1])})
    fdf$fcn <- sapply(V(fakeg), function(x) {allneighbors = neighbors(fakeg, x) ; length(allneighbors[allneighbors$fake_treat == 0])})
    
    # remove the treatment
    fdf$fy0 <- round(fdf$cases_1month - fdf$fake_treat*B1 - fdf$ftn*B2)

    # And now we calculate our test statistic for the randomized uniformity trial and return it
    fy0_rss <- lm(fy0 ~ fake_treat + ftn, data = fdf) %>% 
                        residuals %>% 
                        `^`(2) %>% 
                        sum
    
    random_trials <- append(random_trials, fy0_rss)
}

# Calculate the p value
num <- sum(y0_rss < random_trials)
p <- num/100
p
