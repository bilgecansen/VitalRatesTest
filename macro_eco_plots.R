
# Load packages and data --------------------------------------------------

library(tidyverse)
library(MCMCvis)
library(doSNOW)
library(foreach)
library(rjags)
library(MCMCvis)
library(ggthemes)
library(adehabitatHR)
library(geosphere)
library(grid)
library(BayesFactor)
library(ggrepel)

# Install with devtools::install_github("thomasp85/patchwork")
library(patchwork)

results_dem <- readRDS("results_dem.rds")
station_coord <- readRDS("station_coord.rds")

N <- results_dem$N
nstat_sp <- results_dem$nstat_sp
spcode <- names(results_dem$R)
npop <- results_dem$npop
effort_year <- results_dem$effort_year
mcmc_sum <- results_dem$mcmc_sum

# Calculate average intrinsic growth rate of each population
# using only years with capture effort
R <- results_dem$R
popR <- list()
for (i in 1:17) {
  
  tempR <- c()
  for (h in 1:nrow(R[[i]]$R)) {
    
    index <- which(effort_year[[i]][h,]>0)
    y <- R[[i]]$R[h,index]
    tempR[h] <- prod(y)^(1/length(y))
    
  }
  
  popR[[i]] <- tempR
  
}
names(popR) <- spcode

# Load jags functions
source("models_jags.R")


# Preliminary explorations ------------------------------------------------

vR <- unlist(popR)
vN <- map(N, function(x) x$popN) %>% unlist()

theme_set(theme_bw())
pg1 <- ggplot(mapping = aes(x = vN)) + 
  geom_histogram(bins = 40, col = "orange", fill = "grey", alpha = 0.9) +
  labs(x = bquote(bar("N")),
       y = "Frequency") +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(0,1200,200))

pg2 <- ggplot(mapping = aes(x = log(vR))) + 
  geom_histogram(bins = 40, col = "orange", fill = "grey", alpha = 0.9) +
  labs(x = bquote(bar("r")),
       y = "Frequency") +
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) +
  scale_x_continuous(breaks = seq(-1,1,0.25))

theme_set(theme_bw())
pg1 / pg2 +
  plot_annotation(tag_levels = 'a')
if(!"paper_results" %in% list.files()) dir.create("paper_results")
ggsave("paper_results/macro_figS1.jpeg", width = 4, height = 6, units = "in")



# Abundance-occupancy relationship ----------------------------------------

# Estimate range as minimum convex polygon
coord <- list()
poly <- list()
range <- c()
for (i in 1:length(spcode)) {
  
  index <- which(as.character(station_coord$spec) == as.character(spcode[i]) &
                  station_coord$lat < 99)
  coord[[i]] <- station_coord[index, 4:3]
  coordinates(coord[[i]]) <- c("long", "lat")
  poly[[i]] <- mcp(coord[[i]], percent = 100)
  range[i] <- areaPolygon(poly[[i]]@polygons[[1]]@Polygons[[1]]@coords)/1000000
  
}

# N vs Range
spN <- map_dbl(N, function(x) median(x$spN))

inits <- list(
  list(mu = c(1,10), sigma = c(0.5,0.5), rho = -0.4),
  list(mu = c(0,0), sigma = c(1,5), rho = 0),
  list(mu = c(-10,-1), sigma = c(5,1), rho = 0.4)
)

NvsRange <- cor.jags(data = list(y = cbind(log(spN), log(range)), 
                                 n = length(spN)),
                     inits = inits,
                     n.chains = 3,
                     n.adapt = 1000,
                     n.update = 10000,
                     n.iter = 5000,
                     n.thin = 5,
                     seed_no = 19)

NvsRange$mcmc_sum
rho <- MCMCchains(NvsRange$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)
correlationBF(log(spN), log(range), rscale = "wide")

## Plot
my_text1a <- bquote(rho == 0.06 ~ "(-0.41, 0.52)")
my_grob1a <-  grid.text(my_text1a, x = 0.05, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text1b <- bquote("P("*rho*">0)" ~ "=" ~ "0.52")
my_grob1b <-  grid.text(my_text1b, x=0.05,  y=0.80, gp=gpar(col="orange", fontsize = 10), just = "left")
my_text1c <- my_text1c <- bquote(B[0] == 0.40)
my_grob1c <-  grid.text(my_text1c, x=0.05,  y=0.70, gp=gpar(col="orange", fontsize = 10), just = "left")

g1 <- ggplot(mapping = aes(x = log(range), log(spN))) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) + 
  labs(x = "Geographic Extent",
       y = bquote(bar(bar("N")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1a) +
  annotation_custom(my_grob1b) +
  annotation_custom(my_grob1c)
  

# N vs Number of Populations
NvsNpop <- cor.jags(data = list(y = cbind(log(spN), log(npop)), 
                                n = length(spN)),
                    inits = inits,
                    n.chains = 3,
                    n.adapt = 1000,
                    n.update = 10000,
                    n.iter = 5000,
                    n.thin = 5,
                    seed_no = 19)

NvsNpop$mcmc_sum
rho <- MCMCchains(NvsNpop$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)
correlationBF(log(spN), log(npop), rscale = "wide")

## Plot
my_text2a <- bquote(rho == 0.26 ~ "(-0.20, 0.68)")
my_grob2a <-  grid.text(my_text2a, x = 0.05, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text2b <- bquote("P("*rho*">0)" ~ "=" ~ "0.87")
my_grob2b <-  grid.text(my_text2b, x=0.05,  y = 0.80, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text2c <- my_text2c <- bquote(B[0] == 0.71)
my_grob2c <-  grid.text(my_text2c, x=0.05,  y=0.70, gp=gpar(col="orange", fontsize = 10), just = "left")

g2 <- ggplot(mapping = aes(x = log(npop), log(spN))) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) + 
  labs(x = "Number of Populations",
       y = bquote(bar(bar("N")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob2a) + 
  annotation_custom(my_grob2b) +
  annotation_custom(my_grob2c)

g1 / g2 + 
  plot_annotation(tag_levels = 'a')

ggsave("paper_results/macro_fig2.jpeg", width = 8, height = 8, units = "in")


# Relationship between r and occupancy ------------------------------------

# R vs Range
spR <- map_dbl(popR, function(x) median(log(x)))

RvsRange <- cor.jags(data = list(y = cbind(spR, log(range)), 
                                 n = length(spR)),
                     inits = inits,
                     n.chains = 3,
                     n.adapt = 1000,
                     n.update = 10000,
                     n.iter = 5000,
                     n.thin = 5,
                     seed_no = 19)

RvsRange$mcmc_sum
rho <- MCMCchains(RvsRange$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)
correlationBF(spR, log(range), rscale = "wide")

## Plot
my_text3a <- bquote(rho == 0.13 ~ "(-0.34, 0.59)")
my_grob3a <-  grid.text(my_text3a, x = 0.05, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text3b <- bquote("P("*rho*">0)" ~ "=" ~ "0.71")
my_grob3b <-  grid.text(my_text3b, x=0.05,  y = 0.80, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text3c <- my_text3c <- bquote(B[0] == 0.46)
my_grob3c <-  grid.text(my_text3c, x=0.05,  y = 0.70, gp=gpar(col="orange", fontsize = 10), just = "left")

g3 <- ggplot(mapping = aes(x = log(range), y = spR)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) + 
  labs(x = "Geographic Extent",
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob3a) +
  annotation_custom(my_grob3b) +
  annotation_custom(my_grob3c)

# R vs Number of Populations
RvsNpop <- cor.jags(data = list(y = cbind(spR, log(npop)), 
                                n = length(spR)),
                    inits = inits,
                    n.chains = 3,
                    n.adapt = 1000,
                    n.update = 10000,
                    n.iter = 5000,
                    n.thin = 5,
                    seed_no = 19)

RvsNpop$mcmc_sum
rho <- MCMCchains(RvsNpop$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)
correlationBF(spR, log(npop), rscale = "wide")

## Plot
my_text4a <- bquote(rho == -0.10 ~ "(-0.54, 0.38)")
my_grob4a <-  grid.text(my_text4a, x = 0.05, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text4b <- bquote("P("*rho*">0)" ~ "=" ~ "0.34")
my_grob4b <-  grid.text(my_text3b, x = 0.05,  y = 0.8, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text4c <- my_text4c <- bquote(B[0] == 0.43)
my_grob4c <-  grid.text(my_text4c, x=0.05,  y = 0.70, gp=gpar(col="orange", fontsize = 10), just = "left")

g4 <- ggplot(mapping = aes(x = log(npop), y = spR)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) +  
  labs(x = "Number of Populations",
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob4a) +
  annotation_custom(my_grob4b) +
  annotation_custom(my_grob4c)

g3 / g4 
ggsave("paper_results/macro_fig3.jpeg", width = 8, height = 8, units = "in")


# Relationship between r and N --------------------------------------------

popN <- map(N, function(x) x$popN)

RvsN_pop <- foreach(i=1:17) %do% {
  
  inits <- list(
    list(mu = c(1,10), sigma = c(0.5,0.5), rho = -0.4),
    list(mu = c(0,0), sigma = c(1,5), rho = 0),
    list(mu = c(-10,-1), sigma = c(5,1), rho = 0.4)
  )
  
  cor.jags(data = list(y = cbind(log(popR[[i]]), log(popN[[i]])), 
                       n = length(popR[[i]])),
           inits = inits,
           n.chains = 3,
           n.adapt = 1000,
           n.update = 10000,
           n.iter = 5000,
           n.thin = 5,
           seed_no = 19)
  
}

## Check convergence
map_dbl(RvsN_pop, function(x) max(x$mcmc_sum[,"Rhat"]))
map_dbl(RvsN_pop, function(x) min(x$mcmc_sum[,"n.eff"]))

## Plot
rho_chains <- map(RvsN_pop, function(x) MCMCchains(x$mcmc_samples, params = "rho"))
rho_median <- data.frame(rho_med = map_dbl(rho_chains, median),
                         species = spcode)
rho_50 <- data.frame(rho_min = map_dbl(rho_chains, function(x) quantile(x, 0.25)),
                     rho_max = map_dbl(rho_chains, function(x) quantile(x, 0.75)),
                     species = spcode)
rho_90 <- data.frame(rho_min = map_dbl(rho_chains, function(x) quantile(x, 0.025)),
                     rho_max = map_dbl(rho_chains, function(x) quantile(x, 0.975)),
                     species = spcode)

bf <- map(cor_bf, function(x) round(exp(x@bayesFactor$bf), 1))
bf[[13]] <- ">100"
  
g7 <- ggplot() +
  geom_linerange(data = rho_50, aes(x = spcode, ymin = rho_min, ymax = rho_max), size = 0.9, alpha = 0.8) +
  geom_errorbar(data = rho_90, aes(x = spcode, ymin = rho_min, ymax = rho_max), alpha = 0.5) +
  geom_point(data = rho_median, aes(x = species, y = rho_med), col = "orange", size = 1) +
  geom_text(mapping = aes(x = spcode, y = -0.75, label = bf), size = 3) +
  labs(y = bquote("Correlation" * "  \u2013  " * bar("r") * " vs " * bar("N"))) +
  geom_hline(yintercept = 0, col = "orange", alpha = 0.4) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8)) +
  scale_y_continuous(breaks = seq(-0.5,75, 0.25)) +
  coord_flip()

g7
ggsave("paper_results/macro_fig4.jpeg", width = 6, height = 6, units = "in")


# Variability in DD strength ----------------------------------------------

# species level r vs species level N
RvsN_sp <- cor.jags(data = list(y = cbind(spR, log(spN)), 
                                     n = length(spR)),
                        inits = inits,
                        n.chains = 3,
                        n.adapt = 1000,
                        n.update = 10000,
                        n.iter = 5000,
                        n.thin = 5,
                        seed_no = 19)

RvsN_sp$mcmc_sum
rho <- MCMCchains(RvsN_sp$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)
correlationBF(log(spN), spR, rscale = "wide")
correlationBF(log(spN), spR, rscale = "wide", nullInterval = c(-0.5, 1))

## Plot
my_text5a <- bquote(rho == -0.51 ~ "(-0.84, -0.14)")
my_grob5a <-  grid.text(my_text5a, x = 0.7, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text5b <- bquote("P("*rho*">0)" ~ "=" ~ "0.01")
my_grob5b <-  grid.text(my_text5b, x = 0.7,  y = 0.85, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text5c <- bquote(B[0] == 4.68)
my_grob5c <-  grid.text(my_text5c, x = 0.7, y = 0.80, gp = gpar(col="orange", fontsize = 10), just = "left")
#my_text5d <- bquote(B[i] == 3.1)
#my_grob5d <-  grid.text(my_text5d, x = 0.7, y = 0.75, gp = gpar(col="orange", fontsize = 10), just = "left")

g8 <- ggplot(mapping = aes(x = log(spN), y = spR)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) +   
  labs(x = bquote(bar(bar("N"))),
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob5a) +
  annotation_custom(my_grob5b) +
  annotation_custom(my_grob5c) #+
  #annotation_custom(my_grob5d)

g8
ggsave("paper_results/macro_fig5.jpeg", width = 6, height = 6, units = "in")

dd <- round(spR/spN, 4)
sd(dd)/mean(dd)

# DD plots ----------------------------------------------------------------

beta <- map_dbl(mcmc_sum, function(x) x["beta",1])
zeta <- map_dbl(mcmc_sum, function(x) x["zeta",1])
sad <- map_dbl(mcmc_sum, function(x) x["survival_ad",1])
sjuv <- map_dbl(mcmc_sum, function(x) x["survival_juv",1])
fec <- map_dbl(mcmc_sum, function(x) x["fecundity",1])

# Population size vs DD in Survival
inits <- list(
  list(mu = c(1,10), sigma = c(0.5,0.5), rho = -0.4),
  list(mu = c(0,0), sigma = c(1,5), rho = 0),
  list(mu = c(5,2), sigma = c(5,1), rho = 0.4)
)

NvsB <- cor.jags(data = list(y = cbind(log(spN), -beta), 
                             n = length(spN)),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

NvsB$mcmc_sum
correlationBF(log(spN), -beta, rscale = "wide")
#correlationBF(log(spN), -beta, rscale = "wide", nullInterval = c(-0.3, 0.3))

## Plot
my_text6a <- bquote(rho == -0.57 ~ "(-0.86, -0.21)")
my_grob6a <-  grid.text(my_text6a, x = 0.7, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text6b <- bquote(B[0] == 9.09)
my_grob6b <-  grid.text(my_text6b, x = 0.7, y = 0.85, gp = gpar(col="orange", fontsize = 10), just = "left")
#my_text6c <- bquote(B[i] == 5.42)
#my_grob6c <-  grid.text(my_text6c, x = 0.7, y = 0.80, gp = gpar(col="orange", fontsize = 8), just = "left")

#qNlow <- map_dbl(Nchains, function(x) quantile(log(x$medianN), 0.25))
#qNhigh <- map_dbl(Nchains, function(x) quantile(log(x$medianN), 0.75))
#qBlow <- map_dbl(results_cjspop$mcmc_chains, function(x) quantile(x[,"beta"], 0.25))
#qBhigh <- map_dbl(results_cjspop$mcmc_chains, function(x) quantile(x[,"beta"], 0.75))

g9 <- ggplot(mapping = aes(x = log(spN), y = -beta)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) + 
  #geom_errorbar(mapping = aes(ymin = qBlow, ymax = qBhigh), alpha = 0.5) +
  #geom_errorbarh(mapping = aes(xmin = qNlow, xmax = qNhigh), alpha = 0.5) +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(beta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  scale_y_continuous(breaks = seq(-0.5,0.75, 0.25), limits = c(-0.25,0.75)) +
  annotation_custom(my_grob6a) +
  annotation_custom(my_grob6b) #+
  #annotation_custom(my_grob6c)

# Population size vs DD in Fecundity
NvsZ <- cor.jags(data = list(y = cbind(log(spN), -zeta), 
                             n = length(spN)),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

NvsZ$mcmc_sum
correlationBF(log(spN), -zeta, rscale = "wide")

## Plot
my_text7a <- bquote(rho == 0.30 ~ "(-0.15, -0.71)")
my_grob7a <-  grid.text(my_text7a, x = 0.05, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text7b <- bquote(B[0] == 0.9)
my_grob7b <-  grid.text(my_text7b, x = 0.05, y = 0.85, gp = gpar(col="orange", fontsize = 10), just = "left")

g10 <- ggplot(mapping = aes(x = log(spN), y = -zeta)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) + 
  labs(x = bquote(bar(bar("N"))),
       y = bquote(zeta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob7a) +
  annotation_custom(my_grob7b)

# Growth rate vs DD in Survival
RvsB <- cor.jags(data = list(y = cbind(log(spR), -beta), 
                             n = length(spR)),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

RvsB$mcmc_sum
correlationBF(spR, -beta, rscale = "wide")
#correlationBF(spR, -beta, rscale = "wide", nullInterval = c(-0.3, 0.3))

## Plot
my_text8a <- bquote(rho == -0.53 ~ "(0.17, 0.85)")
my_grob8a <-  grid.text(my_text8a, x = 0.05, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text8b <- bquote(B[0] == 7.18)
my_grob8b <-  grid.text(my_text8b, x = 0.05, y = 0.85, gp = gpar(col="orange", fontsize = 10), just = "left")
#my_text8c <- bquote(B[i] == 4.45)
#my_grob8c <-  grid.text(my_text8c, x = 0.05, y = 0.80, gp = gpar(col="orange", fontsize = 8), just = "left")

g11 <- ggplot(mapping = aes(x = log(spR), y = -beta)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) +
  labs(x = bquote(bar(bar("r"))),
       y = bquote(beta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  scale_y_continuous(breaks = seq(-0.5,0.75, 0.25), limits = c(-0.25,0.75)) +
  annotation_custom(my_grob8a) +
  annotation_custom(my_grob8b) #+
  #annotation_custom(my_grob8c)

# Population size vs DD in Fecundity
RvsZ <- cor.jags(data = list(y = cbind(log(spR), -zeta), 
                             n = length(spR)),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

RvsZ$mcmc_sum
correlationBF(spR, -zeta, rscale = "wide")

## Plot
my_text9a <- bquote(rho == -0.14 ~ "(-0.59, 0.32)")
my_grob9a <-  grid.text(my_text9a, x = 0.5, y = 0.95, gp = gpar(col="orange", fontsize = 10), just = "left")
my_text9b <- bquote(B[0] == 0.46)
my_grob9b <-  grid.text(my_text9b, x = 0.5, y = 0.9, gp = gpar(col="orange", fontsize = 10), just = "left")

g12 <- ggplot(mapping = aes(x = log(spR), y = -zeta)) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) +
  labs(x = bquote(bar(bar("r"))),
       y = bquote(zeta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob9a) +
  annotation_custom(my_grob9b)

(g9 + g10) / (g11 + g12) +
  plot_annotation(tag_levels = 'a')

ggsave("paper_results/macro_fig6.jpeg", width = 10, height = 10, units = "in")


# Explorations statistical patterns ---------------------------------------

pad <- map_dbl(mcmc_sum, function(x) inv.logit(x["gamma[2]",1]))
pjuv <- map_dbl(mcmc_sum, function(x) inv.logit(x["gamma[1]",1]))

# Capture probability vs N and R
tg1 <- ggplot(mapping = aes(x = 1-(1-pad)^4, y = log(spN))) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) +
  #geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Yearly Adult Capture Probability",
       y = bquote(bar(bar("N"))))

tg2 <- ggplot(mapping = aes(x = 1-(1-pjuv)^4, y = log(spN))) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) +
  #geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Yearly Juvenile Capture Probability",
       y = bquote(bar(bar("N"))))

tg3 <- ggplot(mapping = aes(x = 1-(1-pad)^4, y = log(spR))) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) +
  #geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Yearly Adult Capture Probability",
       y = bquote(bar(bar("r"))))

tg4 <- ggplot(mapping = aes(x = 1-(1-pjuv)^4, y = log(spR))) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode), size = 2) +
  #geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Yearly Juvenile Capture Probability",
       y = bquote(bar(bar("r"))))

(tg1 + tg2) / (tg3 + tg4) +
  plot_annotation(tag_levels = 'a')

ggsave("paper_results/macro_figS2.jpeg", width = 8, height = 8, units = "in")

# Patterns when low cp species are removed
index <- which(1-(1-pjuv)^4 >=0.05)

#NvsNpop2 <- cor.jags(data = list(y = cbind(log(spN[index]), log(npop[index])), 
                                #n = length(spN[index])),
                    #inits = inits,
                    #n.chains = 3,
                    #n.adapt = 1000,
                    #n.update = 10000,
                    #n.iter = 5000,
                    #n.thin = 5,
                    #seed_no = 19)

#NvsNpop2$mcmc_sum
#rho <- MCMCchains(NvsNpop2$mcmc_samples, params = "rho")
#length(which(rho>0))/length(rho)

## Plot
#my_text1 <- bquote(rho ~ "=" ~ "0.16" ~ "(-0.42, 0.67)")
#my_grob1 <-  grid.text(my_text1, x=0.7,  y=0.25, gp=gpar(col="orange", fontsize=10))
#my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.72")
#my_grob2 <-  grid.text(my_text2, x=0.7,  y=0.15, gp=gpar(col="orange", fontsize=10))

#tg5 <- ggplot(mapping = aes(x = log(npop[index]), log(spN[index]))) +
  #geom_point() +
  #geom_smooth(method="lm", se=F, color = "orange") +
  #labs(x = "Number of Populations",
       #y = bquote(bar(bar("N")))) + 
  #theme(axis.title.x = element_text(size = 10),
        #axis.title.y = element_text(size = 12)) +
#annotation_custom(my_grob1) +
#annotation_custom(my_grob2)

RvsN_sp2 <- cor.jags(data = list(y = cbind(spR[index], log(spN[index])), 
                                n = length(spR[index])),
                    inits = inits,
                    n.chains = 3,
                    n.adapt = 1000,
                    n.update = 10000,
                    n.iter = 5000,
                    n.thin = 5,
                    seed_no = 19)

RvsN_sp2$mcmc_sum
rho <- MCMCchains(RvsN_sp2$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)
correlationBF(log(spN[index]), spR[index], rscale = "wide")

## Plot
my_text10a <- bquote(rho ~ "=" ~ "-0.34" ~ "(-0.76, 0.26)")
my_grob10a <-  grid.text(my_text10a, x = 0.10,  y = 0.6, gp = gpar(col="orange", fontsize=10), just = "left")
my_text10b <- bquote("P("*rho*">0)" ~ "=" ~ "0.11")
my_grob10b <-  grid.text(my_text10b, x = 0.10,  y = 0.5, gp = gpar(col="orange", fontsize=10), just = "left")
my_text10c <- bquote(B[0] == 0.96)
my_grob10c <-  grid.text(my_text10c, x = 0.10, y = 0.4, gp = gpar(col="orange", fontsize = 10), just = "left")

tg6 <- ggplot(mapping = aes(x = log(spN)[index], y = log(spR)[index])) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode[index]), size = 2) +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob10a) +
  annotation_custom(my_grob10b) +
  annotation_custom(my_grob10c)

NvsB2 <- cor.jags(data = list(y = cbind(log(spN[index]), -beta[index]), 
                             n = length(spN[index])),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

NvsB2$mcmc_sum
correlationBF(log(spN[index]), -beta[index], rscale = "wide")

## Plot
my_text11a <- bquote(rho ~ "=" ~ "-0.45" ~ "(-0.82, -0.11)")
my_grob11a <-  grid.text(my_text11a, x = 0.10,  y = 0.5, gp=gpar(col="orange", fontsize=10), just = "left")
my_text11b <- bquote(B[0] == 1.8)
my_grob11b <-  grid.text(my_text11b, x = 0.10, y = 0.4, gp = gpar(col="orange", fontsize = 10), just = "left")

tg7 <- ggplot(mapping = aes(x = log(spN)[index], y = -beta[index])) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode[index]), size = 2) +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(beta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob11a) +
  annotation_custom(my_grob11b)

RvsB2 <- cor.jags(data = list(y = cbind(log(spR[index]), -beta[index]), 
                             n = length(spR[index])),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

RvsB2$mcmc_sum
correlationBF(spR[index], -beta[index], rscale = "wide")


## Plot
my_text12a <- bquote(rho ~ "=" ~ "0.45" ~ "(-0.12, 0.82)")
my_grob12a <-  grid.text(my_text12a, x = 0.10,  y = 0.5, gp=gpar(col="orange", fontsize=10), just = "left")
my_text12b <- bquote(B[0] == 1.7)
my_grob12b <-  grid.text(my_text12b, x = 0.10, y = 0.4, gp = gpar(col="orange", fontsize = 10), just = "left")

tg8 <- ggplot(mapping = aes(x = log(spR)[index], y = -beta[index])) +
  geom_point(size = 3) +
  geom_label_repel(aes(label = spcode[index]), size = 2) +
  labs(x = bquote(bar(bar("r"))),
       y = bquote(beta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob12a) +
  annotation_custom(my_grob12b)

(tg6 / tg7 / tg8) +
  plot_annotation(tag_levels = 'a')

ggsave("paper_results/macro_figS3.jpeg", width = 6, height = 10, units = "in")

(sd(spR/spN)/mean(spR/spN))*100
(sd(spR[index]/spN[index])/mean(spR[index]/spN[index]))*100

# Additional analysis -----------------------------------------------------

index2 <- which(rho_median$rho_med>0)

NvsNpop3 <- cor.jags(data = list(y = cbind(log(spN[index2]), log(npop[index2])), 
                                 n = length(spN[index2])),
                     inits = inits,
                     n.chains = 3,
                     n.adapt = 1000,
                     n.update = 10000,
                     n.iter = 5000,
                     n.thin = 5,
                     seed_no = 19)

NvsNpop3$mcmc_sum
rho <- MCMCchains(NvsNpop3$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)

RvsNpop2 <- cor.jags(data = list(y = cbind(log(spR[index2]), log(npop[index2])), 
                                 n = length(spR[index2])),
                     inits = inits,
                     n.chains = 3,
                     n.adapt = 1000,
                     n.update = 10000,
                     n.iter = 5000,
                     n.thin = 5,
                     seed_no = 19)

RvsNpop2$mcmc_sum
rho <- MCMCchains(RvsNpop2$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)

NvsRD <- cor.jags(data = list(y = cbind(log(spN), log(npop)/log(range)), 
                                 n = length(spN)),
                     inits = inits,
                     n.chains = 3,
                     n.adapt = 1000,
                     n.update = 10000,
                     n.iter = 5000,
                     n.thin = 5,
                     seed_no = 19)

NvsRD$mcmc_sum


# Tables ------------------------------------------------------------------

beta_low <- map_dbl(mcmc_sum, function(x) x["beta",3])
zeta_low <- map_dbl(mcmc_sum, function(x) x["zeta",3])
sad_low <- map_dbl(mcmc_sum, function(x) x["survival_ad",3])
sjuv_low <- map_dbl(mcmc_sum, function(x) x["survival_juv",3])
fec_low <- map_dbl(mcmc_sum, function(x) x["fecundity",3])
pad_low <- map_dbl(mcmc_sum, function(x) inv.logit(x["gamma[2]",3]))
pjuv_low <- map_dbl(mcmc_sum, function(x) inv.logit(x["gamma[1]",3]))

beta_high <- map_dbl(mcmc_sum, function(x) x["beta",5])
zeta_high <- map_dbl(mcmc_sum, function(x) x["zeta",5])
sad_high <- map_dbl(mcmc_sum, function(x) x["survival_ad",5])
sjuv_high <- map_dbl(mcmc_sum, function(x) x["survival_juv",5])
fec_high <- map_dbl(mcmc_sum, function(x) x["fecundity",5])
pad_high <- map_dbl(mcmc_sum, function(x) inv.logit(x["gamma[2]",5]))
pjuv_high <- map_dbl(mcmc_sum, function(x) inv.logit(x["gamma[1]",5]))

fun_paste <- function(x,low,high) {
  
  paste(round(x,3), "(", sep = " ") %>%
    paste0(., round(low,3)) %>%
    paste0(., ", ", round(high,3)) %>%
    paste0(., ")")
  
}

tableS1 <- data.frame(spcode = spcode,
                      sad = fun_paste(sad, sad_low, sad_high),
                      sjuv = fun_paste(sjuv, sjuv_low, sjuv_high),
                      fec = fun_paste(fec, fec_low, fec_high),
                      beta = fun_paste(beta, beta_low, beta_high),
                      zeta = fun_paste(zeta, zeta_low, zeta_high),
                      pad =  fun_paste(1-(1-pad)^4, 1-(1-pad_low)^4, 1-(1-pad_high)^4),
                      pjuv = fun_paste(1-(1-pjuv)^4, 1-(1-pjuv_low)^4, 1-(1-pjuv_high)^4))

write.csv(tableS1, "tableS1.csv")

tableS2 <- data.frame(spcode = spcode,
                      r = round(log(spR),2),
                      N = round(spN,2),
                      npop = npop,
                      nstat = sum_nstat,
                      range = round(log(range),2))

write.csv(tableS2, "tableS2.csv")

