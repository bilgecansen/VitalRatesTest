
# Load packages and data --------------------------------------------------

library(tidyverse)
library(MCMCvis)
library(doSNOW)
library(foreach)
library(rjags)
library(MCMCvis)
library(ggthemes)
library(adehabitatHR)



# Install with devtools::install_github("thomasp85/patchwork")
library(patchwork)

results_dem <- readRDS("results_dem.rds")

N <- results_dem$N
nstat_sp <- results_dem$nstat_sp
spcode <- names(results_dem$R)

R <- results_dem$R
popR <- list()
for (i in 1:17) {
  
  tempR <- c()
  for (h in 1:nrow(R[[i]]$R)) {
    
    index <- which(results_cjspop$chdata[[i]]$effort_year[h,]>0)
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
ggsave("paper_results/macro_figS1.tiff", width = 4, height = 6, units = "in")


# Hypothesis 1 and 2: Range is positively correlated witk K and R ---------

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

## Check prior posterior overlap
mu_prior <- rnorm(15000, 0,10)
sigma_prior <- runif(15000, 0, 100)
rho_prior <- runif(15000, -1, 1)
MCMCtrace(NvsRange$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(NvsRange$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(NvsRange$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.012" ~ "(-0.460, 0.473)")
my_grob1 <-  grid.text(my_text1, x=0.3,  y=0.9, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.523")
my_grob2 <-  grid.text(my_text2, x=0.3,  y=0.8, gp=gpar(col="orange", fontsize=8))


g1 <- ggplot(mapping = aes(x = log(range), log(spN))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Geographic Extent",
       y = bquote(bar(bar("N")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1) +
  annotation_custom(my_grob2)

# R vs Range
spR <- map_dbl(popR, median)

RvsRange <- cor.jags(data = list(y = cbind(log(spR), log(range)), 
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

## Check prior posterior overlap
MCMCtrace(RvsRange$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(RvsRange$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(RvsRange$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.126" ~ "(-0.367, 0.569)")
my_grob1 <-  grid.text(my_text1, x=0.25,  y=0.90, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.700")
my_grob2 <-  grid.text(my_text2, x=0.25,  y=0.80, gp=gpar(col="orange", fontsize=8))


g2 <- ggplot(mapping = aes(x = log(range), y = log(spR))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Geographic Extent",
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1) +
  annotation_custom(my_grob2)

# N vs Number of Populations
npop <- map_dbl(results_cjspop$chdata, function(x) x$npop)

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

## Check prior posterior overlap
MCMCtrace(NvsNpop$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(NvsNpop$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(NvsNpop$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.261" ~ "(-0.225, 0.657)")
my_grob1 <-  grid.text(my_text1, x=0.3,  y=0.9, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.865")
my_grob2 <-  grid.text(my_text2, x=0.3,  y=0.80, gp=gpar(col="orange", fontsize=8))


g3 <- ggplot(mapping = aes(x = log(npop), log(spN))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Number of Populations",
       y = bquote(bar(bar("N")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1) +
  annotation_custom(my_grob2)

# R vs Number of Populations
RvsNpop <- cor.jags(data = list(y = cbind(log(spR), log(npop)), 
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

## Check prior posterior overlap
MCMCtrace(RvsNpop$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(RvsNpop$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(RvsNpop$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "-0.099" ~ "(-0.551, 0.381)")
my_grob1 <-  grid.text(my_text1, x=0.76,  y=0.90, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.346")
my_grob2 <-  grid.text(my_text2, x=0.76,  y=0.80, gp=gpar(col="orange", fontsize=8))


g4 <- ggplot(mapping = aes(x = log(npop), y = log(spR))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Number of Populations",
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1) +
  annotation_custom(my_grob2)

# N vs Number of Stations
## Estimate range as number of stations captured
sum_nstat <- map_dbl(nstat_sp, sum)

NvsNstat <- cor.jags(data = list(y = cbind(log(spN), log(sum_nstat)), 
                                 n = length(spN)),
                     inits = inits,
                     n.chains = 3,
                     n.adapt = 1000,
                     n.update = 10000,
                     n.iter = 5000,
                     n.thin = 5,
                     seed_no = 19)

NvsNstat$mcmc_sum
rho <- MCMCchains(NvsNstat$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)

## Check prior posterior overlap
MCMCtrace(NvsNstat$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(NvsNstat$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(NvsNstat$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.056" ~ "(-0.424, 0.510)")
my_grob1 <-  grid.text(my_text1, x=0.3,  y=0.9, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.598")
my_grob2 <-  grid.text(my_text2, x=0.3,  y=0.8, gp=gpar(col="orange", fontsize=8))


g5 <- ggplot(mapping = aes(x = log(sum_nstat), log(spN))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Number of Stations",
       y = bquote(bar(bar("N")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1) +
  annotation_custom(my_grob2)

# R vs Number of Stations
RvsNstat <- cor.jags(data = list(y = cbind(log(spR), log(sum_nstat)), 
                                 n = length(spR)),
                     inits = inits,
                     n.chains = 3,
                     n.adapt = 1000,
                     n.update = 10000,
                     n.iter = 5000,
                     n.thin = 5,
                     seed_no = 19)

RvsNstat$mcmc_sum
rho <- MCMCchains(RvsNstat$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)

## Check prior posterior overlap
MCMCtrace(RvsNstat$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(RvsNstat$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(RvsNstat$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.073" ~ "(-0.415, 0.531)")
my_grob1 <-  grid.text(my_text1, x=0.25,  y=0.90, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.620")
my_grob2 <-  grid.text(my_text2, x=0.25,  y=0.80, gp=gpar(col="orange", fontsize=8))


g6 <- ggplot(mapping = aes(x = log(sum_nstat), y = log(spR))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Number of Stations ",
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1) +
  annotation_custom(my_grob2)

(g1 + g2) / (g3 + g4) / (g5 + g6) +
  plot_annotation(tag_levels = 'a')

ggsave("paper_results/macro_fig1.tiff", width = 8, height = 8, units = "in")

# Hypothesis 3: Population level R and N is positively correlated ---------

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
rho_90 <- data.frame(rho_min = map_dbl(rho_chains, function(x) quantile(x, 0.05)),
                     rho_max = map_dbl(rho_chains, function(x) quantile(x, 0.95)),
                     species = spcode)


g7 <- ggplot() +
  geom_linerange(data = rho_50, aes(x = spcode, ymin = rho_min, ymax = rho_max), size = 0.9, alpha = 0.8) +
  geom_errorbar(data = rho_90, aes(x = spcode, ymin = rho_min, ymax = rho_max), alpha = 0.5) +
  geom_point(data = rho_median, aes(x = species, y = rho_med), col = "orange", size = 1) +
  labs(y = bquote("Correlation" * "  \u2013  " * bar("r") * " vs " * bar("N"))) +
  geom_hline(yintercept = 0, col = "orange", alpha = 0.4) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 8)) +
  scale_y_continuous(breaks = seq(-0.5,75, 0.25)) +
  coord_flip()

g7
ggsave("paper_results/macro_fig2.tiff", width = 6, height = 6, units = "in")

# Hypothesis 4: Is avgR correlated with K? ---------------------------

RvsN_sp <- cor.jags(data = list(y = cbind(log(spR), log(spN)), 
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

## Check prior posterior overlap
MCMCtrace(RvsN_sp$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(RvsN_sp$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(RvsN_sp$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "-0.513" ~ "(-0.809, -0.076)")
my_grob1 <-  grid.text(my_text1, x=0.75,  y=0.9, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.013")
my_grob2 <-  grid.text(my_text2, x=0.75,  y=0.8, gp=gpar(col="orange", fontsize=8))


g8 <- ggplot(mapping = aes(x = log(spN), y = log(spR))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1) +
  annotation_custom(my_grob2)

g8
ggsave("paper_results/macro_fig3.tiff", width = 6, height = 6, units = "in")


# DD plots ----------------------------------------------------------------

beta <- map_dbl(results_cjspop$mcmc_sum, function(x) x["beta",1])
zeta <- map_dbl(results_cjspop$mcmc_sum, function(x) x["zeta",1])
sad <- map_dbl(results_cjspop$mcmc_sum, function(x) x["survival_ad",1])
sjuv <- map_dbl(results_cjspop$mcmc_sum, function(x) x["survival_juv",1])
fec <- map_dbl(results_cjspop$mcmc_sum, function(x) x["fecundity",1])

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
rho <- MCMCchains(NvsB$mcmc_samples, params = "rho")

## Check prior posterior overlap
MCMCtrace(NvsB$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(NvsB$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(NvsB$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "-0.565" ~ "(-0.836, -0.153)")
my_grob1 <-  grid.text(my_text1, x=0.35,  y=0.9, gp=gpar(col="orange", fontsize=8))

#qNlow <- map_dbl(Nchains, function(x) quantile(log(x$medianN), 0.25))
#qNhigh <- map_dbl(Nchains, function(x) quantile(log(x$medianN), 0.75))
#qBlow <- map_dbl(results_cjspop$mcmc_chains, function(x) quantile(x[,"beta"], 0.25))
#qBhigh <- map_dbl(results_cjspop$mcmc_chains, function(x) quantile(x[,"beta"], 0.75))

g9 <- ggplot(mapping = aes(x = log(spN), y = -beta)) +
  geom_point() +
  #geom_errorbar(mapping = aes(ymin = qBlow, ymax = qBhigh), alpha = 0.5) +
  #geom_errorbarh(mapping = aes(xmin = qNlow, xmax = qNhigh), alpha = 0.5) +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(beta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  scale_y_continuous(breaks = seq(-0.5,0.75, 0.25), limits = c(-0.25,0.75)) +
  annotation_custom(my_grob1) 

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
rho <- MCMCchains(NvsZ$mcmc_samples, params = "rho")

## Check prior posterior overlap
MCMCtrace(NvsZ$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(NvsZ$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(NvsZ$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.307" ~ "(-0.176, 0.695)")
my_grob1 <-  grid.text(my_text1, x=0.35,  y=0.9, gp=gpar(col="orange", fontsize=8))

g10 <- ggplot(mapping = aes(x = log(spN), y = -zeta)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(zeta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1)

# Growth rate vs DD in Survival
RvsB <- cor.jags(data = list(y = cbind(log(spR), beta), 
                             n = length(spR)),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

RvsB$mcmc_sum
rho <- MCMCchains(RvsB$mcmc_samples, params = "rho")
length(which(rho<0))/length(rho)

## Check prior posterior overlap
MCMCtrace(RvsB$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(RvsB$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(RvsB$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.542" ~ "(0.120, 0.816)")
my_grob1 <-  grid.text(my_text1, x=0.35,  y=0.9, gp=gpar(col="orange", fontsize=8))

g11 <- ggplot(mapping = aes(x = log(spR), y = -beta)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("r"))),
       y = bquote(beta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  scale_y_continuous(breaks = seq(-0.5,0.75, 0.25), limits = c(-0.25,0.75)) +
  annotation_custom(my_grob1)

# Population size vs DD in Fecundity
RvsZ <- cor.jags(data = list(y = cbind(log(spR), zeta), 
                             n = length(spR)),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

RvsZ$mcmc_sum
rho <- MCMCchains(RvsZ$mcmc_samples, params = "rho")

## Check prior posterior overlap
MCMCtrace(RvsZ$mcmc_samples, params = "mu", prior = mu_prior, pdf = F)
MCMCtrace(RvsZ$mcmc_samples, params = "sigma", prior = sigma_prior, pdf = F)
MCMCtrace(RvsZ$mcmc_samples, params = "rho", prior = rho_prior, pdf = F)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "-0.128" ~ "(-0.577, 0.355)")
my_grob1 <-  grid.text(my_text1, x=0.35,  y=0.9, gp=gpar(col="orange", fontsize=8))

g12 <- ggplot(mapping = aes(x = log(spR), y = -zeta)) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("r"))),
       y = bquote(zeta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1)

(g9 + g10) / (g11 + g12) +
  plot_annotation(tag_levels = 'a')

ggsave("paper_results/macro_fig4.tiff", width = 8, height = 6, units = "in")


# Explorations statistical patterns ---------------------------------------

pad <- map_dbl(results_cjspop$mcmc_sum, function(x) inv.logit(x["gamma[2]",1]))
pjuv <- map_dbl(results_cjspop$mcmc_sum, function(x) inv.logit(x["gamma[1]",1]))

# Capture probability vs N and R
tg1 <- ggplot(mapping = aes(x = 1-(1-pad)^4, y = log(spN))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Yearly Adult Capture Probability",
       y = bquote(bar(bar("N"))))

tg2 <- ggplot(mapping = aes(x = 1-(1-pjuv)^4, y = log(spN))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Yearly Juvenile Capture Probability",
       y = bquote(bar(bar("N"))))

tg3 <- ggplot(mapping = aes(x = 1-(1-pad)^4, y = log(spR))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Yearly Adult Capture Probability",
       y = bquote(bar(bar("r"))))

tg4 <- ggplot(mapping = aes(x = 1-(1-pjuv)^4, y = log(spR))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Yearly Juvenile Capture Probability",
       y = bquote(bar(bar("r"))))

(tg1 + tg2) / (tg3 + tg4) +
  plot_annotation(tag_levels = 'a')

ggsave("paper_results/macro_figureS2.tiff", width = 8, height = 8, units = "in")

# Patterns when low cp species are removed
index <- which(1-(1-pjuv)^4 >=0.05)

NvsNpop2 <- cor.jags(data = list(y = cbind(log(spN[index]), log(npop[index])), 
                                n = length(spN[index])),
                    inits = inits,
                    n.chains = 3,
                    n.adapt = 1000,
                    n.update = 10000,
                    n.iter = 5000,
                    n.thin = 5,
                    seed_no = 19)

NvsNpop2$mcmc_sum
rho <- MCMCchains(NvsNpop2$mcmc_samples, params = "rho")
length(which(rho>0))/length(rho)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.159" ~ "(-0.419, 0.665)")
my_grob1 <-  grid.text(my_text1, x=0.7,  y=0.25, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.716")
my_grob2 <-  grid.text(my_text2, x=0.7,  y=0.15, gp=gpar(col="orange", fontsize=8))

tg5 <- ggplot(mapping = aes(x = log(npop[index]), log(spN[index]))) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = "Number of Populations",
       y = bquote(bar(bar("N")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
annotation_custom(my_grob1) +
annotation_custom(my_grob2)

RvsN_sp2 <- cor.jags(data = list(y = cbind(log(spR[index]), log(spN[index])), 
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

## Plot
my_text1 <- bquote(rho ~ "=" ~ "-0.335" ~ "(-0.760, 0.261)")
my_grob1 <-  grid.text(my_text1, x=0.75,  y=0.9, gp=gpar(col="orange", fontsize=8))
my_text2 <- bquote("P("*rho*">0)" ~ "=" ~ "0.116")
my_grob2 <-  grid.text(my_text2, x=0.75,  y=0.8, gp=gpar(col="orange", fontsize=8))

tg6 <- ggplot(mapping = aes(x = log(spN)[index], y = log(spR)[index])) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(bar(bar("r")))) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1) +
  annotation_custom(my_grob2)

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

## Plot
my_text1 <- bquote(rho ~ "=" ~ "-0.450" ~ "(-0.822, -0.110)")
my_grob1 <-  grid.text(my_text1, x=0.25,  y=0.25, gp=gpar(col="orange", fontsize=8))

tg7 <- ggplot(mapping = aes(x = log(spN)[index], y = -beta[index])) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(beta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1)

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
rho <- MCMCchains(RvsB2$mcmc_samples, params = "rho")
length(which(rho<0))/length(rho)

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.443" ~ "(-0.124, 0.817)")
my_grob1 <-  grid.text(my_text1, x=0.75,  y=0.25, gp=gpar(col="orange", fontsize=8))

tg8 <- ggplot(mapping = aes(x = log(spR)[index], y = -beta[index])) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("r"))),
       y = bquote(beta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1)

NvsZ2 <- cor.jags(data = list(y = cbind(log(spN[index]), -zeta[index]), 
                             n = length(spN[index])),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

NvsZ2$mcmc_sum

## Plot
my_text1 <- bquote(rho ~ "=" ~ "0.512" ~ "(0.018, 0.847)")
my_grob1 <-  grid.text(my_text1, x=0.35,  y=0.9, gp=gpar(col="orange", fontsize=8))

tg9 <- ggplot(mapping = aes(x = log(spN)[index], y = -zeta[index])) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("N"))),
       y = bquote(zeta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1)

RvsZ2 <- cor.jags(data = list(y = cbind(log(spR[index]), -zeta[index]), 
                             n = length(spR[index])),
                 inits = inits,
                 n.chains = 3,
                 n.adapt = 1000,
                 n.update = 10000,
                 n.iter = 5000,
                 n.thin = 5,
                 seed_no = 19)

RvsZ2$mcmc_sum

## Plot
my_text1 <- bquote(rho ~ "=" ~ "-0.135" ~ "(-0.640, 0.430)")
my_grob1 <-  grid.text(my_text1, x=0.35,  y=0.9, gp=gpar(col="orange", fontsize=8))

tg10 <- ggplot(mapping = aes(x = log(spR[index]), y = -zeta[index])) +
  geom_point() +
  geom_smooth(method="lm", se=F, color = "orange") +
  labs(x = bquote(bar(bar("r"))),
       y = bquote(zeta)) + 
  theme(axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 12)) +
  annotation_custom(my_grob1)

(tg5 + tg6) / (tg7 + tg8) / (tg9 + tg10) +
  plot_annotation(tag_levels = 'a')

ggsave("paper_results/macro_figureS3.tiff", width = 8, height = 8, units = "in")


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

beta_low <- map_dbl(results_cjspop$mcmc_sum, function(x) x["beta",3])
zeta_low <- map_dbl(results_cjspop$mcmc_sum, function(x) x["zeta",3])
sad_low <- map_dbl(results_cjspop$mcmc_sum, function(x) x["survival_ad",3])
sjuv_low <- map_dbl(results_cjspop$mcmc_sum, function(x) x["survival_juv",3])
fec_low <- map_dbl(results_cjspop$mcmc_sum, function(x) x["fecundity",3])
pad_low <- map_dbl(results_cjspop$mcmc_sum, function(x) inv.logit(x["gamma[2]",3]))
pjuv_low <- map_dbl(results_cjspop$mcmc_sum, function(x) inv.logit(x["gamma[1]",3]))

beta_high <- map_dbl(results_cjspop$mcmc_sum, function(x) x["beta",5])
zeta_high <- map_dbl(results_cjspop$mcmc_sum, function(x) x["zeta",5])
sad_high <- map_dbl(results_cjspop$mcmc_sum, function(x) x["survival_ad",5])
sjuv_high <- map_dbl(results_cjspop$mcmc_sum, function(x) x["survival_juv",5])
fec_high <- map_dbl(results_cjspop$mcmc_sum, function(x) x["fecundity",5])
pad_high <- map_dbl(results_cjspop$mcmc_sum, function(x) inv.logit(x["gamma[2]",5]))
pjuv_high <- map_dbl(results_cjspop$mcmc_sum, function(x) inv.logit(x["gamma[1]",5]))

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

