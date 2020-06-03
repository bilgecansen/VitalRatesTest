
# This script is for wrangling capture history data of multiple species obtained from 
# Montioring Avian Productivity and Survivorship (MAPS) program. It re-shapes the data 
# into a suitable format for Bayesian mark-recapture analysis with JAGS.

# Loaded datasets (band_data, effort_data, midpoint_data, and weather_data) and functions 
# introduced below were already created in a different study. See Şen and Akçakaya 2020 
# for details about this data and fucntions used below.


# Set up directories, load packages and main datasets ---------------------

band_data <- readRDS("data_band.rds")
effort_data <- readRDS("data_effort.rds")
weather_data <- readRDS("data_weather.rds")
midpoint_data <- readRDS("data_midpoint.rds")

if (!"species" %in% list.files()) dir.create("species")

library(tidyverse)
library(foreign)
library(geosphere)
library(measurements)
library(raster)
library(foreach)
library(factoextra)


# Function 1: read.rawdata ------------------------------------------------

read.rawdata <- function(spcode, band_data, effort_data, weather_data) {
  
  # Subsetting the master band data by species code
  band_data <- filter(band_data, spec == spcode) %>%
    arrange(band, year, month)
  
  # Select only the breeding months (May~August)
  band_data <- filter(band_data, month %in% c(5,6,7,8))
  effort_data <- filter(effort_data, month %in% c(5,6,7,8)) %>%
    # subset effort to the populations species was captured
    filter(pop %in% band_data$pop) %>%
    arrange(pop, year, month)
  
  # Reformat days as May 1-day1 and August 31-day 123
  refdays <- function(x) {
    
    if(x[1] == 5) return(x[2])
    if(x[1] == 6) return(x[2] + 31)
    if(x[1] == 7) return(x[2] + 61)
    if(x[1] == 8) return(x[2] + 92)
    
  }
  band_data$day <- unlist(apply(as.matrix(band_data[,4:5]), 1, refdays))
  
  # Remove populations with less than 5 years of captures
  n_trend <- tapply(band_data$year, band_data$pop, function(x) length(unique(x)))
  pop_trend <- names(which(n_trend>4))
  band_data <- filter(band_data, pop %in% pop_trend)
  effort_data <- filter(effort_data, pop %in% pop_trend)
  
  # Subset weather data to the populations species was captured
  unique_pop <- unique(effort_data$pop)
  index <- which(rownames(weather_data[[1]]) %in% unique_pop)
  weather_data <- map(weather_data, function(x) x[index,])
  
  # Remove populations with no weather values 
  # (Here only AMT and the year 1992 is used to detect these populations)
  if (!is.null(nrow(weather_data[[1]]))) {
    
    index2 <- which(!is.na(weather_data[[1]][,1]))
    pop_names <- row.names(weather_data[[1]])[index2]
    band_data <- filter(band_data, pop %in% pop_names)
    effort_data <- filter(effort_data, pop %in% pop_names)
    weather_data <- map(weather_data, function(x) x[index2,])
  
  }
  
  # Centralize weather data
  weather_cent <- lapply(weather_data, function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T))
  
  # Do PCA on weather data
  data_pca <- map_dfc(weather_data, function(x) c(x))[,1:10]
  res_pca <- prcomp(data_pca, scale. = TRUE)
  
  # Select dimensions that explain 80% of the variance
  eig <- get_eigenvalue(res_pca)
  dim_num <- which(eig$cumulative.variance.percent>80)[1]
  
  # PCA values for each pop X year
  res_ind <- get_pca_ind(res_pca)
  npop <- nrow(weather_data$AMT)
  nyear <- ncol(weather_data$AMT)
  
  weather_pca <- list()
  for (i in 1:dim_num) {
    
    weather_pca[[i]] <- matrix(res_ind$coord[,i], nrow = npop, ncol = nyear)
    
  }
  
  data <- list()
  
  if (nrow(band_data)>0) {
    
    data$bands <- band_data
    data$effort <- effort_data
    data$weather_data <- weather_data
    data$weather_cent <- weather_cent
    data$weather_pca <- weather_pca
    data$res_pca <- res_pca
    
    return(data)
    
  } else {
    
    return(NULL)
  }
}


# Function 2: format.ch ---------------------------------------------------

format.ch <- function(bands) {
  
  unique_bands <- unique(bands$band)
  years <- 1992:2008
  months <- 5:8
  
  nband <- length(unique_bands)
  nyear <- length(years)
  nmonth <- length(months)
  
  # Capture history with sub-occasions
  ch_robust <- array(dim = c(nband, nyear, nmonth),
                     dimnames = list(unique_bands, years, months))
  
  for (i in 1:nband) {
    
    d <- filter(bands, band==unique_bands[i])
    
    for (h in 1:nyear) {
      
      if (years[h] %in% d$year) {
        
        d2 <- filter(d, year==years[h])
        
        for (k in 1:nmonth) {
          
          if (months[k] %in% d2$month) {
            
            ch_robust[i,h,k] <- 1
          
          } else ch_robust[i,h,k] <- 0
        
        }#k
      
      } else ch_robust[i,h,1:nmonth] <- 0
    
    }#h
  }#i
  
  # Capture history with only primary occasions
  ch_year <- apply(ch_robust, c(1,2), function(y) ifelse(any(y==1, na.rm = T), 1, 0))
  
  results <- list(ch_robust = ch_robust,
                  ch_year = ch_year,
                  unique_bands = unique_bands,
                  years = years,
                  months = months,
                  nband = nband,
                  nyear = nyear,
                  nmonth = nmonth)
  
  return(results)
}


# Function 3: format.effort -----------------------------------------------

format.effort <- function(effort, ch) {
  
  nyear <- ch$nyear
  nmonth <- ch$nmonth
  years <- ch$years
  months <- ch$months
  
  unique_pop <- unique(effort$pop)
  unique_pop <- unique_pop[order(unique_pop)]
  npop <- length(unique_pop)
  
  # Effort data with sub occasions
  effort_robust <- array(dim = c(npop, nyear, nmonth),
                         dimnames = list(unique_pop, years, months))
  
  for (i in 1:npop) {
    
    d <- filter(effort, pop==unique_pop[i])
    
    for (h in 1:nyear) {
      
      if (years[h] %in% d$year) {
        
        a <- 1
        d2 <- filter(d, year==years[h])
        
        for (k in 1:nmonth) {
          
          if (months[k] %in% d2$month) {
            
            effort_robust[i,h,k] <- d2$effort[a]
            
            a <- a + 1
          
          } else effort_robust[i,h,k] <- 0
        
        }#k
      
      } else effort_robust[i,h,1:nmonth] <- 0
    
    }#h
  }#i
  
  # Effort data with only primary occasions
  effort_year <- apply(effort_robust, c(1,2), sum)
  
  # First and last years with capture effort in a population
  first_pop <- apply(effort_year, 1, function(x) min(which(x>0)))
  last_pop <- apply(effort_year, 1, function(x) max(which(x>0)))
  
  # Centered effort data
  # Only data between first and last capture in a population is used
  # Data centered to have mean 0 and sd 0.5
  effort_cent <- effort_robust
  
  for (i in 1:npop) {
    
    if (first_pop[i]>1) effort_cent[i,1:(first_pop[i]-1),] <- NA
    if (last_pop[i]<nyear) effort_cent[i,(last_pop[i]+1):nyear,] <- NA
  
  }
  
  mean_effort <- mean(effort_cent, na.rm = T)
  sd_effort <- sd(effort_cent, na.rm = T)
  effort_cent <- (effort_cent - mean_effort)/(2*sd_effort)
  
  results <- list(effort_robust = effort_robust,
                  effort_year = effort_year,
                  effort_cent = effort_cent,
                  unique_pop = unique_pop,
                  first_pop = first_pop,
                  last_pop = last_pop,
                  npop = npop,
                  mean_effort = mean_effort,
                  sd_effort = sd_effort)
  
  return(results)
}


# Function 4: format.index ------------------------------------------------

format.index <- function(bands, ch, effort) {
  
  ch_robust <- ch$ch_robust
  ch_year <- ch$ch_year
  unique_bands <- ch$unique_bands
  nband <- ch$nband
  nyear <- ch$nyear
  years <- ch$years
  npop <- effort$npop
  
  # first and last capture years
  first <- apply(ch_year, 1, FUN = function(x) min(which(x==1)))
  last <- apply(ch_year, 1, FUN = function(x) max(which(x==1)))
  
  # first capture month
  detect_first_sub <- function(x) {
    
    f <- min(which(x==1, arr.ind = T)[,1])
    
    return(min(which(x[f,]==1)))
  
    } 
  first_sub <- apply(ch_robust, 1, detect_first_sub)
  
  # Population index
  pop_array <- array(dim = c(nband, nyear),
                     dimnames = list(unique_bands, years))
  
  for (i in 1:nband) {
    
    # count the number of times an individual was captured in a population
    pop_table <- table(filter(bands, band==unique_bands[i])$pop)
    
    # assign the population that has the max number of captures to the individual
    pop_array[i,] <- as.numeric(names(which(pop_table==max(pop_table))))[1]
  
  }
  
  pop <- as.factor(pop_array[,1])
  levels(pop) <- 1:npop
  
  results <- list(first = first,
                  last = last,
                  first_sub = first_sub,
                  pop = pop,
                  pop_array = pop_array)
  
  return(results)
}


# Function 5: format.stage ------------------------------------------------

format.stage <- function(bands, ch) {
  
  nband <- ch$nband
  nyear <- ch$nyear
  years <- ch$years
  unique_bands <- ch$unique_bands
  
  stage_array <- array(dim = c(nband, nyear),
                       dimnames = list(unique_bands, years))
  
  for (i in 1:nband) {
    
    d <- filter(bands, band==unique_bands[i])
    
    for (h in 1:nyear) {
      
      if (years[h] %in% d$year) {
        
        d2 <- filter(d, year==years[h])
        # all stages are the same in a year so take first one
        stage_array[i,h] <- ifelse(d2$age[1]=="A", 2, 1)
      
      } else stage_array[i,h] <- NA
    
    }#h
  }#i
  
  age.f <- function(x) {
    first <- which(!is.na(x))[1]
    # every year after first capture, it's an adult
    if(first<length(x)) x[(first+1):length(x)] <- 2
    return(x)
  }
  
  stage_array2 <- t(apply(stage_array, 1, age.f))
  
  return(stage_array2)
}


# Function 6: format.residents --------------------------------------------

format.residents <- function(bands, ch, index) {
  
  nband <- ch$nband
  nyear <- ch$nyear
  years <- ch$years
  unique_bands <- ch$unique_bands
  first <- index$first
  
  residents <- vector(mode = "numeric", length = nband)
  for (i in 1:nband) {
    
    d <- filter(bands, band==unique_bands[i] & year==years[first[i]])
    residents[i] <- ifelse((max(d$day)-min(d$day))>=10, 1, 0)
  
  }
  
  return(residents)
}


# Function 7: format.Nobs -------------------------------------------------

format.Nobs <- function(effort, ch, stage, index) {
  
  ch_year <- ch$ch_year
  nyear <- ch$nyear
  years <- ch$years
  unique_pop <- effort$unique_pop
  npop <- effort$npop
  pop_array <- index$pop_array
  
  stage_ad <- stage
  stage_ad[stage_ad==1] <- 0
  stage_ad[stage_ad==2] <- 1
  stage_ad[is.na(stage_ad)] <- 0
  
  stage_juv <- stage
  stage_juv[stage_juv==2] <- 0
  stage_juv[is.na(stage_juv)] <- 0
  
  Nobs_ad <- array(dim = c(npop, nyear),
                   dimnames = list(unique_pop, years))
  
  Nobs_juv <- array(dim = c(npop, nyear),
                    dimnames = list(unique_pop, years))
  
  for (i in 1:npop) {
    
    pop <- pop_array
    pop[pop!=unique_pop[i]] <- 0
    pop[pop==unique_pop[i]] <- 1
    
    ab_adult <- pop*ch_year*stage_ad
    ab_adult <- apply(ab_adult, 2, sum)
    Nobs_ad[i,] <- ab_adult
    
    ab_juv <- pop*ch_year*stage_juv
    ab_juv <- apply(ab_juv, 2, sum)
    Nobs_juv[i,] <- ab_juv
  
  }
  
  results <- list(Nobs_ad = Nobs_ad,
                  Nobs_juv = Nobs_juv)
  
  return(results)
}


# Function 8: format.fec --------------------------------------------------

format.fec <- function(effort, Nobs, rawdata) {
  
  Nobs_ad <- Nobs$Nobs_ad
  Nobs_juv <- Nobs$Nobs_juv
  first_pop <- effort$first_pop
  npop <- effort$npop
  weather_cent <- rawdata$weather_cent
  weather_pca <- rawdata$weather_pca
  
  # Assign 0 population size to first capture year 
  for (i in 1:npop) {
    Nobs_ad[i,first_pop[i]] <- 0
    Nobs_juv[i,first_pop[i]] <- 0
  }
  
  fec <- Nobs_juv/Nobs_ad
  fec2 <- fec[!is.infinite(fec) & !is.na(fec)]
  out <- boxplot(fec2)$out
  index <- !is.infinite(fec) & !is.na(fec) & fec<min(out)
  
  # Population and year indices for fecundity data
  pop_index <- which(index, arr.ind = T)[,1]
  year_index <- which(index, arr.ind = T)[,2]
  
  # Fecundity data used in mark-recapture models
  fec_data <- Nobs_juv[index]
  
  # Weather data used in mark-recapture models
  weather_fec_pca <- lapply(weather_pca, function(x) x[index])
  
  results <- list(fec_data = fec_data,
                  pop_index = pop_index,
                  year_index = year_index,
                  weather_fec_pca = weather_fec_pca,
                  nfec = length(fec_data))
  
  return(results)
}


# Function 9: wrap.chdata -------------------------------------------------

wrap.chdata <- function(spcode, band_data, effort_data, weather_data) {
  
  rawdata <- read.rawdata(spcode, band_data, effort_data, weather_data)
  ch <- format.ch(rawdata$bands)
  effort <- format.effort(rawdata$effort, ch)
  index <- format.index(rawdata$bands, ch, effort)
  stage <- format.stage(rawdata$bands, ch)
  residents <- format.residents(rawdata$bands, ch, index)
  Nobs <- format.Nobs(effort, ch, stage, index)
  fec <- format.fec(effort, Nobs, rawdata)
  
  # Latent state set up
  latent_state <- ch$ch_year  
  for (i in 1:nrow(ch$ch_year)) {
    n1 <- min(which(ch$ch_year[i,]==1))
    n2 <- max(which(ch$ch_year[i,]==1))
    latent_state[i, n1:n2] <- 1
    latent_state[i, n1] <- NA 
  }
  latent_state[latent_state==0] <- NA
  
  results <- list()
  results$ch_robust <- ch$ch_robust
  results$ch_year <- ch$ch_year
  results$latent_state <- latent_state
  results$stage <- stage
  results$effort_cent <- effort$effort_cent
  results$effort_year <- effort$effort_year
  results$residents <- residents
  results$first <- index$first
  results$last <- index$last
  results$first_sub <- index$first_sub
  results$first_pop <- effort$first_pop
  results$last_pop <- effort$last_pop
  results$pop <- index$pop
  results$Nobs_ad <- Nobs$Nobs_ad
  results$Nobs_juv <- Nobs$Nobs_juv
  results$fec_data <- fec$fec_data
  results$pop_index <- fec$pop_index
  results$year_index <- fec$year_index
  results$nyear <- ch$nyear
  results$nind <- ch$nband
  results$npop <- effort$npop
  results$nfec <- fec$nfec
  results$mean_effort <- effort$mean_effort
  results$sd_effort <- effort$sd_effort
  results$weather <- rawdata$weather_data
  results$weather_cent <- rawdata$weather_cent
  results$weather_pca <- rawdata$weather_pca
  results$weather_fec_pca <- fec$weather_fec_pca
  results$res_pca <- rawdata$res_pca
  
  return(results)
}


# Wrangle chdata ----------------------------------------------------------

spcode <- as.character(unique(band_data$spec))
ss <- c()
pb <- txtProgressBar(min = 0, max = length(spcode),  style = 3)
for (i in 1:length(spcode)) {
  
  d_juv <- filter(band_data, spec==spcode[i], age == "J")
  l_juv <- length(unique(d_juv$band))
  d_ad <- filter(band_data, spec==spcode[i], age == "A")
  l_ad <- length(unique(d_ad$band))
  
  # First set of filters
  ## 1. Less than 75 adult captured individuals 
  if (l_ad<75) next
  
  ## 2. No juvenile recaptures
  recaps <- group_by(d_juv, band) %>% 
    summarize(nrecap = length(unique(month)))
  if (all(recaps$nrecap==1)) next
  
  try({
    chdata <- wrap.chdata(spcode[i], band_data, effort_data, weather_data)
  }, silent = T)
  if (!exists("chdata")) next
  if (is.null(chdata)) next
  
  # Second set of filters
  juv_index <- apply(chdata$stage, 1, function(x) any(x==1, na.rm = T)) %>% 
    which(.)
  l_juv2 <- length(juv_index)
  
  ## 1. Less than 75 adult captured individuals 
  l_ad2 <- nrow(chdata$stage) - l_juv2
  if (l_ad2<75) next
  
  ## 2. Less than 15 juvenile recaptures
  juv_ch <- chdata$ch_year[juv_index,]
  juv_recaps <- apply(juv_ch, 1, function(x) sum(x)>1) %>%
    which() %>%
    length()
  if (juv_recaps<15) next
  
  ## 3. Less than 15 adult recaptures
  ad_ch <- chdata$ch_year[-juv_index,]
  ad_recaps <- apply(ad_ch, 1, function(x) sum(x)-1) %>%
    sum()
  if (ad_recaps<15) next
  
  ss[i] <- nrow(chdata$ch_year)
  
  folder <- paste("species", spcode[i], sep = "/")
  if(!dir.exists(folder)) dir.create(folder)
  
  file_name <- paste(spcode[i], "chdata.rds", sep = "_")
  saveRDS(chdata, file = paste(folder, file_name, sep = "/"))
  setTxtProgressBar(pb, i)
}

# Record species sample size and number of pca dimensions selected
ss_ch <- c()
nvar <- c()
ss_fec <- c()
fec_list <- list()
splist <- list.files("species")
for (i in 1:length(splist)) {
  filename <- paste(splist[i], "chdata.rds", sep = "_") %>%
    paste("species", splist[i], ., sep = "/")
  z <- readRDS(filename)
  
  fec <- c()
  for (h in 1:z$nfec) {
    
    Nobs_ad <- z$Nobs_ad[z$pop_index[h], z$year_index[h]]
    fec[h] <- z$fec_data[h]/Nobs_ad
    
  }
  
  fec_list[[i]] <- fec
  ss_fec[i] <- length(fec[fec>0])
  nvar[i] <- length(z$weather_pca)
  ss_ch[i] <- nrow(z$ch_year)
}  

ssdata <- data.frame(spcode = splist, ss_ch = ss_ch, ss_fec = ss_fec, nvar = nvar) %>%
  arrange(ss_ch)

# Remove species with fecundity data with less than 5 non-zero obs per pca dimension
index <- which((ssdata$ss_fec/ssdata$nvar)<5)
ssdata2 <- ssdata[-index,]

saveRDS(ssdata2, file = "data_ss.rds")


# Assign single pop to species recorded in multi-pop ----------------------

# KEWA
chdata <- readRDS("species/KEWA/KEWA_chdata.rds")
chdata$pop[names(chdata$pop)==165171230] <- 37
index <- which(names(chdata$year_index)==438 & chdata$year_index==10)
chdata$year_index <- chdata$year_index[-index]
chdata$pop_index <- chdata$pop_index[-index]
chdata$fec_data <- chdata$fec_data[-index]
chdata$nfec <- chdata$nfec-1
saveRDS(chdata, file = "species/KEWA/KEWA_chdata.rds")

# BEWR
chdata <- readRDS("species/BEWR/BEWR_chdata.rds")
index <- which(names(chdata$year_index)==22 & chdata$year_index==16)
chdata$year_index <- chdata$year_index[-index]
chdata$pop_index <- chdata$pop_index[-index]
chdata$fec_data <- chdata$fec_data[-index]
chdata$nfec <- chdata$nfec-1
saveRDS(chdata, file = "species/BEWR/BEWR_chdata.rds")

# BCCH
chdata <- readRDS("species/BCCH/BCCH_chdata.rds")
index <- which(names(chdata$year_index)==22 & chdata$year_index==16)
chdata$year_index <- chdata$year_index[-index]
chdata$pop_index <- chdata$pop_index[-index]
chdata$fec_data <- chdata$fec_data[-index]
chdata$nfec <- chdata$nfec-1
saveRDS(chdata, file = "species/BCCH/BCCH_chdata.rds")


# HOWR
chdata <- readRDS("species/HOWR/HOWR_chdata.rds")
chdata$pop[names(chdata$pop)==214004360] <- 65
saveRDS(chdata, file = "species/HOWR/HOWR_chdata.rds")
