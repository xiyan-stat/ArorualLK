######################################################
#    Function of LatticeKrig based Auroral Map       #
######################################################

LK.auroral.2 <- function(data1, data2, data3, nmlat=363, nmlt=363, time_pt, ratio1, ratio2, NC=NULL, nlevel=NULL, 
                       type = c("LK", "system", "random"), prob = TRUE, plot=TRUE){
  #####################################################################################################
  # Function: Combine interpolation (or/and empirical data) and satellite data to generate auroral maps
  # Note: data1 and data2 are including flux/energy, x, y (from mlat and mlt), time, 
  #       and source (inter, emp, or sate)
  # Input: data1 = a dataframe of temporal interpolation (or/and empirical data)  
  #        data2 = a dataframe of satellite data
  #        data3 = empirical data  
  #        nmlat = dimension of latitude
  #        nmlt = dimension of local time
  #        time_pt = a vector or scalar for time point of interest
  #        ratio1 = the ratio between interpolation  
  #                and satellite data used in LK estimation
  #        ratio2 = the ratio between empirical and interpolation data
  #        NL and nlevel are aurguments from LatticeKrig function: see LatticeKrig
  #        type = types of downsampling, "LK": latticeKrig, "system": systematical sampling, 
  #               "random": random sampling, etc.
  #        prob = generate probabilistic boundary or not using knn (default=TRUE)
  #        plot = generate auroral fitted maps (default=TRUE); if FLASE, no plot
  # Output: res = LK fitted auroral map data
  #####################################################################################################
  
  ## required packages
  require(LatticeKrig)
  require(tidyverse)
  require(class)
  
  `%notin%` <- Negate(`%in%`)
  if (ncol(data1) != 5 | ncol(data2) != 5 | ncol(data3) != 5) {print("Data are not valid")}
  if(ratio1 >= 1/5 | ratio1 < 1/10) {print("Ratio should be greater than 1/10 and less than 1/5")}
  if(type %notin% c("LK", "system", "random")) {print("type is not in c('LK', 'system', 'random')")}
  
  ## keep the data at time_pt
  # interpolation (empirical)
  inter <- data1 %>%
    filter(time %in% time_pt) %>%
    group_by(x, y, source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, x, y, source)
  
  # satellite
  sate <- data2 %>%
    filter(time %in% time_pt) %>%
    group_by(x, y,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, x, y, source)
  
  # empirical
  emp <- data3 %>%
  filter(time %in% time_pt) %>%
    group_by(x, y,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, x, y, source)
  
  # combine three data
  dat <- bind_rows(inter, sate, emp)
  
  ## dimensions
  n1 = nrow(inter)
  n2 = nrow(sate)
  n3 = nrow(emp)
  
  if (n2 < 30) print("Warning: the satellite data size at this time point is small")
  
  ## create test locations for generated map
  # high resolution
  nt <- ceiling(sqrt(nmlat*nmlt/0.55))  # 490
  x <- seq(-0.8, 0.8, len=nt)
  y <- seq(-0.8, 0.8, len=nt)
  
  test.location <- expand.grid(x, y)
  names(test.location) <- c("x", "y")
  
  test.location <- test.location %>% filter(x^2+y^2<=0.7^2)
  
  # low resolution: ratio of satellite data
  nt1 <- ceiling(sqrt(n2*ratio1/0.6))
  x1 <- seq(-0.8, 0.8, len=nt1)
  y1 <- seq(-0.8, 0.8, len=nt1)
  
  test.location.low <- expand.grid(x1, y1)
  names(test.location.low) <- c("x", "y")
  
  test.location.low <- test.location.low %>% filter(x^2+y^2<=0.7^2)
  
  #####################################################
  # generate probabilistic boundary
  if(prob){
    # for fit1
    if(type=="LK"){
      train1 <- inter %>% select(x,y)  # flux>30 -> flux =30 
      test1 <- test.location.low %>% select(x, y)
      cl1 <- factor(1*(inter$flux > 5))
      
      re.knn1 <- knn(train1, test1, cl1, k = 10, prob=TRUE)
      prob1 <- attributes(re.knn1)$prob
      cl1 <- as.numeric(re.knn1)-1
      predict.prob1 <- cl1*prob1 + (1-cl1)*(1-prob1)
      
      prob.low <- cbind(predict.prob1, test.location.low)
      #weight1 <- prob.low$predict.prob1 * 1.5  
      weight1 <- prob.low$predict.prob1
    }else{weight1 = 1}
    
    # for fit2
    train <- dat %>% select(x,y)  # flux>30 -> flux =30 
    test <- test.location %>% select(x, y)
    cl <- factor(1*(dat$flux > 5))
    
    re.knn <- knn(train, test, cl, k = 10, prob=TRUE)
    prob <- attributes(re.knn)$prob
    cl <- as.numeric(re.knn)-1
    predict.prob <- cl*prob + (1-cl)*(1-prob)
    
    prob.high <- cbind(predict.prob, test.location)
    #weight2 <- prob.high$predict.prob * 1.5  
    weight2 <- prob.high$predict.prob
  }else{
    weight1 = 1
    weight2 = 1
  }
  
  ###################################################
  
  ## Downsampling: three ways
  if (type == "LK"){
    inter0 <- inter %>% filter(flux!=0)
    # two inputs for LatticeKrig function
    location <- inter0 %>% select(x, y)   # or energy
    flux.n <- log(inter0$flux)
    
    if (is.null(nlevel) | is.null(NC)){ 
      fit1 <- LatticeKrig(location, flux.n)
    }else{
      fit1 <- LatticeKrig(location, flux.n, nlevel=nlevel,  NC=NC, findAwght=TRUE)
    }
    # conditional simulation
    simu1 <- LKrig.sim.conditional(fit1, x.grid= test.location.low, seed=1)
    # prediction Mean
    pred1 <- predict(fit1, x = test.location.low)
    
    #dat1 <- tibble(flux=exp(simu1$g.draw[,1])*weight1) 
    dat1 <- tibble(flux=exp(pred1[,1])*weight1)
    dat1 <- cbind(dat1, test.location.low)
    dat1$source <- "inter_LK"
  } else if(type == "system"){
    # systematic sampling from interpolation
    # depends on the lattitude and longitude
    index1 <- seq(1, n1, by= ceiling(n1/(n2*ratio1)))
    
    dat1 <- inter[index1,]
    dat1$source <- "inter_sys"
  }else if(type == "random"){
    index2 <- sample(1:n1, floor(n2*ratio1))
    
    dat1 <- inter[index2,]
    dat1$source <- "inter_random"
  }
  
  # downsampling of empirical data
  n4 <- length(index2)
  
  index3 <- sample(1:n3, floor(n4*ratio2))
  emp1 <- emp[index3,]
  emp1$source <- "emp_random"
  
  ## combine the downsampled interpolation (empirical) with satellite daya at this time point
  all <- rbind(emp1, dat1, sate) %>%
    group_by(x, y, source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup()
  
  all0 <- all %>% filter(flux!=0)  # energy
  
  location <- all0 %>% select(x, y)
  flux.n <- log(all0$flux)
  
  if (is.null(nlevel) | is.null(NC)){ 
    fit2 <- LatticeKrig(location, flux.n)
  }else{
    fit2 <- LatticeKrig(location, flux.n, nlevel=nlevel,  NC=NC, findAwght=TRUE)
  }
  cat(print(fit2))
  
  # one realization
  #simu2 <- LKrig.sim.conditional(fit2, x.grid= test.location, seed=1)
  # prediction mean
  pred2 <- predict(fit2, x = test.location)
  
  #dat2 <- tibble(flux=exp(simu2$g.draw[,1])*weight2) 
  dat2 <- tibble(flux=exp(pred2[,1])*weight2)
  dat2 <- cbind(dat2, test.location)
  dat2$source <- "simu"
  
  comb <- rbind(emp1, dat1, sate, dat2)
  
  ## generate a plot of fitted map, downsampled interpolation (empirical) and satellite observations
  if (plot==TRUE){
    p <- comb %>%
      mutate(
        flux = ifelse(flux>30, 30, flux)
      ) %>%
      ggplot(aes(x = x, y = y)) +
      geom_point(aes(colour = flux), size=0.5)+
      scale_colour_distiller(palette = "Spectral")+
      facet_grid(.~source) +
      theme_bw() +
      ggtitle(paste("Time point:",time_pt, ", Ratio:", ratio1, ratio2))
  }
  
  return(list(p=p, res=comb))
}

