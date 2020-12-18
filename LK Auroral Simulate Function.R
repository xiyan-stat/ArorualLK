######################################################################
# Auroral Data 10mins, 30 min interval
######################################################################

LK.auroral.10m <- function(data1, data2, data3, data4, 
                         nmlat=363, nmlt=363, time_pt, 
                         ratio1=1/6, ratio2=1/10, ratio3=1/2, 
                         NCs=c(5,20,30), nlevels= c(1,1,3), kn = c(60,20,10), bmethod){
  ################################################################################
  # NCs = a vector of NC for low, medium and high scales
  # nlevels = a vector of nlevel for low, medium and high scales
  # kn = a vector of knn for low, medium and high scales
  # bmethod = boundary method for low medium and high scales
  ################################################################################
  
  ## required packages
  require(LatticeKrig)
  require(tidyverse)
  require(class)
  #require(sp)
  
  `%notin%` <- Negate(`%in%`)
  if (ncol(data1) != 5 | ncol(data2) != 5 | ncol(data3) != 5) {print("Data are not valid")}
  if(ratio1 >= 1/5 | ratio1 < 1/10) {print("Ratio should be greater than 1/10 and less than 1/5")}
  
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
  
  grd <- data4 %>%
    filter(time %in% time_pt) %>%
    group_by(x, y,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, x, y, source)
  
  ############################################
  ## create test locations for generated map
  # high resolution
  nt <- ceiling(sqrt(nmlat*nmlt))  # 490
 # nt <- 300
  x <- seq(-0.8, 0.8, len=nt)
  y <- seq(-0.8, 0.8, len=nt)
  
  test.location <- expand.grid(x, y)
  names(test.location) <- c("x", "y")
  
  test.location <- test.location %>% filter(x^2+y^2<=0.7^2)
  
  ####################################################
  # outer boundary
  load("boundary.RData")
  b.data <- data.frame(test.location)
  tt <- time_pt/60
  b.data$index <- ifelse((test.location$x-bndry$am[tt])^2 + (test.location$y-bndry$bm[tt])^2 <= bndry$rmax[tt]^2, 1, 0)
  outer <- b.data$index
  
  ###################################################
  ## find the extroplation
  time_ground <- seq(360, 650, 10)
  
  if(time_pt %in% time_ground){
    obs <- bind_rows(sate, inter, grd)
  }else{
    obs <- bind_rows(sate, inter)
  }
  x.coord <- obs %>% pull(x)
  y.coord <- obs %>% pull(y)
  intro.coords <- cbind(x.coord, y.coord)
  
  X <- unique.matrix(intro.coords)
  library(alphahull)
  inregion <- ahull(X, alpha=0.2)  # adjustable para
  
  # extrapolation locations' indicators
  extrap <- ifelse(inahull(inregion, as.matrix(test.location))==FALSE, 1, 0)
  
  #########################################
  # linear rescaling with sate data
  
  #beta0 <- 0.4
  #beta1 <- 2.7
  hours <- round(time_pt/60,0)
  beta0 <- Res.mean[[(1+5*hours)]]
  beta1 <- Res.mean[[(2+5*hours)]]
  emp <- emp %>%
    mutate(flux = beta0 + beta1*flux)
  
  ##########################################
  # LK for empirical data
  emp0 <- emp %>% filter(flux!=0)  
  
  location.emp <- emp0 %>% select(x, y)
  flux.emp <- log(emp0$flux)
  
  fit1 <- LatticeKrig(location.emp, flux.emp)
  cat(print(fit1), "\n")
  pred1 <- predict(fit1, x = test.location)
  
  emp.all <- tibble(flux=exp(pred1[,1]))
  emp.all <- cbind(emp.all, test.location, extrap)
  emp.all$source <- "emp"
  
  ## dimensions
  n1 = nrow(inter)
  n2 = nrow(sate)
  n3 = nrow(emp.all)
  
  if (n2 < 30) print("Warning: the satellite data size at this time point is small")
  
  ###################################################
  ## Downsampling: random sampling
  index1 <- sample(1:n1, floor(n2*ratio1))
  dat1 <- inter[index1,]
  dat1$source <- "inter_random"
  
  # downsampling of empirical data
  emp.in <- emp.all[emp.all$extrap==0,]
  emp.out <- emp.all[emp.all$extrap==1,]
  n.in <- nrow(emp.in)
  n.out <- nrow(emp.out)
  
  n5 <- length(index1)
  index2 <- sample(1:n.in, floor(n5*ratio2))
  emp.in1 <- emp.in[index2,]
  
  index.out <- sample(1:n.out, floor(n.out*0.1))  #ratio can change
  emp.out1 <- emp.out[index.out,]
  
  emp1 <- bind_rows(emp.in1, emp.out1)
  emp1$source <- "emp_random"
  
  emp1 <- emp1 %>% select(flux, x, y, source)
  
  
  ########################################
  # ground base data
  if (time_pt %in% time_ground){
    n4 = nrow(grd)
    
    # downsampling of ground data
    index3 <- sample(1:n4, floor(n2*ratio3))
    grd1 <- grd[index3,]
    grd1$source <- "ground_random"
    
    # combine all data
    dat <- bind_rows(inter, sate, emp1, grd)
    
    ## combine the downsampled interpolation (empirical) with satellite daya at this time point
    all <- rbind(grd1, emp1, dat1, sate) %>%
      group_by(x, y, source) %>%
      summarize(
        flux = mean(flux)
      ) %>%
      ungroup()
  }else{
    print("No observed ground data for this time")
    
    dat <- bind_rows(inter, sate, emp1)  # change the empirical training
    
    all <- rbind(emp1, dat1, sate) %>%
      group_by(x, y, source) %>%
      summarize(
        flux = mean(flux)
      ) %>%
      ungroup()
  }
  
  #####################################################
  # generate probabilistic boundary
  # for weights
  train <- dat %>% select(x,y)  # flux>30 -> flux =30
  test <- test.location %>% select(x, y)
  cl <- factor(1*(dat$flux > 5))
  
  # low scale
  re.knn.l <- knn(train, test, cl, k = kn[1], prob=TRUE)
  # medium scale
  if(kn[2]==kn[1]){
    re.knn.m <- re.knn.l
  }else{
    re.knn.m <- knn(train, test, cl, k = kn[2], prob=TRUE)
  }
  #high scale (full scale)
  if(kn[3]==kn[2]){
    re.knn.h <- re.knn.m
  }else{
    re.knn.h <- knn(train, test, cl, k = kn[3], prob=TRUE)
  }
  
  #re.knn.h <- knn(train, test, cl, k = kn[3], prob=TRUE)
  
  prob.l <- attributes(re.knn.l)$prob
  prob.m <- attributes(re.knn.m)$prob
  prob.h <- attributes(re.knn.h)$prob
  
  cl.l <- as.numeric(re.knn.l)-1
  cl.m <- as.numeric(re.knn.m)-1
  cl.h <- as.numeric(re.knn.h)-1
  
  predict.prob.l <- cl.l*prob.l + (1-cl.l)*(1-prob.l)
  predict.prob.m <- cl.m*prob.m + (1-cl.m)*(1-prob.m)
  predict.prob.h <- cl.h*prob.h + (1-cl.h)*(1-prob.h)
  
  
  #######################################
  all0 <- all %>% filter(flux!=0)  
  
  location <- all0 %>% select(x, y)
  flux.n <- log(all0$flux)
  ########################################
  # combine test.location with observed location
  # len1 <- nrow(test.location)
  # len0 <- nrow(location)
  # 
  # # test.location1 <- bind_rows(test.location, location)
  # # predict.prob.h1 <- c(predict.prob.h, rep(1, len0))
  # 
  # test.location1 <- test.location
  # predict.prob.h1 <- predict.prob.h
  
  prob.bd.l <- cbind(predict.prob.l, test.location)
  prob.bd.m <- cbind(predict.prob.m, test.location)
  prob.bd.h <- cbind(predict.prob.h, test.location)
  
  weight.l <- prob.bd.l$predict.prob.l
  weight.m <- prob.bd.m$predict.prob.m
  weight.h <- prob.bd.h$predict.prob.h
  
  indicator.l <- ifelse(weight.l>0, 1, 0)
  indicator.m <- ifelse(weight.m>0, 1, 0)
  indicator.h <- ifelse(weight.h>0, 1, 0)
  
  # add boundary
  if(bmethod=="1"){
    weight.b.l <- ifelse(outer==1, weight.l>0, 0)
    weight.b.m <- ifelse(outer==1, weight.m>0, 0)
    weight.b.h <- ifelse(outer==1, weight.h>0, 0)
  } 
  if(bmethod=="2"){
    train <- dat %>% select(x,y)  # flux>30 -> flux =30
    test <- test.location %>% select(x, y)
    cl <- factor(1*(dat$flux > 5))
    
    re.knn.s <- knn(train, test, cl, k = 3, prob=TRUE)
    prob.s <- attributes(re.knn.s)$prob
    cl.s <- as.numeric(re.knn.s)-1
    predict.prob.s <- cl.s*prob.s + (1-cl.s)*(1-prob.s)
    weight2.s <- predict.prob.s
    
    indicator2.s <- ifelse(weight2.s>0, 1, 0)
    
    weight.b.l <- ifelse(outer==1, indicator.l, indicator2.s)
    weight.b.m <- ifelse(outer==1, indicator.m, indicator2.s)
    weight.b.h <- ifelse(outer==1, indicator.h, indicator2.s)
  }
  if(bmethod=="3"){
    weight.b.l <- ifelse(outer==1, indicator.l, weight.l^e)
    weight.b.m <- ifelse(outer==1, indicator.m, weight.m^e)
    weight.b.h <- ifelse(outer==1, indicator.h, weight.h^e)
  }
  
  
  fit.l <- LatticeKrig(location, flux.n, nlevel=nlevels[1],  NC=NCs[1])
  fit.m <- LatticeKrig(location, flux.n, nlevel=nlevels[2],  NC=NCs[2])
  #fit.h <- LatticeKrig(location, flux.n, nlevel=nlevels[3],  NC=NCs[3])
  fit.h <- LatticeKrig(location, flux.n)
  cat(print(fit.l),"\n")
  cat(print(fit.m),"\n")
  cat(print(fit.h),"\n")
  
  
  ############################### 
  # one realization
  # simu2 <- LKrig.sim.conditional(fit.h, x.grid= test.location, seed=123) # 1,2,1234,6789,4135062
  # flux_simu <- tibble(flux=exp(simu2$g.draw[,1])*(weight.h>0))
  # flux_simu <- cbind(flux_simu, test.location)
  # flux_simu$source <- "simulation 1"
  
  ################################# 
  # standard error
  # low resolution: ratio of satellite data
  # nt1 <- ceiling(sqrt(n2*ratio1/0.6))
  # x1 <- seq(-0.8, 0.8, len=nt1)
  # y1 <- seq(-0.8, 0.8, len=nt1)
  # 
  # test.location.low <- expand.grid(x1,   y1)
  # names(test.location.low) <- c("x", "y")
  # 
  # test.location.low <- test.location.low %>% filter(x^2+y^2<=0.7^2)
  # 
  # sim2 <- LKrig.sim.conditional(fit1, x.grid= test.location.low, M=100)
  # # sim2 <- LKrig.sim.conditional(fit2, nx=50, ny=50, M=3)
  # dat3 <- tibble(SE = sim2$SE)
  # # dat3 <- tibble(flux = sim2$SE)
  # dat3 <- cbind(dat3, test.location.low)
  # dat3$source <- "uncertainty"
  
  #############################
  # prediction mean
  pred.l <- predict(fit.l, x = test.location)
  pred.m <- predict(fit.m, x = test.location)
  pred.h <- predict(fit.h, x = test.location)
  
  
  # dat.l <- tibble(flux=exp(pred.l[,1])*((weight.l>0)))
  # dat.l <- cbind(dat.l, test.location)
  # dat.l$source <- "a.Low"
  # dat.m <- tibble(flux=exp(pred.m[,1])*((weight.m>0)))
  # dat.m <- cbind(dat.m, test.location)
  # dat.m$source <- "b.Mid"
  # 
  # dat.h <- tibble(flux=exp(pred.h[,1])*((weight.h>0))) 
  # dat.h <- cbind(dat.h, test.location)
  # dat.h$source <- "c.High"
  
  dat.l <- tibble(flux=exp(pred.l[,1])*((weight.b.l)))
  dat.l <- cbind(dat.l, test.location)
  dat.l$source <- "a.Low"
  dat.m <- tibble(flux=exp(pred.m[,1])*((weight.b.m)))
  dat.m <- cbind(dat.m, test.location)
  dat.m$source <- "b.Mid"
  
  dat.h <- tibble(flux=exp(pred.h[,1])*((weight.b.h))) 
  dat.h <- cbind(dat.h, test.location)
  dat.h$source <- "c.High"
  
  #######################################
  # true obs
  # # observed <- all %>% select(flux, x, y)
  # observed <- sate %>% select(flux, x, y)
  # # dat.h1 <- rbind(observed, dat.h)
  # x.coord1 <- sate %>% pull(x)
  # y.coord1 <- sate %>% pull(y)
  # intro.coords1 <- cbind(x.coord1, y.coord1)
  # 
  # X1 <- unique.matrix(intro.coords1)
  # inregion1 <- ahull(X1, alpha=0.01)  # adjustable para
  # 
  # # extrapolation locations' indicators
  # extrap1 <- ifelse(inahull(inregion1, as.matrix(test.location))==FALSE, 1, 0)
  # 
  # dat.h1 <- cbind(dat.h, extrap1)
  # 
  # dat.h2 <- dat.h1 %>% 
  #           filter(extrap1==1) %>%
  #           select(flux, x, y)
  # dat.h3 <- rbind(observed, dat.h2)
  # dat.h3$source <- "d. High+Sate"
  # 
  #dat.h$source <- "simulated map"
  
  #################
  all$source <- "samples"
  
  # dat.h1 <- dat.h
  # dat.h1 <- dat.h1[!duplicated(dat.h1[,2:3]),]
  # dat.h1$source <- "zzz"
  
  #comb <- rbind(all, dat.h, dat.h1)
  #comb <- rbind(all, dat.h)
  comb <- rbind(dat.l, dat.m, dat.h)
  
  ## generate a plot of fitted map
  p <- comb %>%
    mutate(
      flux = ifelse(flux>30, 30, flux)
    ) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(colour = flux), size=0.01)+
    scale_colour_distiller(palette = "Spectral")+
    facet_grid(.~source) +
    theme_bw() +
    ggtitle(paste("UT = 7:30"))
  
  return(list(p=p, res=comb))
}


