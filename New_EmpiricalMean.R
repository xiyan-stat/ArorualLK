#######################################################
# Use the updated data to generate the auroral map
#######################################################
rm(list=ls())

library(tidyverse)
require(class)

#################################################
## Ground base data
load("flux_north_10_new.RData")
data.flux.grd_new <- data.flux.grd_n

################################################
## Satellite and Interpolated Data
load("flux_north_10.RData")

#1. satellite data: data.flux.sate_n
#2. interpolated data: data.flux.inter_n
# Empirical Input

get_HP <- function(kp) {
  A <- c(16.8244, 0.323365, -4.86128)
  B <- c(1.82336, 0.613192, 26.1798)
  if (kp <= 6.9) {
    HP <- A[1] * exp(A[2] * kp) + A[3]
  } else {
    HP <- B[1] * exp(B[2] * kp) + B[3]
  }
  return (HP)
}

get_fitting_value <- function(t, r, AA, BB) {
  ang <- (0: (nrow(AA)-1)) * t
  cosang <- cos(ang)
  sinang <- sin(ang)
  
  ncoeff <- ncol(AA)
  coeff <- array(dim = ncoeff)
  for (i in 1: ncoeff) coeff[i] <- sum(AA[, i] * cosang + BB[, i] * sinang)
  
  F1 <- exp((r - coeff[2]) / coeff[3])
  F2 <- 1 + exp((r - coeff[2]) / coeff[4])
  value <- coeff[1] * F1 / F2^2
  
  return (value)
}

get_guvi_kp_model <- function(kp, mlat, mlt) {
  library(package = "ncdf4")
  nc <- nc_open(filename = 'Data Updated/guvi_auroral_model.nc')
  kp_list <- ncdf4::ncvar_get(nc = nc, varid = "kp")
  AA_flux <- ncdf4::ncvar_get(nc = nc, varid = "AA_flux")
  BB_flux <- ncdf4::ncvar_get(nc = nc, varid = "BB_flux")
  AA_energy <- ncdf4::ncvar_get(nc = nc, varid = "AA_energy")
  BB_energy <- ncdf4::ncvar_get(nc = nc, varid = "BB_energy")
  ncdf4::nc_close(nc = nc)
  
  nkp <- length(kp_list)
  kp_list_mid <- array(dim = nkp-1)
  for (i in 1: (nkp-1)) kp_list_mid[i] <- (kp_list[i] + kp_list[i+1]) / 2
  
  ida <- which.min(abs(kp_list_mid - kp))
  if (kp_list_mid[ida] > kp) ida <- ida - 1
  idb <- ida + 1
  
  if (idb == nkp-1) {
    idb <- nkp - 2
    ida <- idb - 1
  }
  
  x0 <- (kp_list_mid[idb] - kp) / (kp_list_mid[idb] - kp_list_mid[ida])
  x1 <- (kp - kp_list_mid[ida]) / (kp_list_mid[idb] - kp_list_mid[ida])
  
  HPa <- get_HP(kp_list_mid[ida])
  HPb <- get_HP(kp_list_mid[idb])
  HP <- get_HP(kp)
  
  xa <- (HPb - HP) / (HPb - HPa)
  xb <- (HP - HPa) / (HPb - HPa)
  
  n <- length(mlat)
  flux <- array(dim = n)
  energy <- array(dim = n)
  for (i in 1: n) {
    if (mlat[i] < 50) next
    r <- 90 - mlat[i]
    t <- mlt[i] * pi/12
    Va <- get_fitting_value(t, r, AA_flux[, ida, ], BB_flux[, ida, ])
    Vb <- get_fitting_value(t, r, AA_flux[, idb, ], BB_flux[, idb, ])
    flux[i] <- Va*xa + Vb*xb
    Va <- get_fitting_value(t, r, AA_energy[, ida, ], BB_energy[, ida, ])
    Vb <- get_fitting_value(t, r, AA_energy[, idb, ], BB_energy[, idb, ])
    energy[i] <- Va*x0 + Vb*x1
  }
  
  return (list(flux = flux,energy = energy))
}

Emp_Data <- function(index, time, location){
  kp <- approx(x = seq(from = 0, to = 1440, by = 180), y = c(1.3, 6, 6, 5, 4.7, 3.7, 2, 3, 3), xout = time)$y
  
  x <- location[index, ]$x
  y <- location[index, ]$y
  mlat <- 90 -  sqrt(x^2 + y^2)*180/pi
  mlt <- atan2(y, x) * 12/pi
  
  emp <- get_guvi_kp_model(kp, mlat, mlt)$flux   # flux or energy
  return(emp)
}

## Empirical data
  
  # kpt <- approx(x = seq(from = 0, to = 1440, by = 180), y = c(1.3, 6, 6, 5, 4.7, 3.7, 2, 3, 3), xout = time_pt)$y
  # x_e <- c(sate$x, grd$x)
  # y_e <- c(sate$y, grd$y)
  # flux_e <- c(sate$x, grd$y)
  # valid1 <- flux_e > 0
  # valid1[is.na(valid1)] <- FALSE
  # index.e <- length(flux_e[valid1])
  # loc.e <- data.frame(x=x_e[valid1], y=y_e[valid1])
  # if (sum(valid1) >= 1) {
  #     emp_flux <- Emp_Data(1:index.e, time_pt,loc.e)
  #     ratio_e <- median(flux_e[valid1]) / median(emp_flux)
  #   }


# Funtion using the empirical as mean

LK.Empmean.10m <- function(data1, data2,data4, 
                           nmlat=363, nmlt=363, time_pt, 
                           ratio1=1/6, ratio4=1/2,ratio5=1/10,
                           NCs=c(5,20,30), nlevels= c(1,1,3), kn = c(60,20,10)){
  ################################################################################
  # NCs = a vector of NC for low, medium and high scales
  # nlevels = a vector of nlevel for low, medium and high scales
  # 1. inter, 2. emp.in, 3. emp.out, 4. grd, 5. inter itself
  # kn = a vector of knn for low, medium and high scales
  # bmethod = boundary method for low medium and high scales
  ################################################################################
  
  ## required packages
  require(LatticeKrig)
  require(tidyverse)
  require(class)
  #require(sp)
  
  `%notin%` <- Negate(`%in%`)
  if (ncol(data1) != 5 | ncol(data2) != 5) {print("Data are not valid")}
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
  time_ground <- seq(300, 710, 10)
  
  if(time_pt %in% time_ground){
    obs <- bind_rows(sate, inter, grd)
  }else{
    obs <- bind_rows(sate, inter)
  }


  ## dimensions
  n1 = nrow(inter)
  n2 = nrow(sate)
  
  ##############################################################
  # Notice that there might be small sample size of satellite data or no sate/ground base data
  ## Downsampling: random sampling
  if (n2 < 3000) print("Warning: the satellite data size at this time point is small")
  
  if (n2 > 3000) {
    index1 <- sample(1:n1, floor(n2*ratio1))    # ratio1 of satellite  data 
  } else if (n2 < 3000 & time_pt %in% time_ground) {
    n4 <- nrow(grd)
    index1 <- sample(1:n1, floor(n4*ratio1))     # ratio1 of ground base data
  } else if (n2 <3000 & time_pt %notin% time_ground) {
    index1 <- sample(1:n1, floor(n1*ratio5))    # 10% of the interpolated data itself
    n5 <- length(index1)
  }
  
  
  # inter
  dat1 <- inter[index1,]
  dat1$source <- "inter_random"
  # ground base data
  if (time_pt %in% time_ground){
    n4 = nrow(grd)
    
    # downsampling of ground data
    if (n2 > 3000){
      index3 <- sample(1:n4, floor(n4*ratio4))  # ratio4 of satellite data
    } else if (n2 < 3000){
      index3 <- sample(1:n4, floor(n4*1))   # keep all ground base data
    }
    
    grd1 <- grd[index3,]
    grd1$source <- "ground_random"
    
    # combine all data
    dat <- bind_rows(inter, sate,grd)
    
    ## combine the downsampled interpolation (empirical) with satellite daya at this time point
    all <- rbind(grd1, dat1, sate) %>%
      group_by(x, y, source) %>%
      summarize(
        flux = mean(flux)
      ) %>%
      ungroup()
  }else{
    print("No observed ground data for this time")
    
    dat <- bind_rows(inter, sate)  # change the empirical training
    
    all <- rbind(dat1, sate) %>%
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
  ########################################
  # combine test.location with observed location
  
  prob.bd.l <- cbind(predict.prob.l, test.location)
  prob.bd.m <- cbind(predict.prob.m, test.location)
  prob.bd.h <- cbind(predict.prob.h, test.location)
  
  weight.l <- prob.bd.l$predict.prob.l
  weight.m <- prob.bd.m$predict.prob.m
  weight.h <- prob.bd.h$predict.prob.h
  
  ## knn indicator
  indicator.l <- ifelse(weight.l>0, 1, 0)
  indicator.m <- ifelse(weight.m>0, 1, 0)
  indicator.h <- ifelse(weight.h>0, 1, 0)
  
  # knn and outer boundary indicator
    weight.b.l <- ifelse(outer==1, weight.l>0, 0)
    weight.b.m <- ifelse(outer==1, weight.m>0, 0)
    weight.b.h <- ifelse(outer==1, weight.h>0, 0)
    
  #########################################################################################
    ## create test locations for generated map
    # high resolution
    ne <- 50  # 100
    xe <- seq(-0.8, 0.8, len=ne)
    ye <- seq(-0.8, 0.8, len=ne)
    
    emp.location <- expand.grid(xe, ye)
    names(emp.location) <- c("x", "y")
    
    emp.location <- emp.location %>% filter(x^2+y^2<=0.7^2)
  # Extraplate
    obs <- location
    x.coord <- obs %>% pull(x)
    y.coord <- obs %>% pull(y)
    intro.coords <- cbind(x.coord, y.coord)
    
    X <- unique.matrix(intro.coords)
    library(alphahull)
    inregion <- ahull(X, alpha=0.05)  # adjustable para
    
    # extrapolation locations' indicators
    extrap <- ifelse(inahull(inregion, as.matrix(emp.location))==FALSE, 1, 0)
    
    loc.extrap <- emp.location[extrap==1,]
    location1 <- rbind(location, loc.extrap)
    
  # Empirical Inputs
    loc_ind <- nrow(location1) # empirical 
    testloc_ind <- nrow(test.location)
    locex_ind <- nrow(loc.extrap)
    
    # Calculate the ratio
    kpt <- approx(x = seq(from = 0, to = 1440, by = 180), y = c(1.3, 6, 6, 5, 4.7, 3.7, 2, 3, 3), xout = time_pt)$y
    x_e <- c(sate$x, grd$x)
    y_e <- c(sate$y, grd$y)
    flux_e <- c(sate$flux, grd$flux)
    valid1 <- flux_e > 0
    valid1[is.na(valid1)] <- FALSE
    index.e <- length(flux_e[valid1])
    loc.e <- data.frame(x=x_e[valid1], y=y_e[valid1])
    if (sum(valid1) >= 1) {
      emp_flux <- Emp_Data(1:index.e, time_pt,loc.e)
      ratio_e <- median(flux_e[valid1]) / median(emp_flux)
    }
    
    # Empirical data
    z_emp <- sapply(1:loc_ind, Emp_Data, location=location1, time=time_pt)
    empinput <- sapply(1:locex_ind, Emp_Data, location=loc.extrap, time=time_pt)
    empinput <- ratio_e*empinput
    z_new <- sapply(1:testloc_ind, Emp_Data,location = test.location, time=time_pt)
    
    emp.in <- tibble(flux = empinput)
    emp.in <- cbind(emp.in, loc.extrap)
    emp.in$source <- "emp: data"
    emp.z <- tibble(flux = z_emp)
    emp.z <- cbind(emp.z, location1)
    emp.z$source <- "emp: mean"
    emp.out <- tibble(flux = z_new)
    emp.out <- cbind(emp.out, test.location)
    emp.out$source <- "emp: pred"
    
    comb.e <- rbind(emp.in, emp.z, emp.out)
    
    ## generate a plot of fitted map
    pe <- comb.e %>%
      mutate(
        flux = ifelse(flux>30, 30, flux)
      ) %>%
      ggplot(aes(x = x, y = y)) +
      geom_point(aes(colour = flux), size=0.01)+
      scale_colour_distiller(palette = "Spectral")+
      facet_grid(.~source) +
      theme_bw() +
      ggtitle(paste("Empirical Data, Mean and predition, UT = ", time_pt, "mins"))
    
    ############################
    all1 <- rbind(all0, emp.in) %>%
      group_by(x, y) %>%
      summarize(
        flux = mean(flux)
      ) %>%
      ungroup()

  flux.n <- log(all1$flux)
  
  emp.z1 <- emp.z %>%
    group_by(x, y,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup()
  
  
  location2 <- all1 %>% select(x, y)
  
  z_emp1 <- emp.z1$flux
    
    fit.l <- LatticeKrig(location2, flux.n, nlevel=nlevels[1],  NC=NCs[1], Z=z_emp1)
    fit.m <- LatticeKrig(location2, flux.n, nlevel=nlevels[2],  NC=NCs[2], Z=z_emp1)
    fit.h <- LatticeKrig(location2, flux.n, Z=z_emp1)
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
  nt1 <- ceiling(sqrt(n2*ratio1/0.6))
  x1 <- seq(-0.8, 0.8, len=nt1)
  y1 <- seq(-0.8, 0.8, len=nt1)

  test.location.low <- expand.grid(x1,   y1)
  names(test.location.low) <- c("x", "y")

  test.location.low <- test.location.low %>% filter(x^2+y^2<=0.7^2)
  
  ind_loclow <- nrow(test.location.low)
  
  z_new.low <- sapply(1:ind_loclow, Emp_Data,location = test.location.low, time=time_pt) 

  sim2 <- LKrig.sim.conditional(fit.h, x.grid= test.location.low,Z.grid = matrix(z_new.low, ncol=1), M=100)
  # sim2 <- LKrig.sim.conditional(fit2, nx=50, ny=50, M=3)
  uncertainty3 <- tibble(SE = sim2$SE)
  uncertainty3 <- cbind(uncertainty3, test.location.low)
  uncertainty3$source <- "uncertainty.high"
    
  #z_new[is.na(z_new)] = 0
  z_new <- matrix(z_new, ncol=1)
  
  #############################
  # prediction mean
  pred.l <- predict(fit.l, x = test.location, Z=z_new)
  pred.m <- predict(fit.m, x = test.location, Z=z_new)
  pred.h <- predict(fit.h, x = test.location, Z=z_new)
  
  ## without knn and outer boundary
  dat.l <- tibble(flux=exp(pred.l[,1]))
  dat.l <- cbind(dat.l, test.location)
  dat.l$source <- "a.Low"
  dat.m <- tibble(flux=exp(pred.m[,1]))
  dat.m <- cbind(dat.m, test.location)
  dat.m$source <- "b.Mid"
  dat.h <- tibble(flux=exp(pred.h[,1])) 
  dat.h <- cbind(dat.h, test.location)
  dat.h$source <- "c.High"
  
  ## with knn
  dat.l.1 <- tibble(flux=exp(pred.l[,1])*((weight.l>0)))
  dat.l.1 <- cbind(dat.l.1, test.location)
  dat.l.1$source <- "a.Low"
  dat.m.1 <- tibble(flux=exp(pred.m[,1])*((weight.m>0)))
  dat.m.1 <- cbind(dat.m.1, test.location)
  dat.m.1$source <- "b.Mid"
  dat.h.1 <- tibble(flux=exp(pred.h[,1])*((weight.h>0)))
  dat.h.1 <- cbind(dat.h.1, test.location)
  dat.h.1$source <- "c.High"
  
  ## with knn and outer boundary
  dat.l.2 <- tibble(flux=exp(pred.l[,1])*((weight.b.l)))
  dat.l.2 <- cbind(dat.l.2, test.location)
  dat.l.2$source <- "a.Low"
  dat.m.2 <- tibble(flux=exp(pred.m[,1])*((weight.b.m)))
  dat.m.2 <- cbind(dat.m.2, test.location)
  dat.m.2$source <- "b.Mid"
  dat.h.2 <- tibble(flux=exp(pred.h[,1])*((weight.b.h))) 
  dat.h.2 <- cbind(dat.h.2, test.location)
  dat.h.2$source <- "c.High"
  
  #################
  all$source <- "samples"

  comb <- rbind(dat.l, dat.m, dat.h)
  comb.1 <- rbind(dat.l.1, dat.m.1, dat.h.1)
  comb.2 <- rbind(dat.l.2, dat.m.2, dat.h.2)
  
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
    ggtitle(paste("No knn and Outer Boundary, UT = ", time_pt, "mins"))
  
  p1 <- comb.1 %>%
    mutate(
      flux = ifelse(flux>30, 30, flux)
    ) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(colour = flux), size=0.01)+
    scale_colour_distiller(palette = "Spectral")+
    facet_grid(.~source) +
    theme_bw() +
    ggtitle(paste("With knn, UT = ", time_pt, "mins"))
  
  p2 <- comb.2 %>%
    mutate(
      flux = ifelse(flux>30, 30, flux)
    ) %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(colour = flux), size=0.01)+
    scale_colour_distiller(palette = "Spectral")+
    facet_grid(.~source) +
    theme_bw() +
    ggtitle(paste("With knn and Outer Boundary, UT = ", time_pt, "mins"))
  
  pu <- uncertainty3 %>%
    ggplot(aes(x = x, y = y)) +
    geom_point(aes(colour = SE), size=0.01)+
    scale_colour_distiller(palette = "Spectral")+
    facet_grid(.~source) +
    theme_bw() +
    ggtitle(paste("Uncertainty for high, UT = ", time_pt, "mins"))
  
  return(list(p=p,p1 = p1, p2 = p2, pe = pe, pu= pu, res1 = comb, res2 = comb.1, res3 = comb.2, res4 = comb.e, uncertain = uncertainty3))
}

test7 <- LK.Empmean.10m (data1 = data.flux.inter_n, data2 = data.flux.sate_n, data4 = data.flux.grd_new,
                        nmlat=363, nmlt=363, time_pt=420, ratio1=1/6,  ratio4=1,ratio5=1/10,
                        NCs=c(5,20,25), nlevels= c(1,1,3), kn = c(60,20,10))
test11 <- LK.Empmean.10m (data1 = data.flux.inter_n, data2 = data.flux.sate_n, data4 = data.flux.grd_new,
                         nmlat=363, nmlt=363, time_pt=660, ratio1=1/6, ratio4=1,ratio5=1/10,
                         NCs=c(5,20,25), nlevels= c(1,1,3), kn = c(60,20,10))
save(test11, file="empmean_11.RData")
