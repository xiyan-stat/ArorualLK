################################################################
# data fusion

# Note: Empirical is Area, Satellite data is kind of point

# satellite data
data1 <- dat.flux.sate_n
data2 <- dat.flux.emp_n

# empirical Area: #mlt: 0-24 by 0.25 hrs   # mlat: 50-90 by 0.5 deg

# Within time regression
berrocal <- function(data1=data1, data2=data2, time_pt){
  ##############################################################
  # regression of sate on emp data
  # Input: data1 = sate point data at time t
  #        data2 = emp area data at time t
  #        time_pt = at time t
  # Output: beta0 = intercept 
  #         beta1 = ratio/slope of the model
  ###############################################################
  
  # satellite
  sate <- data1 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, mlat, mlt, source)
  
  # empirical
  emp <- data2 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux,mlat, mlt, source)
  
  # set up area vectors
  latitude <- seq(50, 90, 0.5)   # by 0.5 deg 
  localtime <- seq(0, 24, 0.25)   # by 0.25 hrs
  
  x <- c()
  y <- c()
  
  # generate x, y of each area
  for (i in 1: (length(latitude)-1)){
    for (j in 1: (length(localtime)-1)){
      latmin <- latitude[i]
      latmax <- latitude[(i+1)]
      ltmin <- localtime[j]
      ltmax <- localtime[(j+1)]
      
      sate.m <- mean(unlist(sate[sate$mlat>=latmin & sate$mlat<latmax & sate$mlt >= ltmin & sate$mlt <ltmax,1]))
      emp.m <- mean(unlist(emp[emp$mlat>=latmin & emp$mlat<latmax & emp$mlt >= ltmin & emp$mlt <ltmax,1]))
      
      if (!is.na(sate.m)){
        y <- c(y, sate.m)
        x <- c(x, emp.m)
      }
    }
  }

res <- lm(y~x)
beta0 <- as.numeric(res$coef[1])
beta1 <- as.numeric(res$coef[2])

return(list(beta0=beta0, beta1=beta1, x=x, y=y, fit=summary(res)))
}

# times
ts <- seq(0, 1380, 60)
Res.mean <- sapply(ts, berrocal, data1=data1, data2=data2)

beta0s <- Res.mean[1,]
beta1s <- paras[2,]

plot(ts, beta0s)
plot(ts, beta1s)

plot(ts, beta1s,type="l",col="red")
lines(ts, beta0s,col="blue")
  
b1 = berrocal(data1=data1, data2=data2, time_pt=540)

###########################################################################
# Within time regression
berrocal <- function(data1=data1, data2=data2, time_pt){
  ##############################################################
  # regression of sate on emp data
  # Input: data1 = sate point data at time t
  #        data2 = emp area data at time t
  #        time_pt = at time t
  # Output: beta0 = intercept 
  #         beta1 = ratio/slope of the model
  ###############################################################
  
  # satellite
  sate <- data1 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, mlat, mlt, source)
  
  # empirical
  emp <- data2 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux,mlat, mlt, source)
  
  # set up area vectors
  latitude <- seq(50, 90, 0.5)   # by 0.5 deg 
  localtime <- seq(0, 24, 0.25)   # by 0.25 hrs
  
  x <- c()
  y <- c()
  
  # generate x, y of each area
  for (i in 1: (length(latitude)-1)){
    for (j in 1: (length(localtime)-1)){
      latmin <- latitude[i]
      latmax <- latitude[(i+1)]
      ltmin <- localtime[j]
      ltmax <- localtime[(j+1)]
      
      sate.m <- mean(unlist(sate[sate$mlat>=latmin & sate$mlat<latmax & sate$mlt >= ltmin & sate$mlt <ltmax,1]))
      emp.m <- mean(unlist(emp[emp$mlat>=latmin & emp$mlat<latmax & emp$mlt >= ltmin & emp$mlt <ltmax,1]))
      
      if (!is.na(sate.m)){
        y <- c(y, sate.m)
        x <- c(x, emp.m)
      }
    }
  }
  
  res <- summary(lm(y~x))
  beta0 <- as.numeric(res$coef[1])
  beta1 <- as.numeric(res$coef[2])
  
  
  return(list(beta0=beta0, beta1=beta1, x=x, y=y,fit=res))
}

# times
ts <- seq(0, 1380, 10)
Res.mean <- sapply(ts, berrocal, data1=data1, data2=data2)

save(Res.mean, file="Res.mean.RData")


berrocal.point <- function(data1=data1, data2=data2, time_pt){
  ##############################################################
  # regression of sate on emp data
  # Input: data1 = sate point data at time t
  #        data2 = emp area data at time t
  #        time_pt = at time t
  # Output: beta0 = intercept 
  #         beta1 = ratio/slope of the model
  ###############################################################
  
  # satellite
  sate <- data1 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, mlat, mlt, source)
  
  # empirical
  emp <- data2 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux,mlat, mlt, source)
  
  # set up area vectors
  latitude <- seq(50, 90, 0.5)   # by 0.5 deg 
  localtime <- seq(0, 24, 0.25)   # by 0.25 hrs
  
  x <- c()
  y <- c()
  
  # generate x, y of each area
  for (i in 1: (length(latitude)-1)){
    for (j in 1: (length(localtime)-1)){
      latmin <- latitude[i]
      latmax <- latitude[(i+1)]
      ltmin <- localtime[j]
      ltmax <- localtime[(j+1)]
      
      sate.p <- unlist(sate[sate$mlat>=latmin & sate$mlat<latmax & sate$mlt >= ltmin & sate$mlt <ltmax,1])
      emp.g <- unlist(emp[emp$mlat>=latmin & emp$mlat<latmax & emp$mlt >= ltmin & emp$mlt <ltmax,1])
      
      if (length(sate.p)!=0){
        y <- c(y, sate.p)
        x <- c(x, rep(emp.g,length(sate.p)))
      }
    }
  }
  
  res <- lm(y~x)
  beta0 <- as.numeric(res$coef[1])
  beta1 <- as.numeric(res$coef[2])
  
  return(list(beta0=beta0, beta1=beta1, x=x, y=y,fit=summary(res)))
}
#
Res.point <- sapply(ts, berrocal.point, data1=data1, data2=data2)


#b1.1 = berrocal.point(data1=data1, data2=data2, time_pt=540)$beta1

save(Res.mean, Res.point, file="sate_emp.RData")

for (i in 1:24){
  cat(paste("-----------------Regression summary at time:", i-1,"-------------\n"))
  print(Res.mean[[(5+5*(i-1))]])
}

#######################################################################################
berrocal.1 <- function(data1=data1, data2=data2, time_pt){
  ##############################################################
  # regression of sate on emp data
  # Input: data1 = sate point data at time t
  #        data2 = emp area data at time t
  #        time_pt = at time t
  # Output: beta0 = intercept 
  #         beta1 = ratio/slope of the model
  ###############################################################
  
  # satellite
  sate <- data1 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, mlat, mlt, source)
  
  # empirical
  emp <- data2 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux,mlat, mlt, source)
  
  # set up area vectors
  latitude <- seq(50, 90, 0.5)   # by 0.5 deg 
  localtime <- seq(0, 24, 0.25)   # by 0.25 hrs
  
  x <- c()
  y <- c()
  
  # generate x, y of each area
  for (i in 1: (length(latitude)-1)){
    for (j in 1: (length(localtime)-1)){
      latmin <- latitude[i]
      latmax <- latitude[(i+1)]
      ltmin <- localtime[j]
      ltmax <- localtime[(j+1)]
      
      sate.m <- mean(unlist(sate[sate$mlat>=latmin & sate$mlat<latmax & sate$mlt >= ltmin & sate$mlt <ltmax,1]))
      emp.m <- mean(unlist(emp[emp$mlat>=latmin & emp$mlat<latmax & emp$mlt >= ltmin & emp$mlt <ltmax,1]))
      
      if (!is.na(sate.m)){
        y <- c(y, sate.m)
        x <- c(x, emp.m)
      }
    }
  }
  
  index0 <- which(y!=0)
  y <- y[index0]
  x <- x[index0]
  res <- summary(lm(y~x))
  beta0 <- as.numeric(res$coef[1])
  beta1 <- as.numeric(res$coef[2])
  
  return(list(beta0=beta0, beta1=beta1, x=x, y=y, fit=res))
}
berrocal.point1 <- function(data1=data1, data2=data2, time_pt){
  ##############################################################
  # regression of sate on emp data
  # Input: data1 = sate point data at time t
  #        data2 = emp area data at time t
  #        time_pt = at time t
  # Output: beta0 = intercept 
  #         beta1 = ratio/slope of the model
  ###############################################################
  
  # satellite
  sate <- data1 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, mlat, mlt, source)
  
  # empirical
  emp <- data2 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux,mlat, mlt, source)
  
  # set up area vectors
  latitude <- seq(50, 90, 0.5)   # by 0.5 deg 
  localtime <- seq(0, 24, 0.25)   # by 0.25 hrs
  
  x <- c()
  y <- c()
  
  # generate x, y of each area
  for (i in 1: (length(latitude)-1)){
    for (j in 1: (length(localtime)-1)){
      latmin <- latitude[i]
      latmax <- latitude[(i+1)]
      ltmin <- localtime[j]
      ltmax <- localtime[(j+1)]
      
      sate.p <- unlist(sate[sate$mlat>=latmin & sate$mlat<latmax & sate$mlt >= ltmin & sate$mlt <ltmax,1])
      emp.g <- unlist(emp[emp$mlat>=latmin & emp$mlat<latmax & emp$mlt >= ltmin & emp$mlt <ltmax,1])
      
      if (length(sate.p)!=0){
        y <- c(y, sate.p)
        x <- c(x, rep(emp.g,length(sate.p)))
      }
    }
  }
  
  index0 <- which(y!=0)
  y <- y[index0]
  x <- x[index0]
  res <- lm(y~x)
  beta0 <- as.numeric(res$coef[1])
  beta1 <- as.numeric(res$coef[2])
  
  return(list(beta0=beta0, beta1=beta1, x=x, y=y,fit=summary(res)))
}
ts <- seq(0, 1380, 60)
Res.mean1 <- sapply(ts, berrocal.1, data1=data1, data2=data2)
Res.point1 <- sapply(ts, berrocal.point1, data1=data1, data2=data2)
save(Res.mean1, Res.point1, file="sate_emp1.RData")

for (i in 1:24){
  cat(paste("-----------------Regression summary at time:", i-1,"-------------\n"))
  print(Res.mean1[[(5+5*(i-1))]])
}


##########################################################################################
# Spatial Analysis
# Area: lat: 10 degree, lt: 6 hours
berrocal.spatial <- function(data1=data1, data2=data2, time_pt){
  ##############################################################
  # regression of sate on emp data
  # Input: data1 = sate point data at time t
  #        data2 = emp area data at time t
  #        time_pt = at time t
  # Output: beta0 = intercept 
  #         beta1 = ratio/slope of the model
  ###############################################################
  
  # satellite
  sate <- data1 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux, mlat, mlt, source)
  
  # empirical
  emp <- data2 %>%
    filter(time %in% time_pt) %>%
    group_by(mlat, mlt,source) %>%
    summarize(
      flux = mean(flux)
    ) %>%
    ungroup() %>%
    select(flux,mlat, mlt, source)
  
  
  lats <- seq(50,90,10)
  lts <- seq(0,24, 6)
  
  beta0 <- c()
  beta1 <- c()
  m <- 0
  #res <- list()
  xs <- list()
  ys <- list()
  
  for (k in 1:(length(lats)-1)){
    for (t in 1:(length(lts)-1)){
      m <- m + 1
      lat.min <- lats[k]
      lat.max <- lats[(k+1)]
      lt.min <- lts[t]
      lt.max <- lts[(t+1)]
      
      #cat(k,t,"\n")
      #cat(lat.min, lat.max, lt.min, lt.max,"\n")
      x <- c()
      y <- c()
      
      # set up area vectors
      latitude <- seq(lat.min, lat.max, 0.5)   # by 0.5 deg 
      localtime <- seq(lt.min, lt.max, 0.25)   # by 0.25 hrs
      
      # generate x, y of each area
      for (i in 1: (length(latitude)-1)){
        for (j in 1: (length(localtime)-1)){
          latmin <- latitude[i]
          latmax <- latitude[(i+1)]
          ltmin <- localtime[j]
          ltmax <- localtime[(j+1)]
          
          sate.m <- mean(unlist(sate[sate$mlat>=latmin & sate$mlat<latmax & sate$mlt >= ltmin & sate$mlt <ltmax,1]))
          emp.m <- mean(unlist(emp[emp$mlat>=latmin & emp$mlat<latmax & emp$mlt >= ltmin & emp$mlt <ltmax,1]))
          
          if (!is.na(sate.m)){
            y <- c(y, sate.m)
            x <- c(x, emp.m)
          }
        }
      }
      
      
      xs[[m]] = x
      ys[[m]] = y
      
      # index0 <- which(y!=0)
      # y <- y[index0]
      # x <- x[index0]
      # if(length(y!=0)&length(x!=0)){
      #   #xs[[m]] <- x
      #   #ys[[m]] <- y
      #   #res[[m]] <- summary(lm(y~x))
      #   beta0 <- c(beta0, as.numeric(res$coef[1]))
      #   beta1 <- c(beta1, as.numeric(res$coef[2]))
      # }else{
      #   #xs[[m]] <- NA
      #   #ys[[m]] <- NA
      #   #res[[m]] <- NA
      #   beta0 <- c(beta0, NA)
      #   beta1 <- c(beta1, NA)
      # }

    }
  }
  #return(list(beta0=beta0, beta1=beta1, x = xs, y = ys, fit=res))
  #return(list(beta0=beta0, beta1=beta1))
  return(list(x = xs, y = ys))
}

spatialData <- sapply(ts, berrocal.spatial, data1=data1, data2=data2)

b.s = berrocal.spatial(data1=data1, data2=data2, time_pt=540)

save(spatialData, file="spatialData.RData")
###################################################
#
la <- c("50-60", "60-70","70-80","80-90")
lt <- c("0-6","6-12","12-18","18-24")

loc <- expand.grid(lt,la)

b0 <- b.s$beta0
b1 <- b.s$beta1
fitres <- b.s$fit 

SpatialModel.l <- tibble(
  latitude = as.vector(loc[,2]),
  localtime = as.vector(loc[,1]),
  additive = as.vector(b0),
  multiplicative = as.vector(b1)
)

