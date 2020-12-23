################################################################
# data fusion

# Note: Empirical is Area, Satellite data is kind of point

# satellite data
data1 <- dat.flux.sate_n
data2 <- dat.flux.emp_n

# empirical Area: #mlt: 0-24 by 0.25 hrs   # mlat: 50-90 by 0.5 deg

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


