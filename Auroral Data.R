#################################################################
#            R code for LK aurora project -- Data               #
#################################################################

# 1. Linear interpolation + satellite data using LK validate: compare with linear interpolation

rm(list=ls())

library(tidyverse)
require(class)

###################################################
# Source data: satellite data 20140219-20140221
library(ncdf4)

nc <- nc_open(filename = "ssusi_combine_20140219_20140221.nc")  
nc

# Extract the variables: Satellite auroral
# length of time, lattitude and longtitude
norbit <- nc$dim$orbit$len   # 121
nmlat <- nc$dim$nmlat$len    # 363
nmlt <- nc$dim$nmlt$len      # 363

# variables: lattitude & local time (logitude)
mlat <- ncvar_get(nc = nc, "mlat")
mlt <- ncvar_get(nc = nc, "mlt")

# Southern Hemisphere
ut_s <- (ncvar_get(nc = nc, "ut_s") - 24) * 60
flux_s <- ncvar_get(nc = nc, "flux_s")
energy_s <- ncvar_get(nc = nc, "energy_s")

# Northern Hemisphere
ut_n <- (ncvar_get(nc = nc, "ut_n") - 24) * 60
flux_n <- ncvar_get(nc = nc, "flux_n")
energy_n <- ncvar_get(nc = nc, "energy_n")


# close A NetCDF file
nc_close(nc = nc)

################################################################
# Step 1: generate interpolated and satellite data for northern and southern hemisphere *8:45-9:15, t = 540

nhemi <- 2
t <- seq(from = 0, to = 1440, by = 60)  # every 30 mins: try 10 mins 11:00-11:10 660
nt <- length(t)

##########################################
# 1. satellite data
# generate variabels in satellite data
#south
flux_sate_s <- c()
energy_sate_s <- c()
x_sate_s <- c()
y_sate_s <- c()
time_sate_s <- c()
#north
flux_sate_n <- c()
energy_sate_n <- c()
x_sate_n <- c()
y_sate_n <- c()
time_sate_n <- c()

#covert function
xy <- function(mlat, mlt){
  r <- pi/2 - mlat*pi/180
  tt <- mlt*pi/12
  x_res <- r * cos(tt)
  y_res <- r * sin(tt)
  return(list(x_res, y_res))
}

#########################
# satellite data             # check previous time points
for (i in 1:nt){
  v_s = ut_s > t[i] - 15 & ut_s < t[i] + 15
  v_s[is.na(v_s)] = FALSE
  v_n = ut_n > t[i] - 15 & ut_n < t[i] + 15
  v_n[is.na(v_n)] = FALSE
  
  # generate satellite data
  flux_s0 <- flux_s[v_s]
  flux_n0 <- flux_n[v_n]
  energy_s0 <- energy_s[v_s]
  energy_n0 <- energy_n[v_n]
  # length of flux and energy for each time
  m_s = length(flux_s0)
  m_n = length(flux_n0)
  
  mlat_s0 = replicate(norbit, mlat)[v_s]
  mlt_s0 = replicate(norbit, mlt)[v_s]
  mlat_n0 = replicate(norbit, mlat)[v_n]
  mlt_n0 = replicate(norbit, mlt)[v_n]
  
  #convert coordinates function
  x_s0 <- xy(mlat_s0, mlt_s0)[[1]]
  y_s0 <- xy(mlat_s0, mlt_s0)[[2]]
  x_n0 <- xy(mlat_n0, mlt_n0)[[1]]
  y_n0 <- xy(mlat_n0, mlt_n0)[[2]]
  
  #generate variables
  flux_sate_s <- c(flux_sate_s, flux_s0)
  energy_sate_s <- c(energy_sate_s, energy_s0)
  x_sate_s <- c(x_sate_s, x_s0)
  y_sate_s <- c(y_sate_s, y_s0)
  time_sate_s <- c(time_sate_s, rep(t[i], m_s))
  
  flux_sate_n <- c(flux_sate_n, flux_n0)
  energy_sate_n <- c(energy_sate_n, energy_n0)
  x_sate_n <- c(x_sate_n, x_n0)
  y_sate_n <- c(y_sate_n, y_n0)
  time_sate_n <- c(time_sate_n, rep(t[i], m_n))
}


#########################################
# For example,
# northern sate flux data
data.flux.sate_n <- tibble(
  flux = flux_sate_n,
  x = x_sate_n,
  y = y_sate_n,
  time = time_sate_n,
  source = "sate"
) %>%
  filter(!is.na(flux)) 


##########################################
# 2. interpolation
# generate series of flux and enegry for different location and time
flux_inter <- array(dim = c(nmlt, nmlat, nt, nhemi))
energy_inter <- array(dim = c(nmlt, nmlat, nt, nhemi))

for (imlat in 1: nmlat) {
  for (imlt in 1: nmlt) {
    x <- ut_s[imlt, imlat, ]  # time
    y1 <- flux_s[imlt, imlat, ]   # flux
    y2 <- energy_s[imlt, imlat, ]   # energy
    valid <- !(is.na(x) | is.na(y1) | is.na(y2))
    x <- x[valid]
    y1 <- y1[valid]
    y2 <- y2[valid]
    if (length(x) >= 2) {               # threshold parameter: alpha=2; if alpha increases, more accurate. (eg, alpha=5)
      ix <- sort(x, index.return = TRUE)$ix   # x: x < threshold: 1 hour/ 2 hours    
      x <- x[ix]   # 
      if (min(diff(x, lag=1)) < 120){
        # inter_south data
        flux_inter[imlt, imlat, , 1] <- approx(x = x, y = y1[ix], xout = t)$y  # south      # simulation to see how "approx" perform, if no extrapolation, try other
        energy_inter[imlt, imlat, , 1] <- approx(x = x, y = y2[ix], xout = t)$y 
      }
    }
    
    x <- ut_n[imlt, imlat, ]
    y1 <- flux_n[imlt, imlat, ]
    y2 <- energy_n[imlt, imlat, ]
    valid <- !(is.na(x) | is.na(y1) | is.na(y2))
    x <- x[valid]
    y1 <- y1[valid]
    y2 <- y2[valid]
    if (length(x) >= 2) {
      ix <- sort(x, index.return = TRUE)$ix
      x <- x[ix]
      if (min(diff(x, lag=1)) < 120){
      flux_inter[imlt, imlat, , 2] <- approx(x = x, y = y1[ix], xout = t)$y    # north: linear interpolation
      energy_inter[imlt, imlat, , 2] <- approx(x = x, y = y2[ix], xout = t)$y
      }
    }
  }
}

# north:
# Cartesian coordinate in polar view
x_inter <- matrix(NA, nrow=nmlt, ncol=nmlat)
y_inter <- matrix(NA, nrow=nmlt, ncol=nmlat)

for (i in 1:nmlt){
  for (j in 1:nmlat){
    r <- pi/2 - mlat[i,j]*pi/180
    tt <- mlt[i,j]*pi/12
    x_inter[i,j] <- r * cos(tt)
    y_inter[i,j] <- r * sin(tt)
  }
}

#########################################
# For example,
# northern sate flux data
data.flux.inter_n <- tibble(
  flux = as.vector(flux_inter[,,,2]),
  x = rep(as.vector(x_inter), nt),
  y = rep(as.vector(y_inter), nt),
  time = rep(t, each=nmlat*nmlt),
  source = "inter"
) %>%
  filter(!is.na(flux)) 


############################################
# combine the data flux
flux <- bind_rows(data.flux.inter_n, data.flux.sate_n) # original data

flux_30 <- bind_rows(data.flux.inter_n, data.flux.sate_n) %>%
  mutate(
    flux = ifelse(flux>30, 30, flux)            # visulization
  )

# generate the maps
flux_30 %>%
  filter(time %in% 660) %>%  # 690 11:15 - 11:45/ 660 11:00-11:10
  ggplot(aes(x = x, y = y)) +
  geom_point(aes(colour = flux), size=0.1)+
  facet_grid(.~source)+
  scale_colour_distiller(palette = "Spectral")+
  theme_bw()+
  ggtitle("11:00 - 12:00")

# another way to plot
library(fields)

t_temp = 540 # 8:45-9:15   (9:00)
t_temp = 690    # 11:15 - 11:45 (11:30)

data <- data.flux.inter_n %>%  
  #data.flux.sate_n %>%
  mutate(
    flux = ifelse(flux>30, 30, flux)            # visulization
  )

quilt.plot(data$x[data$time %in% t_temp], data$y[data$time %in% t_temp], data$flux[data$time %in% t_temp],
           bty = "n", zlim = range(data$flux[data$time %in% t_temp], na.rm = T),
           xaxt = "n", yaxt = "n", xlim = (pi / 4) * c(-1.1, 1.1),
           ylim = (pi / 4) * c(-1.1, 1.1), xlab = "", ylab = "",
           col = tim.colors(256), nx = 120, ny = 120)
r = (pi / 4) * (4.25 / 4.5); theta <- seq(0, 2 * pi, 0.01)
lines(r * cos(theta), r * sin(theta), lwd = 1.25)
text(0, 1.01 * (pi / 4), "12")
text(1.01 * (pi / 4), 0, "06")
text(-1.01 * (pi / 4), 0, "18")
text(0, -1.01 * (pi / 4), "00")

###############################################
# empirical data



# save the interpolation, empirical, and auroral data

save(flux,
     flux_30,
     data.flux.inter_n, 
     data.flux.sate_n, 
     data.flux.emp_n,
     file="flux_north.RData")  
