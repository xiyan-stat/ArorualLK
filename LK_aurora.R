#################################################################
#            R code for LK aurora project -- LK                 #
#################################################################

rm(list=ls())
load("flux_north.RData")

# Method 1: Using LK of Interpolation

m1 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/10,  NC=10, nlevel=1, type="LK", prob=TRUE)

# Method 2: Using systematic sampling

m2 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/6, NC=10, nlevel=1, type="system", prob=TRUE)


# Method 3: Using random sampling

m3 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/6, NC=10, nlevel=1, type="random", prob=TRUE)


##########################################################
# try different ratios to compare the estimated maps

r1_6 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/6,type="random", prob=TRUE)
r1_7 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/7,type="random", prob=TRUE)
r1_8 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/8, type="random", prob=TRUE)
r1_9 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/9, type="random", prob=TRUE)
r1_10 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 660, ratio= 1/10,type="random", prob=TRUE)

r1_6_1 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, NC=10, nlevel=1, ratio= 1/6,type="random", prob=TRUE)
r1_6_2 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, NC=10, nlevel=2, ratio= 1/6,type="random", prob=TRUE)
r1_6_3 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, NC=10, nlevel=3, ratio= 1/6,type="random", prob=TRUE)

r1_6_4 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, NC=5, nlevel=1, ratio= 1/6,type="random", prob=TRUE)
r1_6_5 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, NC=20, nlevel=1, ratio= 1/6,type="random", prob=TRUE)

r1_6_6 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, NC=20, nlevel=1, ratio= 1/10,type="random", prob=TRUE)

# lk
l1_6 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/6,type="LK", prob=TRUE)
l1_10 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/10,type="LK", prob=TRUE)



# system
s1_6 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690,ratio= 1/6,type="system", prob=TRUE)
s1_10 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 690, ratio= 1/10,type="system", prob=TRUE)

################################################################

r1_6.10 <- LK.auroral(data1 = data.flux.inter_n, data2 = data.flux.sate_n, time_pt = 660, ratio= 1/6,type="random", prob=TRUE)


###############################################################################################
# Three dataset


r2.10.6 <- LK.auroral.2(data1 = data.flux.inter_n, data2 = data.flux.sate_n, data3=data.flux.emp_n, time_pt = 660, ratio1= 1/10, ratio2=1/6, type="random", prob=TRUE)






