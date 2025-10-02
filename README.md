### GLOWORM-FL model

## flsim

**DESCRIPTION**
This function is used to run a single simulation of the GLOWORM-FL model (Rose et al. 2015) in R (https://www.r-project.org/) for either Haemonchus contortus or Teladorsagia
circumcincta infecting sheep/goats/deer or Ostertagia ostertagi infecting cattle.

# USAGE
flsim(init, start, end, input, temp, precip, species, lat)

# ARGUMENTS
init = vector of initial numbers of eggs (E), pre-infective larvae (L), infective L3 in faeces (L3), and infective L3 on herbage (Herbage) in the format
c(E, L, L3, Herbage)
start =  start date of the simulation in the format "YYYY-MM-DD" e.g. "2015-01-01"
end = end date of the simulation in the format "YYYY-MM-DD" e.g. "2015-02-01"
input = vector of daily numbers of eggs deposited on pasture. Can correspond to individual/group fecs or total eggs deposited per day
temp = vector of daily mean air temperature values in degrees centigrade. Vector length must correspond to the length of the simulation (start date
to end date inclusive)
precip = vector of daily total precipitation values in mm. Vector length must correspond to the length of the simulation (start date to end date inclusive)
species = character string specifying the strongyle species to simulate. Choices are "haemonchus", "teladorsagia" or "ostertagia"
lat = latitude of the study site in decimal degrees

# DETAILS
Units need not be specified but must be consistent e.g. if 'init' values are per hectare then 'input' egg counts must also represent eggs deposited per
hectare per day (e.g. based on group faecal egg counts and stocking density), and the output matrix should be interpreted as individuals per hectare. Egg
counts conducted less frequently than daily (as is usually the case) can be interpolated linearly using the parainterp() function to produce the 'input'
vector. Temperature and precipiation time series must not contain NAs. Find NAs using which(is.na(temp)) and replace any NAs with a numeric value e.g. moving average, dataset mean, 0 etc. using
precip[which(is.na(precip))]=0.

# VALUE
A matrix with one row per day and columns showing the time step (time), numbers of eggs (E), pre-infective larvae (L), infective L3 in faeces (L3), total
infective L3 on pasture (Pasture), L3 in soil (Soil) and L3 on herbage (Herbage)

## parainterp()

# DESCRIPTION
This function is used to interpolate parasitological data such as faecal egg counts between sampling dates, to provide a daily time series for use with flsim().

# USAGE
parainterp(para, dates, method = 'linear')

# ARGUMENTS
para = vector of parasitological data
dates = vector of dates corresponding to the parasitological data in the format "YYYY-MM-DD" e.g. "2015-02-01"
method = method of interpolation ('linear' or 'constant'). Defaults to 'linear'.

# DETAILS
Input must include data and date corresponding to the start and end date input into flsim().

# VALUE
Vector of daily estimates for the parasitological data


## EXAMPLE SIMULATION

temp = c(19, 20, 20.5, 20.6, 18, 17, 19, 20, 21, 21, 21.5, 20, 21, 20)
precip = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 10, 5, 0.5)
#generate some random egg counts and dates of egg counts
fecs = sort(runif(n = 4, min = 50, max = 1000), decreasing = FALSE)
fec.dates = c("2000-06-01", "2000-06-05", "2000-06-10", "2000-06-14")
#inteprolate the egg counts using FECinterp
fec = parainterp(para = fecs, dates = fec.dates)
test = flsim(init = c(100, 0, 0, 0), start = "2000-06-01", end = "2000-06-14", input = fec, temp = temp, precip = precip, species = "head(test)
plot(test[,"E"], type = "l", xlab = "days", ylab = "number of eggs")
plot(test[,"Pasture"], type = "l", xlab = "days", ylab = "number of L3 on pasture")


## REFERENCE

Rose, H., Wang, T., van Dijk, J. & Morgan, E. R. 2015. GLOWORM-FL: A simulation model of the effects of climate and climate change on the free-living stages of gastro-intestinal nematode parasites of ruminants. Ecological Modelling. 297, p.232-245
https://doi.org/10.1016/j.ecolmodel.2014.11.033
