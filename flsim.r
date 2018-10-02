# Author: Hannah Vineer
# Copyright: University of Bristol
# Contact: hannah.vineer@liverpool.ac.uk
# Licence: GNU GPL-3 (https://www.r-project.org/Licenses/AGPL-3) (i.e. users must make their code and output open source)


# GLOWORM-FL model function

flsim = function (init, start, end, input, temp, precip, species, lat) 
{
  date.range = seq(as.Date(start), as.Date(end), "days")
  if (length(date.range) != length(temp)) 
    stop("Error: length of climatic data and dates do not match. Ensure climatic data correspond to the start and end dates (one entry per day)")
  global.t = seq(1, length(date.range))
  photoperiod = daylength(lat = lat, doy = date.range)
  if (species == "haemonchus") {
    Pt = (4.95 * exp(0.062 * temp))/100
    PET = 0.55 * ((photoperiod/12)^2) * Pt * 25.4
    PET4 = NA
    PETm4 = NA
    Precip4 = NA
    Precipm4 = NA
    for (z in 1:length(PET)) {
      min = z
      max = ifelse((z + 4) > length(PET), length(PET), 
                   z + 4)
      PET4[z] = sum(x = PET[min:max])
      Precip4[z] = sum(x = precip[min:max])
      max = z
      min = ifelse((z - 4) < 1, 1, z - 4)
      PETm4[z] = sum(x = PET[min:max])
      Precipm4[z] = sum(x = precip[min:max])
      P.E.1 = Precip4/PET4
      P.E.2 = Precipm4/PETm4
    }
    dev.1 = pmax(0, -0.09746 + 0.01063 * temp)
    dev.1 = pmin(1, dev.1)
    mu.1 = pmin(1, exp(-1.47135 - 0.11444 * temp + 0.00327 * 
                         (temp^2)))
    mu.2 = pmin(1, exp(-1.823 - 0.1418 * temp + 0.00405 * 
                         (temp^2)))
    mu.3 = pmin(1, exp(-2.6308 - 0.14407 * temp + 0.00463 * 
                         (temp^2)))
    mu.4 = pmin(1, exp(-3.68423 - 0.2535 * temp + 0.0074 * 
                         (temp^2)))
    mu.5 = mu.3
    h.mig = ifelse(precip >= 2, 0.25, ifelse(P.E.2 >= 1, 
                                             0.051, 0))
    v.mig = (pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * 
                           (temp^2))))
    dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
    mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
    mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
    mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
    mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
    mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
    hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
    vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
    init[5] = (init[4] * (1/v.mig[1])) - init[4]
    egg.correction = ifelse(P.E.1 < 1, 0.1, 1)
    init[1] = init[1] * egg.correction[1]
    print("Haemonchus contortus parameters loaded successfully")
  }
  if (species == "teladorsagia") {
    Pt = (4.95 * exp(0.062 * temp))/100
    PET = 0.55 * ((photoperiod/12)^2) * Pt * 25.4
    PET7 = NA
    PETm7 = NA
    Precip7 = NA
    Precipm7 = NA
    for (z in 1:length(PET)) {
      min = z
      max = ifelse((z + 7) > length(PET), length(PET), 
                   z + 7)
      PET7[z] = sum(x = PET[min:max])
      Precip7[z] = sum(x = precip[min:max])
      max = z
      min = ifelse((z - 7) < 1, 1, z - 7)
      PETm7[z] = sum(x = PET[min:max])
      Precipm7[z] = sum(x = precip[min:max])
      P.E.1 = Precip7/PET7
      P.E.2 = Precipm7/PETm7
    }
    dev.1 = pmax(0, -0.02085 + 0.00467 * temp)
    dev.1 = pmin(1, dev.1)
    mu.1 = pmin(1, exp(-1.62026 - 0.17771 * temp + 0.00629 * 
                         (temp^2)))
    mu.2 = mu.1
    mu.4 = pmin(1, exp(-4.58817 - 0.13996 * temp + 0.00461 * 
                         (temp^2)))
    mu.3 = pmin(1, 10 * mu.4)
    mu.5 = mu.3
    h.mig = ifelse(precip >= 2, 0.21, ifelse(P.E.2 >= 1, 
                                             0.025, 0))
    v.mig = (pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * 
                           (temp^2))))
    dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
    mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
    mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
    mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
    mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
    mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
    hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
    vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
    init[5] = (init[4] * (1/v.mig[1])) - init[4]
    egg.correction = ifelse(P.E.1 < 1, 0.1, 1)
    init[1] = init[1] * egg.correction[1]
    print("Teladorsagia circumcincta parameters loaded successfully")
  }
  if (species == "ostertagia") {
    dev.1 = pmax(0, -0.07258 + 0.00976 * temp)
    dev.1 = pmin(1, dev.1)
    mu.1 = pmin(1, exp(-4.38278 - 0.1064 * temp + 0.0054 * 
                         (temp^2)))
    mu.2 = pmin(1, exp(-4.38278 - 0.1064 * temp + 0.0054 * 
                         (temp^2)))
    mu.4 = pmin(1, exp(-6.388 - 0.2681 * temp + 0.01633 * 
                         (temp^2) - 0.00016 * (temp^3)))
    mu.3 = pmin(1, mu.4 * 10)
    mu.5 = mu.3
    h.mig = ifelse(precip >= 2, 0.06, 0)
    v.mig = (pmax(0, exp(-5.4824 + 0.45392 * temp - 0.01252 * 
                           (temp^2))))
    dev1rate <- approxfun(x = dev.1, method = "linear", rule = 2)
    mu1rate <- approxfun(x = mu.1, method = "linear", rule = 2)
    mu2rate <- approxfun(x = mu.2, method = "linear", rule = 2)
    mu3rate <- approxfun(x = mu.3, method = "linear", rule = 2)
    mu4rate <- approxfun(x = mu.4, method = "linear", rule = 2)
    mu5rate <- approxfun(x = mu.5, method = "linear", rule = 2)
    hmigrate <- approxfun(x = h.mig, method = "linear", rule = 2)
    vmigrate <- approxfun(x = v.mig, method = "linear", rule = 2)
    init[5] = (init[4] * (1/v.mig[1])) - init[4]
    egg.correction = rep(1, length(global.t))
    print("Ostertagia ostertagi parameters loaded successfully")
  }
  if (species != "haemonchus" & species != "teladorsagia" & 
      species != "ostertagia") 
    stop("Error: 'species' parameter undefined. Please define nematode species or check spelling")
  event = data.frame(var = "E", time = global.t, value = input * 
                       egg.correction, method = "add")
  para.dyn = function(t, para.init, para.par) {
    with(as.list(c(para.init, para.par)), {
      dev1 = dev1rate(t)
      mu1 = mu1rate(t)
      mu2 = mu2rate(t)
      mu3 = mu3rate(t)
      mu4 = mu4rate(t)
      mu5 = mu5rate(t)
      m1 = hmigrate(t)
      m2 = vmigrate(t)
      dE = -(dev1 * 2 + mu1) * E
      dL = -(dev1 * 2 + mu2) * L + (dev1 * 2) * E
      dL3 = -(mu3 + m1) * L3 + (dev1 * 2) * L
      dPasture = -mu4 * (Pasture * (1 - m2)) - mu5 * (Pasture * 
                                                        m2) + m1 * L3
      return(list(c(dE = dE, dL = dL, dL3 = dL3, dPasture = dPasture)))
    })
  }
  print("model function loaded successfully")
  para.init = c(E = init[1], L = init[2], L3 = init[3], Pasture = (init[4] + 
                                                                     init[5]))
  print("initial conditions for state variables loaded successfully")
  para.sol = lsoda(y = para.init, times = global.t, func = para.dyn, 
                   parms = NULL, events = list(data = event))
  print("simulation finished - review any warning messages above this line")
  Soil = para.sol[, "Pasture"] * (1 - (vmigrate(global.t)))
  Herbage = para.sol[, "Pasture"] * vmigrate(global.t)
  para.sol = cbind(para.sol, Soil, Herbage)
  print("post-hoc calculations successful")
  return(para.sol)
}




# Interpolation function
parainterp = function(para, dates, method="linear") {
  indices = which(seq(as.Date(dates[1]), as.Date(rev(dates)[1]), "days") %in% as.Date(dates))
  days = seq(1, length(seq(as.Date(dates[1]), as.Date(rev(dates)[1]), "days")), 1)
  interpfun = approxfun(y = para, x = indices, method = method)
  return(interpfun(days))
}



# # Example simulation
# install.packages("deSolve")
# install.packages("geosphere")
# require(deSolve)
# require(geosphere)
# temp = c(19, 20, 20.5, 20.6, 18, 17, 19, 20, 21, 21, 21.5, 20, 21, 20)
# precip = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 10, 5, 0.5)
# #generate some random egg counts and dates of egg counts
# fecs = sort(runif(n = 4, min = 50, max = 1000), decreasing = FALSE)
# fecs = fecs * 2000 * 15 #FECs are given in eggs per gram so multiple by 2000 to get number of eggs excreted per sheep, and multiply by 15 to get total number of eggs excreted by all hosts, assuming 15 sheep per hectare
# fec.dates = c("2000-06-01", "2000-06-05", "2000-06-10", "2000-06-14")
# #inteprolate the egg counts using the interpolation function
# fec = parainterp(para = fecs, dates = fec.dates)
# # Run the model:
# # 'init' is the initial number of eggs, L1/L2 combined, L3 in faeces and L3 on herbage at the start of the simulation
# # 'start' and 'end' dates must be input in the format shown
# # 'input' is the daily vector of eggs input on pasture
# # the length of 'temp' and 'precip' must match the length of 'input' and must be equal to the number of days indicated by the start and end date inclusive
# # 'species' options are 'haemonchus', 'teladorsagia' and 'ostertagia'
# # 'lat' is the latitude of the study site in decimal degrees
# test = flsim(init = c(100, 0, 0, 0), start = "2000-06-01", end = "2000-06-14", input = fec, temp = temp, precip = precip, species = "haemonchus", lat = 48)
# head(test)
# plot(test[,"E"], type = "l", xlab = "days", ylab = "number of eggs")
# plot(test[,"Pasture"], type = "l", xlab = "days", ylab = "number of L3 on pasture") 
# 
