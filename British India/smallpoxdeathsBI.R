setwd("~/Desktop/SOP/Emergence of oscillations/British_India")

library(dplyr)
library(tidyverse)
library(tidyr) 
library(dygraphs)
library(WaveletComp)
library(Rwave)
library(dplR)
library(ggplot2)

#### Plotting Total Deaths from 1870 #####
sp_totaldeaths_BI <- 
  read.table("~/Desktop/SOP/Emergence of oscillations/British_India/smallpox_totaldeaths_BritishIndia.txt",
             quote="\"", comment.char="")

years <- c(1870:1920)

sp_clean <- sp_totaldeaths_BI %>%
  mutate(years,.before = V1) %>%
  mutate(V1 = V1 / 10000) %>%
  rename('deaths' = V1)

dygraph(sp_clean, main= 'Smallpox Deaths in British India from 1870',
        xlab= 'Year', ylab="Total Deaths (in hundred thousands)") %>% 
  dyOptions(stepPlot=TRUE, fillGraph=FALSE)

#### 2d wavelet plot ####

wav_spectra2D_BI <- read.table("~/Desktop/SOP/Emergence of oscillations/British_India/wav_spectra2D_BritishIndia.txt", quote="\"", comment.char="")

#wav_spectra2D_BI$V4 <- wav_spectra2D_BI$V1 + 1870 

wav_spectra2D_BI <- wav_spectra2D_BI %>%
  rename(years = V1, period = V2, power = V3) %>%
  mutate(power = power^2) %>%
  select(years,period,power)

wav_spectra2D_BI <- wav_spectra2D_BI[order(wav_spectra2D_BI$years,wav_spectra2D_BI$period),]


library(wavScalogram)

dt <- 1
time <- seq(0, 50, dt)
signal <- sp_clean$deaths


#### Modifying CWT_WST Function ##### 

wavspect_plot <- function (signal, dt = 1, scales = NULL, powerscales = TRUE, 
          wname = c("MORLET", "DOG", "PAUL", "HAAR", "HAAR2"), wparam = NULL, 
          waverad = NULL, border_effects = c("BE", "PER", "SYM"), makefigure = TRUE, 
          time_values = NULL, energy_density = FALSE, figureperiod = TRUE, 
          xlab = "Time", ylab = NULL, main = NULL, zlim = NULL) {
  wname <- toupper(wname)
  wname <- match.arg(wname)
  if (is.null(waverad)) {
      waverad <- sqrt(2)
  }
  border_effects <- toupper(border_effects)
  border_effects <- match.arg(border_effects)
  nt <- length(signal)
  fourierfactor <- fourier_factor(wname = wname, wparam = wparam)
  if (is.null(scales)) {
    scmin <- 2/fourierfactor
    scmax <- floor(nt/(2 * waverad))
    if (powerscales) {
      scales <- pow2scales(c(scmin, scmax, ceiling(256/log2(scmax/scmin))))
    }
    else {
      scales <- seq(scmin, scmax, by = (scmax - scmin)/256)
    }
    scalesdt <- scales * dt
  }
  else {
    if (powerscales && length(scales) == 3) {
      scales <- pow2scales(scales)
    }
    else {
      if (is.unsorted(scales)) {
        warning("Scales were not sorted.")
        scales <- sort(scales)
      }
      aux <- diff(log2(scales))
      if (powerscales && ((max(aux) - min(aux))/max(aux) > 
                          0.05)) {
        warning("Scales seem like they are not power 2 scales. Powerscales set to FALSE.")
        powerscales <- FALSE
      }
    }
    scalesdt <- scales
    scales <- scales/dt
  }
  ns <- length(scales)
  if (border_effects == "BE") {
    nt_ini <- 1
    nt_end <- nt
  }
  else if (border_effects == "PER") {
    ndatcat <- ceiling(ceiling(waverad * scales[ns])/nt)
    ntcat <- (ndatcat * 2 + 1) * nt
    datcat <- numeric(ntcat)
    for (itcat in 1:ntcat) {
      datcat[itcat] = signal[itcat - floor((itcat - 1)/nt) * 
                               nt]
    }
    signal <- datcat
    nt_ini <- ndatcat * nt + 1
    nt_end <- nt_ini + nt - 1
  }
  else if (border_effects == "SYM") {
    ndatcat <- ceiling(ceiling(waverad * scales[ns])/nt)
    dat_sym <- signal
    aux <- rev(signal)
    for (icat in 1:ndatcat) {
      dat_sym <- c(aux, dat_sym, aux)
      aux <- rev(aux)
    }
    signal <- dat_sym
    nt_ini <- ndatcat * nt + 1
    nt_end <- nt_ini + nt - 1
  }
  nt <- length(signal)
  coefs <- matrix(0, nrow = nt, ncol = ns)
  
    f <- stats::fft(signal)
    k <- 1:trunc(nt/2)
    k <- k * 2 * pi/nt
    k <- c(0, k, -k[trunc((nt - 1)/2):1])
    n <- length(k)
   if (wname == "DOG") {
      if (is.null(wparam)) {
        wparam <- 2
      }
      preexpnt <- (diag(x = scales, ncol = ns) %*% matrix(rep(k, 
                                                              ns), nrow = ns, byrow = T))
      expnt <- -preexpnt^2/2
      Norm <- sqrt(scales * k[2]/gamma(wparam + 0.5)) * 
        sqrt(n)
      daughter <- diag(x = -Norm * ((0+1i)^wparam), ncol = ns) %*% 
        preexpnt^wparam * exp(expnt)
    }
    coefs <- stats::mvfft(t(sweep(daughter, MARGIN = 2, f, 
                                  `*`)), inverse = TRUE)/length(f)

  coefs <- coefs[nt_ini:nt_end, 1:ns] * sqrt(dt)
  nt <- nt_end - nt_ini + 1
  coi_maxscale <- numeric(nt)
  for (i in 1:nt) {
    coi_maxscale[i] <- dt * min(i - 1, nt - i)/waverad
  }
  if (makefigure) {
    if (is.null(time_values)) {
      X <- seq(0, (nt - 1) * dt, dt)
    }
    else {
      if (length(time_values) != nt) {
        warning("Invalid length of time_values vector. Changing to default.")
        X <- seq(0, (nt - 1) * dt, dt)
      }
      else {
        X <- time_values
      }
    }
    ylab_aux <- ylab
    if (figureperiod) {
      Y <- fourierfactor * scalesdt
      coi <- fourierfactor * coi_maxscale
      if (is.null(ylab)) 
        ylab <- "Period"
    }
    else {
      Y <- scalesdt
      coi <- coi_maxscale
      if (is.null(ylab)) 
        ylab <- "Scale"
    }
    if (energy_density) {
      Z <- t(t(abs(coefs)^2)/scalesdt)
      if (is.null(main)) 
        main <- "Wavelet Power Spectrum / Scales"
    }
    else {
      Z <- abs(coefs)^2
      if (is.null(main)) 
        main <- "Wavelet Power Spectrum"
    }
    if (ns > 1) {
      wavPlot(Z = Z, X = X, Y = Y, Ylog = powerscales, 
              coi = coi, Xname = xlab, Yname = ylab, Zname = main, 
              zlim = zlim, Yrev = FALSE)
    }
  }
  return(list(coefs = coefs, scales = scalesdt, fourierfactor = fourierfactor, 
              coi_maxscale = coi_maxscale, ns = ns, X = X, Y=Y, Z = Z, Ylog = powerscales))
}

#### Plotting the Wavelet ######

cwt <- wavspect_plot(signal = signal, dt = 1, energy_density = FALSE, scales = c(0.5,30.1,8),
                     wname = "DOG", xlab = "years (from 1870)", 
                     ylab = "period (years)")

cwt2 <- wavspect_plot(signal = signal, dt = 1, energy_density = TRUE, scales = c(14,46,8),
                    wname = "DOG", xlab = "years (from 1870)", 
                     ylab = "period (years)")

cwt3 <- wavspect_plot(signal = signal, dt = dt, energy_density = FALSE, 
                      wname = "DOG",xlab = "years (from 1870)", ylab = "period (years)")

cwt4 <- wavspect_plot(signal=signal, dt = 1,wname="DOG",makefigure = TRUE, 
                energy_density = TRUE, figureperiod = TRUE,zlim=c(0,60),
                xlab = "years (from 1870)", ylab = "period (years)")
#### DOG ####
test <- DOG(sp_clean$deaths,4,moments = 5,nvoice = 12,plot=FALSE)

z <- abs((Re(test))*15)
x <- seq(0, 50, length.out = nrow(z))
y <- seq(0, 32, length.out = ncol(z))

filled.contour(x,y,z,
               color.palette = function(n) hcl.colors(n, "plasma"),
               plot.title = title(main = "Wavelet Power Spectrum - India",
                                  xlab = "Years, starting 1870", ylab = "Period (years)"),
               key.title = title(main = "Power"),
               plot.axes = { axis(1, seq(0, 50, by = 10))
                 axis(2, seq(0, 32, by = 4)) })

#### Filled Contour ####
testcontour(x,y,z, zlim = c(0,60),
            color.palette = function(n) hcl.colors(n, "sunset"),
            plot.title = title(main = "Wavelet Power Spectrum (India - Log Scale)",
                               xlab = "Years, starting 1870", ylab = "Period (years)"),
            key.title = title(main = "Power"))


#### Global Wavelet Plot ######

wav_global_BI <- read.table("~/Desktop/SOP/Emergence of oscillations/British_India/wav_global_BI.txt", quote="\"", comment.char="")

wav_global_BI <- wav_global_BI %>%
  rename(period = V1, power = V2, sig = V3) %>%
  mutate(power = power / (10^10), sig = sig/(10^10), period = 2^period)

library(scales)

p <- ggplot() + theme_minimal() +
  # blue plot
  geom_path(data=wav_global_BI, aes(x=power, y=period),color='blue') +
  scale_y_continuous(trans = log2_trans()) + xlim(0,3) +
  # red plot
  geom_path(data=wav_global_BI, aes(x=sig, y=period),color='red') + 
  ggtitle('Global Wavelet Spectrum') + xlab('Power') + ylab('Period (years)')

p