setwd("~/Desktop/SOP/Emergence of oscillations/Sweden")

library(dplyr)
library(tidyverse)
library(tidyr) 
library(dygraphs)
library(WaveletComp)
library(Rwave)
library(dplR)
library(ggplot2)
library("biwavelet")
library(plotly)
library(sp)

#### Wavelet Spectrum Plot ####

wav_spectra2D_Sweden <- 
  read.table("~/Desktop/SOP/Emergence of oscillations/Sweden/wav_spectra2D_Sweden.txt", 
             quote="\"", comment.char="")
View(wav_spectra2D_Sweden)

wavspec_swed <- wav_spectra2D_Sweden %>%
  rename(Years = V1, Period = V2, Power = V3) %>% 
  mutate(Power = Power^2) %>%
  select(Years,Period,Power)


#### Contour Custom Function #####
testcontour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                       length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
          ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
          levels = pretty(zlim, nlevels), nlevels = 20, color.palette = function(n) hcl.colors(n, 
                                                                                               "YlOrRd", rev = TRUE), col = color.palette(length(levels) - 
                                                                                                                                            1), plot.title, plot.axes, key.title, key.axes, asp = NA, 
          xaxs = "i", yaxs = "i", las = 1, axes = TRUE, frame.plot = axes, 
          ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "y", xaxs = xaxs, yaxs = yaxs, asp = asp)
  .filled.contour(x, y, z, levels, col)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}


#### Plotting #####


z <- cwt_swed$Z
x <- seq(1, 25, length.out = nrow(z))
y <- seq(1, 32, length.out = ncol(z))


filled.contour(x,y,z,zlim=c(0,60),
            color.palette = function(n) hcl.colors(n, "sunset"),
            plot.title = title(main = "Wavelet Power Spectrum (Sweden) ",
                               xlab = "Years, starting 1774", ylab = "Period (years)"),
            key.title = title(main = "Power"))

testcontour(x,y,z, zlim = c(0,60),
               color.palette = function(n) hcl.colors(n, "sunset"),
               plot.title = title(main = "Wavelet Power Spectrum (Sweden - Log Scale)",
                                  xlab = "Years, starting 1870", ylab = "Period (years)"),
               key.title = title(main = "Power"))

mtext(paste("filled.contour(.) from", R.version.string),
      side = 1, line = 4, adj = 1, cex = .66)

#### Transforming Manually #####

sptotaldeaths_sweden <- read.table("~/Desktop/SOP/Emergence of oscillations/Sweden/smallpox_totaldeaths_Sweden.txt", quote="\"", comment.char="")

years <- c(1774:1800)

sp_clean2 <- sptotaldeaths_sweden %>%
  mutate(years,.before = V1) %>%
  mutate(V1 = V1 / 10000) %>%
  rename('deaths' = V1)

dygraph(sp_clean2, main= 'Smallpox Deaths in Sweden from 1774',
        xlab= 'Year', ylab="Total Deaths (in hundred thousands)") %>% 
  dyOptions(stepPlot=TRUE, fillGraph=FALSE)
#### DOG####
test2 <- DOG(sp_clean2$deaths,4,moments = 15,nvoice = 12,plot=FALSE)

z <- abs(Re(test2)*10^6*50)
x <- seq(1, 25, length.out = nrow(z))
y <- seq(1, 32, length.out = ncol(z))

test3 <- DOG(sp_clean2$deaths,4,moments = 5,nvoice = 12,plot=FALSE)

z <- abs(Re(test3)*6*60)
x <- seq(1, 25, length.out = nrow(z))
y <- seq(1, 32, length.out = ncol(z))

filled.contour(x,y,z,
               color.palette = function(n) hcl.colors(n, "plasma"),
               plot.title = title(main = "Wavelet Power Spectrum - Sweden",
                                  xlab = "Years, starting 1774", ylab = "Period (years)"),
               key.title = title(main = "Power"),
               plot.axes = { axis(1, seq(0, 25, by = 5))
                 axis(2, seq(0, 32, by = 4)) })


signal <- sp_clean2$deaths
dt <- 1
cwt_swed <- wavspect_plot(signal = signal, dt = dt, energy_density = FALSE, zlim=c(0,70),
                      wname = "DOG",xlab = "years (from 1774)", ylab = "period (years)")




