#!/usr/bin/env Rscript
#
# References:
# https://gist.github.com/mollietaylor/4210660
# http://stackoverflow.com/questions/18136468/is-there-a-way-to-add-a-scale-bar-for-linear-distances-to-ggmap

# convert from tcx:
#   gpsbabel -t -i gtrnctr -f test2.tcx -o unicsv -F old.csv

# layout for the image

sink("/dev/null")

library(OpenStreetMap)
library(ggplot2)
library(wq)
library(extrafont)
library(optparse)

option_list <- list(
       make_option(c("-i", "--infile"), default="in.csv",
                   help="File containing the GPS data (csv) [default %default]"),
       make_option(c("-t", "--imgtype"), default="svg",
                   help="Type of the resulting image file (svg, png or pdf) [default %default]"),
       make_option(c("-o", "--outfile"), default="out.svg",
                   help="Filename for the image [default %default]")
       )

opt <- parse_args(OptionParser(option_list=option_list))

imgtype=opt$imgtype

# Name of the input file
infile=opt$infile

# Name of the output file
outfile=opt$outfile

# Output file type
if (imgtype == "svg") {
    svg(outfile)
} else if (imgtype == "png") {
    png(outfile)
} else if (imgtype == "pdf") {
    pdf(outfile)
} else {
    print("Wrong image type selected!")
    q()
}

# Running average (trailing)
ma <- function(x,n=15){filter(x,rep(1/n,n), sides=1)}

# Calculates the geodesic distance between two points specified by radian
# latitude/longitude using Vincenty inverse formula for ellipsoids (vif)
#
# From: [http://www.r-bloggers.com/great-circle-distance-calculations-in-r/]
geoDist <- function(long1, lat1, long2, lat2) {

  # WGS-84 ellipsoid parameters
  a <- 6378137         # length of major axis of the ellipsoid (radius at equator)
  b <- 6356752.314245  # ength of minor axis of the ellipsoid (radius at the poles)
  f <- 1/298.257223563 # flattening of the ellipsoid

  L <- long2-long1 # difference in longitude
  U1 <- atan((1-f) * tan(lat1)) # reduced latitude
  U2 <- atan((1-f) * tan(lat2)) # reduced latitude
  sinU1 <- sin(U1)
  cosU1 <- cos(U1)
  sinU2 <- sin(U2)
  cosU2 <- cos(U2)

  cosSqAlpha <- NULL
  sinSigma <- NULL
  cosSigma <- NULL
  cos2SigmaM <- NULL
  sigma <- NULL

  lambda <- L
  lambdaP <- 0
  iterLimit <- 100
  while (abs(lambda-lambdaP) > 1e-12 & iterLimit>0) {
    sinLambda <- sin(lambda)
    cosLambda <- cos(lambda)
    sinSigma <- sqrt( (cosU2*sinLambda) * (cosU2*sinLambda) +
                      (cosU1*sinU2-sinU1*cosU2*cosLambda) * (cosU1*sinU2-sinU1*cosU2*cosLambda) )
    if (sinSigma==0) return(0)  # Co-incident points
    cosSigma <- sinU1*sinU2 + cosU1*cosU2*cosLambda
    sigma <- atan2(sinSigma, cosSigma)
    sinAlpha <- cosU1 * cosU2 * sinLambda / sinSigma
    cosSqAlpha <- 1 - sinAlpha*sinAlpha
    cos2SigmaM <- cosSigma - 2*sinU1*sinU2/cosSqAlpha
    if (is.na(cos2SigmaM)) cos2SigmaM <- 0  # Equatorial line: cosSqAlpha=0
    C <- f/16*cosSqAlpha*(4+f*(4-3*cosSqAlpha))
    lambdaP <- lambda
    lambda <- L + (1-C) * f * sinAlpha *
              (sigma + C*sinSigma*(cos2SigmaM+C*cosSigma*(-1+2*cos2SigmaM*cos2SigmaM)))
    iterLimit <- iterLimit - 1
  }
  if (iterLimit==0) return(NA)  # formula failed to converge
  uSq <- cosSqAlpha * (a*a - b*b) / (b*b)
  A <- 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)))
  B <- uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)))
  deltaSigma = B*sinSigma*(cos2SigmaM+B/4*(cosSigma*(-1+2*cos2SigmaM^2) -
                                      B/6*cos2SigmaM*(-3+4*sinSigma^2)*(-3+4*cos2SigmaM^2)))
  s <- b*A*(sigma-deltaSigma) / 1000

  return(s) # Distance in km
}

# Include the altitude difference
geoDist3D <- function(p1, p2) {
    d = geoDist(p1[1]*0.01745329, p1[2]*0.01745329, p2[1]*0.01745329, p2[2]*0.01745329)
    dAlt = (p1[3]-p2[3])/1000.0
#   return(sqrt(d*d + dAlt*dAlt))
    return(d)
}


# Load the GPS file
gps <- read.csv(infile, header = TRUE)

# Include the distance
gps$dist.prev = 0.0
gps$dist.prev[2:nrow(gps)] = 
    unlist(sapply(2:nrow(gps), 
           function(i) geoDist3D(c(gps$Longitude[i], gps$Latitude[i], gps$Altitude[i]),
                                 c(gps$Longitude[i-1], gps$Latitude[i-1], gps$Altitude[i-1]))))
gps$dist = 0.0
gps$dist[2:nrow(gps)] =
    sapply(2:nrow(gps), function(i) sum(gps$dist.prev[1:i]))


# Get the map
map =  openmap(c(max(gps$Latitude)+0.005, min(gps$Longitude)-0.005),
               c(min(gps$Latitude)-0.005, max(gps$Longitude)+0.005),
               type="osm")

# Project the map to Long-Lat
mapLL = openproj(map)


# Plot info
tottime = difftime(strptime(paste(gps$Date[nrow(gps)],gps$Time[nrow(gps)]),
                            "%Y/%m/%d %H:%M:%S"),
                   strptime(paste(gps$Date[1],gps$Time[1]),
                            "%Y/%m/%d %H:%M:%S"),
                   unit="h")

# Calculate times in proper unit for summary
toth = floor(tottime)
totm = round((tottime - toth) * 60)
paces = tottime*60*60/gps$dist[nrow(gps)]
pacem = floor(paces/60.0)
paces = round(paces - pacem*60.0)


# Paces measured point-by-point are inaccurage; use 100 m divisions
maxdist = gps$dist[nrow(gps)]
d = seq(1,floor(maxdist*10))*0.1
p = seq(1,floor(maxdist*10))*0.1
index = 1
for (i in 1:floor(maxdist*10)) {
    indexlast = index
    index = which(min(abs(gps$dist-i*0.1))==abs(gps$dist-i*0.1))
    dt = difftime(strptime(paste(gps$Date[index],gps$Time[index]),"%Y/%m/%d %H:%M:%S"),
                  strptime(paste(gps$Date[indexlast],gps$Time[indexlast]),"%Y/%m/%d %H:%M:%S"),
                  unit="s")
    dx = gps$dist[index] - gps$dist[indexlast]
    p[i] = dt/60.0/dx
}

# Plot summary text
p0 <- ggplot() +
    annotate("text", x = 0, y = 0.85, cex=6.0,  hjust=0, fontface="bold", family="Courier New",
             label = paste(sep = "", strftime(gps$Date[1],'%A run'), " (",
                                strftime(gps$Date[1],'%Y/%m/%d'),
                                    " @ ",gps$Time[1],")")) +
    annotate("text", x = 0, y = 0.30, cex=5.0, hjust=0, fontface="plain", family="Courier New",
             lineheight=0.80,
             label = paste(sep = "", "Total distance: ",sprintf("%.2f", gps$dist[nrow(gps)]), ' km',
                    "\n", "Total time: ", toth, ":", sprintf("%02d", totm), "\n",
                    "Average pace: ", pacem, ":", sprintf("%02d", paces), " min/km")) + 
    theme(line = element_blank(),
          text = element_blank(),
          line = element_blank(),
          panel.background = element_blank(),
          title = element_blank()) +
    xlim(0:1) + ylim (0:1)

# Plot the map
p1 <- autoplot(mapLL) +
      geom_point(aes(x = Longitude,
               y = Latitude),
               data = gps,
               colour = rgb(0,0,0.8,0.5),
               size = 2) +
      geom_point(aes(x = Longitude,
               y = Latitude),
               data = gps[1,],
               colour = rgb(0,0,0.8,0.8),
               size = 5) +
      geom_point(aes(x = Longitude,
               y = Latitude),
               data = gps[nrow(gps),],
               colour = rgb(0,0.8,0,0.8),
               size = 5) +
    theme(line = element_blank(),
          text = element_blank(),
          line = element_blank(),
          panel.background = element_blank(),
          title = element_blank())

# Add markers for each kilometer
for (i in 1:floor(gps$dist[nrow(gps)])) {

    index = which(min(abs(gps$dist-i))==abs(gps$dist-i))
    p1 = p1 + geom_point(aes(x = Longitude,
                y = Latitude),
               data = gps[index,],
               colour = rgb(1,1,1,0.6),
                size = 8) +
    annotate("text", x = gps$Longitude[index], y = gps$Latitude[index],
              label = sprintf("%d",i))


}

# Altitude plot
p2 <- ggplot() +
        geom_line(aes(x = dist,
                 y = Altitude),
           size = 1.6,
           data = gps,
           colour=rgb(0,0.8,0,0.7)) +
           xlab("distance (km)") +
           ylab("altitude (m)") +
           theme(axis.title.x = element_text(family="Droid Sans Mono")) +
           theme(axis.title.y = element_text(family="Droid Sans Mono"))


# Pace plot
maxpace = max(p)
if (maxpace > 10.0) {
    maxpace = 10.0
}
p3 <- ggplot() +
    geom_line(aes(x = d,
                 y = p),
              size=1.6,
           col=rgb(0.2,0.2,0.2,0.7)) +
           xlab("distance (km)") +
           ylab("pace (min/km)") +
           theme(axis.title.x = element_text(family="Droid Sans Mono"),
                 axis.title.y = element_text(family="Droid Sans Mono")) +
           coord_cartesian(ylim = c(min(p), maxpace))

          
# Heart rate plot
p4 <- ggplot() +
    geom_line(aes(x = dist,
                 y = Heartrate),
           data = gps,
           size = 1.6,
           col=rgb(0.8,0,0,0.7)) +
           xlab("distance (km)") +
           ylab("heart rate") +
           theme(axis.title.x = element_text(family="Droid Sans Mono")) +
           theme(axis.title.y = element_text(family="Droid Sans Mono"))


layOut(list(p0, 1, 1:3),
       list(p1, 2:4, 1:3),
       list(p2, 5, 1), 
       list(p3, 5, 2),
       list(p4, 5, 3))

