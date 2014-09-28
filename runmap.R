# References:
# https://gist.github.com/mollietaylor/4210660
# http://stackoverflow.com/questions/18136468/is-there-a-way-to-add-a-scale-bar-for-linear-distances-to-ggmap

# convert from tcx:
#   gpsbabel -t -i gtrnctr -f test2.tcx -o unicsv -F old.csv

# layout for the image

library(OpenStreetMap)
library(ggplot2)
library(wq)

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
gps <- read.csv("out.csv",
                    header = TRUE)

# Include the distance
gps$dist.prev = 0.0
gps$dist.prev[2:nrow(gps)] = 
    unlist(sapply(2:nrow(gps), 
           function(i) geoDist3D(c(gps$Longitude[i], gps$Latitude[i], gps$Altitude[i]),
                                 c(gps$Longitude[i-1], gps$Latitude[i-1], gps$Altitude[i-1]))))
gps$dist = 0.0
gps$dist[2:nrow(gps)] =
    sapply(2:nrow(gps), function(i) sum(gps$dist.prev[1:i]))

gps$pace = 0.0
gps$pace[2:nrow(gps)] =
    unlist(sapply(2:nrow(gps),
                 function(i) difftime(strptime(paste(gps$Date[i],gps$Time[i]),"%Y/%m/%d %H:%M:%S"),
                                      strptime(paste(gps$Date[i-1],gps$Time[i-1]),"%Y/%m/%d %H:%M:%S"),
                                      unit="s")/
                             gps$dist.prev[i]/60.0))

# Get the map
map =  openmap(c(max(gps$Latitude)+0.01, min(gps$Longitude)-0.01),
               c(min(gps$Latitude)-0.01, max(gps$Longitude)+0.01),
               type="osm-bw")

# Project the map to Long-Lat
mapLL = openproj(map)


# Layout for the figure
layout(matrix(c(1,1,1,1,2,2,2,3,2,2,2,4,2,2,2,5), 4, 4, byrow = TRUE), 
       heights=c(1,2,2,2))



# Plot info
tottime = difftime(strptime(paste(gps$Date[nrow(gps)],gps$Time[nrow(gps)]),
                            "%Y/%m/%d %H:%M:%S"),
                   strptime(paste(gps$Date[1],gps$Time[1]),
                            "%Y/%m/%d %H:%M:%S"),
                   unit="h")
toth = floor(tottime)
totm = round((tottime - toth) * 60)

p0 <- ggplot() +
    annotate("text", x = 0, y = 25, hjust=0, family="courier", fontface="bold",
             label = paste(sep = "", strftime(gps$Date[1],'%A run'), " (",
                                strftime(gps$Date[1],'%Y/%m/%d'),
                                    " @ ",gps$Time[1],")")) +
    annotate("text", x = 0, y = 25, hjust=0, family="courier", fontface="normal",
             label = paste(sep = "", "Total distance: ",sprintf("%.2f", gps$dist[nrow(gps)]), ' km',
                    "\n", "Total time: ", toth, ":", totm)) + 
    theme(line = element_blank(),
          text = element_blank(),
          line = element_blank(),
          panel.background = element_blank(),
          title = element_blank())

# Plot the map
p1 <- autoplot(mapLL, expand=FALSE, main=strftime(gps$Date[1],'%A')) +
      geom_point(aes(x = Longitude,
               y = Latitude),
               data = gps,
               colour = rgb(0,0,0,0.2),
               size = 3) +
    theme(line = element_blank(),
          text = element_blank(),
          line = element_blank(),
          panel.background = element_blank(),
          title = element_blank())


# Altitude plot
p2 <- ggplot() +
        geom_line(aes(x = dist,
                 y = Altitude),
           size = 1.6,
           data = gps,
           colour=rgb(0,0.8,0,0.7)) +
           xlab("distance (km)") +
           ylab("altitude (m)")

# Heart rate plot
p3 <- ggplot() +
    geom_line(aes(x = dist,
                 y = Heartrate),
           data = gps,
           size = 1.6,
           col=rgb(0.8,0,0,0.7)) +
           xlab("distance (km)") +
           ylab("heart rate")

# Pace plot
avepace = ma(gps$pace, 20)
maxpace = max(na.omit(avepace))
if (maxpace>10) {
    maxpace = 10
}
minpace = min(na.omit(avepace))
p4 <- ggplot() +
    geom_line(aes(x = dist,
                 y = avepace),
              size=1.6,
           data = gps, 
           col=rgb(0.2,0.2,0.2,0.7)) +
           xlab("distance (km)") +
           ylab("pace (min/km)") +
           ylim(c(minpace,maxpace))

layOut(list(p0, 1, 1:3),   # takes three rows and the first column
       list(p1, 2:3, 1:3),    # next three are on separate rows
       list(p2, 4, 1), 
       list(p3, 4, 2),
       list(p4, 4, 3))

