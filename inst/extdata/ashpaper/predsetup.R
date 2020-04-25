pred=icenolakes[,,701:1]
init.ext=c(dim(rs)[1],dim(rs)[2])
keep.thresh=0.05
corners <- (matrix(c( e[1], e[4],
                         e[1], e[3],
                     e[2], e[3],
                     e[2],e[4]),ncol=2,byrow=T))
colnames(corners) <- c("x","y")
rownames(corners) <- c("ul","ll","lr","ur")
keep.thresh = 0.1
sim.epsg=5070
range.epsg=4326
raster.proj='+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs'

