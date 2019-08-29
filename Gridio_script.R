#script for scoring vernal pools at a local, neigbhourhood and landscape level
library(gridio) # remember to use the 32 bit version of R
library(raster)
library(SDMTools)
require(sp)
require(classInt)
require(RColorBrewer)

#required functions
RESCALE <- function (x, nx1, nx2, minx, maxx) 
{ nx = nx1 + (nx2 - nx1) * (x - minx)/(maxx - minx)
return(nx)}

##local Level
#input relevant kernel bandwith distances for local, neighbourhood and regional kernels
SD_L <- 124
SD_N <- 399.6
SD_R <- 399.6

#bring in as landcover grid
forest <- readasciigrid("//psf/Home/Documents/original_grids/lc_f_2005.asc",as.matrix=F)

#bring in vernal pool raster
pools <- read.asc('W:/Marbled Salamander (1)/asciis/pvp.asc')
pools <- raster(pools)

#make gaussian kernel of proper scale
r <- res(pools)
stopifnot(isTRUE(all.equal(r[1], r[2])))
cellsize <- r[1]
k <-  makegaussiankernel(sd=SD_L, cellsize = cellsize)

#convert pvp to spatial object
pools.p <- as.data.frame(rasterToPoints(pools))
head(pools.p)
pools <- subset(pools.p, layer=="1")
results <- numeric(nrow(pools))
for(i in 1:(nrow(pools))){
  results[i] <- calckernel(forest, kernel = k, x = pools$x[i], y = pools$y[i])}
pools$local <- results

##neighborhood kernel
#bring in resistance surfaces
resist_surf <- readasciigrid('W:/Marbled Salamander/Genetics/GIS data/AMMA.asc')                 
#rescale resist_surf
#resist_surf$m <- RESCALE(resist_surf$m, 1, 10, minx=min(resist_surf$m), maxx=max(resist_surf$m))
#bring in empty results grid
result <-  forest
result$m[result$m > 0] <- 0
#s <- forest
#for each vernal pool
for(i in 1:(nrow(pools))){
  s <- result$m
  s <- spread(resist_surf$m, row=y2r(pools$y[i],resist_surf), col=x2c(pools$x[i], resist_surf), 
              sd=SD_N, cellsize=resist_surf$cellsize)
  s[y2r(pools$y[i], resist_surf), x2c(pools$x[i], resist_surf)] <- 0
  # s <- s/sum(s) #renormalize
  result$m <- result$m + s
}
pools$c <- x2c(pools$x, resist_surf)
pools$r <- y2r(pools$y, resist_surf)
pools$nhood <- result$m[as.matrix(pools[,c('r', 'c')])]
#save neighborhood grid as ascii
plot(result)
gridinit()
setwindow(result)
#writegrid(result,"//psf/Home/Documents/Original_Grids/AMOP_n_400")

##regional
#make cumulative surface
#result$m
#sample at pools
#for each pool
for(i in 1:(nrow(pools))){
  s <- result$m
  s <- spread(resist_surf$m, row=y2r(pools$y[i],resist_surf), col=x2c(pools$x[i], resist_surf), 
              sd=SD_R, cellsize=resist_surf$cellsize)
  #s[y2r(pools$y[i], resist_surf), x2c(pools$x[i], resist_surf)] <- 0
  # s <- s/sum(s) #renormalize
  result$m <- result$m + s}
#gridinit()
#setwindow(result)
#writegrid(result,"//psf/Home/Documents/Original_Grids/kernel")
#plot(result)
#pools$c <- x2c(pools$x, resist_surf)
#pools$r <- y2r(pools$y, resist_surf)
#pools$nhood <- result$m[as.matrix(pools[,c('r', 'c')])]

#convert grid to binary
result$m[result$m > 0] <- 1
#t1 <- raster(result)
#plot(t1)
#run patch scan to identify patch ids
res <- patchscan(result$m)
#res matrix to grid
resg <- as.grid(res$patches, cellsize=resist_surf$cellsize, xll=87175, yll=867500)
writegrid(resg,"//psf/Home/Documents/Original_Grids/AMMA_Reg_2000")
#sample pools in ID grid
#sample at pools
pools$c <- x2c(pools$x, result)
pools$r <- y2r(pools$y, result)
pools$patch <- res$patches[as.matrix(pools[,c('r', 'c')])]
#table patch ids
table(pools$patch)
p <- as.data.frame(table(pools$patch))
p$patch <- p$Var1
#assign count for each patch back to each pool
x <-  merge(x = pools, y = p, by = "patch", all.x = TRUE)                   
#pools$reg_count <- x$Freq    
#head(pools)

#rescale local, nhood and reg_count
x$local_rescale <- RESCALE(x$local, 1, 10, minx=min(x$local), maxx=max(x$local))
x$nhood_rescale <- RESCALE(x$nhood, 1, 10, minx=min(x$nhood), maxx=max(x$nhood))    
x$reg_count_rescale <- RESCALE(x$Freq, 1, 10, minx=min(x$Freq), maxx=max(x$Freq))
x$geomean <- (x$local_rescale*x$nhood_rescale*x$reg_count_rescale)^(1/3)
x$mean <- (x$local_rescale*x$nhood_rescale*x$reg_count_rescale)/3
