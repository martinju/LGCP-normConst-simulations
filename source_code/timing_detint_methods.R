
#### Estimating the integral of the exponential latent field

rm(list=ls())

library(spatstat)
library(RandomFields)
library(INLA)
library(inlabru)
library(rgeos)
library(proxy)
library(stringi)
library(prevR)
library(inlabru) # v (2.1.8 or newer)
library(data.table)
library(foreach)
library(doMC)

source("source_code/help_functions.R")


# # Settings providing quick results
these.k.z.samps <- 1:10
seed.z.sampling.initial <- 123467
xy.range <- c(0,16)
xy.range.approx.mesh <- c(-2,18)
npixel <- 20*16
no.MC.intpoints <- npixel*diff(xy.range)/diff(xy.range.approx.mesh)*10 # OK as long as multiplyer of npixel
parallelize.numCores <- 12

boundary.splitter.vec <- c(1,2)
INLA.max.edge.vec <-  c(4,3)
true.field.variance.vec <- c(0.1,0.2)
beta0.true <- 1
kappa.true.vec <- sqrt(8)/c(1) # kappa = sqrt(8)/range
FEM.approx.mesh <- TRUE
INLA.offset <- c(-0.03,-0.1) 
BarycentricManyPoints.nsub <- 29
folder_string = NULL# Provide a folder string here to use pre-computed prep.mesh.list. Need to have exactly the same input.list provided below
do_debug <- FALSE

#UNCOMMENT THE BELOW FOR FULL RUN
# # Settings producing the full simulaton run (several days of computation time), 
# these.k.z.samps <- 1:10000
# seed.z.sampling.initial <- 123467
# xy.range <- c(0,16)
# xy.range.approx.mesh <- c(-2,18)
# npixel <- 20*16
# no.MC.intpoints <- npixel*diff(xy.range)/diff(xy.range.approx.mesh)*10 # OK as long as multiplyer of npixel
# parallelize.numCores <- 12
# 
# boundary.splitter.vec <- c(1,2,4,8)
# INLA.max.edge.vec <-  c(4,3,2,1,0.5,0.25)
# true.field.variance.vec <- c(0.1,0.2,0.5,1)
# beta0.true <- 1
# kappa.true.vec <- sqrt(8)/c(8,6,4,3,2,1,1.5,0.75,0.5,0.25) # kappa = sqrt(8)/range
# FEM.approx.mesh <- TRUE
# INLA.offset <- c(-0.03,-0.1) 
# BarycentricManyPoints.nsub <- 29
# folder_string = NULL# Provide a folder string here to use pre-computed prep.mesh.list. Need to have exactly the same input.list provided below
# do_debug <- FALSE


results_folder <- "results"

#### Various preparations ####
win <- owin(xy.range,xy.range)
nn.mesh.win <- owin(xy.range.approx.mesh+c(-2,2),xy.range.approx.mesh+c(-2,2)) # Only for use in the ppp-function making the nn.mesh object (which is no longer used)
boundary <- rbind(rep(xy.range[1],2),
                  xy.range,
                  rep(xy.range[2],2),
                  rev(xy.range),
                  rep(xy.range[1],2),deparse.level = 3)
A <- boundary

xyvals <- (seq(xy.range.approx.mesh[1],xy.range.approx.mesh[2],length.out=npixel+1)[-(npixel+1)] + seq(xy.range.approx.mesh[1],xy.range.approx.mesh[2],length.out=npixel+1)[-1])/2
rf.xy.grid <- expand.grid(x=xyvals, y=xyvals)
area.pixel = diff(xyvals)[1]^2

boundary.mesh <- rbind(rep(xy.range.approx.mesh[1],2),
                       xy.range.approx.mesh,
                       rep(xy.range.approx.mesh[2],2),
                       rev(xy.range.approx.mesh),
                       rep(xy.range.approx.mesh[1],2),deparse.level = 3)


#### Preparations which can be done outside the looping ####

prep.A.list <- list()

for (l in 1:length(boundary.splitter.vec)){
  boundary.splitter <- boundary.splitter.vec[l]
  
  A_SpatialPolygons <- boundary.splitter.func(boundary = boundary,boundary.splitter = boundary.splitter)
  A_SpatialPolygons2 <- boundary.splitter.func2(boundary = boundary,boundary.splitter = boundary.splitter)
  
  A_SpatialLines <- as(A_SpatialPolygons,"SpatialLines")
  
  poly1 <- A_SpatialPolygons@polygons
  poly2 <- unlist(sapply(poly1, FUN=function(x) x@Polygons))
  A_coords <- unique(data.table::rbindlist(lapply(poly2, FUN=function(x) data.table(x@coords))))
  A_SpatialPoints <- sp::SpatialPoints(A_coords)
  
  A_List_of_SpatialPolygons <- lapply(poly2,FUN = function(x) sp::SpatialPolygons(list(Polygons(list(x),"0"))))
  
  
  Likelihood.OracleTruth.prepare.list <- Likelihood.OracleTruth.prepare.new(A_SpatialPolygons = A_SpatialPolygons,
                                                                            rf.xy.grid = rf.xy.grid)
  
  # NOTE: There is a small bias here, probably, as the points may hit exactly the boundary of domains. Ok when no.MC.intpoints is a multiplyer of npixel*diff(xy.range)/diff(extended.xy.range)
  MonteCarloIntegrator.prepare.list.0 <- MonteCarloIntegrator.prepare.0(boundary = boundary,
                                                                        boundary.splitter = boundary.splitter,
                                                                        no.MC.intpoints = no.MC.intpoints,
                                                                        A_SpatialPolygons = A_SpatialPolygons)
  OmegaArea <- rgeos::gArea(A_SpatialPolygons)
  
  prep.A.list[[l]] <- list(A_SpatialPolygons = A_SpatialPolygons,
                           A_SpatialPolygons2 = A_SpatialPolygons2,
                           A_SpatialLines = A_SpatialLines,
                           A_SpatialPoints = A_SpatialPoints,
                           A_List_of_SpatialPolygons = A_List_of_SpatialPolygons,
                           Likelihood.OracleTruth.prepare.list = Likelihood.OracleTruth.prepare.list,
                           MonteCarloIntegrator.prepare.list.0 = MonteCarloIntegrator.prepare.list.0,
                           OmegaArea = OmegaArea)
}
#### Preparing all the mesh stuff first ####

Voronoi.time.mat <- DualMesh.time.mat <- DualMeshExtra.time.mat <- BarycentricPointSpread.time.mat <- 
  BarycentricPointSpreadManyPoints.time.mat <- matrix(NA,nrow=length(INLA.max.edge.vec),ncol=length(boundary.splitter.vec))


#uu=proc.time()
prep.mesh.list <- list()
for (i in 1:length(INLA.max.edge.vec)){
  INLA.max.edge <- INLA.max.edge.vec[i]
  
  OmegaMesh <- INLA::inla.mesh.2d(loc.domain=boundary,
                                  max.edge=INLA.max.edge, 
                                  cutoff=0.3,
                                  offset=INLA.offset) # Should try to see the effect of settign this larger
  
  basisfunc.grid.mat.0 = prepare.basisfunc.grid.mat.new(OmegaMesh=OmegaMesh,rf.xy.grid = rf.xy.grid)
  basisfunc.grid.mat <- basisfunc.grid.mat.0$mat
  affected.mesh.points <- basisfunc.grid.mat.0$affected.mesh.points
  
  
  inprodmat = basisfunc.grid.mat%*%t(basisfunc.grid.mat)*area.pixel # Could have used inla.mesh.fem(OmegaMesh)$c1, but tests gave slightly better results using my Monte Carlo version
  
  
  ### Finding the mapping from the field points to the mesh locations
  nn.mesh=nncross(X=spatstat::ppp(OmegaMesh$loc[,1],OmegaMesh$loc[,2],window = nn.mesh.win),
                  Y=spatstat::ppp(rf.xy.grid[,1],rf.xy.grid[,2],window=nn.mesh.win))
  
  for (l in 1:length(boundary.splitter.vec)){
    
    
    A_SpatialPolygons <- prep.A.list[[l]]$A_SpatialPolygons
    A_SpatialPolygons2 <- prep.A.list[[l]]$A_SpatialPolygons2
    A_SpatialLines <- prep.A.list[[l]]$A_SpatialLines
    A_SpatialPoints <- prep.A.list[[l]]$A_SpatialPoints
    A_List_of_SpatialPolygons <- prep.A.list[[l]]$A_List_of_SpatialPolygons
    Likelihood.OracleTruth.prepare.list <- prep.A.list[[l]]$Likelihood.OracleTruth.prepare.list
    MonteCarloIntegrator.prepare.list.0 <- prep.A.list[[l]]$MonteCarloIntegrator.prepare.list.0
    OmegaArea <- prep.A.list[[l]]$OmegaArea
    
    ### Mesh dependent method preparations ####
    MonteCarloIntegrator.mapper <- MonteCarloIntegrator.prepare.1(A.xyval.in.A = MonteCarloIntegrator.prepare.list.0$A.xyval.in.A,
                                                                  OmegaMesh = OmegaMesh)
    
    start <- proc.time()
    
    
    CoxLikVoronoi.OmegaVerticesWeight <- CoxLikVoronoi.prepare(OmegaMesh = OmegaMesh,
                                                               A_SpatialPolygons=A_SpatialPolygons)
    
    int1 <- proc.time()
    
    CoxLikDualMesh.OmegaVerticesWeight <- CoxLikDualMesh.prepare(OmegaMesh = OmegaMesh,
                                                                 A_SpatialPolygons=A_SpatialPolygons) 
    
    int2 <- proc.time()
    
    
    CoxLikDualMeshExtra.OmegaVerticesWeightList <- CoxLikDualMeshExtra.prepare.new(OmegaMesh = OmegaMesh,
                                                                                   A_SpatialPolygons = A_SpatialPolygons,
                                                                                   A_SpatialLines = A_SpatialLines,
                                                                                   A_SpatialPoints = A_SpatialPoints)
    
    int3 <- proc.time()
    
    
    CoxLikBarycentricPointSpread.OmegaVerticesWeight <- CoxLikBarycentricPointSpread.prepare.new3(OmegaMesh = OmegaMesh,
                                                                                                  A_SpatialPolygons = A_SpatialPolygons2)
    
    int4 <- proc.time()
    
    
    CoxLikBarycentricPointSpreadManyPoints.OmegaVerticesWeight <- CoxLikBarycentricPointSpread.prepare.new3.nsub(OmegaMesh = OmegaMesh,
                                                                                                                 A_SpatialPolygons = A_SpatialPolygons2,
                                                                                                                 nsub = BarycentricManyPoints.nsub)
    stop <- proc.time()
    
    
    Voronoi.time.mat[i,l] <- int1[3] - start[3]
    DualMesh.time.mat[i,l]<- int2[3] - int1[3]
    DualMeshExtra.time.mat[i,l] <- int3[3] - int2[3]
    BarycentricPointSpread.time.mat[i,l] <- int4[3] - int3[3]
    BarycentricPointSpreadManyPoints.time.mat[i,l] <- stop[3] - int4[3]
    
    print(c(i,l))
    
  }

}

#Voronoi.time.mat <- DualMesh.time.mat <- DualMeshExtra.time.mat <- BarycentricPointSpread.time.mat <- 
#  BarycentricPointSpreadManyPoints.time.mat <- matrix(NA,nrow=length(INLA.max.edge.vec),ncol=length(boundary.splitter.vec))


Voronoi.time.vec <- unlist(lapply(Voronoi.time.list,function(x)x[3]))
DualMesh.time.vec <- unlist(lapply(DualMesh.time.list,function(x)x[3]))
DualMeshExtra.time.vec <- unlist(lapply(DualMeshExtra.time.mat,function(x)x[3]))
BarycentricPointSpread.time.vec <- unlist(lapply(BarycentricPointSpread.time.mat,function(x)x[3]))
BarycentricPointSpreadManyPoints.time.vec <- unlist(lapply(BarycentricPointSpreadManyPoints.time.mat,function(x)x[3]))

dt <- data.table::data.table(cbind(INLA.max.edge=rep(INLA.max.edge.vec,times=length(boundary.splitter.vec)),
                                   boundary.splitter = rep(boundary.splitter.vec,each=length(INLA.max.edge.vec)),
                                   Voronoi = c(Voronoi.time.mat),
                                   DualMesh = c(DualMesh.time.mat),
                                   DualMeshExtra = c(DualMeshExtra.time.mat),
                                   BarycentricPointSpread = c(BarycentricPointSpread.time.mat),
                                   BarycentricPointSpreadManyPoints = c(BarycentricPointSpreadManyPoints.time.mat)))

set.seed(round(proc.time()[3]*1000))
rand_string <- stringi::stri_rand_strings(1,length=6)
folder_string = paste0(format(Sys.time(), "%Y%m%d_%H%M"),"_",rand_string)

dir.create(path = file.path(results_folder,"timing",folder_string),recursive = T)
fwrite(dt,file.path(results_folder,"timing",folder_string,"timing_results.csv"))




