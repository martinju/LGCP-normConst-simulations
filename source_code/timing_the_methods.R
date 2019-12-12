
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

source("Source_code/help_functions.R")


# # Settings providing quick results
#no.z.samps <- 700
these.k.z.samps <- 1:1000#1:2500#1:7000
seed.z.sampling.initial <- 123467
xy.range <- c(0,16)
xy.range.approx.mesh <- c(-2,18)
npixel <- 20*16
no.MC.intpoints <- npixel*diff(xy.range)/diff(xy.range.approx.mesh)*10 # OK as long as multiplyer of npixel
parallelize.numCores <- 12

boundary.splitter.vec <- c(1,2,4,8)
INLA.max.edge.vec <-  c(4,3,2,1,0.5,0.25)#seq(4,0.25,by=-0.25)#c(4,2,1,0.5,0.25)#c(4,2)#c(0.25,0.5,1,2,4)
true.field.variance.vec <- c(0.1,0.2,0.5,1)#0.5#c(0.25,0.5,1)#c(0.5,1)
beta0.true <- 1
kappa.true.vec <- 1 # kappa = sqrt(8)/range
FEM.approx.mesh <- TRUE
INLA.offset <- c(-0.03,-0.1) 
BarycentricManyPoints.nsub <- 29
folder_string = NULL#"20191023_0945_iHcBtj" #"20191021_1602_5KUF25" #"20191018_1210_ZH8vu8"#"20191018_1057_HnU4zC" # Provide a folder string here to use pre-computed prep.mesh.list. Need to have exactly the same input.list provided below
do_debug <- FALSE

# # Big run!
# no.z.samps <- 5000
# seed.z.sampling.initial <- 123467
# xy.range <- c(0,16)
# npixel <- 16*16
# no.MC.intpoints <- npixel*8 # OK as long as multiplyer of npixel
# parallelize.numCores <- 10
# 
# boundary.splitter.vec <- c(1,2,4,8)
# INLA.max.edge.vec <- c(4,2,1,0.5,0.25)#c(0.25,0.5,1,2,4)
# true.field.variance.vec <- c(0.1,0.2,0.5,1)
# beta0.true <- 0
# kappa.true <- 2
# FEM.approx.mesh <- TRUE
# INLA.offset <- -0.03 # This works
# BarycentricManyPoints.nsub <- 100


input.list <- list(these.k.z.samps = these.k.z.samps,
                   seed.z.sampling.initial = seed.z.sampling.initial,
                   xy.range = xy.range,
                   xy.range.approx.mesh = xy.range.approx.mesh,
                   npixel = npixel,
                   no.MC.intpoints = no.MC.intpoints,
                   parallelize.numCores = parallelize.numCores,
                   boundary.splitter.vec = boundary.splitter.vec,
                   INLA.max.edge.vec = INLA.max.edge.vec,
                   true.field.variance.vec = true.field.variance.vec,
                   beta0.true = beta0.true,
                   kappa.true.vec = kappa.true.vec,
                   FEM.approx.mesh = FEM.approx.mesh,
                   INLA.offset = INLA.offset,
                   BarycentricManyPoints.nsub = BarycentricManyPoints.nsub)

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

fwrite(dt,"../Results_new/timing_results.csv")



