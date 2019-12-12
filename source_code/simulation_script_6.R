
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
kappa.true.vec <- sqrt(8)/c(8,6,4,3,2,1,1.5,0.75,0.5,0.25) # kappa = sqrt(8)/range
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

doMC::registerDoMC(parallelize.numCores)


if(is.null(folder_string)){
  
  #uu=proc.time()
  prep.mesh.list <- list()
  for (i in 1:length(INLA.max.edge.vec)){
    INLA.max.edge <- INLA.max.edge.vec[i]
    
    #### Fixing the mesh ####
    #  OmegaMesh <- INLA::inla.mesh.2d(loc.domain=boundary,
    #                            offset=c(.3),
    #                            max.edge=INLA.max.edge, 
    #                            cutoff=.05)
    OmegaMesh <- INLA::inla.mesh.2d(loc.domain=boundary,
                                    max.edge=INLA.max.edge, 
                                    cutoff=0.3,
                                    offset=INLA.offset) # Should try to see the effect of settign this larger
    
    # interior.bnd = inla.mesh.segment(loc=boundary,is.bnd = F)
    # boundary.bnd = inla.mesh.segment(loc=boundary.mesh[-5,],is.bnd=T)
    # 
    # OmegaMesh <- INLA::inla.mesh.2d(loc.domain=boundary,
    #                                 boundary = boundary.bnd,
    #                                 max.edge=INLA.max.edge, 
    #                                 cutoff=0.3)
    #     
    # par(mar=c(0,0,0,0))
    # plot(OmegaMesh, asp=1,add=F)
    # lines(boundary, col=3)
    
    basisfunc.grid.mat.0 = prepare.basisfunc.grid.mat.new(OmegaMesh=OmegaMesh,rf.xy.grid = rf.xy.grid)
    basisfunc.grid.mat <- basisfunc.grid.mat.0$mat
    affected.mesh.points <- basisfunc.grid.mat.0$affected.mesh.points
    
    
    inprodmat = basisfunc.grid.mat%*%t(basisfunc.grid.mat)*area.pixel # Could have used inla.mesh.fem(OmegaMesh)$c1, but tests gave slightly better results using my Monte Carlo version
    
    
    ### Finding the mapping from the field points to the mesh locations
    nn.mesh=nncross(X=spatstat::ppp(OmegaMesh$loc[,1],OmegaMesh$loc[,2],window = nn.mesh.win),
                    Y=spatstat::ppp(rf.xy.grid[,1],rf.xy.grid[,2],window=nn.mesh.win))
    
    prep.mesh.single.A.list <- foreach::foreach(l=1:length(boundary.splitter.vec),
                                                .verbose=FALSE,.inorder=TRUE) %dopar% {
                                                  
                                                  
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
                                                  
                                                  CoxLikVoronoi.OmegaVerticesWeight <- CoxLikVoronoi.prepare(OmegaMesh = OmegaMesh,
                                                                                                             A_SpatialPolygons=A_SpatialPolygons)
                                                  
                                                  CoxLikDualMesh.OmegaVerticesWeight <- CoxLikDualMesh.prepare(OmegaMesh = OmegaMesh,
                                                                                                               A_SpatialPolygons=A_SpatialPolygons) 
                                                  
                                                  CoxLikDualMeshExtra.OmegaVerticesWeightList <- CoxLikDualMeshExtra.prepare.new(OmegaMesh = OmegaMesh,
                                                                                                                                 A_SpatialPolygons = A_SpatialPolygons,
                                                                                                                                 A_SpatialLines = A_SpatialLines,
                                                                                                                                 A_SpatialPoints = A_SpatialPoints)
                                                  
                                                  CoxLikBarycentricPointSpread.OmegaVerticesWeight <- CoxLikBarycentricPointSpread.prepare.new3(OmegaMesh = OmegaMesh,
                                                                                                                                                A_SpatialPolygons = A_SpatialPolygons2)
                                                  
                                                  CoxLikBarycentricPointSpreadManyPoints.OmegaVerticesWeight <- CoxLikBarycentricPointSpread.prepare.new3.nsub(OmegaMesh = OmegaMesh,
                                                                                                                                                               A_SpatialPolygons = A_SpatialPolygons2,
                                                                                                                                                               nsub = BarycentricManyPoints.nsub)
                                                  
                                                  
                                                  ret.list <- list(OmegaMesh = OmegaMesh,
                                                                   nn.mesh = nn.mesh,
                                                                   MonteCarloIntegrator.mapper = MonteCarloIntegrator.mapper,
                                                                   CoxLikVoronoi.OmegaVerticesWeight = CoxLikVoronoi.OmegaVerticesWeight,
                                                                   CoxLikDualMesh.OmegaVerticesWeight = CoxLikDualMesh.OmegaVerticesWeight,
                                                                   CoxLikDualMeshExtra.OmegaVerticesWeightList = CoxLikDualMeshExtra.OmegaVerticesWeightList,
                                                                   CoxLikBarycentricPointSpread.OmegaVerticesWeight = CoxLikBarycentricPointSpread.OmegaVerticesWeight,
                                                                   CoxLikBarycentricPointSpreadManyPoints.OmegaVerticesWeight = CoxLikBarycentricPointSpreadManyPoints.OmegaVerticesWeight,
                                                                   basisfunc.grid.mat = basisfunc.grid.mat,
                                                                   inprodmat = inprodmat,
                                                                   OmegaArea = OmegaArea,
                                                                   affected.mesh.points = affected.mesh.points)
                                                  return(ret.list)
                                                  
                                                }
    print(i)
    prep.mesh.list[[i]] <- prep.mesh.single.A.list
    
  }
  
  # Generate clock-random string
  set.seed(round(proc.time()[3]*1000))
  rand_string <- stringi::stri_rand_strings(1,length=6)
  folder_string = paste0(format(Sys.time(), "%Y%m%d_%H%M"),"_",rand_string)
  
  dir.create(path = paste0("../Results/temp/",folder_string),recursive = T)
  dir.create(path = paste0("../Results/",folder_string),recursive = T)
  
  save(input.list,file = paste0("../Results/",folder_string,"/input.list_run.RData"))
  save(prep.mesh.list,file = paste0("../Results/",folder_string,"/prep.mesh.list_run.RData"))
} else { # I.e. if folder_string is provided
  input.list0 <- input.list
  
  load(file = paste0("../Results/",folder_string,"/input.list_run.RData"))
  load(file = paste0("../Results/",folder_string,"/prep.mesh.list_run.RData"))
  
  
  # Don't care about differences in these.k.z.samps
  input.list$these.k.z.samps = NULL
  input.list0$these.k.z.samps = NULL
  
  if(!isTRUE(all.equal(input.list,input.list0))){
    print(all.equal(input.list,input.list0))
    stop("folder_string does contain same input as provided")
  }
}




doMC::registerDoMC(parallelize.numCores)

#uu=proc.time()
res.dt.list <- list()
#i = j = k= l = u = 1

for (l in 1:length(boundary.splitter.vec)){
  res.dt.list[[l]] <- list()
  boundary.splitter <- boundary.splitter.vec[l]
  Likelihood.OracleTruth.prepare.list <- prep.A.list[[l]]$Likelihood.OracleTruth.prepare.list
  MonteCarloIntegrator.prepare.list.0 <- prep.A.list[[l]]$MonteCarloIntegrator.prepare.list.0
  OmegaArea <- prep.A.list[[l]]$OmegaArea
  
  for (i in 1:length(INLA.max.edge.vec)){
    res.dt.list[[l]][[i]] <- list()
    INLA.max.edge <- INLA.max.edge.vec[i]
    
    OmegaMesh <- prep.mesh.list[[i]][[l]]$OmegaMesh
    nn.mesh <- prep.mesh.list[[i]][[l]]$nn.mesh
    MonteCarloIntegrator.mapper <- prep.mesh.list[[i]][[l]]$MonteCarloIntegrator.mapper
    CoxLikVoronoi.OmegaVerticesWeight <- prep.mesh.list[[i]][[l]]$CoxLikVoronoi.OmegaVerticesWeight
    CoxLikDualMesh.OmegaVerticesWeight <- prep.mesh.list[[i]][[l]]$CoxLikDualMesh.OmegaVerticesWeight
    CoxLikDualMeshExtra.OmegaVerticesWeightList <- prep.mesh.list[[i]][[l]]$CoxLikDualMeshExtra.OmegaVerticesWeightList
    CoxLikBarycentricPointSpread.OmegaVerticesWeight <- prep.mesh.list[[i]][[l]]$CoxLikBarycentricPointSpread.OmegaVerticesWeight
    CoxLikBarycentricPointSpreadManyPoints.OmegaVerticesWeight <- prep.mesh.list[[i]][[l]]$CoxLikBarycentricPointSpreadManyPoints.OmegaVerticesWeight
    basisfunc.grid.mat = prep.mesh.list[[i]][[l]]$basisfunc.grid.mat
    inprodmat = prep.mesh.list[[i]][[l]]$inprodmat
    affected.mesh.points = prep.mesh.list[[i]][[l]]$affected.mesh.points
    
    
    
    for (j in 1:length(true.field.variance.vec)){
      res.dt.list[[l]][[i]][[j]] <- list()
      sigma2x.true <- true.field.variance.vec[j]
      
      for(u in 1:length(kappa.true.vec)){
        kappa.true = kappa.true.vec[u]
        RandomFields::RFoptions(storing=TRUE)
        
        res.mat <- foreach::foreach(k=these.k.z.samps,
                                    .verbose=FALSE,
                                    .inorder=FALSE,
                                    .combine = 'rbind') %dopar% {
                                      
                                      RandomFields::RFoptions(storing=TRUE)
                                      
                                      
                                      #### Sampling the true field ####
                                      this.seed = seed.z.sampling.initial+k
                                      
                                      set.seed(this.seed)
                                      spatstat.options(npixel=npixel)
                                      
                                      lg.s <- spatstat::rLGCP('matern', beta0.true,
                                                              var=sigma2x.true, scale=1/kappa.true, nu=1, win=win)
                                      Lam <- attr(lg.s, 'Lambda')
                                      rf.s <- log(Lam$v)
                                      
                                      ### True int ####
                                      OracleTruth.intval <- OracleTruthIntegral.evaluator.new(Likelihood.OracleTruth.prepare = Likelihood.OracleTruth.prepare.list,
                                                                                              z = rf.s)
                                      
                                      
                                      #### Differnet approximation methods for the integral #####
                                      
                                      if(!FEM.approx.mesh){
                                        rf.s.mesh <- rf.s[nn.mesh$which]
                                      } else {
                                        # This is actually quite fast and better than pre-compuing what needs to be multiplied with rf.s.vec
                                        rf.s.vec = as.vector(rf.s) ### Yes, this is how it has to be, no transformation here. Might do it before rf.s = t(log(Lam$V)), but does not matter as xval=yval.
                                        inprodvec = basisfunc.grid.mat%*%rf.s.vec*area.pixel
                                        rf.s.mesh.affected = as.vector(solve(inprodmat,inprodvec))
                                        rf.s.mesh <- rep(0,OmegaMesh$n)
                                        rf.s.mesh[affected.mesh.points] <- rf.s.mesh.affected
                                        
                                      }
                                      
                                      if (do_debug){
                                        rf.s.mesh.old = rf.s[nn.mesh$which]
                                        rf.s.vec = as.vector(rf.s) ### Yes, this is how it has to be, no transformation here. Might do it before rf.s = t(log(Lam$V)), but does not matter as xval=yval.
                                        inprodvec = basisfunc.grid.mat%*%rf.s.vec*area.pixel
                                        rf.s.mesh.affected = as.vector(solve(inprodmat,inprodvec))
                                        rf.s.mesh.new <- rep(0,OmegaMesh$n)
                                        rf.s.mesh.new[affected.mesh.points] <- rf.s.mesh.affected
                                        
                                        basisfunc.grid.mat.full <- t(INLA::inla.spde.make.A(mesh = OmegaMesh,loc=as.matrix(rf.xy.grid)))
                                        
                                        
                                        funcval.new = t(basisfunc.grid.mat)%*%rf.s.mesh.affected
                                        #funcval.new0 = t(basisfunc.grid.mat.full)%*%rf.s.mesh.new # same as above
                                        funcval.new.mat <- matrix(funcval.new,ncol=length(xyvals))
                                        
                                        funcval.old = t(basisfunc.grid.mat.full)%*%rf.s.mesh.old
                                        funcval.old.mat <- matrix(funcval.old,ncol=length(xyvals))
                                        
                                        range.plot = range(c(funcval.new.mat,funcval.old.mat,rf.s.vec),na.rm=T)
                                        breaks.plot <- seq(range.plot[1],range.plot[2],length.out = 50)
                                        
                                        
                                        dev.off()
                                        library(fields)
                                        par(mfrow=c(2,2))
                                        image.plot(xyvals,xyvals,funcval.old.mat,breaks = breaks.plot,nlevel = length(breaks.plot)-1)
                                        #plot(OmegaMesh,add=T)
                                        plot(A_SpatialLines,add=T,col="purple",lwd=3)
                                        image.plot(xyvals,xyvals,funcval.new.mat,breaks = breaks.plot,nlevel = length(breaks.plot)-1)
                                        # plot(OmegaMesh,add=T)
                                        plot(A_SpatialLines,add=T,col="purple",lwd=3)
                                        
                                        image.plot(xyvals,xyvals,rf.s,breaks = breaks.plot,nlevel = length(breaks.plot)-1)
                                        #plot(OmegaMesh,add=T)
                                        plot(A_SpatialLines,add=T,col="purple",lwd=3)
                                        
                                        image.plot(xyvals,xyvals,rf.s,breaks = breaks.plot,nlevel = length(breaks.plot)-1)
                                        plot(OmegaMesh,add=T)
                                        plot(A_SpatialLines,add=T,col="purple",lwd=3)
                                        
                                        mean((funcval.old.mat-rf.s)^2)
                                        mean((funcval.new.mat-rf.s)^2)
                                        
                                        mean((exp(funcval.old.mat)-exp(rf.s))^2)
                                        mean((exp(funcval.new.mat)-exp(rf.s))^2)
                                        
                                      }
                                      
                                      ### Mesh MC integrator (used in place of analytical expression for speed) ###
                                      
                                      MeshExact.intval <- MonteCarloIntegrator.evaluator.new(MonteCarloIntegrator.mapper = MonteCarloIntegrator.mapper,
                                                                                             MonteCarloIntegrator.prepare.list.0 = MonteCarloIntegrator.prepare.list.0,
                                                                                             z = rf.s.mesh)
                                      
                                      ### Voronoi approximation ####
                                      
                                      Voronoi.intval <- ApproxIntegral.evaluator(OmegaVerticesWeight = CoxLikVoronoi.OmegaVerticesWeight,
                                                                                 z = rf.s.mesh)
                                      
                                      ### Dual mesh SPDE book approximation ####
                                      
                                      DualMesh.intval <- ApproxIntegral.evaluator(OmegaVerticesWeight = CoxLikDualMesh.OmegaVerticesWeight,
                                                                                  z = rf.s.mesh)
                                      
                                      ### Dual Mesh approach with extra integration points from old inlabru package ####
                                      
                                      DualMeshExtra.intval <- ApproxIntegral.mapping.evaluator(intPointsWeight = CoxLikDualMeshExtra.OmegaVerticesWeightList$intPointsWeight,
                                                                                               intPointsMapping = CoxLikDualMeshExtra.OmegaVerticesWeightList$intPointsAMat,
                                                                                               z = rf.s.mesh)
                                      
                                      ### Dual Mesh approach with extra integration points mapped back to original mesh locations ####
                                      
                                      DualMeshExtraMeshMapped.intval <- ApproxIntegral.evaluator(OmegaVerticesWeight = CoxLikDualMeshExtra.OmegaVerticesWeightList$allmeshpointsWeights,
                                                                                                 z = rf.s.mesh)
                                      
                                      ### Barycentric points spreading method implemented in inlabru, mapping back to the original mesh locations ####
                                      
                                      BarycentricPointSpread.intval <- ApproxIntegral.evaluator(OmegaVerticesWeight = CoxLikBarycentricPointSpread.OmegaVerticesWeight,
                                                                                                z = rf.s.mesh)
                                      
                                      ### Barycentric points spreading method implemented in inlabru, mapping back to the original mesh locations, with increased points being spread ####
                                      
                                      BarycentricPointSpreadManyPoints.intval <- ApproxIntegral.evaluator(OmegaVerticesWeight = CoxLikBarycentricPointSpreadManyPoints.OmegaVerticesWeight,
                                                                                                          z = rf.s.mesh)
                                      
                                      ### Returning all results ####
                                      
                                      
                                      res.vec <- c(MeshExact.intval,Voronoi.intval,DualMesh.intval,DualMeshExtra.intval,DualMeshExtraMeshMapped.intval,BarycentricPointSpread.intval,BarycentricPointSpreadManyPoints.intval,
                                                   OracleTruth.intval,INLA.max.edge,sigma2x.true,boundary.splitter,beta0.true,kappa.true,OmegaArea,this.seed)
                                      
                                      return(res.vec)
                                    }
        
        
        colnames(res.mat) = c("MeshExact","Voronoi","DualMesh","DualMeshExtra","DualMeshExtraMeshMapped","BarycentricPointSpread","BarycentricPointSpreadManyPoints",
                              "OracleTruth","INLA.max.edge","sigma2x.true","boundary.splitter","beta0.true","kappa.true","OmegaArea","seed")
        
        res.dt.tmp <- data.table::as.data.table(res.mat)
        
        res.dt.tmp[,folder_string:=folder_string]
        res.dt.tmp[,FEM.approx.mesh:=FEM.approx.mesh]
        
        data.table::fwrite(res.dt.tmp,paste0("../Results/temp/",folder_string,"/res.dt.tmp.i_",i,"__j_",j,"__l_",l,"__u_",u,"_.csv"),append=T)
        
        res.dt.list[[l]][[i]][[j]][[u]] <- res.dt.tmp
        
      }
    }
  }
}

res.dt <- data.table::rbindlist(unlist(unlist(unlist(res.dt.list,recursive = F),recursive = F),recursive = F))

methods <- c("MeshExact","Voronoi","DualMesh","DualMeshExtra","DualMeshExtraMeshMapped","BarycentricPointSpread","BarycentricPointSpreadManyPoints")
bycols <- c("INLA.max.edge","sigma2x.true","boundary.splitter","beta0.true","kappa.true","folder_string","FEM.approx.mesh","OmegaArea")

res.dt <- res.dt[,lapply(.SD, function(x) x/OmegaArea),.SDcols=c(methods,"OracleTruth"),by=c(bycols,"seed")] # Standardizes the estimates to be for every square unit to easier compare the performance for different boundary.splitters


fwrite(res.dt,paste0("../Results/",folder_string,"/full.res.dt.csv"),append=T)
res.dt <- fread(paste0("../Results/",folder_string,"/full.res.dt.csv"))

summary.res.dt <- res.dt[,lapply(.SD, function(x){sqrt(mean((x-OracleTruth)^2))}),.SDcols=methods,by=bycols]

fwrite(summary.res.dt,paste0("../Results/",folder_string,"/res.dt.csv"))

fwrite(res.dt,paste0("../Results/common.full.res.dt.csv"),append = T)


