#### Estimating the integral of the exponential latent field

rm(list=ls())

library(Matrix)
library(INLA)
library(spatstat)
library(fields)

# Help function

prepare.basisfunc.grid.mat.new <- function(OmegaMesh,rf.xy.grid){
  basisfunc.grid.mat <- Matrix::t(INLA::inla.spde.make.A(mesh = OmegaMesh,loc=as.matrix(rf.xy.grid)))
  affected.mesh.points <- rowSums(basisfunc.grid.mat>0)>1 # Removing basis functions which affects 0 or only 1 grid point 
  # (1 to handle situations where two basis functions handle only the same pixel and thus cannot be differed)
  
  basisfunc.grid.mat <- basisfunc.grid.mat[affected.mesh.points,]  # Removes unaffected basis functions
  
  return(list(mat = basisfunc.grid.mat,affected.mesh.points = affected.mesh.points))
}


# # Settings 
seed.z.sampling <- 123
xy.range <- c(0,16)
xy.range.approx.mesh <- c(-2,18)
npixel <- 20*16 # Used for integration
npixels_field_sampling_inv_multiplier <- 2 # npixel times the inverse of this denotes the refinement on the sampled 
                                          # true field (need to be an integer number which npixel can be devided by

INLA.max.edge <- 4
true.field.variance <- 0.5
beta0.true <- 1
kappa.true <- sqrt(8)/c(4) # kappa = sqrt(8)/range
INLA.offset <- c(-0.03,-0.1) 


### Prepping 
boundary <- rbind(rep(xy.range[1],2),
                  xy.range,
                  rep(xy.range[2],2),
                  rev(xy.range),
                  rep(xy.range[1],2),deparse.level = 3)

xyvals <- (seq(xy.range.approx.mesh[1],xy.range.approx.mesh[2],length.out=npixel+1)[-(npixel+1)] + seq(xy.range.approx.mesh[1],xy.range.approx.mesh[2],length.out=npixel+1)[-1])/2
rf.xy.grid <- expand.grid(x=xyvals, y=xyvals)
area.pixel = diff(xyvals)[1]^2

win <- spatstat::owin(xy.range,xy.range)


### Create the mesh ####

OmegaMesh <- INLA::inla.mesh.2d(loc.domain=boundary,
                                max.edge=INLA.max.edge, 
                                cutoff=0.3,
                                offset=INLA.offset) # Should try to see the effect of settign this larger


# Prepare mesh/grid for approximation with of a field

basisfunc.grid.mat.0 = prepare.basisfunc.grid.mat.new(OmegaMesh=OmegaMesh,rf.xy.grid = rf.xy.grid)
basisfunc.grid.mat <- basisfunc.grid.mat.0$mat
affected.mesh.points <- basisfunc.grid.mat.0$affected.mesh.points

inprodmat = basisfunc.grid.mat%*%t(basisfunc.grid.mat)*area.pixel # Could have used inla.mesh.fem(OmegaMesh)$c1, but tests gave slightly better results using my Monte Carlo version

#### Samling a true field to approximate ###

set.seed(seed.z.sampling)
npixel_sampling <- npixel/npixels_field_sampling_inv_multiplier
spatstat.options(npixel=npixel_sampling)
lg.s <- spatstat::rLGCP('matern', beta0.true,
                        var=true.field.variance, scale=1/kappa.true, nu=1, win=win)
Lam <- attr(lg.s, 'Lambda')

if (npixels_field_sampling_inv_multiplier==1){
  rf.s.true.mat <- log(Lam$v)
} else { # This is just a tedious hack to allow the MC integration be more refiend than the sampled field
  rf.s.true.mat <- matrix(NA,ncol = npixel,nrow = npixel)
  for (i in 1:npixel_sampling){
    ind <- 1:npixels_field_sampling_inv_multiplier+(i-1)*npixels_field_sampling_inv_multiplier
    for(j in ind){
      rf.s.true.mat[,j] <- rep(log(Lam$v)[,i],each=npixels_field_sampling_inv_multiplier)
    }
  }
}



#### Approximating the true field ####

rf.s.true.vec = as.vector(rf.s.true.mat) 
inprodvec = basisfunc.grid.mat%*%rf.s.true.vec*area.pixel

rf.s.mesh.affected = as.vector(solve(inprodmat,inprodvec)) # Solving the least squares problem
rf.s.mesh.approx <- rep(0,OmegaMesh$n)
rf.s.mesh.approx[affected.mesh.points] <- rf.s.mesh.affected

#### Image plot of true field and approximation ####
range.plot = range(c(rf.s.true.vec,rf.s.mesh.approx),na.rm=T)
breaks.plot <- seq(range.plot[1],range.plot[2],length.out = 50)


basisfunc.grid.mat.full <- INLA::inla.spde.make.A(mesh = OmegaMesh,loc=as.matrix(rf.xy.grid))
rf.s.mesh.approx.mat = matrix(basisfunc.grid.mat.full%*%rf.s.mesh.approx,ncol=length(xyvals))

par(mfrow=c(2,2))


fields::image.plot(xyvals,xyvals,rf.s.true.mat,breaks = breaks.plot,nlevel = length(breaks.plot)-1)
lines(boundary[,1],boundary[,2],col="purple",lwd=3)
title("True field")

fields::image.plot(xyvals,xyvals,rf.s.mesh.approx.mat,breaks = breaks.plot,nlevel = length(breaks.plot)-1)
lines(boundary[,1],boundary[,2],col="purple",lwd=3)
title("Approx field")

fields::image.plot(xyvals,xyvals,rf.s.mesh.approx.mat,breaks = breaks.plot,nlevel = length(breaks.plot)-1)
plot(OmegaMesh,add=T)
lines(boundary[,1],boundary[,2],col="purple",lwd=3)
title("Approx field with mesh")


fields::image.plot(xyvals,xyvals,rf.s.true.mat-rf.s.mesh.approx.mat)
plot(OmegaMesh,add=T)
lines(boundary[,1],boundary[,2],col="purple",lwd=3)
title("True - approx field with mesh")


# The minimum mean squared error reached
MinMSE <- mean((rf.s.mesh.approx.mat-rf.s.true.mat)^2)

# Adjusting very basis function a bit up and down to see that the reached mse is indeed the minimum
eps <- 0.01

MSE_plus_vec <- rep(0,OmegaMesh$n)
MSE_minus_vec <- rep(0,OmegaMesh$n)

for(i in 1:OmegaMesh$n){
  rf.s.mesh.approx_plus_adjusted <- rf.s.mesh.approx
  rf.s.mesh.approx_minus_adjusted <- rf.s.mesh.approx
  
  rf.s.mesh.approx_plus_adjusted[i] <-   rf.s.mesh.approx_plus_adjusted[i] + eps
  rf.s.mesh.approx_minus_adjusted[i] <-   rf.s.mesh.approx_minus_adjusted[i] - eps
  
  rf.s.mesh.approx.mat_plus_adjusted = as.vector(basisfunc.grid.mat.full%*%rf.s.mesh.approx_plus_adjusted)
  rf.s.mesh.approx.mat_minus_adjusted = as.vector(basisfunc.grid.mat.full%*%rf.s.mesh.approx_minus_adjusted)
  
  MSE_plus_vec[i] <- mean((rf.s.mesh.approx.mat_plus_adjusted-rf.s.true.vec)^2)
  MSE_minus_vec[i] <- mean((rf.s.mesh.approx.mat_minus_adjusted-rf.s.true.vec)^2)
  
}

# Any negative results here indicates our solution is not a local minima
# Note: Very small negatives may be due to Monte Carlo integration uncertainty
min((MSE_plus_vec-MinMSE))
min((MSE_minus_vec-MinMSE))










#### 1D version using a section to approximate a section of the 2D field above, using the same code ####

OmegaMesh_1D <- inla.mesh.1d(loc=seq(xy.range.approx.mesh[1],xy.range.approx.mesh[2],by=INLA.max.edge))
rf.xy.grid_1D <- xyvals

# Prepare mesh/grid for approximation with of a field

basisfunc.grid.mat.0_1D = prepare.basisfunc.grid.mat.new(OmegaMesh=OmegaMesh_1D,rf.xy.grid = rf.xy.grid_1D)
basisfunc.grid.mat_1D <- basisfunc.grid.mat.0_1D$mat
affected.mesh.points_1D <- basisfunc.grid.mat.0_1D$affected.mesh.points

inprodmat_1D = basisfunc.grid.mat_1D%*%t(basisfunc.grid.mat_1D)*sqrt(area.pixel) # Could have used inla.mesh.fem(OmegaMesh)$c1, but tests gave slightly better results using my Monte Carlo version

#### Samling a true field to approximate ###
this_line <- 10
rf.s.true.vec_1D <- rf.s.true.mat[this_line,]

#### Approximating the true field ####

inprodvec_1D = basisfunc.grid.mat_1D%*%rf.s.true.vec_1D*sqrt(area.pixel)

rf.s.mesh.affected_1D = as.vector(solve(inprodmat_1D,inprodvec_1D)) # Solving the least squares problem
rf.s.mesh.approx_1D <- rep(0,OmegaMesh_1D$n)
rf.s.mesh.approx_1D[affected.mesh.points_1D] <- rf.s.mesh.affected_1D


#### Plot of true field and approximation ####
range.plot_1D = range(c(rf.s.true.vec_1D,rf.s.mesh.approx_1D),na.rm=T)
breaks.plot_1D <- seq(range.plot_1D[1],range.plot_1D[2],length.out = 50)


basisfunc.grid.mat.full_1D <- INLA::inla.spde.make.A(mesh = OmegaMesh_1D,loc=as.matrix(rf.xy.grid_1D))
rf.s.mesh.approx.mat_1D = as.vector(basisfunc.grid.mat.full_1D%*%rf.s.mesh.approx_1D)

dev.off()
plot(xyvals,rf.s.true.vec_1D,type="l")
abline(v=OmegaMesh_1D$loc,col="grey")
lines(xyvals,rf.s.mesh.approx.mat_1D,col=2)
title("True and mesh approximate field")

# The minimum mean squared error reached
MinMSE <- mean((rf.s.mesh.approx.mat_1D-rf.s.true.vec_1D)^2)

# Adjusting very basis function a bit up and down to see that the reached mse is indeed the minimum
eps <- 0.01

MSE_plus_vec <- rep(0,OmegaMesh_1D$n)
MSE_minus_vec <- rep(0,OmegaMesh_1D$n)

for(i in 1:OmegaMesh_1D$n){
  rf.s.mesh.approx_1D_plus_adjusted <- rf.s.mesh.approx_1D
  rf.s.mesh.approx_1D_minus_adjusted <- rf.s.mesh.approx_1D

  rf.s.mesh.approx_1D_plus_adjusted[i] <-   rf.s.mesh.approx_1D_plus_adjusted[i] + eps
  rf.s.mesh.approx_1D_minus_adjusted[i] <-   rf.s.mesh.approx_1D_minus_adjusted[i] - eps
  
  rf.s.mesh.approx.mat_1D_plus_adjusted = as.vector(basisfunc.grid.mat.full_1D%*%rf.s.mesh.approx_1D_plus_adjusted)
  rf.s.mesh.approx.mat_1D_minus_adjusted = as.vector(basisfunc.grid.mat.full_1D%*%rf.s.mesh.approx_1D_minus_adjusted)
  
  MSE_plus_vec[i] <- mean((rf.s.mesh.approx.mat_1D_plus_adjusted-rf.s.true.vec_1D)^2)
  MSE_minus_vec[i] <- mean((rf.s.mesh.approx.mat_1D_minus_adjusted-rf.s.true.vec_1D)^2)
  
  }


# Any negative results here indicates it is not a local minima
# Note: Very small negatives may be due to Monte Carlo integration uncertainty
min((MSE_plus_vec-MinMSE))
min((MSE_minus_vec-MinMSE))
