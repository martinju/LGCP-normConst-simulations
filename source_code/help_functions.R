

# polygonList = list()
# for (i in 1:nrow(OmegaMesh$graph$tv)){
#   triangle <- OmegaMesh$graph$tv[i,]
#   triangle.locs <- OmegaMesh$loc[triangle,1:2]
#   polygonList[[i]] = sp::Polygon(coords = triangle.locs[c(1:3,1),])
# }
# OmegaTriangels <- sp::SpatialPolygons(list(sp::Polygons(polygonList,ID=paste0("_"))))
# 
# rf.xy.grid.SpatialPoints = sp::SpatialPoints(rf.xy.grid)
# 
# gInterse

prepare.basisfunc.grid.mat.new <- function(OmegaMesh,rf.xy.grid){
  basisfunc.grid.mat <- t(INLA::inla.spde.make.A(mesh = OmegaMesh,loc=as.matrix(rf.xy.grid)))
  affected.mesh.points <- rowSums(basisfunc.grid.mat>0)>1 # Removing basis functions which affects 0 or only 1 grid point 
  # (1 to handle situations where two basis functions handle only the same pixel and thus cannot be differed)
  
  basisfunc.grid.mat <- basisfunc.grid.mat[affected.mesh.points,]  # Removes unaffected basis functions
  
  return(list(mat = basisfunc.grid.mat,affected.mesh.points = affected.mesh.points))
}



prepare.basisfunc.grid.mat = function(OmegaMesh,rf.xy.grid){
  
  basisfunc.index <- NULL
  basisfunc.vals <- NULL
  
  for (i in 1:nrow(OmegaMesh$graph$tv)){
    
    triangle <- OmegaMesh$graph$tv[i,]
    triangle.locs <- OmegaMesh$loc[triangle,1:2]
    
    triangle.points <- which(sp::point.in.polygon(point.x = rf.xy.grid[,1],point.y = rf.xy.grid[,2],
                                                  pol.x = triangle.locs[,1],pol.y = triangle.locs[,2])>0)
    
    
    X <- cbind(1,triangle.locs)
    Y1 <- c(1,0,0)
    Y2 <- c(0,1,0)
    Y3 <- c(0,0,1)
    
    beta1 <- solve(X,Y1)
    beta2 <- solve(X,Y2)
    beta3 <- solve(X,Y3)
    
    Xmap <- as.matrix(cbind(1,rf.xy.grid[triangle.points,]))
    
    map1 <- Xmap%*%beta1
    map2 <- Xmap%*%beta2
    map3 <- Xmap%*%beta3
    
    basisfunc.index <- rbind(basisfunc.index,
                             cbind(triangle[1],triangle.points),
                             cbind(triangle[2],triangle.points),
                             cbind(triangle[3],triangle.points))
    basisfunc.vals <- c(basisfunc.vals,
                        map1,
                        map2,
                        map3)
    
  }
  these.duplicated <- which(duplicated(basisfunc.index))
  if (length(these.duplicated)>0){
    basisfunc.index = basisfunc.index[-these.duplicated,]
    basisfunc.vals = basisfunc.vals[-these.duplicated]
    
  }
  
  basisfunc.grid.mat = Matrix::sparseMatrix(i =basisfunc.index[,1], 
                                            j=basisfunc.index[,2],
                                            x=basisfunc.vals, 
                                            dims=c(OmegaMesh$n,nrow(rf.xy.grid)))
  
  return(basisfunc.grid.mat)
}




MonteCarloIntegrator.prepare.0 <- function(boundary,boundary.splitter,no.MC.intpoints,A_SpatialPolygons){
  
  boundary.range <- range(boundary[,1])
  no.points.adjusted <- round(no.MC.intpoints/(diff(boundary.range)/boundary.splitter))*(diff(boundary.range)/boundary.splitter)
  
  A.xval <- seq(boundary.range[1],boundary.range[2],length.out=no.points.adjusted+1)
  A.xval <- A.xval[-1]-diff(A.xval)/2
  A.yval <- seq(boundary.range[1],boundary.range[2],length.out=no.points.adjusted+1)
  A.yval <- A.yval[-1]-diff(A.yval)/2
  
  area.rect <- diff(A.xval)[1]*diff(A.yval)[1]
  
  A.xyval <- expand.grid(A.xval,A.yval)
  
  these.Axyval.in.A <- which(point.in.SpatialPolygons(A.xyval[,1], A.xyval[,2], A_SpatialPolygons))
  
  ret <- list(A.xyval.in.A = as.matrix(A.xyval[these.Axyval.in.A,]),
              area.rect = area.rect)
  
  return(ret)
}

boundary.splitter.func <- function(boundary,boundary.splitter){
  
  sub.boundary <- boundary/boundary.splitter
  max.xy <- max(sub.boundary)
  
  boundary.list <- list()
  for(i in 1:boundary.splitter){
    boundary.list[[i]] <- list()
    for(j in 1:boundary.splitter){
      boundary.list[[i]][[j]] <- sub.boundary
    }
  }
  for(i in 1:boundary.splitter){
    for(j in 1:boundary.splitter){
      boundary.list[[i]][[j]][,1] <- boundary.list[[i]][[j]][,1] + (j-1)*max.xy
      boundary.list[[i]][[j]][,2] <- boundary.list[[i]][[j]][,2] + (i-1)*max.xy
    }
  }
  
  ij.include.vec <- NULL
  for(i in 1:boundary.splitter){
    for(j in 1:boundary.splitter){
      if ((i/2) == round(i/2)){
        if((j/2) == round(j/2)){
          ij.include.vec <- rbind(ij.include.vec,c(i,j))
        }
      }
      if ((i/2) != round(i/2)){
        if((j/2) != round(j/2)){
          ij.include.vec <- rbind(ij.include.vec,c(i,j))
        }
      }
    }
  }
  
  k <- 0
  boundary.list.new <- list()
  for(i in 1:boundary.splitter){
    for(j in 1:boundary.splitter){
      if (any((i==ij.include.vec[,1] & j==ij.include.vec[,2]))){
        k <- k+1
        boundary.list.new[[k]] <- boundary.list[[i]][[j]]
      }
    }
  }
  
  
  
  
  
  Polygon.boundary.list <- list()
  for(i in 1:length(boundary.list.new)){
    Polygon.boundary.list[[i]] <- sp::Polygon(boundary.list.new[[i]])
  }
  
  boundary_polys <- sp::Polygons(Polygon.boundary.list, '0')
  boundary_SpatialPolygons <- sp::SpatialPolygons(list(boundary_polys))
  
  return(boundary_SpatialPolygons)
}

boundary.splitter.func2 <- function(boundary,boundary.splitter){
  
  sub.boundary <- boundary/boundary.splitter
  max.xy <- max(sub.boundary)
  
  boundary.list <- list()
  for(i in 1:boundary.splitter){
    boundary.list[[i]] <- list()
    for(j in 1:boundary.splitter){
      boundary.list[[i]][[j]] <- sub.boundary
    }
  }
  for(i in 1:boundary.splitter){
    for(j in 1:boundary.splitter){
      boundary.list[[i]][[j]][,1] <- boundary.list[[i]][[j]][,1] + (j-1)*max.xy
      boundary.list[[i]][[j]][,2] <- boundary.list[[i]][[j]][,2] + (i-1)*max.xy
    }
  }
  
  ij.include.vec <- NULL
  for(i in 1:boundary.splitter){
    for(j in 1:boundary.splitter){
      if ((i/2) == round(i/2)){
        if((j/2) == round(j/2)){
          ij.include.vec <- rbind(ij.include.vec,c(i,j))
        }
      }
      if ((i/2) != round(i/2)){
        if((j/2) != round(j/2)){
          ij.include.vec <- rbind(ij.include.vec,c(i,j))
        }
      }
    }
  }
  
  k <- 0
  boundary.list.new <- list()
  for(i in 1:boundary.splitter){
    for(j in 1:boundary.splitter){
      if (any((i==ij.include.vec[,1] & j==ij.include.vec[,2]))){
        k <- k+1
        boundary.list.new[[k]] <- boundary.list[[i]][[j]]
      }
    }
  }
  
  
  
  Polygons.boundary.list <- list()
  for(i in 1:length(boundary.list.new)){
    Polygons.boundary.list[[i]] <- sp::Polygons(list(sp::Polygon(boundary.list.new[[i]])),ID = paste0(i))
  }
  
  boundary_SpatialPolygons2 <- sp::SpatialPolygons(Polygons.boundary.list)
  
  return(boundary_SpatialPolygons2)
}


Likelihood.OracleTruth.prepare.new <- function(A_SpatialPolygons,rf.xy.grid){
  
  these.rf.xy.grid.in.A <- which(prevR::point.in.SpatialPolygons(rf.xy.grid$x, rf.xy.grid$y, A_SpatialPolygons))
  
  area.rect <- diff(rf.xy.grid$x)[1]^2 # Assuming equal xy-distance
  
  ret <- list(these.rf.xy.grid.in.A = these.rf.xy.grid.in.A,
              area.rect = area.rect)
}


MonteCarloIntegrator.prepare.1 <- function(A.xyval.in.A,OmegaMesh){
  
  inlaA.integrator <- INLA::inla.spde.make.A(mesh = OmegaMesh,loc = A.xyval.in.A)
  
  return(inlaA.integrator)
}

CoxLikVoronoi.prepare <- function(OmegaMesh,A_SpatialPolygons){
  
  # Extracting from the OmegaMesh object 
  OmegaVertices <- OmegaMesh$loc # The coordinates of the Omega vertices
  
  OmegaVerticesWeight <- CreateVoronoiTessellation_new(locationsCoordX=OmegaVertices[,1],
                                                       locationsCoordY=OmegaVertices[,2],
                                                       A_SpatialPolygons=A_SpatialPolygons,
                                                       parallelize.numCores=1)$tileSize
  #intPointsWeight <- VoroniTesselationArea(points=OmegaVertices[,1:2],
  #                                         A = A)
  
  return(OmegaVerticesWeight)
}

CreateVoronoiTessellation_new <- function(locationsCoordX,
                                          locationsCoordY,
                                          A_SpatialPolygons,
                                          parallelize.numCores){
  ## Builing Voronoi triangulated polygons for the complete mesh
  dd <- deldir::deldir(locationsCoordX, locationsCoordY)
  tiles <- deldir::tile.list(dd)  # These tiles defines all the polygons which we are going to use in the modeling
  
  
  
  # Function for checking the area covered by the new polygons within the region with likelihood contribution
  PFunc <- function(p,A_SpatialPolygons) {
    pl <- coo2sp(cbind(p$x, p$y))
    if (rgeos::gIntersects(pl, A_SpatialPolygons)){
      return(rgeos::gArea(rgeos::gIntersection(pl, A_SpatialPolygons)))
    }
    else {
      return(0)
    }
  }
  tileSize <- parallel::mclapply(X = tiles, FUN = PFunc, A_SpatialPolygons = A_SpatialPolygons, mc.cores = parallelize.numCores)
  
  tileSize <- unlist(tileSize)
  
  retList <- list()
  retList$tiles <- tiles
  retList$tileSize <- tileSize
  return(retList)
}

coo2sp <- function(coo) {
  n <- nrow(coo)
  if (any(coo[1,]!=coo[n,]))
    coo <- coo[c(1:n,1),]
  sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coo)), '0')))
}

CoxLikDualMesh.prepare <- function(OmegaMesh,A_SpatialPolygons){
  
  # From book.mesh.dual copied from the spde.inla.book software
  dmesh <- book.mesh.dual(OmegaMesh)
  
  OmegaVerticesWeight <- sapply(1:length(dmesh), function(i) {
    if (rgeos::gIntersects(dmesh[i, ], A_SpatialPolygons))
      return(rgeos::gArea(rgeos::gIntersection(dmesh[i, ], A_SpatialPolygons)))
    else return(0)
  })
  
  return(OmegaVerticesWeight)
}

# Copy of function from R/spde-book-functions.R in the spde-book-files.zip (dated 9 Nov 2018) found here http://www.r-inla.org/spde-book
book.mesh.dual <- function(mesh) { 
  if (mesh$manifold=='R2') {
    ce <- t(sapply(1:nrow(mesh$graph$tv), function(i)
      colMeans(mesh$loc[mesh$graph$tv[i, ], 1:2])))
    library(parallel)
    pls <- mclapply(1:mesh$n, function(i) {
      p <- unique(Reduce('rbind', lapply(1:3, function(k) {
        j <- which(mesh$graph$tv[,k]==i)
        if (length(j)>0) 
          return(rbind(ce[j, , drop=FALSE],
                       cbind(mesh$loc[mesh$graph$tv[j, k], 1] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 1], 
                             mesh$loc[mesh$graph$tv[j, k], 2] +
                               mesh$loc[mesh$graph$tv[j, c(2:3,1)[k]], 2])/2))
        else return(ce[j, , drop=FALSE])
      })))
      j1 <- which(mesh$segm$bnd$idx[,1]==i)
      j2 <- which(mesh$segm$bnd$idx[,2]==i)
      if ((length(j1)>0) | (length(j2)>0)) {
        p <- unique(rbind(mesh$loc[i, 1:2], p,
                          mesh$loc[mesh$segm$bnd$idx[j1, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j1, 2], 1:2]/2, 
                          mesh$loc[mesh$segm$bnd$idx[j2, 1], 1:2]/2 +
                            mesh$loc[mesh$segm$bnd$idx[j2, 2], 1:2]/2))
        yy <- p[,2]-mean(p[,2])/2-mesh$loc[i, 2]/2
        xx <- p[,1]-mean(p[,1])/2-mesh$loc[i, 1]/2
      } else {
        yy <- p[,2]-mesh$loc[i, 2]
        xx <- p[,1]-mesh$loc[i, 1]
      }
      Polygon(p[order(atan2(yy,xx)), ])
    })
    return(SpatialPolygons(lapply(1:mesh$n, function(i)
      Polygons(list(pls[[i]]), i))))
  }
  else stop("It only works for R2!")
}

CoxLikDualMeshExtra.prepare.new <- function(OmegaMesh,A_SpatialPolygons,A_SpatialLines,A_SpatialPoints){
  
  newPoints <- getIntersectionAndInteriorPoints.new(OmegaMesh = OmegaMesh,
                                                    A_SpatialPolygons = A_SpatialPolygons,
                                                    A_SpatialLines = A_SpatialLines,
                                                    A_SpatialPoints = A_SpatialPoints) # Finding which OmegaTriangles have interior points and intersections with the lines of A and 
  
  AMeshList <- Omega2AMesh.new(newPoints=newPoints,
                               OmegaMesh = OmegaMesh,
                               A_SpatialPolygons = A_SpatialPolygons)
  
  intPoints <- unique(do.call(rbind,unlist(AMeshList,recursive = F)))
  
  ## Computing the area of every new triangle (within the A domain)
  newTriangAreaList <- lapply(X = AMeshList,TriangArea2)
  
  ## Linking the new triangles to the vertices they are made of
  AVerticiesPerTriangList <- list()
  for (i in 1:length(AMeshList)){
    if(is.null(AMeshList[[i]])){
      AVerticiesPerTriangList[[i]] <- NULL
    } else {
      AVerticiesPerTriangList[[i]] <- list()
      for(j in 1:length(AMeshList[[i]])){
        crossdists=proxy::dist(AMeshList[[i]][[j]],intPoints,method="euclidean")
        AVerticiesPerTriangList[[i]][[j]] <- apply(X = crossdists,MARGIN = 1,FUN=which.min)
      }
    }
  }
  
  intPointsWeight <- rep(0,nrow(intPoints))
  ## Assigning 1/3 of the area of each triangle to the vertices they are made of
  for(i in 1:length(AVerticiesPerTriangList)){
    if (!is.null(AVerticiesPerTriangList[[i]])){
      for(j in 1:length(AVerticiesPerTriangList[[i]])){
        area <- newTriangAreaList[[i]][[j]]
        vertices <- AVerticiesPerTriangList[[i]][[j]]
        for (k in 1:3){
          intPointsWeight[vertices[k]] =  intPointsWeight[vertices[k]] +area/3
        }
      }
    }
  }
  
  # Mapping the OmegaMesh to the new integration points
  intPointsAMat <- INLA::inla.spde.make.A(mesh=OmegaMesh, loc = intPoints)
  
  
  # Here we also map the weights back to the mesh locations uisng the inlabru function integration_weight_construction
  inter <- list(loc = intPoints, weight = intPointsWeight)
  allmeshpointsWeights <- inlabru:::integration_weight_construction(mesh = OmegaMesh,inter = inter)$weight
  
  
  return(list(intPointsWeight=intPointsWeight,intPointsAMat=intPointsAMat,allmeshpointsWeights = allmeshpointsWeights))
}




getIntersectionAndInteriorPoints.new <- function(OmegaMesh,A_SpatialPolygons,A_SpatialLines,A_SpatialPoints){
  noOmegaTriangles <- nrow(OmegaMesh$graph$tv)
  
  # Creating Line and Polygon objects of all triangles in the Omegamesh
  polylist0 <- list()
  lineList0 <- list()
  for(i in 1:noOmegaTriangles){
    triangle <- OmegaMesh$loc[OmegaMesh$graph$tv[i,],1:2]
    lineList0[[i]] <- sp::Lines(list(sp::Line(rbind(triangle,triangle[1,]))),ID=paste0("_",i,"_"))
    polylist0[[i]] <- sp::Polygons(list(sp::Polygon(rbind(triangle,triangle[1,]),hole=F)),ID=paste0("_",i,"_"))
  }
  OmegaLinesList <- sp::SpatialLines(lineList0)
  OmegaPolyList <- sp::SpatialPolygons(polylist0)
  
  intersectionPoints <- rgeos::gIntersection(OmegaLinesList,A_SpatialLines,byid=T) # Finding line intersections between A and the mesh lines
  interiorPoints <- rgeos::gIntersection(OmegaPolyList,A_SpatialPoints,byid=T) # Finding "corner points" of A which are within the mesh (all should be)
  
  # Merging all info in one data.table, using the information on which triangles the matches/intersections occurs
  
  dt1 <- data.table::data.table(triangleNo =as.numeric(getstr(rownames(intersectionPoints@coords))),intersectionPoints@coords,interior=F)
  dt2 <- data.table::data.table(triangleNo =as.numeric(getstr(rownames(interiorPoints@coords))),interiorPoints@coords,interior=T)
  
  dt <- unique(data.table::rbindlist(list(dt1,dt2)))
  
  coveredOmegaTriangles <- which(rgeos::gCovers(A_SpatialPolygons,OmegaPolyList,byid=T)) # Fiding which triangles are completely covered by the A domain.
  intersectingOmegaTriangles <- sort(unique(dt$triangleNo))
  
  ret <- list(newPointsdt = dt,
              intersectingOmegaTriangles=intersectingOmegaTriangles,
              coveredOmegaTriangles = coveredOmegaTriangles)
  return(ret)
}

getIntersectionAndInteriorPoints.new <- function(OmegaMesh,A_SpatialPolygons,A_SpatialLines,A_SpatialPoints){
  noOmegaTriangles <- nrow(OmegaMesh$graph$tv)
  
  # Creating Line and Polygon objects of all triangles in the Omegamesh
  polylist0 <- list()
  lineList0 <- list()
  for(i in 1:noOmegaTriangles){
    triangle <- OmegaMesh$loc[OmegaMesh$graph$tv[i,],1:2]
    lineList0[[i]] <- sp::Lines(list(sp::Line(rbind(triangle,triangle[1,]))),ID=paste0("_",i,"_"))
    polylist0[[i]] <- sp::Polygons(list(sp::Polygon(rbind(triangle,triangle[1,]),hole=F)),ID=paste0("_",i,"_"))
  }
  OmegaLinesList <- sp::SpatialLines(lineList0)
  OmegaPolyList <- sp::SpatialPolygons(polylist0)
  
  intersectionPoints <- rgeos::gIntersection(OmegaLinesList,A_SpatialLines,byid=T) # Finding line intersections between A and the mesh lines
  interiorPoints <- rgeos::gIntersection(OmegaPolyList,A_SpatialPoints,byid=T) # Finding "corner points" of A which are within the mesh (all should be)
  
  # Merging all info in one data.table, using the information on which triangles the matches/intersections occurs
  
  dt1 <- data.table::data.table(triangleNo =as.numeric(getstr(rownames(intersectionPoints@coords))),intersectionPoints@coords,interior=F)
  dt2 <- data.table::data.table(triangleNo =as.numeric(getstr(rownames(interiorPoints@coords))),interiorPoints@coords,interior=T)
  
  dt <- unique(data.table::rbindlist(list(dt1,dt2)))
  
  coveredOmegaTriangles <- which(rgeos::gCovers(A_SpatialPolygons,OmegaPolyList,byid=T)) # Fiding which triangles are completely covered by the A domain.
  intersectingOmegaTriangles <- sort(unique(dt$triangleNo))
  
  ret <- list(newPointsdt = dt,
              intersectingOmegaTriangles=intersectingOmegaTriangles,
              coveredOmegaTriangles = coveredOmegaTriangles)
  return(ret)
}


Omega2AMesh.new <- function(newPoints,OmegaMesh,A_SpatialPolygons,minTriangArea=10^(-6),maxDigits = 10){
  
  newTriangleList <- list()
  for(i in newPoints$intersectingOmegaTriangles){
    newTriangleList[[i]] <- list()
    triangle <- OmegaMesh$loc[OmegaMesh$graph$tv[i,],1:2]
    extrapoints <- as.matrix(newPoints$newPointsdt[triangleNo==i,.(x,y)])
    both <- rbind(triangle,extrapoints)
    pslgObj <- RTriangle::pslg(P=unique(round(both,digits = maxDigits)))#,PB=c(rep(1,nrow(triangle)),rep(0,nrow(extrapoints))))
    
    newtriangles <- RTriangle::triangulate(p=pslgObj,S=0)
    tempTriangList <- list()
    k <- 1
    for (j in 1:nrow(newtriangles$T)){
      tempTriang <- newtriangles$P[newtriangles$T[j,],]
      tempPolyList <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(rbind(tempTriang,tempTriang[1,]),hole=F)),ID=paste0("_",j,"_"))))
      triangIntersection <- rgeos::gIntersection(A_SpatialPolygons,tempPolyList)
      intersectArea <- rgeos::gArea(triangIntersection)
      if(intersectArea>minTriangArea){ # Needed to remove very flat or tiny triangles
        newTriangleList[[i]][[k]] <- tempTriang
        k <- k+1
      } else {
        newTriangleList[[i]][[k]] <- NULL
      }
    }
    if(length(newTriangleList[[i]])==0){
      newTriangleList[[i]] <- NULL
    }
  }
  
  # Also adding the OmegaTriangles which are completely covered by the A polygon
  for (i in newPoints$coveredOmegaTriangles){
    triangle <- OmegaMesh$loc[OmegaMesh$graph$tv[i,],1:2]
    newTriangleList[[i]] <- list(triangle)
  }
  
  return(newTriangleList)
}


TriangArea <- function(x){
  rgeos::gArea(coo2sp(x))
}
TriangArea2 <- function(x){
  if(is.list(x)){
    sapply(x,TriangArea)
  } else {
    if(is.null(x)){
      0
    } else {
      TriangArea(x)
    }
  }
}

CoxLikBarycentricPointSpread.prepare.new <- function(OmegaMesh,A_List_of_SpatialPolygons){
  
  # Creating new integration points and associated weights with the 28th june 2019 method implemented in inlabru
  if (length(A_List_of_SpatialPolygons)==1){
    ips <- inlabru::ipoints(region = A_List_of_SpatialPolygons[[1]], domain = OmegaMesh)
  } else {
    ips_list <- lapply(A_List_of_SpatialPolygons,function(x) inlabru::ipoints(region = x, domain = OmegaMesh))
    
    ips <- do.call(rbind,ips_list)
    dt <- data.table::data.table(weight=ips$weight,sp::coordinates(ips))
    dt[,weight := sum(weight),by=.(x,y)]
    
    ips$weight <- dt$weight
    dup <- which(duplicated(sp::coordinates(ips)))
    ips <- ips[-dup,]
  }
  
  # Mapping the OmegaMesh to the new integration points
  intPointsAMat <- INLA::inla.spde.make.A(mesh=OmegaMesh, loc = coordinates(ips))
  
  return(list(intPointsWeight=ips$weight,intPointsAMat=intPointsAMat))
}

CoxLikBarycentricPointSpread.prepare.new.nsub <- function(OmegaMesh,A_List_of_SpatialPolygons,nsub=NULL){
  
  # Creating new integration points and associated weights with the 28th june 2019 method implemented in inlabru
  if (length(A_List_of_SpatialPolygons)==1){
    ips <- inlabru::ipoints(region = A_List_of_SpatialPolygons[[1]], domain = OmegaMesh,nsub=nsub)
  } else {
    ips_list <- lapply(A_List_of_SpatialPolygons,function(x) inlabru::ipoints(region = x, domain = OmegaMesh,nsub=nsub))
    
    ips <- do.call(rbind,ips_list)
    dt <- data.table::data.table(weight=ips$weight,sp::coordinates(ips))
    dt[,weight := sum(weight),by=.(x,y)]
    
    ips$weight <- dt$weight
    dup <- which(duplicated(sp::coordinates(ips)))
    ips <- ips[-dup,]
  }
  
  # Mapping the OmegaMesh to the new integration points
  intPointsAMat <- INLA::inla.spde.make.A(mesh=OmegaMesh, loc = coordinates(ips))
  
  return(list(intPointsWeight=ips$weight,intPointsAMat=intPointsAMat))
}


CoxLikBarycentricPointSpread.prepare.new2 <- function(OmegaMesh,A_List_of_SpatialPolygons){
  
  dt0 <- data.table(x=OmegaMesh$loc[,1],y=OmegaMesh$loc[,2])
  
  # Creating new integration points and associated weights with the 28th june 2019 method implemented in inlabru
  if (length(A_List_of_SpatialPolygons)==1){
    ips <- inlabru::ipoints(region = A_List_of_SpatialPolygons[[1]], domain = OmegaMesh)
  } else {
    ips_list <- lapply(A_List_of_SpatialPolygons,function(x) inlabru::ipoints(region = x, domain = OmegaMesh))
    ips <- do.call(rbind,ips_list)
  }
  
  # Matching the 
  dt <- data.table::data.table(weight=ips$weight,sp::coordinates(ips))
  dt[,weight := sum(weight),by=.(x,y)]
  dt <- unique(dt)
  
  ret <- merge(dt0,dt,by=c("x","y"),all.x=T,all.y=T,sort=F) # Important not to sort the outpout DT here
  ret[is.na(weight),weight:=0]
  
  if(!all.equal(ret[,.(x,y)],dt0)){
    stop("weight vector does not match mesh locations")
  } else {
    return(ret$weight)
  }
}

CoxLikBarycentricPointSpread.prepare.new2.nsub <- function(OmegaMesh,A_List_of_SpatialPolygons,nsub=NULL){
  
  dt0 <- data.table(x=OmegaMesh$loc[,1],y=OmegaMesh$loc[,2])
  
  # Creating new integration points and associated weights with the 28th june 2019 method implemented in inlabru
  if (length(A_List_of_SpatialPolygons)==1){
    ips <- inlabru::ipoints(region = A_List_of_SpatialPolygons[[1]], domain = OmegaMesh,nsub=nsub)
  } else {
    ips_list <- lapply(A_List_of_SpatialPolygons,function(x) inlabru::ipoints(region = x, domain = OmegaMesh,nsub=nsub))
    ips <- do.call(rbind,ips_list)
  }
  
  # Matching the 
  dt <- data.table::data.table(weight=ips$weight,sp::coordinates(ips))
  dt[,weight := sum(weight),by=.(x,y)]
  dt <- unique(dt)
  
  ret <- merge(dt0,dt,by=c("x","y"),all.x=T,all.y=T,sort=F) # Important not to sort the outpout DT here
  ret[is.na(weight),weight:=0]
  
  if(!all.equal(ret[,.(x,y)],dt0)){
    stop("weight vector does not match mesh locations")
  } else {
    return(ret$weight)
  }
}

CoxLikBarycentricPointSpread.prepare.new3 <- function(OmegaMesh,A_SpatialPolygons){
  
  dt0 <- data.table(x=OmegaMesh$loc[,1],y=OmegaMesh$loc[,2])
  
  # Creating new integration points and associated weights with the 28th june 2019 method implemented in inlabru
  ips <- inlabru::ipoints(region = A_SpatialPolygons, domain = OmegaMesh)
  
  # Matching the 
  dt <- data.table::data.table(weight=ips$weight,sp::coordinates(ips))
  dt[,weight := sum(weight),by=.(x,y)]
  dt <- unique(dt)
  
  ret <- merge(dt0,dt,by=c("x","y"),all.x=T,all.y=T,sort=F) # Important not to sort the outpout DT here
  ret[is.na(weight),weight:=0]
  
  if(!all.equal(ret[,.(x,y)],dt0)){
    stop("weight vector does not match mesh locations")
  } else {
    return(ret$weight)
  }
}

CoxLikBarycentricPointSpread.prepare.new3.nsub <- function(OmegaMesh,A_SpatialPolygons,nsub=NULL){
  
  dt0 <- data.table(x=OmegaMesh$loc[,1],y=OmegaMesh$loc[,2])
  
  # Creating new integration points and associated weights with the 28th june 2019 method implemented in inlabru
  ips <- inlabru::ipoints(region = A_SpatialPolygons, domain = OmegaMesh,nsub=nsub)
  
  # Matching the 
  dt <- data.table::data.table(weight=ips$weight,sp::coordinates(ips))
  dt[,weight := sum(weight),by=.(x,y)]
  dt <- unique(dt)
  
  ret <- merge(dt0,dt,by=c("x","y"),all.x=T,all.y=T,sort=F) # Important not to sort the outpout DT here
  ret[is.na(weight),weight:=0]
  
  if(!all.equal(ret[,.(x,y)],dt0)){
    stop("weight vector does not match mesh locations")
  } else {
    return(ret$weight)
  }
}



OracleTruthIntegral.evaluator.new <- function(Likelihood.OracleTruth.prepare.list,z){
  Likelihood.OracleTruth.prepare.list$area.rect*sum(exp(z[Likelihood.OracleTruth.prepare.list$these.rf.xy.grid.in.A]))  
}

MonteCarloIntegrator.evaluator.new <- function(MonteCarloIntegrator.mapper,MonteCarloIntegrator.prepare.list.0,z){
  
  ret <- sum(exp(MonteCarloIntegrator.mapper%*%z)*MonteCarloIntegrator.prepare.list.0$area.rect)
  
  return(ret)
}

ApproxIntegral.evaluator <- function(OmegaVerticesWeight,z){
  ret <- sum(OmegaVerticesWeight*exp(z),na.rm=T)
  return(ret)
}


ApproxIntegral.mapping.evaluator <- function(intPointsWeight,intPointsMapping,z){
  ret <- sum(as.vector(intPointsWeight*exp(intPointsMapping%*%z)),na.rm=T)
  return(ret)
}


getstr = function(mystring, initial.character="_", final.character="_")
{
  # check that all 3 inputs are character variables
  if (!is.character(mystring))
  {
    stop('The parent string must be a character variable.')
  }
  
  if (!is.character(initial.character))
  {
    stop('The initial character must be a character variable.')
  }
  
  
  if (!is.character(final.character))
  {
    stop('The final character must be a character variable.')
  }
  
  add=0
  if(initial.character==final.character){add=1}
  
  # pre-allocate a vector to store the extracted strings
  snippet = rep(0, length(mystring))
  
  for (i in 1:length(mystring))
  {
    # extract the initial position
    initial.position = gregexpr(initial.character, mystring[i])[[1]][1] + 1
    
    # extract the final position
    final.position = gregexpr(final.character, mystring[i])[[1]][1+add] - 1
    
    # extract the substring between the initial and final positions, inclusively
    snippet[i] = substr(mystring[i], initial.position, final.position)
  }
  return(snippet)
}

