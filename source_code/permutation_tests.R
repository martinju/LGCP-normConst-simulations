

### Permutation test
library(data.table)

results_folder <- "results"
folder_string <- "20191023_0945_iHcBtj"

DT =    fread(file.path(results_folder,"full",folder_string,"full.res.dt.csv"))


DT1 <- data.table(boundary.splitter = c(1,2,4,8),
                  subDomains = c(1,2,8,32))

DT2 <- data.table(INLA.max.edge = c(4,3, 2,1, 0.50, 0.25),
                  meshNodes = c(126,   183,   359,  1179,  4324, 16429)) # How many mesh nodes that are within the inner mesh region.
# Changing some names
DT[,field_variance:=sigma2x.true]
DT[,sigma2x.true:=NULL]
DT[,field_range:=round(sqrt(8)/kappa.true,5)]
DT[,kappa.true:=NULL]

DT <- merge(DT,DT1,by = "boundary.splitter")
DT <- merge(DT,DT2,by = "INLA.max.edge")
DT[,boundary.splitter:=NULL]
DT[,INLA.max.edge:=NULL]
DT[,meshNodes:=as.factor(meshNodes)]

DT[,folder_string :=NULL]
DT[,FEM.approx.mesh :=NULL]
DT[,OmegaArea :=NULL]
#DT_list[[i]][,OracleTruth :=NULL]
methods <- c("MeshExact","Voronoi","DualMesh","DualMeshExtra","DualMeshExtraMeshMapped","BarycentricPointSpread","BarycentricPointSpreadManyPoints")#,"AvgDetMethods")
bycols <- c("beta0.true","field_variance","field_range","subDomains","meshNodes")

DT_squared_error <- DT[,lapply(.SD, function(x){(x-OracleTruth)^2}),.SDcols=methods,by = bycols]

DT_squared_error[,grp := .GRP,by=bycols]

no_permut <- 100
no_samp <- DT_squared_error[,.N,by=grp]$N[1] # 10000
no_grp <- max(DT_squared_error[,grp])

samp_vec <- c(rep(-1,no_samp/2),rep(1,no_samp/2))

signmat <- matrix(NA,ncol=no_permut,nrow=no_samp)
set.seed(123)
for(i in 1:no_permut){
  signmat[,i] <- sample(samp_vec,size = no_samp,replace = F)
}

method_comb_dt <- as.data.table(t(combn(methods,2)))
colnames(method_comb_dt) <- c("method1","method2")
method_comb_dt[,comb:=.I]

bycols2 <- c(bycols,"grp")
res_DT <- unique(DT_squared_error[,..bycols2])
res_DT <- data.table(merge.data.frame(res_DT,method_comb_dt,all = T))
res_DT[,p_val:=1]

# Just brute-forcing this ()
for (j in 1:no_grp){
  for (i in 1:nrow(method_comb_dt)){
    met1 <- method_comb_dt$method1[i]
    met2 <- method_comb_dt$method2[i]
    comb <- method_comb_dt$comb[i]
    
    vec <- unlist(DT_squared_error[grp==j,..met1] -   DT_squared_error[grp==j,..met2],use.names = F)
    obs_val <- mean(vec)
    permut_val <- as.vector(vec %*% signmat)/no_permut
    p_val <- mean(abs(permut_val) > abs(obs_val))
    #DT_squared_error[,by=grp]
    res_DT[comb == i & grp == j,p_value:=..p_val]
  }
  print(j)
}

res_DT[,comb:=NULL]
res_DT[,grp:=NULL]
# 
# res_DT2 <- copy(res_DT)
# res_DT2[,method1temp:=method1]
# res_DT2[,method1:=method2]
# res_DT2[,method2:=method1temp]
# res_DT2[,method1temp:=NULL]

res_DT <- rbind(res_DT,res_DT2)


mean(res_DT$p_value<0.05)
#dir.create(file.path(results_folder,"permutation_test",folder_string),recursive = T)

fwrite(x = res_DT,file.path(results_folder,"permutation_test",folder_string,"permutation.test.results.csv"))

