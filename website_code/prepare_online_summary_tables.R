
#### Preparing summary tables for online publishing

library(data.table)

res_folder <- "../Results_new/20191023_0945_iHcBtj"

files = c("summary.RMSE.res.dt_10000.csv",
                  "summary.bias.res.dt_10000.csv",
                  "summary.sd.res.dt_10000.csv")
# Need to add skillscore here as well!
names_files = c("RMSE","bias","sd")

methods = c("MeshExact", "Voronoi", "DualMesh", "DualMeshExtra", "DualMeshExtraMeshMapped",
            "BarycentricPointSpread", "BarycentricPointSpreadManyPoints","AvgApproxMethods")

bycols = c("field_variance","field_range","subDomains","meshNodes","beta0.true")

DT1 <- data.table(boundary.splitter = c(1,2,4,8),
                  subDomains = c(1,2,8,32))

DT2 <- data.table(INLA.max.edge = c(4,3, 2,1, 0.50, 0.25),
                  meshNodes = c(126,   183,   359,  1179,  4324, 16429)) # How many mesh nodes that are within the inner mesh region.



DT_list = list()
for (i in 1:length(files)){
  DT_list[[i]] =  fread(file.path(res_folder,files[i]))                  
  
  
  DT_list[[i]][,field_variance:=sigma2x.true]
  DT_list[[i]][,sigma2x.true:=NULL]
  DT_list[[i]][,field_range:=round(sqrt(8)/kappa.true,5)]
  DT_list[[i]][,kappa.true:=NULL]
  
  DT_list[[i]] <- merge(DT_list[[i]],DT1,by = "boundary.splitter")
  DT_list[[i]] <- merge(DT_list[[i]],DT2,by = "INLA.max.edge")
  DT_list[[i]][,boundary.splitter:=NULL]
  DT_list[[i]][,INLA.max.edge:=NULL]
  DT_list[[i]][,meshNodes:=as.factor(meshNodes)]
  
  DT_list[[i]][,folder_string :=NULL]
  DT_list[[i]][,FEM.approx.mesh :=NULL]
  DT_list[[i]][,OmegaArea :=NULL]
  DT_list[[i]][,OracleTruth :=NULL]
  
  DT_list[[i]][,result_type:=..names_files[i]]
  setcolorder(DT_list[[i]],neworder = c("result_type",bycols,methods))

  #colnames(DT_list[[i]])[-(1:length(bycols))] = paste0(names_files[i],"_",methods)
  
}

DT = rbindlist(DT_list)


#DT = DT_list[[1]]
#for (i in 2:length(DT_list)){
#  DT = merge(DT,DT_list[[i]])
#}


setkeyv(DT,c("result_type",bycols))

#fwrite(DT,file.path("../Results_new/online_summary_table.csv"))
fwrite(DT,file.path("../Results_new/online_summary_table_rbinded.csv"))

