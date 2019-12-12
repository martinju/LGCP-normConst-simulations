 # Processing raw inputs to save different types of summaries for visualization etc.

library(data.table)

res_folder <- "../Results_new/"

dirs <- list.dirs(res_folder,full.names = F,recursive = F)
exclude.dirs <- c("old","temp")

data.dirs <- dirs[which(!(dirs %in%exclude.dirs))]

model_names <- c( "MeshMC", "ApproxVoronoi", "ApproxDualMesh", "ApproxFabian", "ApproxInlabru","OracleTruth",
                  "MeshExact", "Voronoi", "DualMesh", "DualMeshExtra", "DualMeshExtraMeshMapped",
                  "BarycentricPointSpread", "BarycentricPointSpreadManyPoints")
detint_methods <- c("Voronoi","DualMesh","DualMeshExtraMeshMapped","BarycentricPointSpread","BarycentricPointSpreadManyPoints")
not_bycols <- "seed"


for (this in data.dirs){
  res.dt <- fread(file.path(res_folder,this,"full.res.dt.csv"))
  
  
  
  methods <- names(res.dt)[names(res.dt) %in% model_names]
  if (all(detint_methods %in% methods)){
    res.dt[,AvgDetMethods:=rowMeans(.SD),.SDcols=detint_methods]
    methods <- c(methods,"AvgDetMethods")
  }
  bycols <- names(res.dt)[!(names(res.dt) %in% c(methods,model_names,not_bycols))]
  

  summary.RMSE.res.dt <- res.dt[,lapply(.SD, function(x){sqrt(mean((x-OracleTruth)^2))}),.SDcols=methods,by=bycols]
  summary.MSE.res.dt <- res.dt[,lapply(.SD, function(x){(mean((x-OracleTruth)^2))}),.SDcols=methods,by=bycols]
  
  summary.sqbias.res.dt <- res.dt[,lapply(.SD, function(x){(mean(x-OracleTruth))^2}),.SDcols=methods,by=bycols]
  summary.variance.res.dt <- res.dt[,lapply(.SD, function(x){var(x)}),.SDcols=methods,by=bycols]
  summary.sd.res.dt <- res.dt[,lapply(.SD, function(x){sd(x)}),.SDcols=methods,by=bycols]
  
  summary.bias.res.dt <- res.dt[,lapply(.SD, function(x){mean(x-OracleTruth)}),.SDcols=methods,by=bycols]
  
  fwrite(x = summary.RMSE.res.dt,file.path(res_folder,this,"summary.RMSE.res.dt.csv"))
  fwrite(x = summary.MSE.res.dt,file.path(res_folder,this,"summary.MSE.res.dt.csv"))
  fwrite(x = summary.sqbias.res.dt,file.path(res_folder,this,"summary.sqbias.res.dt.csv"))
  fwrite(x = summary.variance.res.dt,file.path(res_folder,this,"summary.variance.res.dt.csv"))
  fwrite(x = summary.sd.res.dt,file.path(res_folder,this,"summary.sd.res.dt.csv"))
  fwrite(x = summary.bias.res.dt,file.path(res_folder,this,"summary.bias.res.dt.csv"))
  
  if (this == "20191023_0945_iHcBtj"){ # Additional summary table based in only 10000 first samples for this model
    res.dt <- res.dt[seed < (min(res.dt$seed)+10000),]
    summary.RMSE.res.dt <- res.dt[,lapply(.SD, function(x){sqrt(mean((x-OracleTruth)^2))}),.SDcols=methods,by=bycols]
    summary.MSE.res.dt <- res.dt[,lapply(.SD, function(x){(mean((x-OracleTruth)^2))}),.SDcols=methods,by=bycols]
    
    summary.sqbias.res.dt <- res.dt[,lapply(.SD, function(x){(mean(x-OracleTruth))^2}),.SDcols=methods,by=bycols]
    summary.variance.res.dt <- res.dt[,lapply(.SD, function(x){var(x)}),.SDcols=methods,by=bycols]
    summary.sd.res.dt <- res.dt[,lapply(.SD, function(x){sd(x)}),.SDcols=methods,by=bycols]
    
    summary.bias.res.dt <- res.dt[,lapply(.SD, function(x){mean(x-OracleTruth)}),.SDcols=methods,by=bycols]
    
    fwrite(x = summary.RMSE.res.dt,file.path(res_folder,this,"summary.RMSE.res.dt_10000.csv"))
    fwrite(x = summary.MSE.res.dt,file.path(res_folder,this,"summary.MSE.res.dt_10000.csv"))
    fwrite(x = summary.sqbias.res.dt,file.path(res_folder,this,"summary.sqbias.res.dt_10000.csv"))
    fwrite(x = summary.variance.res.dt,file.path(res_folder,this,"summary.variance.res.dt_10000.csv"))
    fwrite(x = summary.sd.res.dt,file.path(res_folder,this,"summary.sd.res.dt_10000.csv"))
    fwrite(x = summary.bias.res.dt,file.path(res_folder,this,"summary.bias.res.dt_10000.csv"))
  }
  
  print(this)
}


