library(DT) # For producing html table
library(htmltools)
library(htmlwidgets)
library(data.table)

results_folder <- "results"
folder_string <- "20191023_0945_iHcBtj"

names_files = c("RMSE","bias","sd")



#### Changing some terminology for the results table
DT1 <- data.table(boundary.splitter = c(1,2,4,8),
                  subDomains = c(1,2,8,32))

DT2 <- data.table(INLA.max.edge = c(4,3, 2,1, 0.50, 0.25),
                  meshNodes = c(126,   183,   359,  1179,  4324, 16429)) # How many mesh nodes that are within the inner mesh region.

bycols = c("result_type","field_variance","field_range","subDomains","meshNodes","beta0.true")

DT_list = list()
for (i in 1:length(names_files)){
  DT_list[[i]] =    fread(file.path(results_folder,"summary",folder_string,paste0("summary.",names_files[i],".res.dt.csv")))
  
  # Changing some names
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
  #DT_list[[i]][,OracleTruth :=NULL]
  
  DT_list[[i]][,result_type:=..names_files[i]]
  setcolorder(DT_list[[i]],neworder = bycols)
  
  #colnames(DT_list[[i]])[-(1:length(bycols))] = paste0(names_files[i],"_",methods)
  
}

dat = rbindlist(DT_list)

#str(dat)
dat[,(bycols):=lapply(.SD, as.factor),.SDcols=bycols]

methods = colnames(dat)[!(colnames(dat) %in% bycols)]

dat[,(methods):=lapply(.SD, round,digits=9),.SDcols=methods]


methods_names = c("MeshExact", "Voronoi", "DualMesh", "DualMeshExtra", "DualMeshExtra-MeshMapped",
                  "Barycentric-PointSpread", "Barycentric-PointSpread-ManyPoints","AverageDetMesh")

bycols_names = c("field variance","field range", "#sub domains","#mesh nodes","beta0-true")


sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 2, 'Result type'),
      th(colspan = 5, 'Parameters'),
      th(colspan = 8, 'Methods')
    ),
    tr(
      lapply(c(bycols_names,methods_names), th)
    )
  )
))
#print(sketch)

htmltable = datatable(dat, filter = 'top',options = list(
  pageLength = 100, autoWidth = TRUE),container = sketch, rownames = FALSE) %>% 
  formatStyle("result_type",  color = 'white', backgroundColor = '#b399fd', fontWeight = 'bold') %>%
  formatStyle(bycols[-1],  color = 'white', backgroundColor = '#ec8ef5', fontWeight = 'bold') %>%
  formatStyle(methods,  color = 'black', backgroundColor = '#c4f68d', fontWeight = 'bold')

# Adding title as well
header_html = "<h2>LGCP simulation results</h2>
  <p>Full summary table for simulations performed relevant for the paper
Investigating mesh based approximation methods for the normalization constant in the Cox process likelihood by Martin Jullum. 
The source code for the simulations is available
<a href='https://github.com/martinju/LGCP-normConst-simulations'>here</a>.</p>"

online_table = htmlwidgets::prependContent(htmltable,htmltools::HTML(header_html))

DT::saveWidget(online_table, "sim_res.html",background = "#E7ECED",
               title = "Simulation results")


##### Web page for the permutations tests


dat <- fread(file.path(results_folder,"permutation_test",folder_string,"permutation.test.results_10000.csv"))

bycols2 <- c("field_variance","field_range","subDomains","meshNodes","beta0.true")
factorcols <- c(bycols2,"method1","method2")


dat[,(factorcols):=lapply(.SD, as.factor),.SDcols=factorcols]

htmltable = datatable(dat, filter = 'top',options = list(
  pageLength = 100, autoWidth = TRUE), rownames = FALSE) %>% 
  formatStyle(c("method1","method2"),  color = 'white', backgroundColor = '#b399fd', fontWeight = 'bold') %>%
  formatStyle(bycols2,  color = 'white', backgroundColor = '#ec8ef5', fontWeight = 'bold') %>%
  formatStyle("p_value",  color = 'black', backgroundColor = '#c4f68d', fontWeight = 'bold')

# Adding title as well
header_html = "<h2>LGCP simulation results</h2>
  <p>Permutation test for difference between the RMSE of any pair of approximation methods relevant for the paper
Investigating mesh based approximation methods for the normalization constant in the Cox process likelihood by Martin Jullum. 
Every single permutation test is carried out independently for every combination of parameters, using 
10 000 random permutations.
The source code for the simulations and permutation tests is available
<a href='https://github.com/martinju/LGCP-normConst-simulations'>here</a>.</p>"

online_table = htmlwidgets::prependContent(htmltable,htmltools::HTML(header_html))

DT::saveWidget(online_table, "permut_tests.html",background = "#E7ECED",
               title = "Permutation test results")

    



