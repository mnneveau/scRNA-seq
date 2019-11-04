#Single Cell Experiment Objects
sce.2I <- readRDS(file=file.path("C:/Users/meg.neveau/Downloads", "2I.sce.rds"))
sce.2I.imputed <- readRDS(file=file.path("C:/Users/meg.neveau/Downloads", "2I.imputed.sce.rds"))
sce.NoPP.imputed <- readRDS(file=file.path("C:/Users/meg.neveau/Downloads", "No-PP.imputed.sce.rds"))
sce.NoPP <- readRDS(file=file.path("C:/Users/meg.neveau/Downloads", "No-PP.sce.rds"))

library(iSEE)
app <- iSEE(sce.2I)

shiny::runApp(app)