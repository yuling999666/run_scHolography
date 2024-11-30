library(scHolography)
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(Gmisc)
library(dplyr)
library(plotly)
#data_list=c('151507','151508','151509','151510','151669','151670','151671','151672','151673','151674','151675','151676')
############# read data and create seurat object
data_list=c('151507','151508','151509','151510')
create_obj<-function(data_name,slice_order){
  matrix_dir=pathJoin('/home/data3/yuling/xenium_region2_5599/Dataset/DLPFC12',data_name,'filtered_feature_bc_matrix/')
  counts = Seurat::Read10X(data.dir = matrix_dir) 
  data = Seurat::CreateSeuratObject(
    counts = counts , 
    project = 'test', 
    assay = 'Spatial')
  data$slice = slice_order
  data$region = 'test' 
  imgpath = pathJoin('/home/data3/yuling/xenium_region2_5599/Dataset/DLPFC12',data_name,'spatial')
  img = Seurat::Read10X_Image(image.dir = imgpath)  
  Seurat::DefaultAssay(object = img) <- 'Spatial'  
  img = img[colnames(x = data)]  
  data[['image']] = img
 
  meta_1<-read.csv(pathJoin('/home/data3/yuling/xenium_region2_5599/Dataset/DLPFC12', data_name,'gt/tissue_positions_list_GTs.txt'),header=FALSE)
  colnames(meta_1)[7] <- "orig.cluster"
  colnames(meta_1)[1] <- "cell_id"
  meta_1<-meta_1[order(meta_1$cell_id),]
  rownames(meta_1)<-meta_1$cell_id
  data <- AddMetaData(data, metadata =meta_1['orig.cluster'],col.name = "orig.cluster")
  seurat_obj<- RenameCells(data, add.cell.id = paste0("Obj",slice_order))
  return (seurat_obj)
}
###################################

slice_num=1
data_1<-create_obj(data_name = data_list[slice_num],slice_order = slice_num)
data_2<-create_obj(data_name = data_list[slice_num+1],slice_order = slice_num+1)
data_3<-create_obj(data_name = data_list[slice_num+2],slice_order = slice_num+2)
data_4<-create_obj(data_name = data_list[slice_num+3],slice_order = slice_num+3)

data_1<-subset(data_1,cells=which(is.na(data_1$orig.cluster)==F))
data_2<-subset(data_2,cells=which(is.na(data_2$orig.cluster)==F))
data_3<-subset(data_3,cells=which(is.na(data_3$orig.cluster)==F))
data_4<-subset(data_4,cells=which(is.na(data_4$orig.cluster)==F))
####### align data to reference space
#data_1@meta.data[,c("CentroidX","CentroidY")]<-data_1@images$image@coordinates[,c("imagerow","imagecol")]
data_2@images$image@coordinates[,c("imagerow","imagecol")] <- t(t(data_2@images$image@coordinates[,c("imagerow","imagecol")])-c(median(data_2@images$image@coordinates[,c("imagerow")]),median(data_2@images$image@coordinates[,c("imagecol")]))+c(median(data_1@images$image@coordinates[,c("imagerow")]),median(data_1@images$image@coordinates[,c("imagecol")])))
data_3@images$image@coordinates[,c("imagerow","imagecol")]<- t(t(data_3@images$image@coordinates[,c("imagerow","imagecol")])-c(median(data_3@images$image@coordinates[,c("imagerow")]),median(data_3@images$image@coordinates[,c("imagecol")]))+c(median(data_1@images$image@coordinates[,c("imagerow")]),median(data_1@images$image@coordinates[,c("imagecol")])))
data_4@images$image@coordinates[,c("imagerow","imagecol")]<- t(t(data_4@images$image@coordinates[,c("imagerow","imagecol")])-c(median(data_4@images$image@coordinates[,c("imagerow")]),median(data_4@images$image@coordinates[,c("imagecol")]))+c(median(data_1@images$image@coordinates[,c("imagerow")]),median(data_1@images$image@coordinates[,c("imagecol")])))
combined_obj <- merge(
  x = data_2,
  y = data_3,
  project = "CombinedProject")
combined_obj <- merge(
  x = combined_obj,
  y = data_4,
  project = "CombinedProject")

sp.integrated <- dataAlign(data_1,combined_obj,whichReference=1,nPCtoUse = 32,scProcessed = F,future.size = 8000)
scHolography.obj<-trainHolography(sp.integrated,n.slot = 30,n.pcUse = 32,n.pcOut = 32,n.repeat = 30)

### To orient the plot in a more intuitive way
scHolography.obj$scHolography.sc$z3d_sp <- -scHolography.obj$scHolography.sc$z3d_sp
scHolography.obj$scHolography.sc$x3d_sp <- -scHolography.obj$scHolography.sc$x3d_sp
scHolography.obj$scHolography.sc$y3d_sp <- -scHolography.obj$scHolography.sc$y3d_sp
scene = list(camera = list(eye = list(x = -1, z = 0, y = 2)))

### To visualize the 3D structure colored by celltype

fig<-scHolographyPlot(scHolography.obj,color.by = "original_cluster",color=color)%>% plotly::layout(scene = scene)
htmlwidgets::saveWidget(fig, "plot_3d_by_clusters.html")
  
clusterDistanceBoxplot(scHolography.obj,annotationToUse = "celltype",reference.cluster = "orig.cluster",query.cluster.list = c("1","2","3","4","5","6"))
  
  
  
  
  
  
 