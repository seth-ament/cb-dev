# conda activate /local/projects-t3/idea/bherb/software/conda_env/r-4.1.0
  
library( Seurat )
# library( monocle3 )
library( sva )
source('/home/sament/knn-smoothing/knn_smooth.R')

obj = readRDS('Purkinje_cells_combined.seurat.rds')


# metadata

group = rep( 'ctrl' , ncol(obj) )
group[ grep('infl',obj$donor,ignore.case=T) ] = 'infl'
obj$group = group

Idents(obj) = obj$donor

age = rep( NA , ncol(obj) )
age[ obj$age == '9 PCW' ] = (9-40)/52
age[ obj$age == '10 PCW' ] = (10-40)/52
age[ obj$age == '11 PCW' ] = (11-40)/52
age[ obj$age == '12 PCW' ] = (12-40)/52
age[ obj$age == '14 PCW' ] = (14-40)/52
age[ obj$age == '16 PCW' ] = (16-40)/52
age[ obj$age == '17 PCW' ] = (17-40)/52
age[ obj$age == '18 PCW' ] = (18-40)/52
age[ obj$age == '20 PCW' ] = (20-40)/52
age[ grep( '10X' , obj$donor ) ] = 60
age[ obj$donor == 'NIH_HBCC_1218_Hu_Inflammation_4_2yold_Cerebellum' ] = 1533/365
age[ obj$donor == 'UMB1599_WF_Control_2a_hcerebellumVermis_PMI44h_nuclei' ] = 991/365
age[ obj$donor == 'UMB1798_AA_F_INFL_1a_hcerebellumVermis_PMI24h_nuclei' ] = 654/365
age[ obj$donor == 'UMB4327_AA_F_Control_5a_hcerebellumVermis_PMI24h_nuclei' ] = 2093/365
age[ obj$donor == 'UMB4332_AA_M_INFL_5a_hcerebellumVermis_PMI18h_nuclei' ] = 2068/365
age[ obj$donor == 'UMB5180_WM_Control_1a_hcerebellumVermis_PMI25h_nuclei' ] = 628/365
age[ obj$donor == 'UMB5564_AA_M_Control_2a_hcerebellumVermis_PMI13h_nuclei' ] = 940/365
age[ obj$donor == 'UMB_BTB_1284_Hu_Inflammation_1218daysold_Cerebellum' ] = 1218/365
age[ obj$donor == 'UMB_BTB1488_Hu_Control_502daysold_Cerebellum_Vermis_Nuclei_T' ] = 502/365
age[ obj$donor == 'UMB_BTB_1906_Hu_Inflammation_801daysold_Cerebellum' ] = 801/365
age[ obj$donor == 'UMB_BTB_4321_Hu_Inflammation_734daysold_Cerebellum' ] = 734/365
age[ obj$donor == 'UMB_BTB_510_Hu_Control_901daysold_Cerebellum' ] = 901/365
age[ obj$donor == 'UMB_BTB_5282_Hu_Control_1035daysold_Cerebellum' ] = 1035/365
age[ obj$donor == 'UMB_BTB_5576_Hu_Control_1017daysold_Cerebellum' ] = 1017/365
age[ obj$donor == 'UMB_BTB_M3828M_Hu_Control_671daysold_Cerebellum_Vermis_Nucle' ] = 671/365
age[ obj$donor == 'MBC524_WM_Control_47a_hcerebellumL9sec_PMI6h_nuclei_10000' ] = 47
age[ obj$donor == 'MBC524_WM_Control_47a_hcerebellumL9sec_PMI6h_nuclei' ] = 32
age[ obj$donor == 'MBC554_WF_Control_32a_hcerebellumL9sec_PMI7h_nuclei' ] = 32

obj$age = age
obj$log_yrs = log10( age + 0.75 )
Idents(obj) = paste( age , 'yrs' , sep = '_' ) 

donors = unique( obj$donor )
samples = unique(obj$sample_name)
batch = rep(1,ncol(obj))
batch[ which( obj$donor %in% donors[5:10] ) ] = 2
obj$batch = batch

obj.down = subset( obj , downsample = 100 )

expr = GetAssayData(obj.down,slot='counts')
meta = obj.down@meta.data

  mat = as.matrix(expr)
  pct = rowSums(mat>0)/ncol(mat)
  mat = mat[ pct >= 0.03 , ]
  m = meta
  colnames(mat) = rownames(m)

  combat = ComBat_seq(
    counts = mat ,
    batch = factor(m$batch) ,
    group = factor(m$group) ,
    covar_mod = model.matrix( ~ m$log_yrs )
  )

  smoothed = knn_smoothing( mat = combat , k=3 , d=30 ) 
  smoothed = CreateSeuratObject( counts = smoothed , meta.data = m )
  smoothed = NormalizeData(smoothed)

  combat.corr = CreateSeuratObject( counts = combat , meta.data = m )

  saveRDS( smoothed , file = 'smoothed.rds' )
  saveRDS( combat.corr , file = 'combat.corr.rds' )
  saveRDS( obj.down , file = 'merged.downsampled.rds' )

















expr = GetAssayData(obj.down,slot='counts')
cell_meta = obj.down@meta.data[ , c('donor','group','age','log_yrs') ]
gene_meta = as.matrix(rownames(obj.down))
rownames(gene_meta) = rownames(obj.down)
colnames(gene_meta) = 'gene_short_name'
cds <- new_cell_data_set( expr , cell_meta , gene_meta )

cds <- preprocess_cds(cds, num_dim = 10)
# cds <- align_cds(cds, alignment_group = 'group')
  
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds,use_partition=F)

start_cells = colnames(obj.down)[ which( obj.down$age < -0.59 ) ]
cds = order_cells(cds,root_cells=start_cells)

pdf('monocle.umap.pdf')
plot_cells(cds, label_groups_by_cluster=F,  color_cells_by = "donor")
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "pseudotime")
dev.off()






