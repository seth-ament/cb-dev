# conda activate /local/projects-t3/idea/bherb/software/conda_env/r-4.1.0

library( Seurat )

children = readRDS('/local/projects/idea/sament2/cb/2020-10/integration_and_clustering/HsCB.AllTypes.Clean.SeuratObj.rds')
DefaultAssay(children) = 'RNA'
children = subset( children , idents = 'Purkinje' )

fetal = readRDS('/local/projects/idea/sament2/cb/2020-10/fetal/seurat.rds')
fetal = subset( fetal , idents = "01-PC" ) 
DefaultAssay(fetal) = 'RNA'

g = intersect( rownames(children) , rownames(fetal) )

children$donor = children$sample_name
fetal$donor = fetal$sample_id

obj = merge( children , fetal )
obj = subset( obj , features = g )

saveRDS( obj , file = 'Purkinje_cells_combined.seurat.rds' )







