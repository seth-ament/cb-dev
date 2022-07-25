library( Seurat )
library( Signac )
library( GENIE3 )
set.seed(123)

obj = readRDS('/local/projects/idea/sament2/cb/2020-10/age_trajectory/Purkinje/smoothed.rds')

atac.obj = readRDS('/local/projects/idea/sament2/cb/atac/peak_calling/HsCB.macs.integrated.clustered.rds')
DefaultAssay( atac.obj ) = 'peaks'

accessible = AccessiblePeaks( atac.obj , idents = 'InhNeu' , 
			     min.cells = 10 )

base.grn = readRDS('base.grn.rds')
filt.grn = base.grn[ which( base.grn$peak %in% accessible ) , ]

obj = NormalizeData( obj )
expr = as.matrix( GetAssayData( obj ))
expr = log2( expr + 1 )
 
targets = rownames(expr)
n = length( targets )
regulatorList <- list()
for( i in 1:n ) {
  base.regs = filt.grn$tf[ which( filt.grn$symbol == targets[i] ) ]
  regulatorList[[i]] = intersect( base.regs , targets )
}
names(regulatorList) = targets

weightList = list()
for( i in 1:n ) {
  cat( targets[i] , '\n' )
  regs = regulatorList[[i]]
  if( length(regs) < 2 ) next
  weights = GENIE3( expr , 
		      nCores=1,
		      targets= names(regulatorList)[i] , 
		      regulators=regulatorList[[i]] 
		    )
  weights = data.frame( tf = rownames(weights) ,
		        target = targets[i] ,
			weight = weights[,1] )
  weightList[[i]] = weights
}

weightList = do.call( rbind , weightList )

saveRDS( weightList , file = 'GENIE3.linkList.rds' )



