library(JASPAR2020)
library(TFBSTools)
library( MotifDb )

motif.table = readRDS('motif_instances_in_peaks.rds')

human_tfs = readLines('TF_names_v_1.01.txt')

conns = readRDS('conns.downsampled.250k.rds')

peakAnno = read.delim('../annotate_peaks/peak_annotation.ChIPseeker.txt')

motifs = colnames(motif.table)
pfmList = getMatrixByID(JASPAR2020, ID=motifs )
n = length(pfmList)
remap_tf_name = symbol = motif.name = rep( NA , n )
for( i in 1:n ) {
  pfm = pfmList[[i]]
  remap = pfm@tags$remap_tf_name[1]
  if( length(remap) > 0 ) remap_tf_name[i] = remap
  sym = pfm@tags$symbol[1]
  if( length(sym) > 0 ) symbol[i] = sym
  motif.name[i] = pfm@name
}

df = data.frame( motif = motifs , motif.name , remap_tf_name , symbol )
df$tf_try = toupper( gsub( '\\((.*)' , '' , df$motif.name ) )
df$tf = df$tf_try
df$tf[ which( df$tf %in% human_tfs == F ) ] = NA

motifGenes = motifToGene( MotifDb , 
			  motifs = motifs , 
			  c('MotifDb') )

df2 = data.frame( motif = df$motif , 
		  geneSymbol = df$tf , pubmedID = NA ,
		  organism = NA , source = 'name' )

motifGenes = rbind( motifGenes , df2 )
motifGenes = motifGenes[ which( motifGenes$geneSymbol %in% human_tfs ) , ]


all.tfs = sort( unique( motifGenes$geneSymbol ) )
tf.peak = matrix( nrow = nrow(motif.table) , ncol = length(all.tfs) )
for( i in 1:length(all.tfs) ) {
  tf = all.tfs[i]
  motif.match = unique( motifGenes$motif[ 
		        motifGenes$geneSymbol == tf ] )
  idx = which( colnames(motif.table) %in% motif.match )
  sub = motif.table[ , idx ]
  if( length(idx) == 1 ) tf.peak[,i] = sub
  if( length(idx) > 1 ) tf.peak[,i] = as.numeric( rowSums(sub) > 0 )
}
rownames(tf.peak) = rownames(motif.table)
colnames(tf.peak) = all.tfs


summary( colSums( tf.peak ))
summary( rowSums( tf.peak ))
 
# saveRDS( tf.peak , 'predicted_tfbs_per_peak.rds' )
saveRDS( tf.peak , 'predicted_tfbs_per_peak.direct_matches.rds' )

