library( Seurat )
library( Signac )
library( TFBSTools )
library(motifmatchr)
library(JASPAR2020)
library(BSgenome.Hsapiens.UCSC.hg38)


obj = readRDS('../peak_calling/HsCB.macs.integrated.clustered.rds')
DefaultAssay(obj) = 'peaks'

keep.peaks = rownames(obj)[
        grep( 'chr' , as.character(seqnames(granges(obj)) )) ]

obj = subset( obj , features = keep.peaks )

# extract position frequency matrices for the motifs

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", 
	      tax_group = 'vertebrates', 
	      all_versions = FALSE)
)

# add motif information

obj <- AddMotifs(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38 ,
  pfm = pfm
)

mtx = as.matrix(GetMotifData(obj))

motif.table = saveRDS( mtx , file = 'motif_instances_in_peaks.rds' )


