# this script is used to take results from k-means clustering and
# perform filtering and merging of the results to produce a final network

library( Seurat )
library( WGCNA )

# minimum size of a gene co-expression cluster
minModSize = 10
# minimum correlation between a gene and the cluster centroid
kAEtoStay = 0.2
# threshold for merging clusters with strongly correlated patterns
# the "height" corresponds to 1 - Pearson's r
cutHeight = 0.1

km = readRDS('kmeans.k25.rds')


smooth = readRDS('smoothed.rds')

smooth = NormalizeData(smooth)

datExpr0 = as.matrix( GetAssayData(smooth) )

# Z score
datExpr=sweep(datExpr0,1,apply(datExpr0,1,mean),"-")
indx.sd=(apply(datExpr0,1,sd))==0 # these will produce NAs
datExpr=sweep(datExpr,1,apply(datExpr,1,sd),"/")
datExpr[indx.sd,]=0
if(sum(is.na(datExpr))!=0){print("NAs in exprsMTX.Z Zscores!")}
#kmeans

kAE = cor( t(km$centers) , t(datExpr) )

colors = km$cluster

for( i in 1:length(colors) ) {
  r.i = kAE[ colors[i] , i ] 
  if( r.i < kAEtoStay ) colors[i] = 0
}

size = table( colors )
too.small = as.numeric(names(which(size<minModSize))) 
colors[ colors %in% too.small ] = 0


centers = sapply( sort(unique(colors)) , function(i) 
  colMeans(datExpr[ colors == i , ]) )

colnames(centers) = paste( 'AE' , sort(unique(colors)) , sep = '' )

r = cor( centers )

d = as.dist( 1 - r )

cl = cutree( hclust( d , method = 'average' ) , h = cutHeight )

mergeColors = rep(NA,length(colors))
for( i in 1:max(cl) ) {
  idx = as.numeric( gsub( 'AE','', names(cl)[ cl == i ] ))
  mergeColors[ colors %in% idx ] = i
}
mergeColors = mergeColors - 1
names(mergeColors) = names(colors)

MEs = moduleEigengenes( t(datExpr) , mergeColors )

saveRDS( mergeColors , file = 'k25.merged.clusters.rds' )

saveRDS( MEs , file = 'k25.MEs.rds' )







