library( Seurat )

obj = readRDS('/local/projects/idea/sament2/cb/2020-10/age_trajectory/Purkinje/smoothed.rds')

obj = NormalizeData( obj )
expr = as.matrix( GetAssayData( obj ))
expr = log2( expr + 1 )

base.grn = readRDS('base.grn.rds')

links = readRDS('GENIE3.linkList.rds')

# relationship between links and activation / repression

datExpr = t(expr)

r = rep(NA,nrow(links))
for( i in 1:nrow(links) ) {
  r[i] = cor( datExpr[,links$tf[i]] , datExpr[,links$target[i]] )
}

links$r = r

# relationship between links and TFBSs

links$edge = paste( links$tf , links$target )

base.grn$edge = paste( base.grn$tf , base.grn$symbol )

nPeak.prox = table( base.grn$edge[ base.grn$coaccess == Inf ] )
nPeak.dist = table( base.grn$edge[ base.grn$coaccess != Inf ] )

prox = data.frame( edge = names(nPeak.prox) , freq = nPeak.prox )
colnames(prox) = c('edge','name','nProx')

dist = data.frame( edge = names(nPeak.dist) , freq = nPeak.dist )
colnames(dist) = c('edge','name','nDist')


links = merge( links , prox[,-2] , 
	        by = 'edge' , all.x = T ) 

links = merge( links , dist[,-2] ,  
                by = 'edge' , all.x = T )             

links[ is.na(links) ] = 0
links = links[ order( links$weight , decreasing = T ) , ]

strong = head( links , 50000 )
 
t.test( links$weight , links$nProx > 0 )
 
t = table( links$weight > 0.0347 , links$nProx > 0 )

t.test( links$weight , links$nDist > 0 )

kout = table( strong$tf )
kin = table( strong$target )

summary( as.numeric( kin ))

sort( kout )

# activators vs. repressors

tfs = unique( strong$tf )
r.thresh = 0.1
nPos = nNeg = nNS = pBinom = k = est = rep( NA , length(tfs) )
for( i in 1:length(tfs) ) {
  r.i = strong$r[ which( strong$tf == tfs[i] ) ]
  nPos[i] = length(r.i[ r.i > r.thresh ])
  nNeg[i] = length(r.i[ r.i < -1 * r.thresh ])
  nNS[i] = length(r.i[ abs(r.i) < r.thresh ])
  if( sum(nPos[i],nNeg[i]) == 0 ) next
  test = binom.test( x = nPos[i] , n = sum(nPos[i],nNeg[i]) )
  est[i] = test$estimate
  pBinom[i] = test$p.value
  k[i] = length(r.i)
}

tf.func = data.frame(
  tfs , k , nPos , nNeg , nNS , proportion = est , pBinom )

write.table( tf.func , row.names=F , quote=F , sep='\t' ,
	     file = 'tf.activation_vs_repression.txt' )

write.table( strong , row.names=F , quote=F , sep='\t' ,
	     file = 'Purkinje.grn.top_50k_links.txt' )

kout.tfs =
as.matrix( table( strong$tf[ strong$target %in% tfs ] ))


kin.tfs =
as.matrix( table( strong$target[ strong$target %in% tfs ] ))

df.k = merge( kout.tfs , kin.tfs , by = 0 )
colnames(df.k) = c('tf','kOut','kIn')



