
tf.peak = readRDS('predicted_tfbs_per_peak.direct_matches.rds' )

conns0 = readRDS('../cicero/conns.downsampled.250k.rds')
conns = conns0[ which( conns0$coaccess > 0.5 ) , ]

peakAnno = read.delim('../annotate_peaks/peak_annotation.ChIPseeker.txt')

prox = peakAnno[ grep('Promoter',peakAnno$annotation) , 
		 c('peak','SYMBOL') ]
prox$coaccess = Inf

a = merge( conns , peakAnno , by = 1 )
dist = a[ , c('Peak2','SYMBOL','coaccess') ]

colnames(prox) = colnames(dist) = c('peak','symbol','coaccess')

peak.target = rbind( prox , dist )
peak.target = peak.target[ 
  order( peak.target$coaccess , decreasing = T ) , ]
peak.target = peak.target[ is.na( peak.target$symbol ) == F , ]
peak.target = peak.target[ duplicated( peak.target[,1:2] ) == F , ]

peak.target.anno = merge( peak.target , peakAnno , by = 1 )

##

peaks = rownames(tf.peak)
tfs = colnames(tf.peak)
tf.peak.list = list()
for( i in 1:ncol(tf.peak) ) {
  tmp = tf.peak[,i]
  df = data.frame( tf = tfs[i] ,
		   peak = peaks[ which(tmp>0) ] )
  tf.peak.list[[i]] = df
}

tf.peak.links = do.call( rbind , tf.peak.list )

base.grn = merge( tf.peak.links , peak.target , by.x = 2 , by.y = 1 )
base.uni = base.grn[ duplicated( base.grn[,2:3] ) == F , 2:3 ]
k.in = table( base.uni$symbol )
k.out = table( base.uni$tf )
summary( as.numeric( k.in ) )
summary( as.numeric( k.out ))

saveRDS( base.grn , file = 'base.grn.rds' )




