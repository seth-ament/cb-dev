
colors = readRDS('k25.merged.clusters.rds')
 
sets = readRDS('../../genelist/Combined_CB_Genelists_2022-03-12.rds')

mods = 1:22
n = length(sets)
m = length(mods)

u = names(colors)
setname = module = p = o = or = setsize = modsize = matrix( NA , nrow = n , ncol = m )
for( i in 1:n ) {
  for ( j in 1:m ) {
    mod = names(colors)[ colors == j ]
    set = sets[[i]]
    t = table( u %in% mod , u %in% set )
    test = fisher.test( t , alternative = 'greater' )
    p[i,j] = test$p.value
    or[i,j] = test$estimate
    o[i,j] = t[2,2]
    setsize[i,j] = length(set)
    modsize[i,j] = length(mod)
    setname[i,j] = names(sets)[i]
    module[i,j] = paste( 'M' , j , sep = '' )
  }
}

fdr = apply( p , 2 , p.adjust , method = 'BH' )

df = data.frame( 
  module = as.vector(module) ,
  geneset = as.vector(setname) ,
  nMod = as.vector(modsize) ,
  nSet = as.vector(setsize) ,
  nOverlap = as.vector(o) ,
  OR = as.vector(or) ,
  P = as.vector(p) ,
  FDR = as.vector( fdr )
)
df = df[ order( df$P ) , ]
df$bonf = p.adjust( df$P , method = 'bonferroni' )

enrich = df

write.table( enrich , sep='\t' , quote=F , row.names=F ,
	     file = 'Curated_gene_sets.enrichment_for_Purkinje_modules.txt' )






