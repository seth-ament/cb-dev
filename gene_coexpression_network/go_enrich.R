options(stringsAsFactors=F)
library( GO.db )

colors = readRDS('k25.merged.clusters.rds')

sets = readRDS('og.Hs.egGO2ALLEGS.hgnc_symbol.2022-04-11.rds')

setsize = sapply( 1:length(sets) ,
  function(i) length(sets[[i]]) )
sets = sets[ setsize >= 5 ]
setsize = setsize[ setsize >= 5 ]

goterms <- Term(GOTERM)

u = names(colors)
 
terms = names(sets)
nTerm = length(terms)
mods = sort( unique( colors ) )
nMod = length(unique(colors))
term = mod = o = or = p = matrix( NA , ncol = nMod , nrow = nTerm )
for( i in 1:nTerm ) {
  for( j in 1:nMod ) {
    set1 = names(colors)[ colors == mods[j] ]
    set2 = sets[[i]]
    t = table( u %in% set1 , u %in% set2 )
    if( nrow(t) != 2 | ncol(t) != 2 ) next
    test = fisher.test( t , alternative = 'greater' )
    o[i,j] = t[2,2]
    or[i,j] = test$estimate
    p[i,j] = test$p.value
    mod[i,j] = mods[j]
    term[i,j] = terms[i]
  }
}
colnames(o) = colnames(or) = colnames(p) = paste('M',mods,sep='')
rownames(o) = rownames(or) = rownames(p) = terms

q = apply( p , 2 , p.adjust , method = 'BH' )

df = data.frame(
  term = as.vector( term ) ,
  mod = as.vector( mod ) , 
  o = as.vector( o ) ,
  or = as.vector( or ) , 
  p = as.vector( p ) ,
  fdr = as.vector( q )
)
df = df[ is.na( df$mod ) == F , ]

modsize = rep(NA,nMod)
for( i in 1:nMod ) {
  set1 = intersect( u , names(colors)[ colors == mods[i] ] )
  modsize[i] = length(set1)
}
names(modsize) = mods

df$nMod = modsize[ as.character(df$mod) ]

termsize = rep(NA,nTerm)
for( i in 1:nTerm ) {
   set2 = intersect( u , sets[[i]] )
   termsize[i] = length(set2)
}
names(termsize) = names(sets)

df$nTerm = termsize[ df$term ]

df$description = goterms[ df$term ]

df = df[ order( df$p ) , ]

df$bonf = p.adjust( df$p , method = 'bonferroni' )

sig = df[ df$p < 0.05 , ]
 
write.table( sig , row.names=F , quote=F , sep='\t' ,
  file = 'Purkinje.merged_clusters.go_enrich.txt' )






