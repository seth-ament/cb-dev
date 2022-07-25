library( Seurat )
library( limma )
type = 'Purkinje'

dirs = Sys.glob('*/')
types = gsub('/','',dirs)
types = setdiff( types , c('nobatchcorr','Choroid_Plexus','Lymphoid_T_cells','1000umi') )

# linear modeling of eigengenes

for( i in 1:length(types) ) {
  type = types[i]
  alt_type = gsub('_',' ',type)
  cat(type,'\n')
  setwd( type )

  obj = readRDS( paste( alt_type , '.ComBat.Smoothed.rds' , sep = '' ))
  net = readRDS( paste( type , '.WGCNA.net.rds' , sep='' ))
  MEs = net$MEs
  obj@meta.data = cbind( obj@meta.data , MEs )

  cells = which( obj$Age_days < 365*10 )
  noadult = subset(obj,cells = cells )

  m = noadult@meta.data
 
  p = beta = matrix( ncol = 4 , nrow = max(net$colors) )
  for( i in 1:max(net$colors) ) {
    y = MEs[ cells , paste('ME',i,sep='') ]
    fit = lm( y ~ Group + Age_days + Sex , data = m )
    p[i,] = drop1( fit , ~. , test='F' )$'Pr(>F)'
    beta[i,] = coef(fit)
  }
  rownames(p) = rownames(beta) = paste('M',1:max(net$colors),sep='')
  colnames(p) = colnames(beta) = c('Intercept','Group','Age','Sex')

  m.all = obj@meta.data
  m.all$RIN[ is.na(m.all$RIN) ] = mean( m.all$RIN , na.rm = T )
  p.all = beta.all = matrix( ncol = 4 , nrow = max(net$colors) )
  for( i in 1:max(net$colors) ) {
    y = MEs[ , paste('ME',i,sep='') ]
    fit = lm( y ~ Group + Age_days + Sex , data = m.all )
    p.all[i,] = drop1( fit , ~. , test='F' )$'Pr(>F)'
    beta.all[i,] = coef(fit)
  }
  rownames(p.all) = rownames(beta.all) = paste('M',1:max(net$colors),sep='')
  colnames(p.all) = colnames(beta.all) = c('Intercept','Group','Age','Sex')

  df = data.frame(
    module = rownames(p) ,
    size = as.vector(table( net$colors ))[-1] ,
    beta.infl = beta[,2] , 
    p.infl = p[,2] ,
    beta.age_noadult = beta[,3] ,
    p.age_noadult = p[,3] ,
    beta.sex = beta[,4] ,
    p.sex = p[,4] ,
    beta.age_all = beta.all[,3] ,
    p.age_all = p.all[,3] )

  write.table( df , quote=F , row.names=F , sep='\t' ,
    file = 'ME.lm.pvals.txt' )

  setwd('..')

}


# relate modules to traits with FRY
library( limma )
for( i in 1:length(types) ) {
  type = types[i]
  alt_type = gsub('_',' ',type)
  cat(type,'\n')
  setwd( type )
  obj = readRDS( paste( alt_type , '.ComBat.Smoothed.rds' , sep = '' ))
  net = readRDS( paste( type , '.WGCNA.net.rds' , sep='' ))
  MEs = net$MEs
  obj@meta.data = cbind( obj@meta.data , MEs )
  cells = which( obj$Age_days < 365*10 )
  noadult = subset(obj,cells = cells )
  y = as.matrix( GetAssayData( noadult ) )
  idx = list()
  for( k in 1:max(net$colors) ) {
    idx[[k]] = which( net$colors == k )
  }
  names(idx) = paste( 'M' , 1:max(net$colors) , sep = '' )
  m = noadult@meta.data
  design = model.matrix(~ m$Group + m$Age_days +m$Sex +m$batch +m$PMI +m$RIN)
  colnames(design) = gsub('m\\$','',colnames(design))
  fry.infl = fry( y = y , index = idx , design = design , contrast = 2 )
  fry.age = fry( y = y , index = idx , design = design , contrast = 3 )
  m2 = obj@meta.data
  m2$RIN[ is.na(m2$RIN) ] = mean(m2$RIN,na.rm=T)
  y2 = as.matrix( GetAssayData( obj))
  d2 = model.matrix(~ m2$Age_days +m2$Group+m2$Sex +m2$batch+m2$PMI +m2$RIN)
  colnames(d2) = gsub('m\\$','',colnames(d2))
  fry.age_all = fry( y = y2 , index = idx , design = d2 , contrast = 2 )
  fry = merge( fry.infl[,1:3] , fry.age[,2:3] , by = 0 )
  fry = merge( fry , fry.age_all[,2:3] , by.x = 1 , by.y = 0 )
  colnames(fry) = c('Module','Size','Direction.Infl','PValue.Infl',
    'Direction.Age_noadult','PValue.Age_noadult',
    'Direction.Age_all','PValue.Age_all')
  write.table( fry , sep='\t' , row.names = F , quote=F ,
    file = 'fry.wgcna_modules.txt' )
  setwd('..')
}






meta = noadult@meta.data
pvals = coefs = matrix( NA , nrow = max(net$colors) , ncol = 6 )
for( i in 1:max(net$colors) ) {
  yvar = paste( 'ME' , i , sep = '' )
  fit = lm( meta[,yvar] ~ 
    Group + Age_days + Sex + RIN + PMI , data = meta )
  drop = drop1( fit , ~. , test='F' )
  pvals[i,] = drop[,'Pr(>F)']
  coefs[i,] = coef(fit)
}
colnames(pvals) = colnames(coefs) = names(coef(fit))
rownames(pvals) = rownames(coefs) = paste('ME',1:max(net$colors),sep='')




