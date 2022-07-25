library( Seurat )

obj = readRDS( 'PCs.SeuratObj.withMEs.rds' )

cells = colnames(obj)[ obj$group == 'ctrl' ]
ctrl = subset( obj , cells = cells )

pdf('violins_by_age_and_subtype.ctrl_only.top_tfs.pdf',width=7,height=7)
par( mfrow = c(1,6) )
VlnPlot( ctrl , group.by = 'age.subtype' ,
        features = c('FOXP1','ZNF423','EBF3','TCF12','TEAD1','EBF1') ,
        combine = T )
dev.off()



####


subtypes = readRDS('../../fig1/PC.subtypes.seuratObj.rds')
ids = c('ALDOC_1','Non-ALDOC_1','ALDOC_2','Non-ALDOC_2')
names(ids) = 0:3
subtypes = RenameIdents( subtypes , ids )
subtypes = Idents(subtypes)
subtypes = subtypes[ names(subtypes) %in% colnames(obj) ]

subtype = rep( 'Fetal' , ncol(obj) )
names(subtype) = colnames(obj)
subtype[ names(subtypes) ] = as.character(subtypes)
obj$subtype = subtype

Idents(obj) = factor( obj$age , levels = sort( unique( obj$age )) )

obj@meta.data = cbind( obj@meta.data , MEs$eigengenes )

obj$simple_subtype = gsub( '_(.*)','',obj$subtype )
obj$simple_subtype = factor( obj$simple_subtype , levels = c('Fetal','Non-ALDOC','ALDOC') )

cells = colnames(obj)[ obj$age == (16-40)/52 ]
obj = subset( obj , cells  = cells , invert = T )
 
age.bin = rep(NA,ncol(obj))
age.bin[ obj$age >= (9-40)/52 & obj$age <= (10-40)/52 ] = '9-10 pcw'
# age.bin[ obj$age == (10-40)/52 ] = '10 pcw'
age.bin[ obj$age >= (11-40)/52 & obj$age <= (12-40)/52 ] = '11-12 pcw'
# age.bin[ obj$age == (12-40)/52 ] = '12 pcw'
age.bin[ obj$age >= (14-40)/52 & obj$age <= (17-40)/52 ] = '14-17 pcw'
# age.bin[ obj$age == (16-40)/52 ] = '16 pcw'
age.bin[ obj$age >= (18-40)/52 & obj$age <= (-20/52) ] = '18-20 pcw'
# age.bin[ obj$age == (18-40)/52 ] = '18 pcw'
# age.bin[ obj$age == (20-40)/52 ] = '20 pcw'
age.bin[ obj$age > 1 & obj$age < 2.1 ] = '1 yr'
age.bin[ obj$age > 2.1 & obj$age < 7 ] = '2-5 yrs'
age.bin[ obj$age > 7 ] = 'Adult'

obj$age.bin = factor( age.bin , levels = c(
  '9-10 pcw','11-12 pcw','14-17 pcw',
  '18-20 pcw','1 yr','2-5 yrs','Adult' ) )

obj$age.subtype = paste( obj$age.bin , obj$simple_subtype , sep = ' ' )
obj$age.subtype = gsub( ' Fetal' , '' , obj$age.subtype )
obj$age.subtype = factor( obj$age.subtype , levels = c(
  '9-10 pcw','11-12 pcw','14-17 pcw',
  '18-20 pcw',
  '1 yr Non-ALDOC','1 yr ALDOC',
  '2-5 yrs Non-ALDOC','2-5 yrs ALDOC',
  'Adult Non-ALDOC','Adult ALDOC') )
						     
obj$sub.group = paste( obj$simple_subtype , obj$group , sep = '_' )
obj$sub.group[ obj$sub.group == 'Fetal_ctrl' ] = 'Fetal'
obj$sub.group = factor( obj$sub.group , levels = c(
  'Fetal','ALDOC_ctrl','ALDOC_infl','Non-ALDOC_ctrl','Non-ALDOC_infl') )

saveRDS( obj , file = 'PCs.SeuratObj.withMEs.rds' )

cells = colnames(obj)[ obj$group == 'ctrl' ]
ctrl = subset( obj , cells = cells )

pdf('violins_by_age_and_subtype.ctrl_only.top_tfs.pdf',width=7,height=3)
VlnPlot( ctrl , group.by = 'age.subtype' , 
	features = c('FOXP1','FOXP2','RORA','KLF9','EBF1','EBF3','TCF12','KLF12') ,
	combine = F )	
dev.off()

cells = colnames(obj)[ obj$age > 0 & obj$age < 10 ]
child = subset( obj , cells = cells )

pdf('violin_by_infl.pdf',height=2.25,width=2.25)
VlnPlot( child , group.by = 'group' , 
	features = paste( 'ME' , 1:max(MEs$validColors) , sep = '' ) , combine = F )
dev.off()

## 

x = child@meta.data

fit = lm( ME7 ~ group * simple_subtype + log_yrs , data = x )
drop1( fit , ~. , test = 'F' )

t = table( x$group , as.character(x$simple_subtype) )
fisher.test( t )

subtype.counts = as.matrix( table( child$donor , child$simple_subtype ) )[,-1]

samples = unique( child@meta.data[,c('sample_name','log_yrs','group')] )
rownames(samples) = samples[,1]
samples = samples[ rownames(subtype.counts) , ]

dat = data.frame( samples , NonALDOC = subtype.counts[,1] ,
		ALDOC_1 = subtype.counts[,2] , ALDOC_2 = subtype.counts[,3] )

prop = cbind( dat$NonALDOC , rowSums(dat[,c('ALDOC_1','ALDOC_2')]) )
da.fit = glm( formula =  prop ~  dat$log_yrs + dat$group , family = binomial )
drop1( da.fit , ~. , test = 'Chisq' )

n = max( MEs$validColors )
infl.de.pvals = matrix( NA , nrow = n , ncol = 4 )
for( i in 1:n ) {
  y = x[,paste('ME',i,sep='')]
  fit = lm( y ~ group * simple_subtype + log_yrs , data = x )
  infl.de.pvals[i,] = drop1( fit ,  ~. , test = 'F' )$'Pr(>F)'[-1]
}
infl.de.pvals = as.data.frame( infl.de.pvals )
colnames(infl.de.pvals) = c('group','subtype','log_yrs','groupXtype')
rownames(infl.de.pvals) = paste( 'ME' , 1:n , sep = '' )
infl.de.pvals$FDR.group = p.adjust( infl.de.pvals$group )
infl.de.pvals$FDR.groupXtype = p.adjust( infl.de.pvals$groupXtype )
infl.de.pvals[ infl.de.pvals$FDR.groupXtype < 0.1 | infl.de.pvals$FDR.group < 0.1 , ]
 
##

options(stringsAsFactors=F)
  
anno = read.delim('/local/projects/idea/sament2/cb/2020-10/fig2/BioMart_symbol_biotype_2022-03-12.txt')
protein_coding = unique( anno$Gene.name[ anno$Source.of.gene.name == 'HGNC Symbol' &
                                   anno$Gene.type == 'protein_coding' ] )

 
degs = read.delim('/local/projects/idea/sament2/cb/2020-10/inflammation_degs/HuCB.AllTypes.Inflammation.DGE.2021-06-28.txt')
types = sort( unique( degs$celltype ) )

filt = degs[ which( degs$FDR.infection < 0.01 &
             abs(degs$logFC.infection) > log(1.5) &
             abs(degs$logFC.sexM) < log(2) &
             ( degs$pct.infection > 0.2 | degs$pct.control > 0.2 )  &
             degs$gene %in% protein_coding )
        , ]

down = filt$gene[ filt$celltype == 'Purkinje' & filt$logFC.infection < 0 ]
up = filt$gene[ filt$celltype == 'Purkinje' & filt$logFC.infection > 0 ]

colors = MEs$validColors
o = p = or = rep(NA,max(colors))
u = names(colors)
for( i in 1:max(colors) ) {
  mod = names(colors)[ colors == i ]
  t = table( u %in% mod , u %in% down )
  test = fisher.test( t , alternative = 'greater' )
  or[i] = test$estimate
  p[i] = test$p.value
  o[i] = t[2,2]
}
df.down = data.frame( Module = paste( 'ME' , 1:22 , sep = '' ) , 
		      Direction = 'Down' , o , or , p ) 

o = p = or = rep(NA,max(colors))
u = names(colors)
for( i in 1:max(colors) ) {
  mod = names(colors)[ colors == i ]
  t = table( u %in% mod , u %in% up )
  test = fisher.test( t , alternative = 'greater' )
  or[i] = test$estimate
  p[i] = test$p.value
  o[i] = t[2,2]
} 
df.up = data.frame( Module = paste( 'ME' , 1:22 , sep = '' ) , 
		    Direction = 'Up' , o , or , p )

df = rbind( df.down , df.up )
df = df[ order( df$p ) , ]
df$FDR = p.adjust( df$p )

# Module Direction   o        or            p          FDR
#   ME18        Up  59 26.224945 9.723473e-51 4.278328e-49
#   ME15        Up  58  9.973911 1.166024e-29 5.013904e-28
#    ME7      Down 160  2.823723 2.685987e-23 1.128115e-21
#    ME2      Down  52  3.417163 4.348481e-12 1.782877e-10
#   ME12      Down 144  2.055128 5.658833e-12 2.263533e-10
#   ME14      Down  52  2.821556 1.649457e-09 6.432884e-08
#   ME20      Down  59  2.006517 5.024068e-06 1.909146e-04




