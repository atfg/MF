###############################################################################
# Angela's reminders for starting jobs, ignore this comment block: 
# bsub -R "rusage[mem=20GB]" -Is bash 
# module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.0.0; module switch R R/4.0.0; R
###############################################################################

# this loads up libraries, sets up some global varibles and functions 
# and loads up the seurat objects created by Luca so takes a few minutes to run
source("/icgc/dkfzlsdf/analysis/B210/angela/atfg_github/MF/initialise_env.R")
source("/icgc/dkfzlsdf/analysis/B210/angela/atfg_github/MF/functions.R")
source("/icgc/dkfzlsdf/analysis/B210/angela/atfg_github/MF/vagina_functions.R")

###############################################################################
# analyses
###############################################################################

# only the normal samples
samples = paste("N",1:5,sep="")

a = get_single_sample_vagina_seurats(c("N1","N2","N3","N4","N5"))
seu = a[[1]]
markers = a[[2]]
rm(a); gc()

a = data.frame( Embeddings(seu[["umap"]]), seu@meta.data ) 
pdf(file, width=2*6, height=6 )
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.1) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=type) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.1) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + facet_wrap( ~ type)
ggplot( a, aes(UMAP_1,UMAP_2, color=cell_type) ) + geom_point(shape=20) + theme_bw() + facet_wrap( ~ type)
ggplot( a, aes(UMAP_1,UMAP_2, color=li_celltype) ) + geom_point(shape=20) + theme_bw() + facet_wrap( ~ type)
dev.off()


###############################################################################
# integrate endometrium, vagina and mf
###############################################################################
seu.integrated = get_integrated_endometrium_vagina_mf()
markers = seu.integrated[[2]]
seu.integrated = seu.integrated[[1]]
seu.integrated = FindClusters(seu.integrated, verbose = TRUE, resolution=1)

a = data.frame( seu.integrated@meta.data ) 
a$type = "MF"
a$type[a$donor %in% c("57","19","63")] = "endometrium"
a$type[a$donor %in% paste("N",1:5,sep="")] = "vagina"
a$cell_type[a$li_celltype == "Fibroblasts"] = "Stromal fibroblasts"
a$cell_type[a$li_celltype %in% c("B cells","Plasma B cells","T cells")] = "Lymphocytes"
a$cell_type[a$li_celltype == "Macrophage"] = "Macrophages"
a$cell_type[a$li_celltype == "Epithelial"] = "Unciliated epithelia 1"
a$cell_type[a$li_celltype == "Endothelial cells"] = "Endothelia"

majority = table( a[,c("cell_type","integrated_snn_res.1")] )
majority = apply( majority, 2, function(x) { rownames(majority)[which(x == max(x))] } )
a$quake_majority = majority[as.character(a$integrated_snn_res.1)]

majority = table( a[,c("li_celltype","integrated_snn_res.1")] )
majority = apply( majority, 2, function(x) { rownames(majority)[which(x == max(x))] } )
a$li_majority = majority[as.character(a$integrated_snn_res.1)]

table(a$cell_type)
a$cell_type = sub(" cells","",a$cell_type)
a$cell_type = sub("Ciliated","Ciliated epithelial",a$cell_type)
a$cell_type = sub("Endothelia","Endothelial",a$cell_type)
a$cell_type = sub("Stromal fibroblasts","Stromal",a$cell_type)

seu.integrated@meta.data = a


a = data.frame( Embeddings(seu.integrated[["umap"]]), seu.integrated@meta.data ) 

file = paste("scrna_umaps_integrated_endometrium_vagina_mf.pdf", sep="" )
pdf(file, width=2*6, height=6 )
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.1) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=type) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.1) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + facet_wrap( ~ type)
ggplot( a, aes(UMAP_1,UMAP_2, color=cell_type) ) + geom_point(shape=20) + theme_bw() + facet_wrap( ~ type)
ggplot( a, aes(UMAP_1,UMAP_2, color=li_celltype) ) + geom_point(shape=20) + theme_bw() + facet_wrap( ~ type)
dev.off()


###############################################################################
# SVM_reject from sctransformed joined endometrium and vagina matrices
###############################################################################

library("e1071")
th = 0.7

qseu = get_single_sample_quake_seurats()
vseu = get_single_sample_vagina_seurats(c("N1","N2","N3"))[[1]]

qsamples=c("57")
vsamples=c("N1")

################################################
# training set
################################################
seu = get_joined_quake_vagina(qsamples,vsamples)
trdat = log(seu[["RNA"]]@data+1)
trdat = t(trdat)
trdat = as.data.frame(trdat)
trmeta = seu@meta.data

################################################
# test set (1 vagina, 1 quake, 1 vento-tormo )
################################################
tseu = get_single_sample_vagina_seurats("N2") 

tedat = tseu[[1]][["N2"]][["RNA"]]@data
tedat = log(tedat+1)
tedat = t(tedat)
tedat = as.data.frame(tedat)
tedat = tedat[,colnames(tedat)%in% colnames(trdat)] #genes
pred = predict(svmfit,tedat,probability=T) # gets it quite right!
table( tseu[[1]][["N2"]]@meta.data$cell_type, pred )

tseu = get_single_sample_quake_seurats()[["19"]] 
tedat = tseu[["RNA"]]@data
tedat = log(tedat+1)
tedat = t(tedat)
tedat = as.data.frame(tedat)
tedat = tedat[,colnames(tedat)%in% colnames(trdat)]
pred = predict(svmfit,tedat,probability=T) 
table( tseu@meta.data$cell_type, pred ) # doesn't work very well

tseu = get_single_sample_tormo_seurats()
tedat = tseu[["RNA"]]@data
tedat = log(tedat+1)
tedat = t(tedat)
tedat = as.data.frame(tedat)
tedat = tedat[,colnames(tedat)%in% colnames(trdat)]

################################################
# prediction set
################################################
samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020afM","MFCON007dfM") 
iseu = get_integrated_cellbender_seurat( samples ) # for the truth


sample = "MFCON007dcM"
tseu = get_single_sample_seurats( sample )[[1]][[1]]
tedat = tseu[["RNA"]]@data
tedat = log(tedat+1)
tedat = t(tedat)
tedat = as.data.frame(tedat)
tedat = tedat[,colnames(tedat)%in% colnames(trdat)] #genes
for( gene in colnames(trdat)[!(colnames(trdat) %in% colnames(tedat))] ) {
	tedat = cbind(tedat,0)
	colnames(tedat)[ncol(tedat)] = gene
}
tedat = tedat[,colnames(trdat)]
a = predict(svmfit,tedat,probability=T) 
pred = as.character(a)
pred[apply(attr(a, "probabilities"),1,max) < 0.7] = "unknown"
table(pred)

iseu[[1]]@meta.data$cell_type_level2[iseu[[1]]@meta.data$sample == sample]

################################################
# common set of genes between training, test and prediction sets
################################################
genes = rownames(trdat)[rownames(trdat) %in% rownames(mseu[[sample]][["RNA"]]@data)]

################################################
# train the model
################################################
svmfit = svm( trdat, factor(trmeta$cell_type), kernel = "linear", cost = 10, scale = FALSE, probability=T ) # [genes,]
save(svmfit,file=paste("svmfit_linear_endo",paste(qsamples,collapse="_"),"vagina", paste(vsamples,collapse="_"), ".RData", sep="_"))

################################################

# test data

# prediction dataset

#subset to common genes



save(svmfit,file="svmfit_linear_quake_19.RData")



head(attr(pred, "probabilities"))

a = as.character(pred)
a[apply(attr(pred, "probabilities"),1,max) < 0.7] = "unknown"
names(a) = names(pred)
table(a)


seu.integrated = get_integrated_mf_quake()
tab = seu.integrated@meta.data
tab$svm_end = NA
tab$svm_end = a[sapply(strsplit(rownames(tab),"_"),"[[",1)]



###########
seu = get_sctransform_endometrium_vagina_mf()
tab = data.frame( Embeddings(seu[["umap"]]), seu@meta.data ) 

pdf("sctransform_endometrium_vagina.pdf", width=6, height=6 )
ggplot(tab, aes(UMAP_1,UMAP_2, color=SCT_snn_res.0.2) ) + geom_point(shape=20) + theme_bw()
ggplot(tab, aes(UMAP_1,UMAP_2, color=type) ) + geom_point(shape=20) + theme_bw() 
ggplot(tab, aes(UMAP_1,UMAP_2, color=cell_type) ) + geom_point(shape=20) + theme_bw() 
ggplot(tab, aes(UMAP_1,UMAP_2, color=cell_type) ) + geom_point(shape=20) + theme_bw() + facet_wrap( ~ type)
dev.off()

###########

###############################################################################
# are there any tissue-specific clusters?
###############################################################################

table( a[,c("type","integrated_snn_res.1") ] )
unique( a[,c("integrated_snn_res.1","li_majority")])
# MF and endometrium do not have Plasma B cells and only small numbers of B cells than vagina
# MF and endometrium have much lower numbers of smooth muscle cells than vagina
# vagina has much fewer T cells than MF and endometrium

###############################################################################
# DE between fibroblasts in endometrium vs vagina
###############################################################################

smarkers = FindMarkers(seu.integrated, ident.1="endometrium", ident.2="vagina", group.by="type", subset.ident = c("0","2","3","4","8","10","12"), verbose = TRUE)
# % of cells the gene is expressed in
ind = a$cell_type %in% c("Ciliated epithelial","Unciliated epithelia 1","Unciliated epithelia 2","Epithelial") & a$type %in% "vagina"
b = rowSums(seu.integrated@assays$RNA@counts[,ind] > 2)
b1 = b/sum(ind)
ind = a$cell_type %in% c("Ciliated epithelial","Unciliated epithelia 1","Unciliated epithelia 2","Epithelial") & a$type %in% "endometrium"
b = rowSums(seu.integrated@assays$RNA@counts[,ind] > 2)
b2 = b/sum(ind)
pp = data.frame(vagina=b1,endometrium=b2)
pp = pp[rowSums(pp) > 0,]
pp = pp[order(apply(pp,1,max),decreasing=T),]
pp[ pp$endometrium == 0 & pp$vagina > 0.5,]

a = data.frame( Embeddings(seu.integrated[["umap"]]), seu.integrated@meta.data ) 
genes = rownames(pp[pp$endometrium == 0 & pp$vagina > 0.5,])
for( gene in genes ) {
	a[,gene] = seu.integrated@assays$RNA@counts[gene,]
}


file = paste("scrna_umaps_integrated_endometrium_vagina_mf_markers.pdf", sep="" )
pdf(file, width=2*6, height=6 )
for( gene in genes ) {
	p = ggplot( a[ a$integrated_snn_res.1 %in% c("0","2","3","4","8","10","12"),], aes(UMAP_1,UMAP_2, color=get(gene)) ) 
	p = p + geom_point(alpha=0.7, shape=20) + theme_bw() + facet_wrap( ~ type )  + scale_color_gradient(low = "blue", high = "red") + ggtitle(gene)
	print(p)
}
dev.off()

	
genes = rownames(pp[ pp$endometrium == 0 & pp$vagina > 0.5,])
ind = a$quake_majority %in% c("Ciliated epithelial","Unciliated epithelia 1","Unciliated epithelia 2","Epithelial") | a$li_celltype %in% c("Epithelial cells")
b = log10(as.matrix(seu.integrated@assays$RNA@counts[genes,ind])+1)
pc = prcomp(t(b))
pc = data.frame(pc$x)
pc = cbind(pc,a[rownames(pc),])

library("gplots")
#b = log10(as.matrix(seu.integrated@assays$RNA@counts[genes,a$cell_type %in% c("Ciliated epithelial","Unciliated epithelia 1","Unciliated epithelia 2","Epithelial")])+1)
pdf("scrnaintegrated_endometrium_vagina_mf_pca.pdf", width=6, height=6 )
ggplot(pc,aes(PC1,PC2,color=type)) + geom_point()
dev.off()


