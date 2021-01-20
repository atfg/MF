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
seu.integrated = FindClusters(seu.integrated, verbose = TRUE, resolution=1)

file = paste("scrna_umaps_integrated_endometrium_vagina_mf.pdf", sep="" )
a = data.frame( Embeddings(seu.integrated[["umap"]]), seu.integrated@meta.data ) 
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
a$cell_type_majority = 
		
b = factor(as.character(a$integrated_snn_res.1))
levels(b) = list(majority)
		
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

################################################
# training set
################################################
seu = get_joined_quake_vagina( qsamples=c("57"), vsamples=c("N1"))

################################################
# test set (1 vagina, 1 quake, 1 vento-tormo )
################################################


################################################
# prediction set
################################################
sample = "MFCON007dcM"
mseu = get_single_sample_cellbender_seurats( sample )[[1]]

################################################
# common set of genes between training, test and prediction sets
################################################
genes = rownames(trdat)[rownames(trdat) %in% rownames(mseu[[sample]][["RNA"]]@data)]

################################################
# train the model
################################################
trdat = log(seu[["RNA"]]@data[genes,]+1)
trdat = t(trdat)
trdat = as.data.frame(trdat)
trmeta = seu@meta.data

svmfit = svm( trdat, trmeta, kernel = "linear", cost = 10, scale = FALSE, probability=T)


################################################

# test data

# prediction dataset
tedat = mseu[["N1"]][["RNA"]]@data
tedat = log(tedat+1)
tedat = t(tedat)
tedat = as.data.frame(tedat)
#subset to common genes
tedat = tedat[,genes]


save(svmfit,file="svmfit_linear_quake_19.RData")


pred = predict(svmfit,tedat,probability=T)

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
