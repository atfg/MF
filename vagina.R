###############################################################################
# Angela's reminders for starting jobs, ignore this comment block: 
# bsub -R "rusage[mem=20GB]" -Is bash 
# module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.0.0; module switch R R/4.0.0; R
###############################################################################

# this loads up libraries, sets up some global varibles and functions 
# and loads up the seurat objects created by Luca so takes a few minutes to run
source("/icgc/dkfzlsdf/analysis/B210/angela/atfg_github/MF/initialise_env.R")
source("/icgc/dkfzlsdf/analysis/B210/angela/atfg_github/MF/functions.R")

###############################################################################
# analyses
###############################################################################

# only the normal samples
samples = paste("N",1:5,sep="")

adat = list()
markers = list()
for( sample in samples ) {
	dat = Read10X(data.dir = paste("/icgc/dkfzlsdf/analysis/B210/references_data/2020_li/", sample, "/", sep="") )	
	adat[[sample]] = CreateSeuratObject(counts = dat)
}

a = get_single_sample_vagina_seurats()
adat = a[[1]]
markers = a[[2]]
rm(a); gc()




############

th = 0.7

library("e1071")

qsamples = c("57","19","63")
qseu = get_single_sample_human_data_quake()
qseu = qseu[qsamples]

samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON020afM","MFCON007dfM") 
mseu = get_single_sample_cellbender_seurats( samples )[[1]]

trdat = qseu[["19"]][["RNA"]]@data
trdat = log(trdat+1)
trdat = t(trdat)
trdat = as.data.frame(trdat)

tedat = mseu[[1]][["RNA"]]@data
tedat = log(tedat+1)
tedat = t(tedat)
tedat = as.data.frame(tedat)

#subset to common genes
genes = colnames(trdat)[colnames(trdat) %in% colnames(tedat)]
trdat = trdat[,genes]
tedat = tedat[,genes]

svmfit = svm( trdat, factor(qseu[["19"]]@meta.data$cell_type), kernel = "linear", cost = 10, scale = FALSE,probability=T)
save(svmfit,file="svmfit_linear_quake_19.RData")

print(svmfit)


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





