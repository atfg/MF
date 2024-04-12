###############################################################################
# Angela's reminders for starting jobs, ignore this comment block: 
# bsub -R "rusage[mem=20GB]" -Is bash 
# module unload R; module unload python; module unload gcc; module load anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 gcc/11.1.0 R/4.3.1 cmake/3.21.0; R
## module unload R; module unload python; module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.1.0 gcc/7.2.0; R
###############################################################################

source("/omics/groups/OE0433/internal//angela/atfg_github/MF/initialise_env.R")


###############################################################################
# Deconvolve EPFL data using our samples or Quake's
###############################################################################

library("DeconRNASeq")
library("org.Hs.eg.db")
library("DESeq2")

dat = get_endoCounts(normalized=TRUE)
tab = dat$counts
meta = dat$meta
rm(dat)
gc()

# read in reference counts
load("mf_counts_deconvolute.RData")
ref = sdat[,grepl("MFCON007dcM",colnames(sdat)) & !grepl("endothelial",colnames(sdat))]


sref = ref[rownames(stab),] 

de = DESeqDataSetFromMatrix(countData = cbind(stab,sref), colData = data.frame(exp=c(rep("endo",ncol(stab)),rep("mf",ncol(sref)))), design= ~ exp)
de = estimateSizeFactors(de)
ntab = counts(de, normalized=TRUE)


res = DeconRNASeq(as.data.frame(ntab[rownames(ntab)%in%unique(unlist(goi)),colnames(stab)]), as.data.frame(ntab[rownames(ntab)%in%unique(unlist(goi)),colnames(sref)]))
a = t(res$out.all)
colnames(a) = colnames(stab)

b = ntab[,colnames(sref)]
b = data.frame(GeneSymbol=rownames(b),b)
write.table(b,file="temp_ref.txt",sep="\t",quote=F,row.names=F)      

b = ntab[,colnames(stab)]
b = data.frame(GeneSymbol=rownames(b),b)
write.table(b,file="temp_mix.txt",sep="\t",quote=F,row.names=F)      


##

samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020afM","MFCON007dfM") 
a = get_integrated_cellbender_seurat( samples, update=F )
seu.integrated = a[[1]]
imarkers = a[[2]]
rm(a)
gc()

a = data.frame( Embeddings(seu.integrated[["umap"]]), seu.integrated@meta.data ) 
a$type = "biopsy"
a$type[is.na( a$donor)] = "mf"
b = rownames(seu.integrated@meta.data)
nn = b
nn[grepl("_6",b)] = sub("_6","_1",b[grepl("_6",b)])
nn[grepl("_7",b)] = sub("_7","_2",b[grepl("_7",b)])
nn[grepl("_8",b)] = sub("_8","_4",b[grepl("_8",b)])
nn[grepl("_9",b)] = sub("_9","_5",b[grepl("_9",b)])
nn[grepl("_4",b)] = sub("_4","_6",b[grepl("_4",b)])
nn[grepl("_5",b)] = sub("_5","_7",b[grepl("_5",b)])
rownames(a) = nn
mfm = iseu@meta.data
mfm = mfm[ mfm$sample != "MFCON018bfM",]
a$ct_qu = a$cell_type
a$ct_mf[match(rownames(mfm),rownames(a))] = mfm$cell_type_level2
#a$cell_type[match(rownames(mfm),rownames(a))] = mfm$cell_type_level1
a$cell_type[match(rownames(mfm),rownames(a))] = mfm$cell_type_level2
a$donor[match(rownames(mfm),rownames(a))] = mfm$sample
a$cell_type_level1[grepl("Ciliated|cilliated|epithelial|epithelia|glands",a$cell_type)] = "epithelial"
a$cell_type_level1[grepl("Stroma|stroma",a$cell_type)] = "stromal"
a$cell_type_level1[grepl("Lymphocytes|Macrophages|T-cells",a$cell_type)] = "leukocytes"
a$cell_type_level1[grepl("Endothelia|muscle",a$cell_type)] = "endothelial"

seu.integrated@meta.data = a


dat = seu.integrated[["RNA"]]@counts
ndat = c()
a = seu.integrated@meta.data
for( ct in unique(a$cell_type_level1) ) {
	for(donor in unique(a$donor) ) {
		ndat = cbind(ndat, rowSums(dat[,a$cell_type_level1 == ct & a$donor == donor]) )
		colnames(ndat)[ncol(ndat)] = paste(donor,ct,sep="_")
	}
}
sdat = ndat[,grepl("MF",colnames(ndat))]
sdat = sdat[rowSums(sdat)>0,]


##############

save(sdat,file="mf_counts_deconvolute.RData")

load("mf_counts_deconvolute.RData")
de = DESeqDataSetFromMatrix(countData = sdat, colData = data.frame(cell_type=sapply(strsplit(colnames(sdat),"_"),"[[",2)), design= ~ cell_type)
de = estimateSizeFactors(de)
sdat = counts(de, normalized=TRUE)


res = DeconRNASeq(as.data.frame(ndat), as.data.frame(ref[,grepl("MFCON007dcM",colnames(ref))]))
a = t(res$out.all)
head(a)


###############################################################################
# Deconvolve EPFL data using Gregersen data
#
# module unload R; module unload python; module load anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.5.1 geos/3.11.0 proj/6.2.1 R/4.3.1 gcc/7.2.0 sqlite/3.38.5; R
# module load gcc/11.1.0 sqlite/3.38.5 gdal/3.5.1 geos/3.11.0 proj/6.2.1 R/4.3.1; R
.libPaths(paste("/omics/groups/OE0433/internal/software/R-library-",paste(R.version$major,R.version$minor,sep="."),sep=""))
install.packages("sf")
packageUrl <- "https://cran.r-project.org/src/contrib/Archive/sf/sf_0.9-0.tar.gz" 
install.packages(packageUrl, repos=NULL, type = "source")

###############################################################################

library("DeconRNASeq")

# make reference from gregersen control samples
meta = get_metadata_gregersen()
seu = get_single_sample_gregersen_seurats( update=F )

cts = unique(unlist(sapply( seu, function(x) { unique(x@meta.data$clusterID) } )))
cts = cts[!is.na(cts)]
cts = cts[cts != "Unknown"]
scounts = list()
for( sample in names(seu) ) {
	cat("at",sample,"\n")
	a = seu[[sample]]@meta.data
	scounts[[sample]] = c()
	for( ct in cts ) {
		cat("\tat",ct,"\n")
		ind = a$clusterID %in% ct & a$pheno %in% "Control"
		if( sum(ind) > 1 ) {
			scounts[[sample]] = cbind(scounts[[sample]], rowSums(seu[[sample]][["RNA"]]@counts[,ind]) )
		} else if( sum(ind) == 1 ) {
			scounts[[sample]] = cbind(scounts[[sample]], seu[[sample]][["RNA"]]@counts[,ind] )
		} else {
			scounts[[sample]] = cbind(scounts[[sample]], rep(0,nrow(seu[[sample]][["RNA"]]@counts)) )
		}
		colnames(scounts[[sample]])[ncol(scounts[[sample]])] = ct
	}
	rownames(scounts[[sample]]) = rownames(seu[[sample]][["RNA"]]@counts)
}

acounts = scounts[[1]]
for( sample in names(scounts)[2:length(names(scounts))] ) {
	acounts = acounts + scounts[[sample]]
}
meta = get_metadata_gregersen()
a = table(meta$clusterID)
a = a[colnames(acounts)]
#for( i in 1:ncol(acounts) ) {
#	acounts[,i] = acounts[,i]/a[i]
#}
de = DESeqDataSetFromMatrix(countData = acounts, colData = data.frame(cell_type=colnames(acounts)), design= ~ cell_type)
de = estimateSizeFactors(de)
sdat = counts(de, normalized=TRUE)


# get endo data
dat = get_endoCounts(normalized=TRUE)
tab = dat$counts
meta = dat$meta
rm(dat)
gc()

target = as.data.frame(tab)
reference = as.data.frame(sdat)
reference = reference[rowSums(reference) > 0,]
genes = rownames(target)[rownames(target) %in% rownames(reference)]
target = target[genes,]
reference = reference[genes,]

write.table(reference, file="endo_deconvolute_reference.tsv", sep="\t")

# deconvolute
res = DeconRNASeq(target, reference)
a = t(res$out.all)
colnames(a) = colnames(tab)
a[,1:3]

########### normalising together

acounts = scounts[[1]]
for( sample in names(scounts)[2:length(names(scounts))] ) {
	acounts = acounts + scounts[[sample]]
}

dat = get_endoCounts(normalized=FALSE)
tab = dat$counts

genes = rownames(acounts)[rownames(acounts) %in% rownames(reference)]
jtab = cbind(acounts[genes,],tab[genes,])
de = DESeqDataSetFromMatrix(countData = jtab, colData = data.frame(type=c(rep("ref",ncol(acounts)),rep("targ",ncol(tab)))), design= ~ type)
de = estimateSizeFactors(de)
sdat = counts(de, normalized=TRUE)

target = as.data.frame(sdat[,colnames(tab)])
reference = as.data.frame(sdat[,colnames(acounts)])

res = DeconRNASeq(target, reference)
a = t(res$out.all)
colnames(a) = colnames(tab)
a[,1:7]

###############################################################################
# try music
###############################################################################

library("MuSiC")
library("SingleCellExperiment")

# load bulk data
dat = get_endoCounts(normalized=FALSE)
bulk = dat$counts
emeta = dat$meta
bulk = ExpressionSet(assayData=as.matrix(bulk))
rm(dat); gc()

meta = get_metadata_gregersen()
seu = get_single_sample_gregersen_seurats( update=F )
gc()

acounts = c()
for( sample in names(seu) ) {
	cat("at",sample,"\n")
	a = seu[[sample]][["RNA"]]@counts
	colnames(a) = paste(colnames(a),".",sample,sep="")
	acounts = cbind(acounts,a)
}
gc()
smeta = meta[colnames(acounts),]
ind = !is.na(smeta$clusterID)
acounts = acounts[,ind]
smeta = smeta[ind,]
ref = SingleCellExperiment(list(counts=acounts), colData=DataFrame(smeta))
rm(acounts)
gc()

save(ref,bulk,file="music_input.RData")

props = music_prop(bulk.mtx = exprs(bulk), sc.sce = ref, clusters = 'clusterID', samples = 'subjectID', verbose = T)
save(props,file="music_output.RData")
a = props$Est.prop.weighted
a = a[,colSums(a) > 0]
a = as.data.frame(a)
semeta = emeta[rownames(a),]
aprops = sapply( split(a, paste(semeta$cell_type,semeta$type)), colMeans )
aprops = round(aprops * 100)
a$sample = rownames(a)
a = melt(a)
a = cbind(a,emeta[a$sample,])

pdf("music_endo_with_gregersen_proportions.pdf", width=4*4, height=4)
ggplot(a,aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45p" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45p" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45m" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45m" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "all" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "all" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
dev.off()

# level 2
props_level_2 = music_prop(bulk.mtx = exprs(bulk), sc.sce = ref, clusters = 'cell_type_level_2', samples = 'subjectID', verbose = T)
save(props_level_2,file="music_output_level_2.RData")
a = props_level_2$Est.prop.weighted
a = a[,colSums(a) > 0]
a = as.data.frame(a)
semeta = emeta[rownames(a),]
aprops = sapply( split(a, paste(semeta$cell_type,semeta$type)), colMeans )
aprops = round(aprops * 100)
a$sample = rownames(a)
a = melt(a)
a = cbind(a,emeta[a$sample,])
a = a[order(a$cell_type,a$type,a$endo),]
a$sample = factor(a$sample,levels=unique(a$sample),ordered=T)
a$variable = factor( a$variable, levels=rev(c("myeloid","pdc", "granulocyte","b","T","uNK","stromal","EC","epithelial","unknown")), ordered=T )

pdf("music_endo_with_gregersen_proportions_level_2.pdf", width=4*4, height=4)
ggplot(a,aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45p" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45p" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45m" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45m" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "all" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "all" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
dev.off()

pdf("music_endo_with_gregersen_proportions_level_2b.pdf", width=4, height=4*4)
ggplot(a,aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic() + coord_flip()
dev.off()

# level 0
props_level_0 = music_prop(bulk.mtx = exprs(bulk), sc.sce = ref, clusters = 'cell_type', samples = 'subjectID', verbose = T)
save(props_level_0,file="music_output_level_0.RData")

a = props_level_0$Est.prop.weighted
a = a[,colSums(a) > 0]
a = as.data.frame(a)
semeta = emeta[rownames(a),]
aprops = sapply( split(a, paste(semeta$cell_type,semeta$type)), colMeans )
aprops = round(aprops * 100)
a$sample = rownames(a)
a = melt(a)
a = cbind(a,emeta[a$sample,])

pdf("music_endo_with_gregersen_proportions_level_0.pdf", width=4*4, height=4)
ggplot(a,aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45p" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45p" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45m" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "Cd45m" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "all" & a$type == "ex_vivo",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
ggplot(a[a$cell_type == "all" & a$type == "in_vitro",],aes(sample,value,fill=variable)) +  geom_bar(position="stack", stat="identity") + theme_classic()
dev.off()