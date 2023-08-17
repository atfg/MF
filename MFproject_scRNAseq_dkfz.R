###############################################################################
# Angela's reminders for starting jobs, ignore this comment block: 
# bsub -R "rusage[mem=20GB]" -Is bash 
# module unload R; module unload python; module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.1.0 gcc/7.2.0; R
###############################################################################

# this loads up libraries, sets up some global varibles and functions 
# and loads up the seurat objects created by Luca so takes a few minutes to run
source("/omics/groups/OE0433/internal//angela/atfg_github/MF/initialise_env.R")

samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM") 
#"MFCON020afM","MFCON007dfM" are failed samples


###############################################################################
# plots and clustering sample by sample
###############################################################################
samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020afM","MFCON007dfM") 
a = get_single_sample_seurats( samples, update=F )
seu = a[[1]]
markers = a[[2]]
rm(a)
gc()
aqc = get_cellranger_qc_output( samples )

for( name in names(seu) ) {
	cat("\tat",name,"\n")
	
	file = paste("scrna_umaps_", name, ".pdf", sep="" )
	
	pdf(file, width=6, height=6 )
	a = data.frame( Embeddings(seu[[name]][["umap"]]), seu[[name]]@meta.data ) 
	p = ggplot( a, aes(UMAP_1,UMAP_2, color=SCT_snn_res.0.2) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
	print(p)
	p = ggplot( a, aes(UMAP_1,UMAP_2, color=SCT_snn_res.1) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
	print(p)
	p = ggplot( a, aes(UMAP_1,UMAP_2, color=cell_type_single_cellbender) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
	print(p)
	mks = unique(unlist(goi))
	mks = mks[ mks %in% rownames(seu[[name]])]
	for( i in 1:length(mks) ) {
		mk = mks[i]
		p = FeaturePlot(seu[[name]], features = mk)
		print(p)
	}
	for( cluster in unique(markers[[name]]$cluster) ) {
		cat("\tat",cluster,"\n")
		mks = markers[[name]][ markers[[name]]$cluster == cluster, ]
		for( i in 1:min(10,nrow(mks)) ) {
			mk = mks[i,"gene"]
			p = FeaturePlot(seu[[name]], features = mk)
			print(p)
		}
	}
	dev.off()
}

coi = seu[["MFCON007dcM"]]@assays$RNA@counts["OVGP1",] > 0
head(sort(seu[["MFCON007dcM"]]@assays$SCT[,coi][,1],decreasing=T))

epi = Embeddings(seu[["MFCON007dcM"]][["umap"]])[,1] <= 1

dat = as.matrix(seu[["MFCON007dcM"]]@assays$SCT[,epi])
dim(dat)
dat = dat[rowSums(dat > 0) >= 5,]
dim(dat)
coi = dat["OVGP1",] > 0

res = c()
for( i in 1:nrow(dat) ) {
	progress(i,nrow(dat))
	res = rbind(res,coef(summary(lm(exp ~ coi, data=data.frame(exp=dat[i,],coi) )))["coiTRUE",c("Estimate","Pr(>|t|)")])
}
rownames(res) = rownames(dat)
res = res[order(abs(res[,"Estimate"]),decreasing=T),]

pdf("temp.pdf")
for( i in 1:30 ) {
	tab = data.frame(exp=dat[rownames(res)[i],], coi )
	p = ggplot(tab,aes(coi,exp)) + geom_jitter() + theme_bw() + ggtitle(rownames(res)[i])
	print(p)
}
dev.off()

stem( seu[["MFCON007dcM"]]@assays$RNA@counts["KRT14",] )

###############################################################################
# integrated
###############################################################################
samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020afM","MFCON007dfM") 
a = get_integrated_seurat( samples, update=F )
iseu = a[[1]]
imarkers = a[[2]]
rm(a)
gc()

file = paste("scrna_umaps_integrated_clusters_", paste(samples,collapse="_"), ".pdf", sep="" )
pdf(file, width=6, height=6 )
a = data.frame( Embeddings(iseu[["umap"]]), iseu@meta.data ) 
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.0.2) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
ggplot( a, aes(UMAP_1,UMAP_2, color=sample) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=cell_type_single_cellbender) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
mks = unique(unlist(goi))
mks = mks[ mks %in% rownames(iseu)]
for( i in 1:length(mks) ) {
	mk = mks[i]
	p = FeaturePlot(iseu, features = mk)
	print(p)
}
for( cluster in unique(imarkers$cluster) ) {
	cat("\tat",cluster,"\n")
	mks = imarkers[ imarkers$cluster == cluster & imarkers$p_val < 0.1 & imarkers$avg_logFC > 0 & imarkers$in_goi, ]
	for( i in 1:nrow(mks) ) {
		mk = mks[i,"gene"]
		p = FeaturePlot(iseu, features = mk)
		print(p)
	}
}
dev.off()

###############################################################################
# cellbender, plots sample by sample
###############################################################################
samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020afM","MFCON007dfM") 
samples = c("MFCON007dcM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020acM")  #"MFCON020afM","MFCON007dfM"
samples = c("MFCON007dcM","MFCON007dfM") 
a = get_single_sample_cellbender_seurats( samples, update=F)
seu = a[[1]]
markers = a[[2]]
rm(a)
gc()
aqc = get_cellranger_qc_output( samples )
aqc$cellbender_number_cells = sapply( seu, ncol )

# plot UMAPs colored by cluster and GOI
for( name in names(seu) ) {
	cat("\tat",name,"\n")
	file = paste("scrna_umaps_cellbender_", name, ".pdf", sep="" )
	
	pdf(file, width=6, height=6 )
	a = data.frame( Embeddings(seu[[name]][["umap"]]), seu[[name]]@meta.data ) 
	p = ggplot( a, aes(UMAP_1,UMAP_2, color=seurat_clusters) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
	print(p)
	mks = unique(unlist(goi))
	mks = mks[ mks %in% rownames(seu[[name]])]
	for( i in 1:length(mks) ) {
		mk = mks[i]
		p = FeaturePlot(seu[[name]], features = mk)
		print(p)
	}
	for( cluster in unique(markers[[name]]$cluster) ) {
		cat("\tat cluster",cluster,"\n")
		mks = markers[[name]][ markers[[name]]$cluster == cluster, ]
		for( i in 1:min(10,nrow(mks)) ) {
			mk = mks[i,"gene"]
			p = FeaturePlot(seu[[name]], features = mk)
			print(p)
		}
	}
	dev.off()
}




###############################################################################
# cellbender, integrated
###############################################################################
samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020afM","MFCON007dfM") 
#samples = c("MFCON007dcM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020acM") 
#samples = c("MFCON007dcM","MFCON010dfM") 
#samples = c("MFCON007dcM","MFCON007dfM") 
a = get_integrated_cellbender_seurat( samples, update=F )
iseu = a[[1]]
imarkers = a[[2]]
rm(a)
gc()
imarkers$in_goi = imarkers$gene %in% unlist(goi)

file = paste("scrna_umaps_cellbender_integrated_clusters_", paste(samples,collapse="_"), ".pdf", sep="" )
pdf(file, width=6, height=6 )
a = data.frame( Embeddings(iseu[["umap"]]), iseu@meta.data ) 
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.0.1) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.0.2) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.0.4) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.0.6) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=integrated_snn_res.1) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=sample) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw()

mks = unique(unlist(goi))
mks = mks[ mks %in% rownames(iseu)]
for( i in 1:length(mks) ) {
	mk = mks[i]
	p = FeaturePlot(iseu, features = mk)
	print(p)
}

for( cluster in unique(imarkers$cluster) ) {
	cat("\tat",cluster,"\n")
	mks = imarkers[ imarkers$cluster == cluster & imarkers$p_val < 0.1 & imarkers$avg_logFC > 0 & imarkers$in_goi, ]
	for( i in 1:nrow(mks) ) {
		mk = mks[i,"gene"]
		p = FeaturePlot(iseu, features = mk)
		print(p)
	}
}
dev.off()

file = paste("scrna_umaps_cellbender_integrated_clusters_", paste(samples,collapse="_"), "_markers.pdf", sep="" )
pdf(file, width=6, height=5 )
Idents(iseu) = iseu@meta.data$cell_type_level2
features = c("PAEP","EZR","NEAT1","VMP1","CD69","GNLY","NKG7","IGFBP1","LEFTY2","CLDN10","SPP1","MGP","VIM","COL12A1")
VlnPlot(iseu, features = features, stack=T, pt.size=0)
DotPlot(iseu, features = features) + RotatedAxis()
features = c("PAEP","EZR","NEAT1","CD69","IGFBP1","CLDN10","SPP1","MGP","VIM")
DotPlot(iseu, features = features) + RotatedAxis()
dev.off()

#table( iseu@meta.data[,c("cell_type_level2","integrated_snn_res.0.2")] )
imarkers[ imarkers$cluster == 6 & imarkers$p_val_adj < 0.05,][1:20,]

pdf("temp.pdf", width=6, height=6 )
a = data.frame( Embeddings(iseu[["umap"]]), iseu@meta.data ) 
ggplot( a, aes(UMAP_1,UMAP_2, color=cell_type_level1) ) + geom_point(alpha=0.4, shape=20) + theme_bw() + cols1_scale
ggplot( a, aes(UMAP_1,UMAP_2, color=cell_type_level2) ) + geom_point(alpha=0.4, shape=20) + theme_bw() + cols2_scale
library("reshape")
b = sapply( split( factor(a$cell_type_level2), a$sample ), table )
b = t(apply(b,2,function(x){ x/sum(x)}))
sample_order = names(sort(rowSums(b[,grepl("glands",colnames(b))]),decreasing=T))
sample_order = substring(sample_order,7,10)
b = melt(b)
b$X1 = substring(b$X1,7,10)
b$X1 = factor(b$X1,levels=sample_order,ordered=T)
b$X2 = factor(b$X2,levels=rev(types_level2),ordered=T)
ggplot( b, aes(X1,value,fill=X2) ) + geom_bar(stat="identity") + theme_bw() + cols2_fill + xlab("sample") + ylab("fraction")
dev.off()


###############################################################################
# garnett
###############################################################################

cds = to_cds( seu[["MFCON007dcM"]])

marker_check = check_markers(cds, "/icgc/dkfzlsdf/analysis/B210/angela/atfg_github/MF/classification_markers.txt", 
		db=org.Hs.eg.db, cds_gene_id_type = "SYMBOL", marker_file_gene_id_type = "SYMBOL")

file = paste("garnett_check_markers.pdf", sep="" )
pdf(file, width=6, height=12 )
plot_markers(marker_check)
dev.off()


###############################################################################
# load Quake data
###############################################################################
samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM","MFCON020afM","MFCON007dfM") 
a = get_integrated_cellbender_seurat( samples, update=F )
iseu = a[[1]]
rm(a)
gc()
#iseu@meta.data$donor = names(seu)[as.numeric(sapply(strsplit(colnames(iseu),"_"),"[[",2))]
#iseu@meta.data$type = "flow"
#iseu@meta.data$type[substr( iseu@meta.data$sample, 10, 10) == "c"] = "clumps"
#iseu@meta.data$cell_type_level1 = NA
#ind = as.character(iseu@meta.data$integrated_snn_res.0.2)
#iseu@meta.data$cell_type_level1[ind %in% c("0","1","5")] = "epithelial"
#iseu@meta.data$cell_type_level1[ind %in% c("2","3","4")] = "stromal"
#iseu@meta.data$cell_type_level1[ind %in% c("6")] = "leukocytes"
#
#iseu@meta.data$cell_type_level2 = NA
#ind = as.character(iseu@meta.data$integrated_snn_res.0.2)
#iseu@meta.data$cell_type_level2[ind == "0"] = "glands"
#iseu@meta.data$cell_type_level2[ind == "1"] = "exhausted glands"
#iseu@meta.data$cell_type_level2[ind == "5"] = "cilliated glands"
#iseu@meta.data$cell_type_level2[ind == "2"] = "decidualised stroma"
#iseu@meta.data$cell_type_level2[ind == "3"] = "stroma/mesenchyme"
#iseu@meta.data$cell_type_level2[ind == "4"] = "smooth muscle"
#iseu@meta.data$cell_type_level2[ind == "6"] = "T-cells"

a = get_integrated_mf_quake()
seu.integrated = a[[1]]
markers = a[[2]]
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
a$cell_type_level1[grepl("Endothelia|muscle",a$cell_type)] = "smooth muscle/endothelial"
table(a$cell_type)
table(a$cell_type_level1)

table( a$seurat_clusters, a$ct_qu )
table( a$seurat_clusters, a$ct_mf )


library("reshape")
b = table( a$cell_type_level1, a$type )
b = apply( b, 2, function(x) { x/sum(x)*100 } )
b = melt(b)
pdf("comparison_quake_proportions.pdf")
ggplot(b,aes(Var.2,value,fill=Var.1)) + geom_col() + theme_classic()
dev.off()

pdf("temp.pdf", width=2*6, height=6 )
ggplot(a,aes(UMAP_1,UMAP_2,color=integrated_snn_res.0.2)) + geom_point() + theme_bw()
ggplot(a,aes(UMAP_1,UMAP_2,color=cell_type_level1)) + geom_point() + theme_bw() + cols1_scale
ggplot(a,aes(UMAP_1,UMAP_2,color=cell_type_level1)) + geom_point() + theme_bw() + facet_wrap( ~ type ) + cols1_scale
ggplot(a,aes(UMAP_1,UMAP_2,color=donor)) + geom_point() + theme_bw()
dev.off()

ct = unique(a$cell_type_level1)[1]
tt = data.frame(rowSums(seu.integrated[["RNA"]]@counts[,a$cell_type_level1 == ct & a$type == "biopsy"]),
		rowSums(seu.integrated[["RNA"]]@counts[,a$cell_type_level1 == ct & a$type == "mf"]))
colnames(tt) = c(paste(ct,"_biopsy",sep=""),paste(ct,"_mf",sep=""))
dat = tt
for( ct in unique(a$cell_type_level1)[-1] ) {
	tt = data.frame(rowSums(seu.integrated[["RNA"]]@counts[,a$cell_type_level1 == ct & a$type == "biopsy"]),
			rowSums(seu.integrated[["RNA"]]@counts[,a$cell_type_level1 == ct & a$type == "mf"]))
	colnames(tt) = c(paste(ct,"_biopsy",sep=""),paste(ct,"_mf",sep=""))
	dat = cbind(dat,tt)
}
dat = dat[,colSums(dat)>0]
dat = dat[(rowSums(dat[,grepl("mf",colnames(dat))]) > 0) & (rowSums(dat[,grepl("biopsy",colnames(dat))]) > 0),]
dat = log10(dat+1)
cm = cor(dat,method="spearman")

tab = cm[c("stromal_mf","epithelial_mf","leukocytes_mf"),c("stromal_biopsy","epithelial_biopsy","leukocytes_biopsy")] #,"smooth muscle_biopsy","endothelial_biopsy")]
#tab = apply(tab,1,function(x) { (x-mean(x))/sd(x) } )
#tab = t(tab)
tab = melt(tab)
tab = as.data.frame(tab)
tab$X1 = factor(tab$X1,levels=c("epithelial_mf","stromal_mf","leukocytes_mf"),ordered=T)
tab$X2 = factor(tab$X2,levels=c("epithelial_biopsy","stromal_biopsy","leukocytes_biopsy","endothelial_biopsy","smooth muscle_biopsy"),ordered=T)

pdf("temp.pdf", width=2*6, height=6 )
ggplot(tab,aes(X2,X1,fill=value)) + geom_tile() + theme_bw() + scale_fill_distiller()
dev.off()

ggplot(dat,aes(stromal_biopsy,stromal_mf)) + geom_point() + theme_bw()
ggplot(dat,aes(stromal_biopsy,epithelial_mf)) + geom_point() + theme_bw()
ggplot(dat,aes(stromal_biopsy,leukocytes_mf)) + geom_point() + theme_bw()
dev.off()

#cp /icgc/dkfzlsdf/analysis/B210/angela/mf/human_endometrium_quake.RData /icgc/dkfzlsdf/analysis/OE0538_projects/DO-0002/human_endometrium_quake_seurat_sct_by_donor.RData

###############################################################################
# old 
###############################################################################

maxk = 20
a = get_integrated_data( samples, k=20, update=F )
seurat = a[[1]]
seurat.markers = a[[2]]
rm(a)

s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes
seurat = CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

ids = c("stroma_0", "epithelial_1", "epithelial_2", "stroma_3", "stroma_4", "epithelial_5",
		"epithelial_6","stroma_7","epithelial_8","stroma_9","leukocytes_10","leukocytes_11",
		"endothelial_12","epithelial_13")
names(ids) = levels(seurat@meta.data$seurat_clusters)
seurat = SetIdent(seurat,value="seurat_clusters")
seurat = RenameIdents(seurat, ids)

pdf(paste("plots/integrated_umap_maxk", maxk, ".pdf",sep=""), width=4, height=4)
p = DimPlot(seurat, reduction = "umap", group.by = "donor")
print(p)
#DimPlot(seurat, reduction = "umap", group.by = "percent.mt")
p = DimPlot(seurat, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend()
print(p)
p = DimPlot(seurat, reduction = "umap", group.by = "Phase", label = FALSE, repel = TRUE) 
print(p)
p = RidgePlot(seurat, features = s.genes[s.genes %in% rownames(seurat)], ncol = 2)
print(p)
p = RidgePlot(seurat, features = g2m.genes[g2m.genes %in% rownames(seurat)], ncol = 2)
print(p)
dev.off()


filename = paste("plots/integrated_seurat_maxk", maxk, "_goi.pdf",sep="")
pdf(filename,width=10, height=10)
for( g in unique(unlist(goi)) ) {
	cat("\tat",g,"\n")
	p = try( FeaturePlot(object = seurat, features = g, cols = c("lightgrey", "blue") ) )
	print(p)
}
dev.off()

# slingshot
# select only epithelial
#dat = t(as.matrix(seurat@assays$integrated@data))
dat = seurat$pca@cell.embeddings
cls = as.character(Idents(seurat))
ind = grepl("epithelial",cls)
sg = slingshot(dat[ind,])

#colors = colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
#plotcol = colors[cut(sim$slingPseudotime_1, breaks=100)]



t = slingPseudotime(sg)[,"curve1"]
Y = seurat@assays$integrated@data[,ind]
var100 = names(sort(apply(Y,1,var),decreasing = TRUE))[1:100]
Y = Y[var100,]
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
			d <- data.frame(z=z, t=t)
			suppressWarnings({
						tmp <- suppressWarnings(gam(z ~ lo(t), data=d))
					})
			p <- summary(tmp)[3][[1]][2,3]
			p
		})
# gam.pval
topgenes = names(sort(gam.pval, decreasing = FALSE))[1:100]
heatdata = as.matrix(seurat@assays$integrated@data)[,ind][topgenes, order(t, na.last = NA)]
heatclus = cls[ind][order(t, na.last = NA)]


pdf(paste("plots/integrated_seurat_maxk", maxk, "_slignshot_epithelial.pdf",sep=""),width=6, height=6)
cols = brewer.pal(length(unique(cls[ind])),'Spectral')
cols = cols[as.numeric(factor(cls[ind]))]
names(cols) = cls[ind]
colkey = cols[!duplicated(cols)]
plot(reducedDims(sg), pch=16, col=cols, asp = 1)
lines(sg, lwd=2, col='black')
legend("bottomright",legend=names(colkey),fill=colkey)
heatmap(heatdata, Colv = NA, ColSideColors = brewer.pal(length(unique(cls[ind])),"Spectral")[as.numeric(factor(heatclus))])
dev.off()

# stroma
dat = seurat$pca@cell.embeddings
cls = as.character(Idents(seurat))
ind = grepl("stroma",cls)
sg_stroma = slingshot(dat[ind,])

pdf(paste("plots/integrated_seurat_maxk", maxk, "_slignshot_stroma.pdf",sep=""),width=6, height=6)
cols = brewer.pal(length(unique(cls[ind])),'Spectral')
cols = cols[as.numeric(factor(cls[ind]))]
names(cols) = cls[ind]
colkey = cols[!duplicated(cols)]
plot(reducedDims(sg), pch=16, col=cols, asp = 1)
lines(sg, lwd=2, col='black')
legend("bottomright",legend=names(colkey),fill=colkey)
dev.off()

###############################################################################
# cell bender, failed samples
###############################################################################
sample = "MFCON007dfM"
i = which(datamap$sample_full == sample)
file = paste(datamap$base_dir[i],"outs/raw_gene_bc_matrices/GRCh38/cellbender_matrix_filtered.h5", sep="")
dat = Read10X_h5(filename = file, use.names = TRUE)
bseu = CreateSeuratObject(counts = dat)
rm(dat)
gc()

bseu = SCTransform(bseu, verbose = TRUE)
bseu = RunPCA(bseu, verbose = TRUE)
bseu = RunUMAP(bseu, dims = 1:30, verbose = TRUE)
bseu = FindNeighbors(bseu, dims = 1:30, verbose = TRUE)
bseu = FindClusters(bseu, verbose = TRUE, resolution=0.2)

bseu@meta.data$donor = datamap$sample_full[i]
bseu = PercentageFeatureSet(bseu, pattern = "^MT-", col.name="percent.mt" )

markers = FindAllMarkers(bseu, verbose = TRUE)
seu = bseu
rm(bseu)
gc()

file = paste("scrna_umaps_cellbender_", sample, ".pdf", sep="" )
pdf(file, width=6, height=6 )
a = data.frame( Embeddings(seu[["umap"]]), seu@meta.data ) 
p = ggplot( a, aes(UMAP_1,UMAP_2, color=seurat_clusters) ) + geom_point(alpha=0.3, shape=20) + theme_bw() 
print(p)
mks = unique(unlist(goi))
mks = mks[ mks %in% rownames(seu)]
for( i in 1:length(mks) ) {
	mk = mks[i]
	p = FeaturePlot(seu, features = mk)
	print(p)
}
dev.off()


###############################################################################
# DE between quake and MF
###############################################################################

a = get_integrated_mf_quake()
seu.integrated = a[[1]]
markers = a[[2]]
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

ind = rowSums(ndat >= 3) >=3
sdat = log10(ndat[ind,]+1)
sdevs = apply(sdat,1,sd)
sdat = sdat[order(sdevs,decreasing=T),][1:2000,]
pc = prcomp( t(sdat) )

a = as.data.frame(pc$x[,1:4])
a$type = c("biopsy","MF")[grepl("MF",rownames(a))+1]
a$ct = sapply(strsplit(rownames(a),"_"),"[[",2)
a$sample = sapply(strsplit(rownames(a),"_"),"[[",1)

pdf("comparison_quake_pseudocount_PCA.pdf", width=4, height=4)
ggplot(a,aes(PC1,PC2,color=type)) + geom_point() + theme_classic()

ggplot(a,aes(PC2,PC3,color=type)) + geom_point() + theme_classic()
ggplot(a,aes(PC2,PC3,color=ct)) + geom_point() + theme_classic()

ggplot(a,aes(PC3,PC4,color=type)) + geom_point() + theme_classic()
ggplot(a,aes(PC3,PC4,color=ct)) + geom_point() + theme_classic()
dev.off()

library("DESeq2")
library("ggrepel")
library("apeglm")

ares = list()
for( ct in c("stromal","epithelial","leukocytes") ) {
	cat("at",ct,"\n")
	ind = grepl(ct,rownames(a))
	dat = ndat[,ind]
	# filtering is quite important
	dat = dat[rowSums(dat > 0) >= 3,]
	de = DESeqDataSetFromMatrix(countData = dat, colData = a[ind,], design= ~ type)
	de = DESeq(de)
	
	res = lfcShrink(de, coef="type_MF_vs_biopsy", type="apeglm")
	#res = as.data.frame(results(de))	
	res = res[!is.na(res$padj),]
	res = res[order(res$log2FoldChange,decreasing=T),]
	res$gene = rownames(res)
	ares[[ct]] = res
	write.table(res[,c("gene","log2FoldChange")],file=paste0("comparison_quake_DE_",ct,".rnk"),quote=F,row.names=F,col.names=F,sep="\t")
}

pdf("comparison_quake_DE.pdf")
for( ct in c("stromal","epithelial","leukocytes") ) {
	tab = as.data.frame(ares[[ct]])
	p = ggplot(tab,aes(log10(baseMean),log2FoldChange,label=gene)) + geom_point(alpha=0.5) + geom_text_repel(data=tab[c(1:15,(nrow(tab)-15):nrow(tab)),], max.overlaps = Inf) + theme_classic() + ggtitle(ct)
	print(p)
}
dev.off()


# there's a error - GNB2L1 and RACK1 are the exact same gene



###############################################################################
# Deconvolve EPFL data using our samples or Quake's
#
# module unload R; module unload python; module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.2.0 gcc/7.2.0; R
#module unload R; module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.0.0; R
###############################################################################

dir = "/omics/groups/OE0433/internal/angela/mf/"
options(width=230)
setwd(dir)
.libPaths(paste("/omics/groups/OE0433/internal/software/R-library-",paste(R.version$major,R.version$minor,sep="."),sep=""))

# read in count table
tab = read.table("STAR-HTSeqTags.-Jun18.counts.txt", check.names=F)

# read in metadata
meta = data.frame( id=sapply(strsplit(colnames(tab),".",fixed=T),"[[",1),
		endo_stage=as.numeric(sapply(strsplit(colnames(tab),".",fixed=T),"[[",2)),
		contraceptive=sapply(strsplit(colnames(tab),".",fixed=T),"[[",3),
		cell_type=sapply(strsplit(colnames(tab),".",fixed=T),"[[",4),
		unknown=as.numeric(sapply(strsplit(colnames(tab),".",fixed=T),"[[",5)) )
meta$endo_stage[is.na(meta$endo_stage)] = 0
meta$endo_stage = factor(meta$endo_stage)

# which method to use?

library("DeconRNASeq")
library("DEseq2")

de = DESeqDataSetFromMatrix(countData = tab, colData = meta, design= ~ endo_stage)
de = estimateSizeFactors(de)
ntab = counts(de, normalized=TRUE)


res = DeconRNASeq(as.data.frame(ndat), as.data.frame(ref))

a = t(apply( res$out.all, 1, function(x) {  } ))


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


