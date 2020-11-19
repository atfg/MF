###############################################################################
# Angela's reminders for starting jobs, ignore this comment block: 
# bsub -R "rusage[mem=20GB]" -Is bash 
# module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.0.0; R
###############################################################################

# this loads up libraries, sets up some global varibles and functions 
# and loads up the seurat objects created by Luca so takes a few minutes to run
source("/icgc/dkfzlsdf/analysis/B210/angela/atfg_github/MF/initialise_env.R")

samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON018bfM") 
#"MFCON020afM","MFCON007dfM" are failed samples

###############################################################################
# some analysis
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
# cell bender, sample by sample
###############################################################################
seu = get_single_sample_cellbender_seurats( samples, update=F)[[1]]
aqc = get_cellranger_qc_output()
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
	dev.off()
}

# QC plots
file = paste("scrna_cellbender_qc.pdf", sep="" )
pdf(file, width=6, height=6 )
for( name in names(seu) ) {
	cat("\tat",name,"\n")
	a = data.frame( seu[[name]]@meta.data ) 
	p = ggplot( a, aes(nFeature_RNA,nCount_RNA, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title=name) 
	print(p)
	a = a[ a$nCount_RNA < quantile(a$nCount_RNA,probs=seq(0,1,0.01))[99] & 
					a$nFeature_RNA < quantile(a$nFeature_RNA,probs=seq(0,1,0.01))[99] &
					a$percent.mt < 30 &
					a$nFeature_RNA > 50,]
	p = ggplot( a, aes(nFeature_RNA,nCount_RNA, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title=name) 
	print(p)
	p = ggplot( a, aes(nFeature_SCT,nCount_SCT, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title=name) 
	print(p)
}
dev.off()


###############################################################################
# cell bender, integrated
###############################################################################
a = get_integrated_cellbender_seurat()
iseu = a[[1]]
imarkers = a[[2]]
rm(a)
gc()

a = c()
for( name in names(seu) ) {
	a = rbind(a,seu[[name]]@meta.data)
}

iseu@meta.data$donor = names(seu)[as.numeric(sapply(strsplit(colnames(iseu),"_"),"[[",2))]


file = paste("scrna_umaps_cellbender_integrated_clusters.pdf", sep="" )
pdf(file, width=6, height=6 )
a = data.frame( Embeddings(iseu[["umap"]]), iseu@meta.data ) 
ggplot( a, aes(UMAP_1,UMAP_2, color=seurat_clusters) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
ggplot( a, aes(UMAP_1,UMAP_2, color=donor) ) + geom_point(alpha=0.3, shape=20) + theme_bw()
for( cluster in unique(markers$cluster) ) {
	cat("\tat",cluster,"\n")
	mks = markers[ markers$cluster == cluster, ]
	for( i in 1:min(10,nrow(mks)) ) {
		mk = mks[i,"gene"]
		p = FeaturePlot(iseu, features = mk)
		print(p)
	}
}
dev.off()

file = paste("scrna_umaps_cellbender_integrated.pdf", sep="" )
pdf(file, width=6, height=6 )
mks = unique(unlist(goi))
mks = mks[ mks %in% rownames(iseu)]
for( i in 1:length(mks) ) {
	mk = mks[i]
	p = FeaturePlot(iseu, features = mk)
	print(p)
}
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
# load Roser's data
###############################################################################


