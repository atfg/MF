# takes cellranger's filtered matrices and creates an integrated seurat object
get_integrated_data <- function( samples, k, update=F ) {
	kfile = paste("integrated_seurat_maxk",k,".RData",sep="")
	mfile = paste("integrated_seurat_maxk",k,"_markers.RData",sep="")
	
	if( !file.exists(kfile) | update ) {	
		cat( "reading in matrices\n")
		adat = list()
		for( sample in samples ) {
			cat("\tat", sample, "\n")
			dat = Read10X(data.dir = paste("matrices/", sample, sep="") )	
			adat[[sample]] = CreateSeuratObject(counts = dat)
		}
		
		cat( "normalising, finding varibale features and adding QC\n")
		for( sample in samples ) {
			cat("\tat", sample, "\n")
			adat[[sample]] = NormalizeData(adat[[sample]], verbose = TRUE)
			adat[[sample]] = FindVariableFeatures(adat[[sample]], selection.method = "vst", nfeatures = 2000, verbose = TRUE)
			adat[[sample]]@meta.data$donor = sample
			adat[[sample]] = PercentageFeatureSet(adat[[sample]], pattern = "^MT-", col.name="percent.mt" )
			adat[[sample]] = adat[[sample]][,adat[[sample]]@meta.data$percent.mt < mt_th]
		}
		
		# k.filter has to be smaller than the smalest number of cells
		cat( "finding anchors and integrating over a range of k filters\n")
		for( maxk in c(20,30,40,50) ) {
			cat("\tat", maxk, "\n")
			filename = paste("integrated_seurat_maxk",maxk,".RData",sep="")
			markerfile = paste("integrated_seurat_maxk",maxk,"_markers.RData",sep="")
			embedfile = paste("integrated_seurat_maxk",maxk,"_forsleepwalker.RData",sep="")
			
			anchors = FindIntegrationAnchors(object.list = adat[c("MFCON020acM","MFCON007efM","MFCON010dfM","MFCON007dcM","MFCON018bfM")], dims = 1:maxk, k.filter=110) 
			seurat = IntegrateData(anchorset = anchors, dims = 1:maxk)
			DefaultAssay(seurat) = "integrated"
			
			seurat = ScaleData(seurat, verbose = FALSE)
			seurat = RunPCA(seurat, npcs = maxk, verbose = FALSE)
			seurat = RunUMAP(seurat, reduction = "pca", dims = 1:maxk)
			
			seurat = FindNeighbors(object = seurat)
			seurat = FindClusters(object = seurat)
			
			save(seurat,file=filename)
			cat("\tsaved to", filename,"\n")
			
			seurat_umap = seurat$umap@cell.embeddings
			seurat_pca = seurat$pca@cell.embeddings
			save(seurat_umap,seurat_pca,file=embedfile)
			cat("\tsaved to", embedfile,"\n")
			
			cat( "finding all cluster markers\n")
			seurat.markers = FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
			save(seurat.markers, file=markerfile)
			cat("\tsaved to", markerfile,"\n")
		}
	} 
	cat("\tloading from", kfile, "\n")
	load(kfile)
	load(mfile)
	
	return(list(seurat,seurat.markers))
}

# run cell bender on samples in the directory
run_cellbender <- function() {
	for( sample in samples ) {
		path = datamap[ match(sample, datamap$sample_full), "base_dir"]
		cmd = paste( "bsub -R 'rusage[mem=20GB]' '/icgc/dkfzlsdf/analysis/B210/lab_github/Umay/CellBender.sh ", path,"'", sep="" )
		cat(cmd,"\n")
	}	
}




