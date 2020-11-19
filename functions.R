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

get_single_sample_cellbender_seurats <- function(samples, update=F) {
	dm = datamap[ datamap$sample_full %in% samples, ]
	seu = list()
	markers = list()
	for( i in 1:nrow(dm) ) {
		cat("\tat",i,"\n")
		
		filename = paste(dm$base_dir[i],"outs/raw_gene_bc_matrices/GRCh38/cellbender_seurat.RData",sep="" )
		if( !file.exists(filename) | update ) {	
			file = paste(dm$base_dir[i],"outs/raw_gene_bc_matrices/GRCh38/cellbender_matrix_filtered.h5", sep="")
			dat = Read10X_h5(filename = file, use.names = TRUE)
			bseu = CreateSeuratObject(counts = dat)
			rm(dat)
			gc()
			
			bseu@meta.data$sample = dm$sample_full[i]
			bseu = PercentageFeatureSet(bseu, pattern = "^MT-", col.name="percent.mt" )
			
			a = data.frame( bseu@meta.data ) 
			ind = a$nCount_RNA < quantile(a$nCount_RNA,probs=seq(0,1,0.01))[99] & 
					a$nFeature_RNA < quantile(a$nFeature_RNA,probs=seq(0,1,0.01))[99] &
					a$percent.mt < 30 &
					a$nFeature_RNA > 50
			bseu = bseu[,ind]
			
			bseu = SCTransform(bseu, verbose = TRUE)
			bseu = RunPCA(bseu, verbose = TRUE)
			bseu = RunUMAP(bseu, dims = 1:30, verbose = TRUE)
			bseu = FindNeighbors(bseu, dims = 1:30, verbose = TRUE)
			bseu = FindClusters(bseu, verbose = TRUE, resolution=0.2)
			
			markers[[dm$sample_full[i]]] = FindAllMarkers(bseu, verbose = TRUE)
			seu[[dm$sample_full[i]]] = bseu
				
			rm(bseu)
			gc()
			save(seu,markers,file=filename)
			cat("\tsaved to",filename,"\n")
		} else{
			cat("\tloading from", filename, "\n")
			load(filename)
		}
	}
	
	# manual annotation of cell types
#	name = "MFCON007efM"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual = NA
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters %in% c(0,4)] = "epithelial"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters %in% c(1,2)] = "stromal"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters == 3] = "leukocytes"
#	
#	name = "MFCON010dfM"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual = NA
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters %in% c(4)] = "epithelial"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters %in% c(0,1,2,3)] = "stromal"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters == 5] = "erythrocytes"
#	
#	name = "MFCON018bfM"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual = "unknown"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters %in% c(0)] = "epithelial"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters %in% c(3)] = "stromal"
#	seu[[name]]@meta.data$cell_type_seurat_clusters_0.2_manual[seu[[name]]@meta.data$seurat_clusters == 4] = "leukocytes"
	
	return(list(seu,markers))
}

to_cds <- function(seud) {
	pd = new("AnnotatedDataFrame", data = seud@meta.data)
	a = as(seud[["RNA"]]@counts, "dgCMatrix")
	fd = new("AnnotatedDataFrame", data = data.frame(row.names=rownames(a), gene_short_name=rownames(a)))
	cds = newCellDataSet(a, phenoData = pd, featureData = fd)
	cds = estimateSizeFactors(cds)
}

get_human_data <- function() {
	file = "endometrium_all_seurat.RData"
	
	if( !file.exists(file) ) {
		filename = "/icgc/dkfzlsdf/analysis/B210/data/sanger_human_endometrium/endometrium_all.h5ad"
		dat = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
		
		meta = h5read( filename, "/obs" )
		rownames(meta) = meta$index
		var = h5read( filename, "/var" )
		counts = h5read( filename, "/X" )
		colnames(counts) = meta$index
		rownames(counts) = var$index
		rvt_umap = h5read( filename, "/obsm" )
		
		map = h5read( filename, "/uns" )
		for( name in c("batch","clinical","clusters","day","individual","leiden","location","phase","sample","treatment","type") ) {
			meta[,name] = map[[paste(name,"categories",sep="_")]][as.factor(meta[,name])]
		}
		meta$rvt_UMAP1 = rvt_umap$X_umap[1,]
		meta$rvt_UMAP2 = rvt_umap$X_umap[2,]
		rm(rvt_umap)
		
		seu = CreateSeuratObject( counts, project = "HumanBiopsies", assay = "RNA", min.cells = 0, 
				min.features = 0, names.field = 1, names.delim = "_", meta.data = meta )
		save(seu,file=file)
	} else {
		load(file)
	}
	return(seu)
}


# CCA
get_integrated_cellbender_seurat <- function() {
	filename = "integrated_cellbender_seurat.RData"
	
	if( !file.exists(filename) | update ) {
		seu.list = get_single_sample_cellbender_seurats()[[1]]
		
		seu.features = SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)
		seu.list = PrepSCTIntegration(object.list = seu.list, anchor.features = seu.features, verbose = TRUE)
		seu.anchors = FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", anchor.features = seu.features, verbose = TRUE)
		seu.integrated = IntegrateData(anchorset = seu.anchors, normalization.method = "SCT", verbose = TRUE)
		
		seu.integrated = RunPCA(seu.integrated, verbose = TRUE)
		seu.integrated = RunUMAP(seu.integrated, dims = 1:30)
		
		seu.integrated = FindNeighbors(seu.integrated, dims = 1:30, verbose = TRUE)
		seu.integrated = FindClusters(seu.integrated, verbose = TRUE, resolution=0.2)
		markers = FindAllMarkers(seu.integrated, verbose = TRUE)
		
		save(seu.integrated, markers, file=filename)
		cat("\tsaved to",filename,"\n")
	} else{
		cat("\tloading from", filename, "\n")
		load(filename)
	}
	return(seu.integrated)
}

# run cell bender on samples in the directory
run_cellbender <- function() {
	
	for( sample in samples ) {
		path = datamap[ match(sample, datamap$sample_full), "base_dir"]
		cmd = paste( "bsub -R 'rusage[mem=20GB]' 'chmod u+x /icgc/dkfzlsdf/analysis/B210/lab_github/Umay/CellBender.sh; ", 
				"/icgc/dkfzlsdf/analysis/B210/angela/atfg_github/MF/CellBender.sh ", path,"'", sep="" )
		cat(cmd,"\n")
	}	
	
	
}

get_cellranger_qc_output <- function() {
	dm = datamap[ datamap$sample_full %in% samples, ]
	
	aqc = c()
	for( i in 1:nrow(dm) ) {
		file = paste( dm$base_dir[i],"outs/metrics_summary.csv", sep="")
		aqc = rbind(aqc, read.csv(file) )
	}
	rownames(aqc) = dm$sample_full
	
	return(aqc)
}




