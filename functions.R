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
		cat("\tat",dm$sample_full[i],"\n")
		
		filename = paste(dm$base_dir[i],"outs/raw_gene_bc_matrices/GRCh38/cellbender_seurat.RData",sep="" )
		if( !file.exists(filename) | update ) {	
			file = paste(dm$base_dir[i],"outs/raw_gene_bc_matrices/GRCh38/cellbender_matrix_filtered.h5", sep="")
			dat = Read10X_h5(filename = file, use.names = TRUE)
			bseu = CreateSeuratObject(counts = dat)
			rm(dat)
			gc()
			
			bseu@meta.data$sample = dm$sample_full[i]
			bseu = PercentageFeatureSet(bseu, pattern = "^MT-", col.name="percent.mt" )
			
			file = paste("scrna_cellbender_qc_", dm$sample_full[i], ".pdf", sep="" )
			pdf(file, width=6, height=6 )
			
			# filter cells 
			a = data.frame( bseu@meta.data ) 
			p = ggplot( a, aes(nFeature_RNA,nCount_RNA, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title=paste("before filtering",name)) 
			print(p)
			
			ind = a$nCount_RNA < quantile(a$nCount_RNA,probs=seq(0,1,0.01))[99] & 
					a$nFeature_RNA < quantile(a$nFeature_RNA,probs=seq(0,1,0.01))[99] &
					a$percent.mt < 10 &
					a$nFeature_RNA > 100
			bseu = bseu[,ind]
			
			bseu = SCTransform(bseu, verbose = TRUE)
			bseu = RunPCA(bseu, verbose = TRUE)
			bseu = RunUMAP(bseu, dims = 1:30, verbose = TRUE)
			bseu = FindNeighbors(bseu, dims = 1:30, verbose = TRUE)
			bseu = FindClusters(bseu, verbose = TRUE, resolution=0.2)
			
			bmarkers = FindAllMarkers(bseu, verbose = TRUE)
			
			a = data.frame( bseu@meta.data )
			p = ggplot( a, aes(nFeature_RNA,nCount_RNA, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title=paste("after filtering",name)) 
			print(p)
			dev.off()
			
			save(bseu,bmarkers,file=filename)
			cat("\tsaved to",filename,"\n")
		} else{
			cat("\tloading from", filename, "\n")
			load(filename)
		}
		markers[[dm$sample_full[i]]] = bmarkers
		seu[[dm$sample_full[i]]] = bseu
		rm(bseu)
		rm(bmarkers)
		gc()
	}
	
	# refine some 
	name = "MFCON007efM"
	seu[[name]] = FindClusters(seu[[name]], verbose = TRUE, resolution=1)
	seu[[name]]@meta.data$cell_type_cellbender_single = "stromal"
	seu[[name]]@meta.data$cell_type_cellbender_single[seu[[name]]@meta.data$SCT_snn_res.1 %in% c(3)] = "epithelial"
	seu[[name]]@meta.data$cell_type_cellbender_single[seu[[name]]@meta.data$SCT_snn_res.1 %in% c(1)] = "leukocytes"
	
# manual annotation of cell types
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

get_single_sample_seurats <- function(samples, update=F) {
	dm = datamap[ datamap$sample_full %in% samples, ]
	seu = list()
	markers = list()
	for( i in 1:nrow(dm) ) {
		cat("\tat",dm$sample_full[i],"\n")
		
		filename = paste(dm$base_dir[i],"outs/seurat.RData",sep="" )
		if( !file.exists(filename) | update ) {	
			file = paste(dm$base_dir[i],"outs/filtered_gene_bc_matrices_h5.h5", sep="")
			dat = Read10X_h5(filename = file, use.names = TRUE)
			bseu = CreateSeuratObject(counts = dat)
			rm(dat)
			gc()
			
			bseu@meta.data$sample = dm$sample_full[i]
			bseu = PercentageFeatureSet(bseu, pattern = "^MT-", col.name="percent.mt" )
			
			file = paste("scrna_qc_", dm$sample_full[i], ".pdf", sep="" )
			pdf(file, width=6, height=6 )
			
			# filter cells 
			a = data.frame( bseu@meta.data ) 
			p = ggplot( a, aes(nFeature_RNA,nCount_RNA, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title=paste("before filtering",name)) 
			print(p)
			
			ind = a$nCount_RNA < quantile(a$nCount_RNA,probs=seq(0,1,0.01))[99] & 
					a$nFeature_RNA < quantile(a$nFeature_RNA,probs=seq(0,1,0.01))[99] &
					a$percent.mt < 10 &
					a$nFeature_RNA > 100
			bseu = bseu[,ind]
			
			bseu = SCTransform(bseu, verbose = TRUE)
			bseu = RunPCA(bseu, verbose = TRUE)
			bseu = RunUMAP(bseu, dims = 1:30, verbose = TRUE)
			bseu = FindNeighbors(bseu, dims = 1:30, verbose = TRUE)
			bseu = FindClusters(bseu, verbose = TRUE, resolution=1)
			
			bmarkers = FindAllMarkers(bseu, verbose = TRUE)
			
			a = data.frame( bseu@meta.data )
			p = ggplot( a, aes(nFeature_RNA,nCount_RNA, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title=paste("after filtering",name)) 
			print(p)
			dev.off()
			
			save(bseu,bmarkers,file=filename)
			cat("\tsaved to",filename,"\n")
		} else{
			cat("\tloading from", filename, "\n")
			load(filename)
		}
		markers[[dm$sample_full[i]]] = bmarkers
		seu[[dm$sample_full[i]]] = bseu
		rm(bseu)
		rm(bmarkers)
		gc()
	}
	
	for( name in names(seu) ) {
		seu[[name]] = FindClusters(seu[[name]], verbose = TRUE, resolution=1)
	}
	
	cat( "\tsetting cell types\n")
	for( name in names(seu) ) {
		if( name == "MFCON007efM" ) {
			seu[[name]]@meta.data$cell_type_single_cellbender = "stromal"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.1 %in% c("0")] = "epithelial"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.1 %in% c("4")] = "leukocytes"
		} else if( name == "MFCON010dfM" ) {
			seu[[name]]@meta.data$cell_type_single_cellbender = "stromal"
		} else if( name == "MFCON018bfM" ) {
			seu[[name]]@meta.data$cell_type_single_cellbender = "stromal"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.1 %in% c("1")] = "epithelial"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.1 %in% c("2")] = "leukocytes"	
		} else if( name == "MFCON020afM" ) {
			seu[[name]]@meta.data$cell_type_single_cellbender = "unknown"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.1 %in% c("0","2")] = "epithelial"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.1 %in% c("1","3","6")] = "stromal"
		} else if( name == "MFCON007dfM" ) { # failed?
			seu[[name]]@meta.data$cell_type_single_cellbender = "unknown"
		} else if( name == "MFCON007dcM" ) {
			seu[[name]]@meta.data$cell_type_single_cellbender = "unknown"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.0.2 %in% c("2")] = "leukocytes"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.0.2 %in% c("0")] = "epithelial"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.0.2 %in% c("1")] = "stromal"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.0.2 %in% c("3")] = "endothelial"
		} else if( name == "MFCON020acM" ) {
			seu[[name]]@meta.data$cell_type_single_cellbender = "unknown"
			seu[[name]]@meta.data$cell_type_single_cellbender[seu[[name]]@meta.data$SCT_snn_res.1 %in% c("0","2")] = "epithelial"	
		}
	}
	
	return(list(seu,markers))
}


to_cds <- function(seud) {
	pd = new("AnnotatedDataFrame", data = seud@meta.data)
	a = as(seud[["RNA"]]@counts, "dgCMatrix")
	fd = new("AnnotatedDataFrame", data = data.frame(row.names=rownames(a), gene_short_name=rownames(a)))
	cds = newCellDataSet(a, phenoData = pd, featureData = fd)
	cds = estimateSizeFactors(cds)
}

get_human_data <- function( sample ) {
	file = paste("human_endometrium_normalised_individual_",sample,".RData",sep="")
	
	if( !file.exists(file) ) {
		bseu = get_raw_human_data( sample ) 
		bseu = PercentageFeatureSet(bseu, pattern = "^MT-", col.name="percent.mt" )
		bseu = SCTransform(bseu, verbose = TRUE)
		bseu = RunPCA(bseu, verbose = TRUE)
		bseu = RunUMAP(bseu, dims = 1:30, verbose = TRUE)
		bseu = FindNeighbors(bseu, dims = 1:30, verbose = TRUE)
		bseu = FindClusters(bseu, verbose = TRUE, resolution=0.2)
		bmarkers = FindAllMarkers(bseu, verbose = TRUE)
		
		save(bseu, bmarkers, file=file)
		cat("\tsaved to",file,"\n")
	} else {
		cat("\tloading from",file,"\n")
		load(file)
	}
	
}

get_raw_human_data <- function( samples=c("A10","A13") ) {
	file = "endometrium_all_seurat.RData"
	
	for( sample in samples ) {
		filename = paste("human_endometrium_individual_",sample,".RData",sep="")
		load(filename)
		
	}
	
	if( !file.exists(file) ) {
		cat("reading 10x data\n")
		filename = "/icgc/dkfzlsdf/analysis/B210/data/sanger_human_endometrium/endometrium_all.h5ad"
		
		#dat = Read10X_h5(filename, use.names = TRUE, unique.features = TRUE)
		
		# doesn't work on R 4.0.0 but works on R 3.5.2, problem with the rhdf5 library
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
		
		ind = meta$treatment == "C" & meta$clinical == "U" & meta$location == "ENMY"
		meta = meta[ind,]
		counts = counts[,ind]
		gc()
		
		for( individual in unique(meta$individual) ) {
			filename = paste("human_endometrium_individual_",individual,".RData",sep="")
			
			cat("creating seurat object\n")
			ind = meta$individual == individual
			seu = CreateSeuratObject( counts[,ind], project = "HumanBiopsies", assay = "RNA", min.cells = 0, 
					min.features = 0, names.field = 1, names.delim = "_", meta.data = meta[ind,] )
			cat("saving to", filename,"\n")
			save(seu,file=filename)	
			rm(seu)
		}
		
		cat("creating seurat object\n")
		seu = CreateSeuratObject( counts, project = "HumanBiopsies", assay = "RNA", min.cells = 0, 
				min.features = 0, names.field = 1, names.delim = "_", meta.data = meta )
		cat("saving to", file,"\n")
		save(seu,file=file)
	} else {
		cat("loading from", file,"\n")
		load(file)
	}
	return(seu)
}

# only grch38 
get_gene_name_from_synonym <- function( gnames ) {
	syn = read.table("gene_name_synonyms_grch38.tsv", stringsAsFactors=F, sep="\t", header=T, comment.char="", fill=T, quote="")
	a = data.frame(orig=gnames,name=syn$Gene.name[match(gnames,syn$Gene.name)],syn=syn$Gene.name[match(gnames,syn$Gene.Synonym)])
	nnames = a$name
	nnames[is.na(nnames)] = a$syn[is.na(nnames)]
	
#	syn = read.table("gene_name_synonyms_grch37.tsv", stringsAsFactors=F, sep="\t", header=T, comment.char="", fill=T, quote="")
#	a = cbind(a,data.frame(name37=syn$Gene.name[match(a$orig,syn$Gene.name)],syn37=syn$Gene.name[match(a$orig,syn$Gene.Synonym)]))
	
	return(nnames)
}

get_single_sample_gregersen_seurats <- function( update=F ) {
	meta = read.table("/omics/groups/OE0433/internal/references_data/GSE203191/GSE203191_Shih_endo_meta_frame.tsv", stringsAsFactors=F, header=T, sep="\t")
	
	basedir = "/omics/groups/OE0433/internal/references_data/GSE203191"
	files = dir(basedir, pattern=".h5", full.names=T )
	names(files) = sapply( strsplit( sapply(strsplit( sub("_filtered_feature_bc_matrix.h5","", files ),"/",fixed=T),"[[",8), "[0-9]_" ), "[[", 2 )
	
	gseu = list()
	for( i in 1:length(files) ) {
		file = files[i]
		cat( "\tat", file, "\n" )
		filename = sub(".h5", "_seurat.RData", file, fixed=T )
		
		if( !file.exists(filename) | update ) {	
			cat("\treading in h5\n")
			dat = Read10X_h5(filename = file, use.names = TRUE)
			bseu = CreateSeuratObject(counts = dat)
			rm(dat)
			gc()
			
			cat("\tadding metadata\n")
			smeta = meta[meta$run == names(files)[i],]
			bseu@meta.data = cbind(bseu@meta.data, smeta[match(rownames(bseu@meta.data), smeta$barcode),] )
			bseu@meta.data$run = names(files)[i]
			
			cat("\tpercent mitochondria\n")
			bseu = PercentageFeatureSet(bseu, pattern = "^MT-", col.name="percent.mt" )
			
			ifile = paste("scrna_qc_gregersen_", names(files)[i], ".pdf", sep="" )
			cat("\tplotting QC to", ifile, "\n")
			pdf(ifile, width=6, height=6 )
			
			a = data.frame( bseu@meta.data ) 
			p = ggplot( a, aes(nFeature_RNA,nCount_RNA, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title="before filtering") 
			print(p)
			
			cat("\tfilter cells\n")
			ind = a$nCount_RNA < quantile(a$nCount_RNA,probs=seq(0,1,0.01))[99] & 
					a$nFeature_RNA < quantile(a$nFeature_RNA,probs=seq(0,1,0.01))[99] &
					a$percent.mt < 10 &
					a$nFeature_RNA > 100
			bseu = bseu[,ind]
			
			cat("\tSCTransform, PCA, UMAP, Neighbours, Clusters\n")
			bseu = SCTransform(bseu, verbose = TRUE)
			bseu = RunPCA(bseu, verbose = TRUE)
			bseu = RunUMAP(bseu, dims = 1:30, verbose = TRUE)
			bseu = FindNeighbors(bseu, dims = 1:30, verbose = TRUE)
			bseu = FindClusters(bseu, verbose = TRUE, resolution=1)
			
			a = data.frame( bseu@meta.data )
			p = ggplot( a, aes(nFeature_RNA,nCount_RNA, color=percent.mt) ) + geom_point(alpha=0.3, shape=20) + theme_bw() + labs(title="after filtering") 
			print(p)
			dev.off()
			
			cat("\tsaving\n")
			save(bseu,file=filename)
			cat("\tsaved to",filename,"\n")
		} else{
			cat("\tloading from", filename, "\n")
			load(filename)
		}	
		gseu[[names(files)[i]]] = bseu
		rm(bseu)
		gc()
	}
	return(gseu)
}

gregersen_clusterID_to_cell_type <- function( ids ) {
	cell_type = rep(NA,length(ids))
	cell_type[ids %in% c("B","CD4T","CD8T1","CD8T2","CD8T3","CD8T4","Granulocyte","Myeloid1","Myeloid2","Myeloid3","uNK1","uNK2","pDC")] = "immune"
	cell_type[ids %in% c("EC-like")] = "EC"
	cell_type[ids %in% c("Epithelial1","Epithelial2","Epithelial3")] = "epithelial"
	cell_type[ids %in% c("Stromal")] = "stromal"
	
	return(cell_type)
}

gregersen_clusterID_to_cell_type_level_2 <- function( ids ) {
	cell_type = tolower(ids)
	cell_type[ids %in% c("CD4T","CD8T1","CD8T2","CD8T3","CD8T4")] = "T"
	cell_type[ids %in% c("Myeloid1","Myeloid2","Myeloid3")] = "myeloid"
	cell_type[ids %in% c("uNK1","uNK2")] = "uNK"
	cell_type[ids %in% c("EC-like")] = "EC"
	cell_type[ids %in% c("Epithelial1","Epithelial2","Epithelial3")] = "epithelial"
	cell_type[ids %in% c("Stromal")] = "stromal"
	
	return(cell_type)
}

get_pseudocounts_gregersen <- function( seu ) {
	dat = c()
	for( sample in names(seu) ) {
		seu[[sample]]@meta.data$cell_type = gregersen_clusterID_to_cell_type( seu[[sample]]@meta.data$clusterID )
		cts = unique(seu[[sample]]@meta.data$cell_type)
		cts = cts[!is.na(cts)]
		for( ct in cts ) {
			donors = unique(seu[[sample]]@meta.data$subjectID)
			donors = donors[!is.na(donors)]
			for( donor in donors ) {
				ind = seu[[sample]]@meta.data$cell_type %in% ct & seu[[sample]]@meta.data$subjectID %in% donor
				if( sum(ind) > 1 ) {
					a = rowSums(seu[[sample]][["RNA"]]@counts[,ind])
					dat = cbind(dat,a)
					colnames(dat)[ncol(dat)] = paste(sample,donor,ct,sep="-")		
				}
			}
		}
	}
	dat = dat[rowSums(dat) > 0,]
	
	return(dat)
}



get_single_sample_quake_seurats <- function() {
	filename = "human_endometrium_quake.RData"
	
	if( !file.exists(filename) ) {
		qseu = readRDS("/omics/groups/OE0433/internal/references_data/GSE111976_ct_endo_10x.rds")
		a = rowSums(qseu)
		qseu = qseu[ a > 0, ]
		
		gnames = get_gene_name_from_synonym(rownames(qseu))
		ind = !is.na(gnames)
		
		gnames = gnames[ind]
		qseu = qseu[ind,]
		rownames(qseu) = gnames
		
		ind = which( rownames(qseu) %in% rownames(qseu)[duplicated(rownames(qseu))] )
		a = qseu[ind,]
		ind = !(rownames(qseu) %in% rownames(a))
		qseu = qseu[ind,]
		
		a = a[order(rowSums(a),decreasing=T),]
		a = a[!duplicated(rownames(a)),]
		qseu = rbind(qseu,a)
		
		meta = read.table("/omics/groups/OE0433/internal/references_data/GSE111976_summary_10x_day_donor_ctype.csv",stringsAsFactors=F,fill=T,sep=",",header=T)
		rownames(meta) = meta$X
		meta = meta[,-1]
		phase = read.table("/omics/groups/OE0433/internal/references_data/GSE111976_summary_10x_donor_phase.csv",stringsAsFactors=F,fill=T,sep=",",header=T)
		meta$phase = phase[match(meta$donor,phase$donor),"phase_canonical"]
		meta$donor = as.character(meta$donor)
		
		seu = list()
		for( donor in unique(meta$donor) ) {
			cat( "\tat", donor, "\n" )
			seu[[donor]] = CreateSeuratObject( qseu[,meta$donor == donor], 
					project = "HumanBiopsies", assay = "RNA", min.cells = 0, min.features = 0,
					names.field = 1, names.delim = "_", 
					meta.data = meta[meta$donor == donor,] )
			seu[[donor]] = PercentageFeatureSet(seu[[donor]], pattern = "^MT-", col.name="percent.mt" )
			seu[[donor]] = SCTransform(seu[[donor]], verbose = TRUE)
		}
		
		for( donor in unique(meta$donor) ) {
			seu[[donor]] = RunPCA(seu[[donor]], verbose = TRUE)
			seu[[donor]] = RunUMAP(seu[[donor]], dims = 1:30, verbose = TRUE)
		}
		
		for( donor in unique(meta$donor) ) {
			bseu = seu[[donor]]
			save(bseu,file=paste("human_endometrium_quake_",donor,".RData",sep=""))
		}
		
		cat("\tsaving to", filename,"\n")
		save( seu, file=filename )
	} else {
		cat("\tloading from", filename,"\n")
		load(filename)
	}
	
	return(seu)
}

get_integrated_mf_quake <- function() {
	filename = "integrated_mf_quake.RData"
	filename2 = "integrated_mf_quake_markers.RData"
	
	if( !file.exists(filename) ) {
		# select quake samples
		qsamples = c("57","19","63")
		qseu = get_single_sample_quake_seurats()
		qseu = qseu[qsamples]
		
		samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON020afM","MFCON007dfM") 
		mseu = get_single_sample_cellbender_seurats( samples )[[1]]
		
		seu = c(qseu,mseu)
		seu = seu[c(qsamples,samples)]
		rm(qseu)
		rm(mseu)
		
		seu.features = SelectIntegrationFeatures(object.list = seu, nfeatures = 3000)
		
		cat("\tpreping integration\n")
		seu = PrepSCTIntegration(object.list = seu, anchor.features = seu.features, verbose = TRUE)
		
		cat("\tfinding anchors\n")
		seu.anchors = FindIntegrationAnchors(object.list = seu, normalization.method = "SCT", anchor.features = seu.features, 
				reference = which(names(seu) == "19"),verbose = TRUE)
		
		cat("\tintegrating\n")
		seu.integrated = IntegrateData(anchorset = seu.anchors, normalization.method = "SCT", verbose = TRUE)
		
		seu.integrated = RunPCA(seu.integrated, verbose = TRUE)
		seu.integrated = RunUMAP(seu.integrated, dims = 1:30)
		
		seu.integrated = FindNeighbors(seu.integrated, dims = 1:30, verbose = TRUE)
		seu.integrated = FindClusters(seu.integrated, verbose = TRUE, resolution=0.2)
		
		cat("\tsaving to", filename,"\n")
		save( seu.integrated, file=filename )
		
		markers = FindAllMarkers(seu.integrated, verbose = TRUE)
		
		save( markers, file=filename2 )
	} else {
		load(filename)
		load(filename2)
	}
	return(list(seu.integrated,markers))
}

# read in EPFL data
get_endoCounts <- function( normalized=FALSE ) {
	# read in endometriosis count table
	cat("\treading counts\n")
	tab = read.table("STAR-HTSeqTags.-Jun18.counts.txt", check.names=F)
	tab = tab[rowSums(tab) > 0,]
	
	# create metadata from column names
	cat("\tcreating metadata\n")
	meta = data.frame( id=sapply(strsplit(colnames(tab),".",fixed=T),"[[",1),
			endo_stage=as.numeric(sapply(strsplit(colnames(tab),".",fixed=T),"[[",2)),
			contraceptive=sapply(strsplit(colnames(tab),".",fixed=T),"[[",3),
			cell_type=sapply(strsplit(colnames(tab),".",fixed=T),"[[",4),
			index=as.numeric(sapply(strsplit(colnames(tab),".",fixed=T),"[[",5)) )
	meta$endo_stage[is.na(meta$endo_stage)] = 0
	meta$endo = meta$endo_stage != 0
	meta$endo_stage = factor(meta$endo_stage)
	rownames(meta) = colnames(tab)
	#a = read.table("ex_vivo_samples.txt",stringsAsFactors=F)[,1]
	meta$type = "ex_vivo"
	
	a = read.table("in_vitro_samples.txt",stringsAsFactors=F)[,1]
	meta[a,"type"] = "in_vitro"
	
	cat("\tgetting gene symbols\n")
	annot = getAnnotation()
	
	a = annot[match(rownames(tab),annot$hs_ensembl_gene_id),"gene"]
	ind = !is.na(a)
	tab = tab[ind,]
	a = a[ind]
	
	ind = !duplicated(a)
	tab = tab[ind,]
	a = a[ind]
	
	rownames(tab) = a
	
	if( normalized ) {
		cat("\tnormalizing\n")
		library("DESeq2")
		
		de = DESeqDataSetFromMatrix(countData = tab, colData = meta, design= ~ endo )
		de = estimateSizeFactors(de)
		tab = counts(de, normalized=TRUE)
	}
	
	return(list(counts=tab,meta=meta))
}

getAnnotation <- function( update=F ) {
	filename = "2023_annotation.RData"
	
	if( file.exists(filename) & !update ) {
		cat("loading annotation from", filename, "\n")
		load( filename )
	} else {
		library("biomaRt")
		# human 
		ensembl = useMart(host="http://jul2018.archive.ensembl.org/", "ENSEMBL_MART_ENSEMBL")
		ensembl = useDataset(mart=ensembl, "hsapiens_gene_ensembl")
		
		cat("downloading human annotation\n")
		annot = getBM(mart=ensembl, attributes = c("chromosome_name", "ensembl_gene_id", "external_gene_name", "strand", 
						"start_position", "end_position", "gene_biotype"))
		annot$chromosome_name = paste("chr",annot$chromosome_name,sep="")
		annot = annot[annot$chromosome_name %in% paste("chr",c(1:22,"X","Y"),sep=""),]
		colnames(annot) = paste("hs_",colnames(annot),sep="")
		annot$gene = annot$hs_external_gene_name
		
		cat("\tsaving to", filename,"\n")
		save(annot,file=filename)
	}
	
	return(annot)
}

get_single_sample_tormo_seurats <- function() {
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
		load( file )
	}
	return(seu)
	
}

# get quake data
get_integrated_human_data_quake <- function() {
	filename = "human_endometrium_integrated_quake.RData"
	
	if( !file.exists(filename) ) {
		seu = get_single_sample_quake_seurats()
		
		seu.features = SelectIntegrationFeatures(object.list = seu, nfeatures = 3000)
		
		cat("\tpreping integration\n")
		seu = PrepSCTIntegration(object.list = seu, anchor.features = seu.features, verbose = TRUE)
		
		cat("\tfinding anchors\n")
		#seu.anchors = FindIntegrationAnchors(object.list = seu, normalization.method = "SCT", anchor.features = seu.features, 
		#		reference = which(names(seu) == ),verbose = TRUE)
		
		cat("\tintegrating\n")
		seu.integrated = IntegrateData(anchorset = seu.anchors, normalization.method = "SCT", verbose = TRUE)
		
		seu.integrated = RunPCA(seu.integrated, verbose = TRUE)
		seu.integrated = RunUMAP(seu.integrated, dims = 1:30)
		
		seu.integrated = FindNeighbors(seu.integrated, dims = 1:30, verbose = TRUE)
		seu.integrated = FindClusters(seu.integrated, verbose = TRUE, resolution=0.2)
		markers = FindAllMarkers(seu.integrated, verbose = TRUE)
		cat("\tsaving to", filename,"\n")
		save( seu.integrated, markers, file=filename )
	} else {
		cat("\tloading from", filename,"\n")
		load(filename)
	}
	
	return(list(seu.integrated ,markers))
}



# reference based CCA
get_integrated_cellbender_seurat <- function( samples, update=F ) {
	samples = sort(samples)
	filename = paste( "integrated_cellbender_seurat_", paste(samples,collapse="_"),".RData", sep="" )
	
	if( !file.exists(filename) | update ) {
		cat("\tloading objects\n")
		seu.list = get_single_sample_cellbender_seurats( samples )[[1]]
		cat("\tselecting integration features\n")
		seu.features = SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)
		cat("\tpreping integration\n")
		seu.list = PrepSCTIntegration(object.list = seu.list, anchor.features = seu.features, verbose = TRUE)
		cat("\tfinding anchors\n")
		seu.anchors = FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", anchor.features = seu.features, 
				reference = which(names(seu.list) == "MFCON007dcM"),verbose = TRUE)
		cat("\tintegrating\n")
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
	
	seu.integrated@meta.data$donor = substr(seu.integrated@meta.data$sample,1,9)
	seu.integrated@meta.data$type = "flow"
	seu.integrated@meta.data$type[substr( seu.integrated@meta.data$sample, 10, 10) == "c"] = "clumps"
	for( res in c(0.1,0.4,0.6,0.8,1) ) {
		seu.integrated = FindClusters(seu.integrated, verbose = TRUE, resolution=res)
	}
	seu.integrated@meta.data$cell_type_level1 = NA
	ind = as.character(seu.integrated@meta.data$integrated_snn_res.0.2)
	seu.integrated@meta.data$cell_type_level1[ind %in% c("0","1","5")] = "epithelial"
	seu.integrated@meta.data$cell_type_level1[ind %in% c("2","3","4")] = "stromal"
	seu.integrated@meta.data$cell_type_level1[ind %in% c("6")] = "leukocytes"
	
	seu.integrated@meta.data$cell_type_level2 = NA
	ind = as.character(seu.integrated@meta.data$integrated_snn_res.0.2)
	seu.integrated@meta.data$cell_type_level2[ind == "0"] = "glands"
	seu.integrated@meta.data$cell_type_level2[ind == "1"] = "exhausted glands"
	seu.integrated@meta.data$cell_type_level2[ind == "5"] = "cilliated glands"
	seu.integrated@meta.data$cell_type_level2[ind == "2"] = "decidualised stroma"
	seu.integrated@meta.data$cell_type_level2[ind == "3"] = "stroma/mesenchyme"
	seu.integrated@meta.data$cell_type_level2[ind == "4"] = "smooth muscle"
	seu.integrated@meta.data$cell_type_level2[ind == "6"] = "T-cells"
	
	return(list(seu.integrated,markers))
}

get_integrated_seurat <- function( samples, update=F ) {
	samples = sort(samples)
	filename = paste( "integrated_seurat_", paste(samples,collapse="_"),".RData", sep="" )
	
	if( !file.exists(filename) | update ) {
		cat("\tloading objects\n")
		seu.list = get_single_sample_seurats( samples )[[1]]
		cat("\tselecting integration features\n")
		seu.features = SelectIntegrationFeatures(object.list = seu.list, nfeatures = 3000)
		cat("\tpreping integration\n")
		seu.list = PrepSCTIntegration(object.list = seu.list, anchor.features = seu.features, verbose = TRUE)
		cat("\tfinding anchors\n")
		seu.anchors = FindIntegrationAnchors(object.list = seu.list, normalization.method = "SCT", anchor.features = seu.features, 
				reference = which(names(seu.list) == "MFCON007dcM"),verbose = TRUE, k.filter=103)
		cat("\tintegrating\n")
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
	return(list(seu.integrated,markers))
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

get_cellranger_qc_output <- function( samples ) {
	dm = datamap[ datamap$sample_full %in% samples, ]
	
	aqc = c()
	for( i in 1:nrow(dm) ) {
		file = paste( dm$base_dir[i],"outs/metrics_summary.csv", sep="")
		aqc = rbind(aqc, read.csv(file) )
	}
	rownames(aqc) = dm$sample_full
	
	return(aqc)
}
