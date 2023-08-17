
get_joined_quake_vagina <- function( qsamples=c("57"), vsamples=c("N1"), update=F ) {
	samples = sort(c(qsamples,vsamples))
	filename = paste( "joined_quake_vagina_seurat_", paste(samples,collapse="_"),".RData", sep="" )
	
	if( !file.exists(filename) | update ) {
		vdat = c()
		vmeta = c()
		ct = read.table( "/icgc/dkfzlsdf/analysis/B210/references_data/2020_li/celltypes.tsv", sep="\t", stringsAsFactors=F )	
		for( sample in vsamples ) {
			cat("\tat", sample, "\n")
			dat = Read10X(data.dir = paste("/icgc/dkfzlsdf/analysis/B210/references_data/2020_li/", sample, "/", sep="") )
			vdat = cbind(vdat,dat)
			
			sct = ct[ ct$sample == sample, ]
			sct = sct[ match(colnames(dat), sct$barcode), ]
			meta = data.frame( day=NA, donor=sample, cell_type=sct$cell.types, cell_name=sct$barcode, phase=NA )
			vmeta = rbind(vmeta, meta)
		}
		vmeta$cell_name = colnames(vdat)
		
		qseu = readRDS("/icgc/dkfzlsdf/analysis/B210/references_data/GSE111976_ct_endo_10x.rds")
		meta = read.table("/icgc/dkfzlsdf/analysis/B210/references_data/GSE111976_summary_10x_day_donor_ctype.csv",stringsAsFactors=F,fill=T,sep=",",header=T)
		rownames(meta) = meta$X
		meta = meta[,-1]
		phase = read.table("/icgc/dkfzlsdf/analysis/B210/references_data/GSE111976_summary_10x_donor_phase.csv",stringsAsFactors=F,fill=T,sep=",",header=T)
		meta$phase = phase[match(meta$donor,phase$donor),"phase_canonical"]
		meta$donor = as.character(meta$donor)
		
		qdat = qseu[,meta$donor %in% qsamples]
		qmeta = meta[meta$donor %in% qsamples,]
		
		vmeta$type = "vagina"
		qmeta$type = "endometrium"
		
		colnames(vdat) = rownames(vmeta) = paste(vmeta$cell_name,vmeta$type,vmeta$donor,sep="_")
		colnames(qdat) = rownames(qmeta) = paste(qmeta$cell_name,qmeta$type,qmeta$donor,sep="_")
		meta = rbind(vmeta,qmeta)
		
		meta[ meta$cell_type %in% c("B cells","T cells","Lymphocytes","Plasma B cells","Lymphocytes"), "cell_type"] = "lymphocyte" 
		meta[ meta$cell_type %in% c("Endothelia","Endothelial cells"), "cell_type" ] = "endothelial"
		meta[ meta$cell_type %in% c("Fibroblasts"), "cell_type" ] = "vaginal fibroblasts"
		meta[ meta$cell_type %in% c("Stromal fibroblasts"), "cell_type" ] = "endometrial fibroblasts"
		meta[ meta$cell_type %in% c("Ciliated"), "cell_type" ] = "endometrial ciliated epithelial"
		meta[ meta$cell_type %in% c("Unciliated epithelia 1","Unciliated epithelia 2"), "cell_type" ] = "endometrial unciliated epithelial"
		meta[ meta$cell_type %in% c("Smooth muscle cells"), "cell_type" ] = "smooth muscle"
		meta[ meta$cell_type %in% c("Macrophage","Macrophages"), "cell_type" ] = "macrophage"
		meta[ meta$cell_type %in% c("Epithelial cells"), "cell_type" ] = "vaginal epithelial"
		table(meta[,c("cell_type","type")])
		
		genes = rownames(vdat)[rownames(vdat) %in% rownames(qdat)]
		dat = cbind(vdat[genes,],qdat[genes,])
		rm(qdat);rm(vdat)
		gc()
		
		seu = CreateSeuratObject(counts = dat, meta.data=meta )
		
		seu = PercentageFeatureSet(seu, pattern = "^MT-", col.name="percent.mt" )
		seu = SCTransform(seu, verbose = TRUE)
		seu = RunPCA(seu, verbose = TRUE)
		seu = RunUMAP(seu, dims = 1:30, verbose = TRUE)
		seu = FindNeighbors(seu, dims = 1:30, verbose = TRUE)
		seu = FindClusters(seu, verbose = TRUE, resolution=0.2)
		seu = FindClusters(seu, verbose = TRUE, resolution=1)
		
		save( seu, file=filename )
	} else{
		cat("\tloading from", filename, "\n")
		load(filename)
	}
	return(seu)
}

get_single_sample_vagina_seurats <- function(samples, update=F) {
	filename ="li_single_sample_seurats.RData"
	
	if( !file.exists(filename) | update ) {	
		# SCT transform
		adat = list()
		markers = list()
		for( sample in samples ) {
			cat("\tat", sample, "\n")
			dat = Read10X(data.dir = paste("/icgc/dkfzlsdf/analysis/B210/references_data/2020_li/", sample, "/", sep="") )	
			adat[[sample]] = CreateSeuratObject(counts = dat)
			
			bseu = CreateSeuratObject(counts = dat)
			
			bseu = adat[[sample]]
			bseu@meta.data$donor = sample
			bseu = PercentageFeatureSet(bseu, pattern = "^MT-", col.name="percent.mt" )
			
			a = data.frame( bseu@meta.data ) 
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
			
			markers[[sample]] = bmarkers
			adat[[sample]] = bseu
		}
		
		ct = read.table( "/icgc/dkfzlsdf/analysis/B210/references_data/2020_li/celltypes.tsv", sep="\t", stringsAsFactors=F )
		for( sample in samples ) {
			sct = ct[ ct$sample == sample, ]
			sct = sct[ match(rownames(adat[[sample]]@meta.data), sct$barcode), ]
			adat[[sample]]@meta.data$li_celltype = sct$cell.types
		}	
		save(adat,markers,file=filename)	
	} else {
		cat("\tloading from", filename, "\n")
		load(filename)
	}
	
	for( sample in samples ) {
		adat[[sample]]@meta.data$cell_type = adat[[sample]]@meta.data$li_celltype
		adat[[sample]]@meta.data$cell_type[ adat[[sample]]@meta.data$li_celltype %in% c("Fibroblasts")] = "Stromal fibroblasts"
		adat[[sample]]@meta.data$cell_type[ adat[[sample]]@meta.data$li_celltype %in% c("B cells","Plasma B cells","T cells")] = "Lymphocytes"
		adat[[sample]]@meta.data$cell_type[ adat[[sample]]@meta.data$li_celltype == "Macrophage"] = "Macrophages"
		adat[[sample]]@meta.data$cell_type[ adat[[sample]]@meta.data$li_celltype == "Endothelial cells"] = "Endothelia"
	}
	
	return(list(adat,markers))
}

# sctransform of joined matrices for endometrium and vagina
get_joined_endometrium_vagina <- function( update=F ) {
	filename = "sctransform_endometrium_vagina_mf.RData" 
	
	if( !file.exists(filename) | update ) {
		cat("\tloading quake data\n")
		qseu = readRDS("/icgc/dkfzlsdf/analysis/B210/references_data/GSE111976_ct_endo_10x.rds")
		meta = read.table("/icgc/dkfzlsdf/analysis/B210/references_data/GSE111976_summary_10x_day_donor_ctype.csv",stringsAsFactors=F,fill=T,sep=",",header=T)
		rownames(meta) = meta$X
		meta = meta[,-1]
		phase = read.table("/icgc/dkfzlsdf/analysis/B210/references_data/GSE111976_summary_10x_donor_phase.csv",stringsAsFactors=F,fill=T,sep=",",header=T)
		meta$phase = phase[match(meta$donor,phase$donor),"phase_canonical"]
		meta$donor = as.character(meta$donor)
		
		donor = "57"
		qdat = qseu[,meta$donor == donor]
		qmeta = meta[meta$donor == donor,]
		qmeta$type = "endometrium"
		
		cat("\tloading vagina data\n")
		sample = "N1"
		vdat = Read10X(data.dir = paste("/icgc/dkfzlsdf/analysis/B210/references_data/2020_li/", sample, "/", sep="") )	
		ct = read.table( "/icgc/dkfzlsdf/analysis/B210/references_data/2020_li/celltypes.tsv", sep="\t", stringsAsFactors=F )
		ct = ct[ct$sample == sample,]
		vmeta = data.frame(day=NA,donor=rep("N1",ncol(vdat)),cell_type=ct[match(colnames(vdat),ct$barcode),"cell.types"], cell_name=colnames(vdat), phase=NA)
		vmeta$type = "vagina"
		rownames(vmeta) = colnames(vdat)
		
		genes = rownames(qdat)[rownames(qdat) %in% rownames(vdat)]
		qdat = qdat[genes,]
		vdat = vdat[genes,]
		meta = rbind(qmeta,vmeta)
		
		cat("\tsctransform\n")
		seu = CreateSeuratObject(counts = cbind(qdat,vdat), min.cells = 0, min.features = 0, meta.data=meta )
		seu = PercentageFeatureSet(seu, pattern = "^MT-", col.name="percent.mt" )
		seu = SCTransform(seu, verbose = TRUE)
		seu = RunPCA(seu, verbose = TRUE)
		seu = RunUMAP(seu, dims = 1:50, verbose = TRUE)
		
		seu = FindNeighbors(seu, dims = 1:50, verbose = TRUE)
		seu = FindClusters(seu, verbose = TRUE, resolution=0.2)
	} else {
		load(filename)
	}
	
	return( seu )
}

# integrate 3 endometriums, all vaginas, all MFs; using one endometrium and one vagina as reference for CCA
# return a list with a seurat object and the list of cluster markers
get_integrated_endometrium_vagina_mf <- function(update=F) {
	filename = "integrated_endometrium_vagina_mf.RData"
	filename2 = "integrated_endometrium_vagina_mf_markers.RData"
	
	if( !file.exists(filename) | update ) {	
		qsamples = c("57","19","63")
		qseu = get_single_sample_quake_seurats()
		qseu = qseu[qsamples]
		
		samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON020afM","MFCON007dfM") 
		#mseu = get_single_sample_cellbender_seurats( samples )[[1]]
		mseu = get_single_sample_cellbender_seurats( samples )[[1]]
		
		vseu = get_single_sample_vagina_seurats(paste("N",1:5,sep=""))[[1]]
		vsamples = names(vseu)
		
		seu = c(qseu,mseu,vseu)
		seu = seu[c(qsamples,samples,vsamples)]
		rm(qseu)
		rm(mseu)
		rm(vseu)
		gc()
		
		seu.features = SelectIntegrationFeatures(object.list = seu, nfeatures = 3000)
		
		cat("\tpreping integration\n")
		seu = PrepSCTIntegration(object.list = seu, anchor.features = seu.features, verbose = TRUE)
		
		# reference based integration
		cat("\tfinding anchors\n")
		seu.anchors = FindIntegrationAnchors(object.list = seu, normalization.method = "SCT", anchor.features = seu.features, 
				reference = which(names(seu) %in% c("57","N1")),verbose = TRUE)
		save(seu.anchors,file="integrated_endometrium_vagina_mf_anchors.RData")
		
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
		cat("\tloading from", filename, "\n")
		load(filename)
		load(filename2)
	}
	return(list(seu.integrated,markers))
}
