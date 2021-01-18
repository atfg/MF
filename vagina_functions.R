
get_single_sample_vagina_seurats <- function(update=F) {
	filename ="li_single_sample_seurats.RData"
	
	if( !file.exists(filename) | update ) {	
		# SCT transform
		for( sample in samples ) {
			cat("\tat", sample, "\n")
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
	return(list(adat,markers))
}

get_integrated_endometrium_vagina_mf <- function(update=F) {
	filename = "integrated_endometrium_vagina_mf.RData"
	filename2 = "integrated_endometrium_vagina_mf_markers.RData"
	
	if( !file.exists(filename) | update ) {	
		qsamples = c("57","19","63")
		qseu = get_single_sample_human_data_quake()
		qseu = qseu[qsamples]
		
		samples = c("MFCON007dcM","MFCON020acM","MFCON007efM","MFCON010dfM","MFCON020afM","MFCON007dfM") 
		mseu = get_single_sample_cellbender_seurats( samples )[[1]]
		
		vseu = get_single_sample_vagina_seurats()[[1]]
		vsamples = names(vseu)
		
		seu = c(qseu,mseu,vseu)
		seu = seu[c(qsamples,samples,vsamples)]
		rm(qseu)
		rm(mseu)
		rm(vseu)
		
		seu.features = SelectIntegrationFeatures(object.list = seu, nfeatures = 3000)
		
		cat("\tpreping integration\n")
		seu = PrepSCTIntegration(object.list = seu, anchor.features = seu.features, verbose = TRUE)
		
		# reference based integration
		cat("\tfinding anchors\n")
		seu.anchors = FindIntegrationAnchors(object.list = seu, normalization.method = "SCT", anchor.features = seu.features, 
				reference = which(names(seu) == "19"),verbose = TRUE)
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
