###############################################################################
# Angela's reminders for starting jobs, ignore this comment block: 
# bsub -R "rusage[mem=20GB]" -Is bash 
# module unload R; module unload python; module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.1.0 gcc/7.2.0; R
###############################################################################

options(width=230)
.libPaths(paste("/omics/groups/OE0433/internal/software/R-library-",paste(R.version$major,R.version$minor,sep="."),sep=""))

setwd("/omics/groups/OE0433/internal/angela/mf/")
library("Seurat")
library("ggplot2")
library("DESeq2")

meta = read.table("/omics/groups/OE0433/internal/references_data/GSE203191/GSE203191_Shih_endo_meta_frame.tsv", stringsAsFactors=F, header=T, sep="\t")
meta$cell_type = gregersen_clusterID_to_cell_type( meta$clusterID )
meta$cell_type_level_2 = gregersen_clusterID_to_cell_type_level_2( meta$clusterID )

seu = get_single_sample_gregersen_seurats( update=F )

dat = get_pseudocounts_gregersen( update=F )

# plot cell abundances


a = table( meta[,c("cell_type","subjectID")] )
a = as.data.frame(matrix(a, ncol = ncol(a), dimnames = dimnames(a)))
rowMeans(a)
a = apply(a,2,function(x) { x/sum(x) } )
rowMeans(a)
a = as.data.frame(as.table(a))
a$type = meta$pheno[match(a$Var2,meta$subjectID)]

pdf("gregersen_cell_proportions.pdf", width=3*4, height=4)
ggplot(a,aes(Var2,Freq,fill=Var1)) + geom_bar(position="stack", stat="identity") + facet_wrap( ~ type, scales = "free_x" ) + theme_classic()
dev.off()


###############################################################################
# Gregersen DE
# 26 DE genes in stroma at 10% FDR, 16 at 5%, enriched for markers in NK cells
# 0 DE in epithelial - but we don't have epithelial either
# 802 DE in leukocytes 
###############################################################################
ndat = get_pseudocounts_gregersen(seu)

smeta = data.frame( run=sapply(strsplit(colnames(ndat),"-"),"[[",1), donor=sapply(strsplit(colnames(ndat),"-"),"[[",2), cell_type=sapply(strsplit(colnames(ndat),"-"),"[[",3) ) 
smeta = cbind(smeta,type=meta[match(smeta$donor,meta$subjectID),c("pheno")])

ares = list()
for( ct in c("stromal","epithelial","immune") ) {
	cat("at",ct,"\n")
	
	ind = smeta$cell_type == ct & smeta$type %in% c("Control","Diagnosed")
 	dat = ndat[,ind]
	ssmeta = smeta[ind,]
	# filtering is quite important
	dat = dat[rowSums( dat[,ssmeta$type == "Control"] > 0 ) >= 3 | rowSums( dat[,ssmeta$type == "Diagnosed"] > 0 ) >= 3,]
	#dat = dat[rowSums(dat > 0) >= 3,]
	de = DESeqDataSetFromMatrix(countData = dat, colData = ssmeta, design= ~ type)
	de = DESeq(de)
	
	res = lfcShrink(de, coef="type_Diagnosed_vs_Control", type="apeglm")
	#res = as.data.frame(results(de))	
	res = res[!is.na(res$padj),]
	res = res[order(res$log2FoldChange,decreasing=T),]
	res$gene = rownames(res)

	ares[[ct]] = res
	write.table(res[,c("gene","log2FoldChange")],file=paste0("gregersen_DE_",ct,".rnk"),quote=F,row.names=F,col.names=F,sep="\t")
}
save(ares,file="gregersen_DE_genes.RData")

sapply( ares, function(x) { table( x$padj < 0.1 ) } )

###############################################################################
# Gregersen cell type proportions
###############################################################################

a = unlist( sapply( split( meta$clusterID, meta$subjectID ), function(x) { table(x)/length(x) } ))
a = data.frame( donor=sapply(strsplit(names(a),".",fixed=T),"[[",1), 
		cell_type = sapply(strsplit(names(a),".",fixed=T),"[[",2), freq=a )

pdf( "gregersen_cell_proportions_level2.pdf" )
ggplot(a,aes(donor,freq,fill=cell_type)) + geom_bar(position="stack", stat="identity") + theme_classic()
dev.off()

# immune only
smeta = meta[ meta$cell_type == "immune", ]
a = unlist( sapply( split( smeta$clusterID, smeta$subjectID ), function(x) { table(x)/length(x) } ))
a = data.frame( donor=sapply(strsplit(names(a),".",fixed=T),"[[",1), 
		cell_type = sapply(strsplit(names(a),".",fixed=T),"[[",2), freq=a )

pdf( "gregersen_cell_proportions_level2_immune.pdf" )
ggplot(a,aes(donor,freq,fill=cell_type)) + geom_bar(position="stack", stat="identity") + theme_classic()
dev.off()

sort(table( smeta$cell_type_level_2 )/nrow(smeta))

###############################################################################
# compare gregersen DE to our DE
###############################################################################

# scatter plot?
load("gregersen_DE_genes.RData")
load("DE_results_leukocytes_ex_vivo.RData")
load("DE_results_stromal_in_vitro.RData")

resg = as.data.frame(ares$immune)
colnames(resg) = paste0("g_",colnames(resg))
reso = res_leukocytes_ex_vivo
colnames(reso) = paste0("o_",colnames(reso))

resa = cbind( resg, reso[rownames(resg),] )
resa = resa[!is.na(resa$g_pvalue) & !is.na(resa$o_pvalue),]
# only keep significant in one or the other analysis
resa = resa[ resa$g_pvalue < 0.1 | resa$o_pvalue < 0.1,]
resa$sig = "not"
resa$sig[resa$g_pvalue < 0.1 & resa$o_pvalue > 0.1] = "greg only"
resa$sig[resa$g_pvalue > 0.1 & resa$o_pvalue < 0.1] = "ours only"
resa$sig[resa$g_pvalue < 0.1 & resa$o_pvalue < 0.1] = "both"

pdf( "DE_comparison_gregersen_ours.pdf")
ggplot(resa,aes(o_log2FoldChange,g_log2FoldChange,color=sig)) + geom_point()+ theme_classic()
ggplot(resa,aes(o_pvalue,g_pvalue,color=sig)) + geom_point()+ theme_classic()
dev.off()

resa_immune = resa
ind = resa_immune$sig %in% "both" & ((resa_immune$g_log2FoldChange > 0 & resa_immune$o_log2FoldChange > 0) | (resa_immune$g_log2FoldChange < 0 & resa_immune$o_log2FoldChange < 0))
a = rownames(resa_immune[ind,])
save( a, file="DE_comparison_gregersen_ours_overlapping_immune.RData" )

# IL6/JAK/STAT3 pathway enriched in ours:
a = resa_immune[c("CRLF2","IL3RA","IL1R1","TNFRSF1A","IFNGR1","IL18R1","HMOX1",
				"CXCL10","LTBR","IL10RB","MAP3K8","CSF1","PIK3R5","SOCS3","SOCS1",
				"CSF2RA","JUN","CSF3R","A2M","PIM1","IFNGR2","CXCL1","ACVR1B","STAT3",
				"MYD88","IL1R2","CCR1","IFNAR1","TNFRSF1B"),]
a = a[ !is.na(a$sig),]
a = a[order(a$g_pvalue),]
pdf( "DE_comparison_gregersen_ours_il6path.pdf")
ggplot(a,aes(o_log2FoldChange,g_log2FoldChange)) + geom_point() + xlim(c(-3,3)) + ylim(c(-1,1)) + theme_classic()
dev.off()

######

resg = as.data.frame(ares$stromal)
colnames(resg) = paste0("g_",colnames(resg))
reso = res_stromal_in_vitro
colnames(reso) = paste0("o_",colnames(reso))

resa = cbind( resg, reso[rownames(resg),] )
resa = resa[!is.na(resa$g_pvalue) & !is.na(resa$o_pvalue),]
resa = resa[ resa$g_pvalue < 0.1 | resa$o_pvalue < 0.1,]
resa$sig = "not"
resa$sig[resa$g_pvalue < 0.1 & resa$o_pvalue > 0.1] = "greg only"
resa$sig[resa$g_pvalue > 0.1 & resa$o_pvalue < 0.1] = "ours only"
resa$sig[resa$g_pvalue < 0.1 & resa$o_pvalue < 0.1] = "both"

resa_stromal = resa
