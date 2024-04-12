###############################################################################
# Angela's reminders for starting jobs, ignore this comment block: 
# bsub -R "rusage[mem=20GB]" -Is bash 
# module unload R; module unload python; module load anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.1.0 gcc/7.2.0; R
###############################################################################

source("/omics/groups/OE0433/internal//angela/atfg_github/MF/initialise_env.R")

# Gregersen has:
# ~4000 epithelial cells / 872 normal 
# ~32000 immune cells / 11000 normal
# ~7000 stromal cells / 2000 normal

# we have 2000 epithelial, 91 immune, 2000 stromal 
meta = get_metadata_gregersen()

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
# Gregersen DE level 0
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

load("gregersen_DE_genes.RData")
sapply( ares, function(x) { table( x$padj < 0.1 ) } )

###############################################################################
# Gregersen DE level 3 
###############################################################################
ndat = get_pseudocounts_gregersen_level_3(seu)
save(ndat,file="gregersen_DE_level_3_counts.RData")

smeta = data.frame( run=sapply(strsplit(colnames(ndat),"-"),"[[",1), donor=sapply(strsplit(colnames(ndat),"-"),"[[",2), cell_type=sapply(strsplit(colnames(ndat),"-"),"[[",3) ) 
smeta = cbind(smeta,type=meta[match(smeta$donor,meta$subjectID),c("pheno")])

ares = list()
for( ct in unique(smeta$cell_type) ) {
	cat("at",ct,"\n")
	filename = paste0("gregersen_DE_",ct,".RData")
	
	if( !file.exists(filename) ) {
		ind = smeta$cell_type == ct & smeta$type %in% c("Control","Diagnosed")
		dat = ndat[,ind]
		ssmeta = smeta[ind,]
		if( length(unique(ssmeta$type)) > 1 & all(table(ssmeta$type) > 1) ) {
			# filtering is quite important
			dat = dat[rowSums( dat[,ssmeta$type == "Control"] > 0 ) >= 3 | rowSums( dat[,ssmeta$type == "Diagnosed"] > 0 ) >= 3,]
			de = DESeqDataSetFromMatrix(countData = dat, colData = ssmeta, design= ~ type)
			de = DESeq(de)
			
			res = lfcShrink(de, coef="type_Diagnosed_vs_Control", type="apeglm")
			res = res[!is.na(res$padj),]
			res = res[order(res$log2FoldChange,decreasing=T),]
			res$gene = rownames(res)
			
			save( res, file=filename )
			write.table(res[,c("gene","log2FoldChange")],file=sub(".RData",".rnk",filename),quote=F,row.names=F,col.names=F,sep="\t")
		}
	} else { 
		load(filename)
	}
	ares[[ct]] = res
}
save(ares,ndat,file="gregersen_DE_genes_level_3.RData")

load("gregersen_DE_genes_level_3.RData")

sres = sapply( ares, function(x) { x[order(x$log2FoldChange),] } )
sres = sapply( ares, function(x) { x[x$padj < 0.1,] } )

sort(sapply( ares, function(x) { sum( x$padj < 0.1 ) } ))

sapply( ares, function(x) { x["IL6R",] } )

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
# compare gregersen DE to our DE - level 1
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

###############################################################################
# compare gregersen DE to our DE - level 3
###############################################################################

load("gregersen_DE_genes_level_3.RData")
load("DE_results_leukocytes_ex_vivo.RData")

gene = "CXCL8"
gene = "HIF1A"
gene = "LCP2"
gene = "IL1R1"
gene = "IL6R"
a = as.data.frame(t(sapply(ares, function(x) { as.data.frame(x[rownames(x) %in% gene,c("log2FoldChange","pvalue")]) } )))
a$sig = a$pvalue < 0.05
a


###############################################################################
# score Gregersen by IL6 
###############################################################################
meta = get_metadata_gregersen()
seu = get_single_sample_gregersen_seurats( update=F )

# get expression quantiles
for( i in 1:length(seu) ) {
	cat("at", names(seu)[i], "\n")
	
}

# load IL6 pathway genes
genesets = list(il6=read.table( "HALLMARK_IL6_JAK_STAT3_SIGNALING.gmt", fill=T, skip=2 )[,1],
		inf=read.table( "HALLMARK_INFLAMMATORY_RESPONSE.txt", fill=T, skip=2 )[,1] )

library("AUCell")
maxrank = 800

scores = c()
for( i in 1:length(seu) ) {
	cat("at", names(seu)[i], "\n")
	
	filename = paste("gregersen_aucell_cell_rankings_", names(seu)[i], ".RData", sep="")
	if( !file.exists(filename) ) {
		cell_rankings = AUCell_buildRankings(seu[[i]][["SCT"]]@data, plotStats=F )
		cat("saving to", filename, "\n")
		save(cell_rankings,file=filename)
	} else {
		cat("loading from", filename,"\n")
		load(filename)
	}
	cells_AUC = AUCell_calcAUC(genesets, cell_rankings, aucMaxRank=maxrank)
	scores = rbind(scores, data.frame( seu[[i]]@meta.data, score_il6=getAUC(cells_AUC)["il6",], score_inf=getAUC(cells_AUC)["inf",] ) )
}
scores$level_1 = gregersen_clusterID_to_cell_type(scores$clusterID)
scores$level_2 = gregersen_clusterID_to_cell_type_level_2(scores$clusterID)
scores$endo = !(scores$pheno %in% c("Control"))
scores = scores[!is.na(scores$level_2),]

save(scores,file="gregersen_aucell_il6.RData")

###############################################################################
# IL6 scores in immune cells
###############################################################################

load("gregersen_aucell_il6.RData")
scores$level_3 = gregersen_get_level3( scores )
scores$score = scores$score_il6
save(scores,file="gregersen_temp.RData")

load("~/Downloads/gregersen_temp.RData")

sscores = scores[ scores$level_1 == "immune",]
sscores = sscores[!is.na(sscores$level_2),]
sscores$level_2 = factor(sscores$level_2, levels=c("myeloid","pdc","granulocyte","uNK","T","b"),ordered=T)

pdf( "gregersen_aucell_il6.pdf", width=4*1.5, height=4)
ggplot(sscores,aes(level_2,score,fill=pheno)) + geom_violin() + theme_classic() 
ggplot(sscores,aes(level_2,score,fill=pheno)) + geom_boxplot() + theme_classic() 
dev.off()

pdf( "gregersen_aucell_il6_clusterID.pdf", width=4, height=4*1.5)
ggplot(sscores[ sscores$pheno %in% c("Control","Diagnosed"),],aes(clusterID,score,fill=pheno)) + geom_boxplot() + theme_classic() + coord_flip()
ggplot(sscores,aes(clusterID,score,fill=pheno)) + geom_boxplot() + theme_classic() + coord_flip()
ggplot(sscores,aes(clusterID,score,fill=endo)) + geom_boxplot() + theme_classic() + coord_flip()
dev.off()

pdf( "gregersen_aucell_il6_level_3.pdf", width=4, height=4*1.5)
ggplot(scores[ scores$pheno %in% c("Control","Diagnosed"),],aes(level_3,score,fill=pheno)) + geom_boxplot() + theme_classic() + coord_flip()
ggplot(scores,aes(level_3,score,fill=pheno)) + geom_boxplot() + theme_classic() + coord_flip()
ggplot(scores,aes(level_3,score,fill=endo)) + geom_boxplot() + theme_classic() + coord_flip()

ggplot(scores[ scores$pheno %in% c("Control","Diagnosed"),],aes(level_3,score_inf,fill=pheno)) + geom_boxplot() + theme_classic() + coord_flip()
ggplot(scores,aes(level_3,score_inf,fill=pheno)) + geom_boxplot() + theme_classic() + coord_flip()
ggplot(scores,aes(level_3,score_inf,fill=endo)) + geom_boxplot() + theme_classic() + coord_flip()
dev.off()

summary(lm( score ~ pheno, data=sscores[sscores$level_2 == "myeloid",] ))

library("lmerTest")
summary( lmer( score ~ pheno + (1 | clusterID) + (1 | subjectID), data=sscores[sscores$level_2 == "myeloid",] ) )

summary( lmer( score ~ pheno + (1 | clusterID) + (1 | subjectID), data=sscores[sscores$level_2 == "myeloid",] ) )


sscores$clusterID = relevel(factor(sscores$clusterID),ref="Myeloid2")
summary( lmer( score ~ endo * clusterID + (1 | subjectID), data=sscores ) )

for( ct in unique(as.character(sscores$clusterID)) ) {
	cat(ct, coef(summary( lmer( score ~ endo + (1 | subjectID), data=sscores[sscores$clusterID == ct,] ) ))["endoTRUE","Pr(>|t|)"], "\n" )
}

ce = c()
for( ct in unique(as.character(sscores$clusterID)) ) {
	ce = rbind(ce, coef(summary( lm( log2(score+0.0000001) ~ endo, data=sscores[sscores$clusterID == ct,] ) ))["endoTRUE",c("Estimate","Pr(>|t|)")] )
}
rownames(ce) = unique(as.character(sscores$clusterID))
ce = ce[order(ce[,2]),]
ce = as.data.frame(ce)
ce$padj = p.adjust(ce[,2],method="BH")
ce$sig = ce$padj < 0.05
ce$cell_type = rownames(ce)
ce = ce[order(ce[,1],decreasing=T),]
ce$cell_type = factor(ce$cell_type,levels=ce$cell_type,ordered=T)
pdf("~/Downloads/gregersen_aucell_il6_lolipops.pdf", width=4, height=3)
ggplot( ce, aes( Estimate, cell_type,color=-log10(padj) ) ) + geom_point() + theme_classic() + geom_vline( xintercept=0 ) + geom_segment(aes(x=0,y=cell_type,xend=Estimate,yend=cell_type))
dev.off()


sscores$clusterID = relevel(factor(sscores$clusterID),ref="Myeloid1")

###############################################################################
# Inflammation score lolipops
###############################################################################

load("~/Downloads/gregersen_temp.RData")

sscores = scores

ce = c()
for( ct in unique(as.character(sscores$level_3)) ) {
	dat = sscores[sscores$level_3 == ct,]
	ss = table( dat$endo )
	if( length(ss) == 2 & all(ss > 1) ) {
		ce = rbind(ce, coef(summary( lm( log2(score_inf+0.0000001) ~ endo, data=dat ) ))["endoTRUE",c("Estimate","Std. Error","Pr(>|t|)")] )	
		rownames(ce)[nrow(ce)] = ct
	}
}
ce = ce[order(ce[,2]),]
ce = as.data.frame(ce)
ce$padj = p.adjust(ce[,2],method="BH")
ce$sig = ce$padj < 0.05
ce$cell_type = rownames(ce)
ce = ce[order(ce[,1],decreasing=T),]
ce$cell_type = factor(ce$cell_type,levels=ce$cell_type,ordered=T)
colnames(ce) = c("Estimate","se","pval","padj","sig","cell_type")
ce$level_1 = scores$level_1[match(ce$cell_type,scores$level_3)]
ce$level_2 = scores$clusterID[match(ce$cell_type,scores$level_3)]
ce$level_1b = ce$level_1
ce$level_1b[ce$level_2 %in% c("CD8T1","CD8T2","CD8T3","CD8T4","B","uNK1","uNK2","uNK3","CD4T")] = "lymphoid"
ce$level_1b[ce$level_2 %in% c("Granulocyte","Myeloid1","Myeloid2","Myeloid3","pDC")] = "myeloid"
ce = ce[!is.na(ce$level_1) & ce$cell_type != "Unknown",]

pdf("~/Downloads/gregersen_aucell_lolipops.pdf", width=4, height=3)
ggplot( ce, aes( Estimate, cell_type,color=-log10(padj) ) ) + geom_point() + theme_classic() + geom_vline( xintercept=0 ) + geom_segment(aes(x=0,y=cell_type,xend=Estimate,yend=cell_type))
for( ct in unique(ce$level_1b) ) {
	p = ggplot( ce[ ce$level_1b == ct,], aes( Estimate, cell_type, color=-log10(padj) ) ) + geom_point() + theme_classic() + geom_vline( xintercept=0 ) + geom_errorbar(aes(xmin=Estimate-se, xmax=Estimate+se), width=.2) 
	print(p)
	
}
dev.off()


###############################################################################
# IL6 scores in all cell types
###############################################################################

load("gregersen_aucell_il6.RData")
sscores = scores

pdf( "gregersen_aucell_il6_clusterID_all_types.pdf", width=4, height=4*1.5)
ggplot(sscores,aes(clusterID,score_il6,fill=endo)) + geom_boxplot() + theme_classic() + coord_flip()
ggplot(sscores,aes(clusterID,score_inf,fill=endo)) + geom_boxplot() + theme_classic() + coord_flip()
dev.off()


summary( lm( score ~ endo * clusterID, data=sscores ) )