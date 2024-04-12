###############################################################################
# Angela's reminders for starting jobs, ignore this comment block: 
# bsub -R "rusage[mem=20GB]" -Is bash 
# module unload R; module unload python; module load python/3.7.0 anaconda3/2019.07 libpng/1.6.37 hdf5/1.8.18 python/3.6.1 gdal/3.0.2 R/4.1.0 gcc/7.2.0; R
###############################################################################

# this loads up libraries, sets up some global varibles and functions 
# and loads up the seurat objects created by Luca so takes a few minutes to run
source("/omics/groups/OE0433/internal//angela/atfg_github/MF/initialise_env.R")


###############################################################################
# Volume of EPFL data
###############################################################################

meta2 = get_metadata_epfl()

# number of low versus high volume samples
table(sapply( split( meta2$cells, paste( meta2$Sample.date, meta2$ID ) ), function(x) { all(x == "all") } ))

# number of individual MF samples
nrow(unique( meta2[,c("ID","Sample.date")] ))

a = meta2[ !is.na(meta2$Volume), ]
a = a[!duplicated(paste(a$ID,a$Sample.date)), ]

pdf("endo_volume.pdf", width=4.5, height=4)
p = ggplot(a,aes(ID,Volume,color=factor(Endo.stage))) + geom_jitter(shape=21,alpha=0.8,size=4,width = 0.25) + theme_classic()
print(p)
p = ggplot(a,aes(factor(Endo.stage),Volume,fill=factor(ID))) + geom_jitter(shape=21,color="black",size=3, width = 0.25) + theme_classic()
print(p)
p = ggplot(a,aes(factor(Endo.stage),Volume,fill=factor(ID))) + geom_jitter(shape=21,color="black",size=3, width = 0.25) + theme_classic()
p = p + scale_fill_manual(values = c("5821"="#ffffff", "9834"="#f0f0f0", "9989"="#d9d9d9", "6177"="#bdbdbd", 
				"4760"="#969696", "1114" = "#737373","8096" = "#525252","4859" = "#252525","6652"="#000000",
				"2403"="#fee5d9","2976"="#fcbba1","4602"="#fc9272","27"="#fb6a4a","9018"="#ef3b2c","2674"="#cb181d","4210"="#99000d") ) 
print(p)
p = ggplot(a,aes(factor(Age),Volume,fill=factor(ID))) + geom_jitter(shape=21,color="black",size=3) + theme_classic()
print(p)
dev.off()
   

library("lmerTest")
fit = lmer( Volume ~ Age + (1 | ID), data=a )
summary(fit)

one.way = aov(Volume ~ ID, data = a)
summary(one.way)

one.way = aov(Volume ~ ID, data = a[a$Endo.stage == "0",])
summary(one.way)


###############################################################################
# PCA of EPFL data
###############################################################################

dat = get_endoCounts(normalized=TRUE)
tab = dat$counts
meta = dat$meta
rm(dat)
gc()

pc = prcomp( t(tab), center=TRUE, scale=TRUE )
a = data.frame(meta,pc$x[,1:5])

summary(pc)

pdf("endo_all_PCA.pdf", width=4.5, height=4)
p = ggplot(a,aes(PC1,PC2,col=paste(endo,type))) + geom_point() + theme_classic()
p = p + scale_colour_manual(values = c("TRUE ex_vivo"="red", "FALSE ex_vivo"="orange", "TRUE in_vitro"="darkblue", "FALSE in_vitro"="lightblue") ) 
print(p)
ggplot(a,aes(PC1,PC2,col=cell_type)) + geom_point() + theme_classic()
ggplot(a,aes(PC2,PC3,col=cell_type,shape=type)) + geom_point() + theme_classic()
ggplot(a,aes(PC3,PC4,col=cell_type,shape=type)) + geom_point() + theme_classic()
p = ggplot(a[a$type == "ex_vivo",],aes(PC2,PC4,col=cell_type)) + geom_point(size=4, alpha=0.5) + theme_classic()
p = p + geom_point(aes(fill=cell_type), col="black", shape=21, size=2, data=a[a$type == "in_vitro",])
print(p)
ggplot(a[a$type == "ex_vivo",],aes(PC2,PC4,col=cell_type)) + geom_point(size=4, alpha=0.5) + theme_classic()
dev.off()

a = pc$rotation[,"PC1"]
a = sort(a,decreasing=T)
a = a[c(1:100,(length(a)-100):length(a))]
a = data.frame(gene=factor(names(a),ordered=T), loadings=a)
pdf("endo_all_PCA_PC1_loadings_barchart.pdf", width=3*4, height=4)
ggplot(a,aes(gene,loadings)) + geom_col() + theme_bw()
dev.off()


n = 100
a = pc$rotation[,"PC1"]
a = sort(a,decreasing=T)
a = names(a)[1:n]
write.table(a,file="endo_all_PCA_PC1_loadings_top_in_vitro.txt",quote=F,col.names=F,row.names=F,sep="\t")
a = pc$rotation[,"PC1"]
a = names(sort(a,decreasing=F))[1:n]
write.table(a,file="endo_all_PCA_PC1_loadings_top_ex_vivo.txt",quote=F,col.names=F,row.names=F,sep="\t")

a = pc$rotation[,"PC1"]
a = sort(a,decreasing=T)
write.table(a,file="endo_all_PCA_PC1_loadings.rnk",quote=F,col.names=F,row.names=T,sep="\t")

a = read.table( "/Users/filimon/Documents/projects/MF_noShare/Results/Manuscript/supplementary/PanglaoDB_Augmented_2021_table.txt", header=T, sep="\t" )
a = a[order(a$Adjusted.P.value),][a$Adjusted.P.value < 0.01,]
a = a[1:min(5,nrow(a)),]
pdf("~/Downloads/endo_all_PCA_PC1_loadings_top_ex_vivo_barchart.pdf", width=5, height=4)
ggplot( a, aes( factor(Term, levels=rev(unique(Term)), ordered=T ), -log10(Adjusted.P.value), fill=Odds.Ratio ) ) + geom_col() + coord_flip() 
dev.off()

###############################################################################
# DE of EPFL data
###############################################################################

library("ggrepel")

dat = get_endoCounts(normalized=FALSE)
tab = dat$counts
meta = dat$meta
rm(dat)
gc()

############################
# all samples together modelling cell type and in vivo vitro
# i.e. genes that are commonly disregulated in all cell types and in all conditions
# 2022 DE genes
############################
de = DESeqDataSetFromMatrix(countData = tab, colData = meta, design= ~ type + endo + cell_type )
de = DESeq(de)
resultsNames(de)
res = as.data.frame(results(de, contrast=list("endoTRUE")))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
sres = res[res$padj < 0.05,]
dat = counts(de,normalized=TRUE)
goi = rownames(sres)[1:10]
#dat = data.frame(meta,t(dat[goi,]))
res$gene = rownames(res)
sum(res$padj < 0.05)

res_all = res
dat_all = dat
save( res_all, dat_all, file="DE_results_all.RData" )

load("DE_results_all.RData" )
res = res_all
#dat = dat_all
res$gene = rownames(res)
sres = res[res$padj < 0.05,]
goi = rownames(sres)[1:10]

pdf("endo_DE_all_together.pdf", width=3, height=3)
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.5) + geom_hline(yintercept=-log10(0.05)) + theme_classic()
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
for( g in goi ) {
	p = ggplot(dat_all,aes(endo,log10(get(g)), color=type, shape=cell_type)) + geom_jitter(alpha=0.7,size=3) + theme_classic() + ggtitle(g) 
	print(p)
	p = ggplot(dat_all,aes(endo,log10(get(g)), fill=paste(cell_type,type))) + geom_boxplot() + theme_classic() + ggtitle(g) 
	print(p)
}
dev.off()

a = data.frame(gene=rownames(res),logfc=res$log2FoldChange)
a = a[order(a$logfc,decreasing=T),]
write.table(a,file="endo_DE_all_together.rnk",row.names=FALSE,col.names=FALSE, quote=F,sep="\t")

pdf("endo_DE_all_together_heatmap.pdf", width=6, height=6)
tab = log10(t(dat_all[,8:ncol(dat_all)])+1)
cc = c("black","red")[as.numeric(dat_all$endo)+1]
ind = order(dat_all$endo)
tab = tab[,ind]
cc = cc[ind]
heatmap.2(as.matrix(tab), trace="none",cexRow=1/log10(nrow(dat[goi,]))-0.3, col=colorRampPalette(c("blue","white","red"))(10),
		ColSideColors=cc, Colv=F )
dev.off()

# classifier endometriosis using:
# MF ex vivo all (8 vs 29) - i.e. mostly leukocytes
# MF ex vivo CD45+ (5 vs 9) - leukocytes - not worth it, too few?
# MF ex vivo CD45- (6 vs 7) - stromal, too few?
# MF in vitro all+CD45- (8 vs 19) - stromal 
# not enough samples anywhere for a classifier

# DE comparisons to make:
# to find genes that are systematically found to be DE in either in vitro and ex vivo and in 
# all cell type, we analysed all of them together - all other analysis should be included in this
# one then

############################
# "all" ex vivo: similar to the above
# 29 vs 8
# 358 DE genes
############################

ind = meta$cell_type == "all" & meta$type == "ex_vivo"
table( meta[ind,"endo"] )
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ endo )
de = DESeq(de)
res = as.data.frame(results(de))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
sum(res$padj < 0.05)
sres = res[res$padj < 0.05,]
dat = counts(de,normalized=TRUE)
goi = rownames(sres)[1:15]
dat = data.frame(colData(de),t(dat[goi,]))
res$gene = rownames(res)

res_all_ex_vivo = res
dat_all_ex_vivo = dat
save( res_all_ex_vivo, dat_all_ex_vivo, file="DE_results_all_ex_vivo.RData" )

pdf("endo_DE_all_ex_vivo.pdf", width=5, height=5)
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.5) + geom_hline(yintercept=-log10(0.05)) + theme_classic()
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
for( g in goi ) {
	p = ggplot(dat,aes(endo,log10(get(g)))) + geom_jitter(alpha=0.7,size=3) + theme_classic() + ggtitle(g) 
	print(p)
}
dev.off()


############################
# all mixtures of cell types ex vivo only
# 8 endo vs 29 ctrl
# 358 genes significant
############################

# all vitro 
# 4 endo vs 10 ctrl
# 9 genes significant
ind = meta$cell_type == "all" & meta$type == "in_vitro"
table( meta[ind,"endo"] )
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ endo )
de = DESeq(de)
res = as.data.frame(results(de))
res = res[!is.na(res$padj),]
res1 = res[order(abs(res$log2FoldChange),decreasing=T),]
sum(res1$padj < 0.05)
res1[ res1$padj < 0.05, ]
# does this overlap best with cd45- or leukocytes like in the pca?

############################
# stromal only, reproducible ex vivo and in vitro
# 10 endo vs 16 ctrl
# 9 genes significant
############################
ind = meta$cell_type == "Cd45m"
table( meta[ind,"endo"] )
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ type + endo )
de = DESeq(de)
resultsNames(de)
res = as.data.frame(results(de,contrast=list("endoTRUE")))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
sum(res$padj < 0.05)
res[ res$padj < 0.05, ]
res$gene = rownames(res)
dat = counts(de,normalized=TRUE)

res_stromal = res
dat_stromal = dat
save( res_stromal, dat_stromal, file="DE_results_stromal.RData" )

load("DE_results_stromal.RData" )
dec = unique(c(goi$stromal_decidualisation,goi$decidualisation_warren,goi$roser_decidual_stromal))
a = res_stromal[dec,]
a[ a$padj < 0.1,]


pdf("endo_DE_stromal_all.pdf", width=3, height=3)
goi = rownames(res[res$padj < 0.05,])[1:min(10,sum(res$padj < 0.05))]
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
for( g in goi ) {
	p = try( ggplot(dat,aes(endo,log10(get(g)), color=type, shape=cell_type)) + geom_jitter(alpha=0.7,size=3) + theme_classic() + ggtitle(g) )
	try( print(p) )
	p = try( ggplot(dat,aes(endo,log10(get(g)), fill=type)) + geom_boxplot() + theme_classic() + ggtitle(g) )
	try( print(p) )
}
dev.off()

############################
# stromal ex vivo only - still a mixture of stromal and epithelial
# 6 endo vs 7 ctrl
# 2 genes significant
############################
ind = meta$cell_type == "Cd45m" & meta$type == "ex_vivo"
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ endo )
de = DESeq(de)
res = as.data.frame(results(de))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
sum(res$padj < 0.05)
res$gene = rownames(res)

res_stromal_ex_vivo = res
dat_stromal_ex_vivo = dat
save( res_stromal_ex_vivo, dat_stromal_ex_vivo, file="DE_results_stromal_ex_vivo.RData" )

############################
# stromal vitro - this is the purest stromal population
# 4 endo vs 9 ctrl
# 256 genes significant
############################
ind = meta$cell_type == "Cd45m" & meta$type == "in_vitro"
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ endo )
de = DESeq(de)
res = as.data.frame(results(de,contrast=list("endoTRUE")))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
sum(res$padj < 0.05)
res[ res$padj < 0.05, ]
res$gene = rownames(res)
dat = counts(de,normalized=TRUE)

res_stromal_in_vitro = res
dat_stromal_in_vitro = dat
save( res_stromal_in_vitro, dat_stromal_in_vitro, file="DE_results_stromal_in_vitro.RData" )

load( "DE_results_stromal_in_vitro.RData" )

pdf("endo_DE_stromal_in_vitro.pdf", width=3, height=3)
res = res_stromal_in_vitro
dat = dat_stromal_in_vitro
goi = c("PENK","MLLT11","NME4","HAND2-AS1","S100A4","MMP3","FN1","PSG4","BPGM","PDE4B")
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)

res = res_stromal_in_vitro
dat = dat_stromal_in_vitro
res = res[order(res$padj),]
goi = rownames(res[res$padj < 0.05,])[1:min(10,sum(res$padj < 0.05))]
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)

res = res_stromal_in_vitro
dat = dat_stromal_in_vitro
goi = rownames(res[res$padj < 0.05,])[1:min(10,sum(res$padj < 0.05))]
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
for( g in goi ) {
	p = try( ggplot(dat,aes(endo,log10(get(g)) )) + geom_jitter(alpha=0.7,size=3) + theme_classic() + ggtitle(g) )
	try( print(p) )
	p = try( ggplot(dat,aes(endo,log10(get(g)) )) + geom_boxplot() + theme_classic() + ggtitle(g) )
	try( print(p) )
}
dev.off()

a = data.frame(gene=rownames(res_stromal_in_vitro),logfc=res_stromal_in_vitro$log2FoldChange)
a = a[order(a$logfc,decreasing=T),]
write.table(a,file="endo_DE_stromal_in_vitro.rnk",row.names=FALSE,col.names=FALSE, quote=F,sep="\t")

a = rownames(res_stromal_in_vitro)[res_stromal_in_vitro$padj < 0.05 & res_stromal_in_vitro$log2FoldChange > 0]
write.table(a,file="endo_DE_stromal_in_vitro_sig_up.txt",row.names=FALSE,col.names=FALSE, quote=F,sep="\t")
a = rownames(res_stromal_in_vitro)[res_stromal_in_vitro$padj < 0.05 & res_stromal_in_vitro$log2FoldChange < 0]
write.table(a,file="endo_DE_stromal_in_vitro_sig_down.txt",row.names=FALSE,col.names=FALSE, quote=F,sep="\t")

dec = unique(c(goi$stromal_decidualisation,goi$decidualisation_warren,goi$roser_decidual_stromal,goi$quake_stromal_decidualisation))
a = res_stromal_in_vitro[dec,]
a[ a$padj < 0.1,]

# no correlation whatsoever
a = cbind(res1$log2FoldChange, res2[rownames(res1),"log2FoldChange"])
a = a[rowSums(is.na(a)) == 0,]

library("gplots")
pdf("endo_DE_stromal_in_vitro_heatmap.pdf", width=6, height=6)
dec = unique(c(goi$stromal_decidualisation,goi$decidualisation_warren,goi$roser_decidual_stromal,goi$quake_stromal_decidualisation))
dec = dec[ res[dec,"pvalue"] < 0.05 ]
cc = c("black","red")[as.numeric(factor(colData(de)$endo))]
ind = order(cc)
dat = dat_stromal_in_vitro
dat = dat[dec,ind]
cc = cc[ind]
heatmap.2(as.matrix(dat), trace="none", scale="row",cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10),
		ColSideColors=cc, Colv=F )
dev.off()

# look at decidualised genes
a = unlist(goi)
a = a[grepl("deci",names(a))]
a = unique(a)
b = res_stromal_in_vitro[a,]
b = b[b$pvalue < 0.05,]

# quake genes peaking in stroma at decidualization
library("gplots")
pdf("endo_DE_stromal_in_vitro_heatmap_quake_decidua.pdf", width=6, height=6)
a = read.table("quake_genes_peaking_phase_4_stromal.txt")[,1]
b = res_stromal_in_vitro[a,]
b = b[b$padj < 0.05,]
b = b[!is.na(b[,1]),]
dec = unique(c(goi$stromal_decidualisation,goi$decidualisation_warren,goi$roser_decidual_stromal,goi$quake_stromal_decidualisation))
dec = c(dec[ res[dec,"pvalue"] < 0.05 ],rownames(b))
dec = unique(dec)
cc = c("black","red")[as.numeric(factor(colData(de)$endo))]
ind = order(cc)
dat = dat_stromal_in_vitro
dat = dat[dec,ind]
cc = cc[ind]
heatmap.2(as.matrix(dat), trace="none", scale="row",cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10),
		ColSideColors=cc, Colv=F )
dev.off()

############################
# leukocytes ex vivo
# 5 endo vs 9 ctrl
# 1607 genes significant
############################

library("DESeq2")

dat = get_endoCounts(normalized=FALSE)
tab = dat$counts
meta = dat$meta
rm(dat)
gc()

# remove one cd45p ex vivo endometriosis sample that has bad contamination from stromal cells
# determined by deconvolution "2976.2.0.Cd45p.37"
bad = "2976.2.0.Cd45p.37"
a = props_level_2$Est.prop.weighted
a[bad,]

ind = meta$cell_type == "Cd45p" & meta$type == "ex_vivo" & (rownames(meta) != bad)
table( meta[ind,"endo"] )
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ endo )
de = DESeq(de)
res = as.data.frame(results(de))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
sum(res$padj < 0.05)
dat = counts(de,normalized=TRUE)

res_leukocytes_ex_vivo = res
dat_leukocytes_ex_vivo = dat
save( res_leukocytes_ex_vivo, dat_leukocytes_ex_vivo, file="DE_results_leukocytes_ex_vivo.RData" )

load( "DE_results_leukocytes_ex_vivo.RData" )
sres = res_leukocytes_ex_vivo[ res_leukocytes_ex_vivo$padj < 0.05,]

a = data.frame(gene=rownames(res_leukocytes_ex_vivo),logfc=res_leukocytes_ex_vivo$log2FoldChange)
a = a[order(a$logfc,decreasing=T),]
write.table(a,file="endo_DE_leukocytes_ex_vivo.rnk",row.names=FALSE,col.names=FALSE, quote=F,sep="\t")

a = rownames(res_leukocytes_ex_vivo)[res_leukocytes_ex_vivo$padj < 0.05 & res_leukocytes_ex_vivo$log2FoldChange > 0]
write.table(a,file="endo_DE_leukocytes_ex_vivo_sig_up.txt",row.names=FALSE,col.names=FALSE, quote=F,sep="\t")
a = rownames(res_leukocytes_ex_vivo)[res_leukocytes_ex_vivo$padj < 0.05 & res_leukocytes_ex_vivo$log2FoldChange < 0]
write.table(a,file="endo_DE_leukocytes_ex_vivo_sig_down.txt",row.names=FALSE,col.names=FALSE, quote=F,sep="\t")


pdf("endo_DE_leukocytes_ex_vivo.pdf", width=3, height=3)
res$gene = rownames(res)
res = res[order(res$padj),]
goi = rownames(res[res$padj < 0.05,])[1:min(10,sum(res$padj < 0.05))]
dat = counts(de,normalized=TRUE)
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)

res = res[order(abs(res$log2FoldChange),decreasing=T),]
goi = rownames(res[res$padj < 0.05,])[1:min(10,sum(res$padj < 0.05))]
dat = counts(de,normalized=TRUE)
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
for( g in goi ) {
	p = try( ggplot(dat,aes(endo,log10(get(g)))) + geom_jitter(alpha=0.7,size=3) + theme_classic() + ggtitle(g) )
	try( print(p) )
	p = try( ggplot(dat,aes(endo,log10(get(g)))) + geom_boxplot() + theme_classic() + ggtitle(g) )
	try( print(p) )
}
dev.off()

# plot GSEA results
tab = read.table("/Users/filimon/gsea_home/output/mar04/endo_DE_leukocytes_ex_vivo.GseaPreranked.1709648817451/gsea_report_for_na_pos_1709648817451.tsv", sep="\t", quote="", fill=T, header=T )
tab$direction = "up"
tab2 = read.table("/Users/filimon/gsea_home/output/mar04/endo_DE_leukocytes_ex_vivo.GseaPreranked.1709648817451/gsea_report_for_na_neg_1709648817451.tsv", sep="\t", quote="", fill=T, header=T )
tab2$direction = "down"
tab = rbind(tab,tab2)
rm(tab2)
stab = tab[tab$FDR.q.val < 0.1,]
stab$NAME = sub("HALLMARK_","",stab$NAME)
stab = stab[order(stab$direction,-log10(stab$FDR.q.val+0.00000001)),]
stab$NAME = str_to_sentence(gsub("_"," ",stab$NAME)) 
stab$NAME = factor(stab$NAME, ordered=T, levels=unique(stab$NAME) )

pdf("~/Downloads/endo_DE_leukocytes_ex_vivo_GSEA.pdf", height=3, width=5)
ggplot(stab,aes(direction,NAME)) + geom_point(aes(size=-log10(FDR.q.val+0.00000001),color=abs(as.numeric(NES)))) + theme_classic()
dev.off()

# plot heatmap of genes overlapping with gregersen
load( "DE_results_leukocytes_ex_vivo.RData" )
res = res_leukocytes_ex_vivo
dat = dat_leukocytes_ex_vivo
load("DE_comparison_gregersen_ours_overlapping_immune.RData")

library("gplots")
pdf("DE_comparison_gregersen_ours_overlapping_immune_heatmap.pdf", width=6, height=6)
cc = c("black","red")[as.numeric(factor(colData(de)$endo))]
ind = order(cc)
dat = dat_leukocytes_ex_vivo
dat = dat[a,ind]
cc = cc[ind]
heatmap.2(as.matrix(dat), trace="none", scale="row",cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10),
		ColSideColors=cc, Colv=F )
dev.off()


# leukocytes vitro - makes no sense
# 3 endo vs 2 ctrl
# 17 genes significant
ind = meta$cell_type == "Cd45p" & meta$type == "in_vitro"
table( meta[ind,"endo"] )
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ endo )
de = DESeq(de)
res = as.data.frame(results(de))
res = res[!is.na(res$padj),]
res2 = res[order(abs(res$log2FoldChange),decreasing=T),]
sum(res2$padj < 0.05)



############################
# DE between CD45- in vitro vs ex vivo
############################
library("ggrepel")

ind = !meta$endo & ((meta$cell_type == "Cd45m" & meta$type == "ex_vivo") | (meta$type == "in_vitro"))
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ type )
de = DESeq(de)
res = as.data.frame(results(de))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
sum(res$padj < 0.05)
res[ res$padj < 0.05, ]
res$gene = rownames(res)

goi = rownames(res[res$padj < 0.05,])[1:min(10,sum(res$padj < 0.05))]
dat = counts(de,normalized=TRUE)
dat = data.frame(colData(de),t(dat[goi,]))

res_stromal_in_vitro_vs_ex_vivo = res
dat_stromal_in_vitro_vs_ex_vivo = dat
save( res_stromal_in_vitro_vs_ex_vivo, dat_stromal_in_vitro_vs_ex_vivo, file="DE_results_stromal_in_vitro_vs_ex_vivo.RData" )

load("DE_results_stromal_in_vitro_vs_ex_vivo.RData")
a = data.frame(gene=rownames(res_stromal_in_vitro_vs_ex_vivo),logfc=res_stromal_in_vitro_vs_ex_vivo$log2FoldChange)
a = a[order(a$logfc,decreasing=T),]
write.table(a,file="DE_results_stromal_in_vitro_vs_ex_vivo.rnk",row.names=FALSE,col.names=FALSE, quote=F,sep="\t")


pdf("endo_DE_stromal_in_vitro_vs_ex_vivo.pdf", width=3, height=3)
goi = rownames(res[res$padj < 0.05 & res$log2FoldChange < 0,])[1:min(10,sum(res$padj < 0.05))]
goi = c(goi, rownames(res[res$padj < 0.05 & res$log2FoldChange > 0,])[1:min(10,sum(res$padj < 0.05))])
dat = counts(de,normalized=TRUE)
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
for( g in goi ) {
	p = try( ggplot(dat,aes(type,log10(get(g)), color=cell_type)) + geom_jitter(alpha=0.7,size=3) + theme_classic() + ggtitle(g) )
	try( print(p) )
	
	p = try( ggplot(dat,aes(type,log10(get(g)), fill=cell_type)) + geom_boxplot() + theme_classic() + ggtitle(g) )
	try( print(p) )
	
	p = try( ggplot(dat,aes(type,log10(get(g)))) + geom_boxplot() + theme_classic() + ggtitle(g) )
	try( print(p) )
}
dev.off()

library("gplots")
pdf("endo_DE_stromal_in_vitro_vs_ex_vivo_heatmap.pdf", width=6, height=6)
n = 5
goi = rownames(res[res$padj < 0.05 & res$log2FoldChange < 0,])[1:min(n,sum(res$padj < 0.05))]
goi = c(goi, rownames(res[res$padj < 0.05 & res$log2FoldChange > 0,])[1:min(n,sum(res$padj < 0.05))])
dat = counts(de,normalized=TRUE)
cc = c("black","red")[as.numeric(factor(colData(de)$type))]
ind = order(cc)
dat = dat[,ind]
cc = cc[ind]
heatmap.2(as.matrix(dat[goi,]), trace="none", scale="row",cexRow=1/log10(nrow(dat[goi,]))-0.3, col=colorRampPalette(c("blue","white","red"))(10),
		ColSideColors=cc, Colv=F )
dev.off()


a = res[,c("gene","log2FoldChange")]
a = a[order(a$log2FoldChange,decreasing=T),]
write.table(a,file="endo_DE_stromal_in_vitro_vs_ex_vivo.rnk",quote=F,col.names=F,row.names=F,sep="\t")
a = res[res$padj < 0.05 & res$log2FoldChange > 0,c("gene")]
write.table(a,file="endo_DE_stromal_in_vitro_vs_ex_vivo_genes_up.txt",quote=F,col.names=F,row.names=F,sep="\t")
a = res[res$padj < 0.05 & res$log2FoldChange < 0,c("gene")]
write.table(a,file="endo_DE_stromal_in_vitro_vs_ex_vivo_genes_down.txt",quote=F,col.names=F,row.names=F,sep="\t")

a = read.table("endo_DE_stromal_in_vitro_vs_ex_vivo_genes_up_enrich_gtex_table.txt", fill=T, header=T, sep="\t")
pdf("endo_DE_stromal_in_vitro_vs_ex_vivo_genes_up_enrich_gtex.pdf", width=4, height=4)
p = ggplot(a,aes(Term,-log10(Adjusted.P.value),label=Term)) + geom_point() + theme_classic() 
p = p + geom_text_repel(color="red",data=a[a$Adjusted.P.value < 0.05,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
dev.off()


############################
# severity vs gene expression, including non-endometriosis
# nothing convincing
############################
ind = meta$cell_type == "Cd45p" & meta$type == "ex_vivo"
table( meta[ind,"endo"] )
meta$endo_stage = as.numeric(meta$endo_stage)
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ endo_stage )
de = DESeq(de)
res = as.data.frame(results(de))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
res$gene = rownames(res)
sum(res$padj < 0.05)

pdf("endo_DE_severity_with_controls.pdf", width=3, height=3)
goi = rownames(res[res$padj < 0.05,])[1:min(10,sum(res$padj < 0.05))]
dat = counts(de,normalized=TRUE)
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
for( g in goi ) {
	p = try( ggplot(dat,aes(endo_stage,log10(get(g)))) + geom_jitter(alpha=0.7,size=3) + theme_classic() + ggtitle(g) )
	try( print(p) )
	p = try( ggplot(dat,aes(factor(endo_stage,ordered=T,levels=1:4),log10(get(g)))) + geom_boxplot() + theme_classic() + ggtitle(g) )
	try( print(p) )
}
dev.off()


############################
# severity vs gene expression, only endometriosis
# nothing convincing
############################

ind = meta$cell_type == "Cd45p" & meta$type == "ex_vivo" & meta$endo_stage > 1
table( meta[ind,"endo"] )
meta$endo_stage = as.numeric(meta$endo_stage)
de = DESeqDataSetFromMatrix(countData = tab[,ind], colData = meta[ind,], design= ~ endo_stage )
de = DESeq(de)
res = as.data.frame(results(de))
res = res[!is.na(res$padj),]
res = res[order(abs(res$log2FoldChange),decreasing=T),]
res$gene = rownames(res)
sum(res$padj < 0.05)

pdf("endo_DE_severity_without_controls.pdf", width=3, height=3)
goi = rownames(res[res$padj < 0.05,])[1:min(10,sum(res$padj < 0.05))]
dat = counts(de,normalized=TRUE)
dat = data.frame(colData(de),t(dat[goi,]))
p = ggplot(res,aes(log2FoldChange,-log10(padj),label=gene)) + geom_point(alpha=0.7, shape=19) + geom_hline(yintercept=-log10(0.05)) 
p = p + theme_classic() + geom_text_repel(color="red",data=res[res$gene %in% goi,],min.segment.length = 0, ylim=c(1.30103, Inf) )
print(p)
for( g in goi ) {
	p = try( ggplot(dat,aes(factor(endo_stage,ordered=T,levels=1:4),log10(get(g)))) + geom_jitter(alpha=0.7,size=3) + theme_classic() + ggtitle(g) )
	try( print(p) )
	p = try( ggplot(dat,aes(factor(endo_stage,ordered=T,levels=1:4),log10(get(g)))) + geom_boxplot() + theme_classic() + ggtitle(g) )
	try( print(p) )
}
dev.off()


###############################################################################
# DE compare the 3 models
###############################################################################

meta = get_metadata_epfl()
load( "DE_results_all.RData" )
load( "DE_results_leukocytes_ex_vivo.RData" )
load( "DE_results_stromal_in_vitro.RData" )

s1 = rownames(res_all)[res_all$padj < 0.05]
s2 = rownames(res_leukocytes_ex_vivo)[res_leukocytes_ex_vivo$padj < 0.05]
s3 = rownames(res_stromal_in_vitro)[res_stromal_in_vitro$padj < 0.05]

sort(table(c(s1,s2,s3)),decreasing=T)[1:10]

gene = "ORM2"
res_all[gene,]
res_leukocytes_ex_vivo[gene,]
res_stromal_in_vitro[gene,]

res_all[c("TNF","MMP9","TIMP1"),]
res_leukocytes_ex_vivo[c("TNF","MMP9","TIMP1"),]
res_stromal_in_vitro[c("TNF","MMP9","TIMP1"),]

ora = c()
for( ct in c("stromal_in_vitro", "leukocytes_ex_vivo") ) {
	for( direction in c("up","down") ) {
		for( gs in c("reactome", "msigdb") ) {
			a = read.table(paste0("endo_DE_", ct, "_sig_", direction, "_enrichr_", gs, ".txt"), sep="\t", header=T)
			a$ct = ct
			a$dir = direction
			a$gs = gs
			ora = rbind(ora,a)
		}
	}
}
ora = ora[order(ora$Adjusted.P.value),]

###############################################################################
# DE heatmap with pathways stromal in vitro
###############################################################################

library("gplots")
meta = get_metadata_epfl()

maxpath = 12-1
maxgenes = 156
sora = ora[ ora$ct == "stromal_in_vitro" & ora$Adjusted.P.value < 0.05 & ora$gs == "msigdb",]
sora = rbind( sora[sora$dir == "up",][1:3,], sora[sora$dir == "down",][1:(maxpath-3),] )
oratab = unlist(apply( sora, 1, function(x) { paste(strsplit(x["Genes"],";")[[1]],x["Term"],sep="_") } ))
oratab = data.frame( gene=sapply(strsplit(oratab,"_"),"[[",1),  path=sapply(strsplit(oratab,"_"),"[[",2) )       

res = res_stromal_in_vitro

dec = read.table("quake_genes_peaking_phase_4_stromal.txt")[,1]
dec = unique(c(dec,goi$stromal_decidualisation,goi$decidualisation_warren,goi$roser_decidual_stromal,goi$quake_stromal_decidualisation))
dec = dec[dec %in% res$gene[res$padj < 0.05]]
oratab = rbind(oratab,data.frame(gene=dec,path="decidualization"))

gg = sample(oratab$gene, min(nrow(oratab),150))
gg = c(gg,rownames(res)[order(res$padj)][1:60])
gg = c(gg,rownames(res)[order(abs(res$log2FoldChange),decreasing=T)][1:60])
gg = unique(gg)
length(gg)

res = res[ rownames(res) %in% gg,]
dat = dat_stromal_in_vitro[rownames(res),]
smeta = meta[match( colnames(dat), meta$allI ),]
rownames(smeta) = colnames(dat)
smeta$cc = c("black","red")[as.numeric(smeta$endo)+1]
ind = order(smeta$cc)
dat = dat[,ind]
smeta = smeta[ind,]

gtab = data.frame( gene=rownames(dat), path=oratab$path[match(rownames(dat),oratab$gene)] )
gtab$col = brewer.pal(maxpath, "Paired")[as.numeric(factor(gtab$path))]
gtab$col[is.na(gtab$col)] = "black"

na.rm <- TRUE
scaleRow <- function(x) {
	rm <- rowMeans(x, na.rm = na.rm)
	x <- sweep(x, 1, rm)
	sx <- apply(x, 1, sd, na.rm = na.rm)
	x <- sweep(x, 1, sx, "/")
	return(round(x, 6))
}
a = scaleRow(dat)
ind = order(gtab$col,-1*apply(a[,smeta$endo],1,median),decreasing=F)
dat = dat[ind,]
gtab = gtab[ind,]

pdf("endo_DE_stromal_in_vitro_heatmap.pdf", width=3, height=2*5)
heatmap.2(as.matrix(dat), trace="none", scale="row", cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10), RowSideColors=gtab$col, ColSideColors=smeta$cc, Colv=F, Rowv=F )
heatmap.2(as.matrix(dat), trace="none", scale="row", cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10), RowSideColors=gtab$col, ColSideColors=smeta$cc, Colv=T, Rowv=T )
dev.off()


###############################################################################
# DE heatmap with pathways leukocytes ex vivo
###############################################################################

library("gplots")
library("RColorBrewer")
meta = get_metadata_epfl()
load( "DE_results_leukocytes_ex_vivo.RData" )

maxpath = 12
sora = ora[ ora$ct == "leukocytes_ex_vivo" & ora$Adjusted.P.value < 0.05 & ora$gs == "msigdb",]
sora = rbind( sora[sora$dir == "up",][1:3,], sora[sora$dir == "down",][1:(maxpath-3),] )
oratab = unlist(apply( sora, 1, function(x) { paste(strsplit(x["Genes"],";")[[1]],x["Term"],sep="_") } ))
oratab = data.frame( gene=sapply(strsplit(oratab,"_"),"[[",1),  path=sapply(strsplit(oratab,"_"),"[[",2) )       

res = res_leukocytes_ex_vivo

gg = rownames(res)[order(res$padj)][1:10]
gg = c(gg,rownames(res)[order(abs(res$log2FoldChange),decreasing=T)][1:10])
gg = unique(c(gg, sample(oratab$gene, 170 ) ))
length(gg)

res = res[ rownames(res) %in% gg,]
#res = res[order(res$padj),]
#res = res[res$padj < 0.05,][1:256,]
dat = dat_leukocytes_ex_vivo[rownames(res),]
smeta = meta[match( colnames(dat), meta$allI ),]
rownames(smeta) = colnames(dat)
smeta["4210.4.0.Cd45p.24","endo"] = TRUE
smeta$cc = c("black","red")[as.numeric(smeta$endo)+1]
ind = order(smeta$cc)
dat = dat[,ind]
smeta = smeta[ind,]

gtab = data.frame( gene=rownames(dat), path=oratab$path[match(rownames(dat),oratab$gene)] )
gtab$col = brewer.pal(maxpath, "Paired")[as.numeric(factor(gtab$path))]
gtab$col[is.na(gtab$col)] = "black"

na.rm <- TRUE
scaleRow <- function(x) {
	rm <- rowMeans(x, na.rm = na.rm)
	x <- sweep(x, 1, rm)
	sx <- apply(x, 1, sd, na.rm = na.rm)
	x <- sweep(x, 1, sx, "/")
	return(round(x, 6))
}
a = scaleRow(dat)
ind = order(gtab$col,-1*apply(a[,smeta$endo],1,median),decreasing=F)
dat = dat[ind,]
gtab = gtab[ind,]

pdf("endo_DE_leukocytes_ex_vivo_heatmap.pdf", width=3, height=2*5)
heatmap.2(as.matrix(dat), trace="none", scale="row", cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10), RowSideColors=gtab$col, ColSideColors=smeta$cc, Colv=F, Rowv=F )
heatmap.2(as.matrix(dat), trace="none", scale="row", cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10), RowSideColors=gtab$col, ColSideColors=smeta$cc, Colv=T, Rowv=T )
dev.off()


###############################################################################
# DE heatmap all samples
###############################################################################

library("gplots")
library("RColorBrewer")

meta = get_metadata_epfl()
a = load( "DE_results_all.RData" )

#maxpath = 12
#sora = ora[ ora$ct == "leukocytes_ex_vivo" & ora$Adjusted.P.value < 0.05 & ora$gs == "msigdb",]
#sora = rbind( sora[sora$dir == "up",][1:3,], sora[sora$dir == "down",][1:(maxpath-3),] )
#oratab = unlist(apply( sora, 1, function(x) { paste(strsplit(x["Genes"],";")[[1]],x["Term"],sep="_") } ))
#oratab = data.frame( gene=sapply(strsplit(oratab,"_"),"[[",1),  path=sapply(strsplit(oratab,"_"),"[[",2) )       

res = res_all

#gg = rownames(res)[order(res$padj)][1:10]
#gg = c(gg,rownames(res)[order(abs(res$log2FoldChange),decreasing=T)][1:10])
#gg = unique(c(gg, sample(oratab$gene, 170 ) ))
#length(gg)

#res = res[ rownames(res) %in% gg,]
#res = res[order(res$padj),]
res = res[res$padj < 0.05,][1:256,]
dat = dat_all[rownames(res),]
smeta = meta[match( colnames(dat), meta$allI ),]
rownames(smeta) = colnames(dat)
smeta$cc = c("black","red")[as.numeric(smeta$endo)+1]
ind = order(smeta$cc)
dat = dat[,ind]
smeta = smeta[ind,]

#gtab = data.frame( gene=rownames(dat), path=oratab$path[match(rownames(dat),oratab$gene)] )
#gtab$col = brewer.pal(maxpath, "Paired")[as.numeric(factor(gtab$path))]
#gtab$col[is.na(gtab$col)] = "black"

na.rm <- TRUE
scaleRow <- function(x) {
	rm <- rowMeans(x, na.rm = na.rm)
	x <- sweep(x, 1, rm)
	sx <- apply(x, 1, sd, na.rm = na.rm)
	x <- sweep(x, 1, sx, "/")
	return(round(x, 6))
}
a = scaleRow(dat)
#ind = order(gtab$col,-1*apply(a[,smeta$endo],1,median),decreasing=F)
ind = order(-1*apply(a[,smeta$endo],1,median),decreasing=F)
dat = dat[ind,]
#gtab = gtab[ind,]

pdf("endo_DE_all_heatmap.pdf", width=3, height=2*5)
heatmap.2(as.matrix(dat), trace="none", scale="row", cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10), ColSideColors=smeta$cc, Colv=F, Rowv=F )
heatmap.2(as.matrix(dat), trace="none", scale="row", cexRow=1/log10(nrow(dat))-0.3, col=colorRampPalette(c("blue","white","red"))(10), ColSideColors=smeta$cc, Colv=T, Rowv=T )
dev.off()

###############################################################################
# plot all together
###############################################################################

dat = get_endoCounts(normalized=TRUE)
tab = dat$counts
meta = dat$meta
rm(dat)
gc()

load("DE_results_all.RData")

a = res_all[res_all$padj < 0.05,]
g = a[order(a$log2FoldChange,decreasing=T),][1:min(10,nrow(a)),]
g = rbind(g,a[order(a$log2FoldChange,decreasing=F),][1:min(10,nrow(a)),])

cc = rainbow(8, start=0, end=.3)[as.numeric(as.factor(paste(meta$cell_type,meta$type)))]
cc = c("black","red")[as.numeric(as.factor(meta$endo))]

library("gplots")
pdf("DE_heatmap.pdf")
heatmap.2(as.matrix(tab[rownames(g),]),trace="none",scale="row",ColSideColors=cc)
dev.off()
