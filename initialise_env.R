# change to your own working directory
dir = "/omics/groups/OE0433/internal/angela/mf/"
# set to your preferred screen width
options(width=230)

setwd(dir)

###############################################################################
# Load libraries 
###############################################################################

.libPaths(paste("/omics/groups/OE0433/internal/software/R-library-",paste(R.version$major,R.version$minor,sep="."),sep=""))

library("Seurat")
#library("pagoda2")
library("RColorBrewer")
library("ggplot2")
#library("velocyto.R")
library("dplyr")
#library("garnett")
library("org.Hs.eg.db")
library("rhdf5")
#library("slingshot")
#library("gam")

#BiocManager::install("")
#library("devtools")
#devtools::install_github("cole-trapnell-lab/garnett")

source("/omics/groups/OE0433/internal/angela/atfg_github/MF/functions.R")

###############################################################################
# Set up some useful global variables 
###############################################################################

datamap = data.frame( batch=c(rep("jun17",3),"all","org","10xDemo1","all","aug17","aug17","aug17","aug17"), #
		sample_full=c("MFCON007efM","MFCON010dfM","MFCON018bfM","MFmg","orga","pbmc4k","MFpg","MFCON020afM","MFCON007dfM","MFCON007dcM","MFCON020acM"), 
		run=c(rep(23156,3),"",21346,"10xPBMC","aggregate","24192-25608","24192-25608","24192-25608","24192-25608"), 
		lane=c(6,7,8,"6-7-8",NA,"4k",NA,"4839STDY7131581","4839STDY7131582","4839STDY7131583","4839STDY7131584"),
		root="/omics/groups/OE0433/internal/data/mf/",
		filter_ngenes = c(5200), filter_percent.mito = c(0.3), nPCs=c(5), stringsAsFactors=F, 
		nPC = c(11,11,5,NA,NA,NA,NA,11,11,20,11),
		maxPC = c(50,50,20,NA,NA,NA,NA,50,50,50,50),
		k=c(5,4,5,3,3,NA,NA,NA,NA,NA,NA) )
datamap$sample = substring(datamap$sample_full,nchar(datamap$sample_full)-4,nchar(datamap$sample_full)-1)
datamap$base_dir=paste(datamap$root,"cellranger201_count_", datamap$run, "_", datamap$lane, "_GRCh38/",sep="")
datamap$matrix_dir=paste(datamap$base_dir,"/outs/filtered_gene_bc_matrices/GRCh38/",sep="")
datamap$irods_dir=paste("/seq/",datamap$run,"/cellranger/", sep="")


goi = list( cilia = c("TUBA1A", "CTH", "EZR","FOXJ1"),
		ciliated_epithelium = c("SLC1748","DRG1","CFAP52","CFAP46","TUBA4B","FAM52B","DNAH3","AKAP14",
				"ARMC4","ENKUR","RSPH4A","TEKT1","SNTN","DNAH10","ZMYND10","DYNLR52","CFAP45","PIFO",
				"CFAP43","HYDIN","DNAH7","MAATS1","DNAH8","DNAAF1","SPAG17","CFAP53","DNAH5","TPPP3",
				"DNAH11","RSPH1","FOXJ1","WDR66"),# from Wang/Quake
		ciliated_epithelium_2 = c("DYDC2","CDHR3","SNTN","DYNLRB2","FAM183A","PIFO","ARMC3","CAPSL",
				"ROPN1L","CCDC173","ZBBX","DNAH12","DNAH7","RSPH1","SPATA18","SPAG17","MNS1","TMEM231",
				"LRRIQ1","IQCG","AGR3","TPPP3","FHAD1"),
		secretory = c("PAEP", "MTHFD1","ESR1"),
		up_glands_vs_stroma=c("CDH1","CLDN10","EPCAM","PAX8","MUC1","PAEP", "KLC11","MUC20","FOXA2","SOX17","KLF5"),
		mucosal_secretory=c("PAX8","MUC1"),
		glandular_products=c("PAEP","KLC11","MUC20"),
		up_stroma_vs_gland=c("THY1","NT5E","IFITM1","COL8A1","COL12A1","COL13A1",
				"LAMA1","MMP11","MMP2","MMP12","MMP27","MMP3","TIMP2","CTGF"),
		endothelial=c("PECAM1","CD34","VWF"), # PECAM1 and CD34 confirmed by Roser
		smooth_muscle = c("MGP"),
		epithelial = c("MUC16"), # AKA CA125, in epithelia of the fallopian tube, endometrium, endocervix and ovaries
		endometrial_epithelium = c("ALCAM"),
		vaginal_epithelium = c("KRT4","KRT14"), # KRT4 - Predominantly non-keratinised squamous epithelia, ectocervix and vagina, KR14 is expressed in vaginal and cervical but not in uterine epithelium
		stromal=c("ANPEP","VIM"),
		krjutskov_endometrial_epithelial = "CD9",
		krjutskov_endometrial_stromal = "ANPEP",
		krjutskov_basal_glandular_epithelial_cells = "FUT4", #SSEA-1 or FUT4
		petra_stem=c("MYC","CCND2","VGLL4","LCP1","LGR5","VDR","WWP1"),
		epithelial_stem=c("PROM1","AXIN2","LRIG1","SOX9"),
		endometrial_stem=c("FUT4","FUT9","SOX9"),
		epithelial_glandular_vs_luminal_marker_glandular = c("ITGA1"),
		epithelial_glandular_vs_luminal_marker_luminal = c("WNT7A","SVIL"),
		epithelial_glandular_vs_luminal_top_de = c("HPGD","SULT1E1","LGR5","VTCN1","ITGA1"),
		cemsel_stem=c("PDGFRB"),
		gargett_mesenchymal_stem=c("ITGB1","CD44","NT5E","THY1","ENG","PDGFRB","MCAM"),
		gargett_menstrual_blood_ercs=c("TERT","POU5F1","SSEA4","NANOG"), # MSCs or stromal fibroblasts, STRO-1 has no gene symbol because no seq is available
		proliferation=c("MKI67"), # not expressed
		differentiation=c("PAEP","SPP1", "17HSDB2", "LIF"),
		stromal_decidualisation = c("IGFBP1","PRL","PRLR","MMP3","MMP9","LEFTY2","FOXO1","IL15","HAND2","HAND2-AS1"), # PRLR is mainly expressed in epithelial cells (IRF1)
		quake_stromal_decidualisation = c("CRYAB","S100A4","DKK1","FOXO1","IL15","FGF7","LMCD1"),
		epithelial_decidualisation = c("PRLR","IRF1","FOXO1"),
		decidualisation_warren = c("IGFBP1", "SST", "PRL", "BCL2L11", "WNT5A", "FOXO1"),
		leukocytes=c("PTPRC","HLA-B"),
		nk=c("NCAM1"),
		nk_mucosal=c("ITGAE"),
		tcells=c("ITGAE","CD69","CXCR4","CD8A", "CD8B", "CD4"),
		apoptosis=c("CASP3"),
		ovary_specific=c("OVGP1","MUM1L1","RBP1","CDH11","NCAM1"),
		endometrium_specific=c("MME"),
		up_fallopian_vs_endometrium=c("FMO3"),
		up_endometrium_vs_fallopian=c("DMBT1"),
		cervix_fallopian_specfic=c("MUC4"),
		blending=c("PTPRC_HLA-B"),
		stromal_odd=c("NEAT1","MALAT1"),
		msc=c("NT5E","THY1","ENG"), # and do not express CD45, CD34, ITGAX, CD14, CD19, CD79A and HLA-DR
		haem=c("PTPRC","CD34","ITGAX","CD14","CD19","CD79A","HLA-DRA","HLA-DRB1"),
		cycling=c("MCM2","CENPF","MK167"),  # Roser
		roser_mf1=c("DKK1","MCAM", "PDGFRB", "REN", "RGS5"), # MCAM (CD146) and PDGFRB are pericyte markers
		roser_mf2=c("DKK1"),
		roser_mf3=c("DKK1","ACTA2"), #spongiosa
		roser_mf4=c(), # decidua compacta
		roser_mf5=c(), # decidua compacta
		roser_nk1 = c("ITGA1","CD9","CD39"), # endometrial resident NK cells
		roser_nk2 = c("ITGA1","CD9","CDHR1"),
		roser_nk3 = c("ITGA1","CD9","CXCR4"),
		roser_pericytes = c("MCAM","PDGFRB","REN","RGS5"),
		roser_glands_active = c("FXYD2","MT1G","S100A1","MT1F"),
		roser_glands_exhausted = c("TMEM49","XIST","CP","VMP1","NEAT1"),
		roser_decidual_stromal = c("DKK1"),
		roser_spongiosa_vs_compacta = c("ACTA2"),
		hla=c("HLA-A")
) 
dissociation_effects=c('ACTG1','ANKRD1','ARID5A','ATF3','ATF4','BAG3','BHLHE40','BRD2',
		'BTG1','BTG2','CCNL1','CCRN4L','CEBPB','CEBPD','CEBPG','CSRNP1','CXCL1','CYR61',
		'DCN','DDX3X','DDX5','DES','DNAJA1','DNAJB1','DNAJB4','DUSP1','DUSP8','EGR1','EGR2',
		'EIF1','EIF5','ERF','ERRFI1','FAM132B','FOS','FOSB','FOSL2','GADD45A','GCC1','GEM',
		'H3F3B','HIPK3','HSP90AA1','HSP90AB1','HSPA1A','HSPA1B','HSPA5','HSPA8','HSPB1',
		'HSPH1','ID3','IDI1','IER2','IER3 ','IFRD1','IL6','IRF1','IRF8','ITPKC','JUN','JUNB',
		'JUND','KLF2','KLF4','KLF6','KLF9','LITAF','LMNA','MAFF','MAFK','MCL1','MIDN',
		'MIR22HG','MT1','MT2','MYADM','MYC','MYD88','NCKAP5L','NCOA7','NFKBIA','NFKBIZ',
		'NOP58','NPPC','NR4A1','ODC1','OSGIN1','OXNAD1','PCF11','PDE4B','PER1','PHLDA1',
		'PNP','PNRC1','PPP1CC','PPP1R15A','PXDC1','RAP1B','RASSF1','RHOB','RHOH','RIPK1',
		'SAT1','SBNO2','SDC4','SERPINE1','SKIL','SLC10A6','SLC38A2','SLC41A1','SOCS3','SQSTM1',
		'SRF','SRSF5','SRSF7','STAT3','TAGLN2','TIPARP','TNFAIP3','TNFAIP6','TPM3','TPPP3',
		'TRA2A','TRA2B','TRIB1','TUBB4B','TUBB6','UBC','USP2','WAC','ZC3H12A','ZFAND5','ZFP36',
		'ZFP36L1','ZFP36L2','ZYX','GADD45G','HSPE1','IER5','KCNE4')


# module load R/3.5.2, don't have to do this, it is the default, doesn't work from inside screen
# module load python/3.6.1 # this is required for Seurat to run
# module load gsl/2.5 

# BiocManager::install("pcaMethods")
# library(devtools)
# install_github("velocyto-team/velocyto.R")

cols1 = c("#961D4E","#4C6085","#FCB97D","#F1DA87","#BEC06B")
names(cols1) = c("epithelial","stromal","leukocytes","endothelial","smooth muscle")
cols1_scale = scale_color_manual(breaks=names(cols1), values=cols1)
cols1_fill = scale_fill_manual(breaks=names(cols1), values=cols1)

cols2 = c("#961D4E","#CE84AD","#661431","#FCB97D","#93B5C6","#4C6085","#2D394E")
names(cols2) = c("glands","exhausted glands","cilliated glands","T-cells","decidualised stroma","stroma/mesenchyme","smooth muscle")
cols2_scale = scale_color_manual(breaks=names(cols2), values=cols2)
cols2_fill = scale_fill_manual(breaks=names(cols2), values=cols2)

# ordered list
types_level2 = c("exhausted glands","glands","cilliated glands","decidualised stroma","stroma/mesenchyme","smooth muscle","T-cells")

