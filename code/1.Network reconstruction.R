## Complete Pipeline for Network Construction from Transcriptome Dataset ##

############### Step 0: Preparation ###############

# Load NetBID2 package
library(NetBID2)
library(tidyverse)

# Define main working directory and project name
project_main_dir <- './c9orf72' # user defined main directory for the project, one main directory could have multiple project folders, distinguished by project name.
current_date <- format(Sys.time(), "%Y-%m-%d") # optional, if user like to add current date to name the project folder.
project_name <- sprintf('project_%s',current_date) # project name for the project folders under main directory.

# Create a hierarchcial working directory and return a list contains the hierarchcial working directory information
# This list object (network.par) is an ESSENTIAL variable in network construction pipeline
network.par  <- NetBID.network.dir.create(project_main_dir=project_main_dir,project_name=project_name)


############### Step 1: Load in gene expression datasets for network construction (exp-load) ###############
# Creat eset 
# 对counts raw进行Deseq2标准化之后再进行log2转换和基因过滤获得最终所需表达矩阵
counts_mat <- read_tsv('data/counts.txt') %>% column_to_rownames(.,'gene_name')

# counts raw数据需要先进行RNASeqCount.normalize.scale()标准化处理
counts_mat <- RNASeqCount.normalize.scale(counts_mat)

# 然后判断是否需要log2转换函数
med_val <- median(apply(counts_mat,2,median)); print(med_val)
if(med_val>16){counts_mat <- log2(counts_mat)}

# 过滤低表达基因
sum(is.na(counts_mat))
choose1 <- apply(counts_mat<= quantile(as.matrix(counts_mat), probs = 0.05), 1, sum)<= ncol(counts_mat) * 0.90 #这里要把dataframe转换成矩阵不然会报错
print(table(choose1))
counts_mat <- counts_mat[choose1,]

phenotype_info <- read_tsv('data/phenotype_info.txt')%>% column_to_rownames(.,'geo_accession')

sum(is.na(counts_mat))

eset <- generate.eset(exp_mat = counts_mat, phenotype_info = phenotype_info,feature_info = NULL, annotation_info = "")

identical(colnames(eset@assayData$exprs),rownames(eset@phenoData))

# Add the variable into network.par. ESSENTIAL STEP.
network.par$net.eset <- eset

# Save Step 1 network.par as RData
NetBID.saveRData(network.par = network.par,step='exp-load')

# QC for eset
draw.eset.QC(network.par$net.eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='after_QC_',
             emb_plot_type='2D.interactive') #Error: no more error handlers available (recursive errors?)


############### Step 2: Check sample cluster analysis, optional (exp-cluster) ###############

# Reload network.par RData from Step 2
network.par <- list()
network.par$out.dir.DATA <- 'c9orf72/project_2022-10-15/DATA/' 
NetBID.loadRData(network.par = network.par,step='exp-load')

# Select the most variable genes across samples
mat <- exprs(network.par$net.eset)
choose1 <- IQR.filter(exp_mat=mat,use_genes=rownames(mat),thre = 0.5)
print(table(choose1))
mat <- mat[choose1,]


# Generate temporary eset
tmp_net_eset <- generate.eset(exp_mat=mat, phenotype_info=pData(network.par$net.eset)[colnames(mat),],
                              feature_info=fData(network.par$net.eset)[rownames(mat),], annotation_info=annotation(network.par$net.eset))
# QC plot for IQR filtered eset
draw.eset.QC(tmp_net_eset,outdir=network.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='Cluster_',
             emb_plot_type='2D.interactive')

# The following scripts provide various ways to visualize and check if the IQR filter selected genes
# can be used to perform good sample cluster analysis (predicted labels vs. real labels). Figures will be displayed instead of saving as files.

# Extract phenotype information data frame from eset
phe <- pData(network.par$net.eset)
# Extract all "cluster-meaningful" phenotype columns
intgroup <- get_int_group(network.par$net.eset)
# Cluster analysis using Kmean and plot result using PCA biplot (pca+kmeans in 2D)
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.emb.kmeans(mat=mat,embedding_method = "pca", all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]))
}
# Cluster analysis using Kmeans and plot result using PCA (pca+kmeans in 3D)
for(i in 1:length(intgroup)){
  print(intgroup[i])
  pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,intgroup[i]),plot_type='3D')
  print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe,intgroup[i]))))
}
# Pick one phenotype column "subgroup" from the demo eset and show various plots NetBID2 can create
use_int <- c('tissue','group')
pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D')
pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.ellipse')
pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='2D.text')
pred_label <- draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),plot_type='3D')
print(table(list(pred_label=pred_label,obs_label=get_obs_label(phe, use_int))))
draw.clustComp(pred_label,obs_label=get_obs_label(phe,use_int),outlier_cex=1,low_K=10)

# 交互式2D作图
draw.emb.kmeans(mat=mat,all_k = NULL,obs_label=get_obs_label(phe,use_int),
                plot_type='2D.interactive')

############### Step 3: Prepare files to run SJARACNe (sjaracne-prep) ###############

# Reload network.par RData from Step 1
network.par <- list()
network.par$out.dir.DATA <- 'c9orf72/project_2022-10-15/DATA/' 
NetBID.loadRData(network.par = network.par,step='exp-load')
# Load database
db.preload(use_level='gene',use_spe='human',update=FALSE)

# Converts gene ID into the corresponding TF/SIG list
use_gene_type <- 'external_gene_name' # user-defined
use_genes <- rownames(fData(network.par$net.eset))
use_list  <- get.TF_SIG.list(use_genes,use_gene_type=use_gene_type)

# Select Frontal Cortex samples for analysis    
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe[phe$tissue=='Frontal Cortex',]) # here is using all samples, users can modify
prj.name <- 'Frontal_Cortex' # if use different samples, need to change the project name
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0.5,IQR.loose_thre = 0.1,
                 SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)

# Select Cerebellum samples for analysis    
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe[phe$tissue=='Cerebellum',]) # here is using all samples, users can modify
prj.name <- 'Cerebellum' # if use different samples, need to change the project name
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0.5,IQR.loose_thre = 0.1,
                 SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)

# Select all samples for analysis    
phe <- pData(network.par$net.eset)
use.samples <- rownames(phe) # here is using all samples, users can modify
prj.name <- 'Cortex_Cerebellum' # if use different samples, need to change the project name
SJAracne.prepare(eset=network.par$net.eset,use.samples=use.samples,
                 TF_list=use_list$tf,SIG_list=use_list$sig,
                 IQR.thre = 0.5,IQR.loose_thre = 0.1,
                 SJAR.project_name=prj.name,SJAR.main_dir=network.par$out.dir.SJAR)
