## Complete Pipeline for Driver Estimation and Master Table Creation ##

############### Step 0: Preparation ###############

# Load NetBID2 package
library(NetBID2)

# Get the constructed network data
network.dir <- 'c9orf72/project_2022-10-15' 
network.project.name <- 'project_2022-10-15' 

TF_Cerebellum <- 'c9orf72/project_2022-10-15/SJAR/Cerebellum/output_tf/consensus_network_ncol_.txt'
Sig_Cerebellum <- 'c9orf72/project_2022-10-15/SJAR/Cerebellum/output_sig/consensus_network_ncol_.txt'


# Define main working directory and project name
project_main_dir <- 'c9orf72/' # user defined main directory for the project, one main directory could

project_Cerebellum <- sprintf('driver_%s','Cerebellum') # project name for the project folders under main directory

# Create a hierarchcial working directory and return a list contains the hierarchcial working directory information
# This list object (analysis.par) is an ESSENTIAL variable in driver estimation pipeline


analysis.par  <- NetBID.analysis.dir.create(project_main_dir=project_main_dir, project_name=project_Cerebellum,
                                            network_dir=network.dir, tf.network.file=TF_Cerebellum, sig.network.file=Sig_Cerebellum)


############### Step 1: Load in gene expression dataset for analysis (exp-load, exp-cluster, exp-QC) ###############

# If use the same expression dataset as in the network construction, just reload it directly
load(sprintf('%s/DATA/network.par.Step.exp-load.RData',network.dir)) # RData saved in the network construction step
analysis.par$cal.eset <- network.par$net.eset

# Save Step 1 network.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='exp-load')

############### Step 2: Read in network files and calcualte driver activity (act-get) ###############

# Reload network.par RData from Step 1
analysis.par <- list()

analysis.par$out.dir.DATA <- 'c9orf72/driver_Cerebellum/DATA'

NetBID.loadRData(analysis.par=analysis.par,step='exp-load')

# Get network information
analysis.par$tf.network  <- get.SJAracne.network(network_file=analysis.par$tf.network.file)
analysis.par$sig.network <- get.SJAracne.network(network_file=analysis.par$sig.network.file)

# Creat QC report for the network
draw.network.QC(analysis.par$tf.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='TF_net_',html_info_limit=FALSE)
draw.network.QC(analysis.par$sig.network$igraph_obj,outdir=analysis.par$out.dir.QC,prefix='SIG_net_',html_info_limit=TRUE)


# Merge network first
analysis.par$merge.network <- merge_TF_SIG.network(TF_network=analysis.par$tf.network,SIG_network=analysis.par$sig.network)

# Get activity matrix
ac_mat <- cal.Activity(target_list=analysis.par$merge.network$target_list,cal_mat=exprs(analysis.par$cal.eset),es.method='weightedmean')

# Create eset using activity matrix
analysis.par$merge.ac.eset <- generate.eset(exp_mat=ac_mat,phenotype_info=pData(analysis.par$cal.eset)[colnames(ac_mat),],
                                            feature_info=NULL,annotation_info='activity in net-dataset')

# QC plot for activity eset
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=NULL,do.logtransform=FALSE,prefix='AC_',emb_plot_type ='2D.interactive')
# intgroup应该是默认前6个有意义的分组，超过可能需要手动设置要展示的组别

intgroup <- get_int_group(analysis.par$merge.ac.eset)
draw.eset.QC(analysis.par$merge.ac.eset,outdir=analysis.par$out.dir.QC,intgroup=intgroup,do.logtransform=FALSE,prefix='AC2_',emb_plot_type ='2D.interactive')
# Too many samples for drawing the correlation plot, will select the first 30 samples !
##             If want to specify the sample list, please choose the subset of eset as input !

# Save Step 2 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-get')

############### Step 3: Get differential expression (DE) / differential activity (DA) for drivers (act-DA) ###############

# Reload network.par RData from Step 2
analysis.par <- list()
analysis.par$out.dir.DATA <- 'c9orf72/driver_Cerebellum/DATA'
NetBID.loadRData(analysis.par=analysis.par,step='act-get')

# Create empty list to store comparison result
analysis.par$DE <- list()
analysis.par$DA <- list()

# First comparison: Cluster2 vs. Cluster4
comp_name <- 'ALS Vs Healthy' # Each comparison must has a name
# Get sample names from each compared group
phe_info <- pData(analysis.par$cal.eset)
G1  <- rownames(phe_info)[which((phe_info$group=='ALS') & (phe_info$tissue=='Cerebellum'))] # Experiment group
G0  <- rownames(phe_info)[which((phe_info$group=='Healthy') & (phe_info$tissue=='Cerebellum'))] # Control group
DE_gene_bid <- getDE.BID.2G(eset=analysis.par$cal.eset,G1=G1,G0=G0,G1_name='ALS',G0_name='Healthy')
DA_driver_bid   <- getDE.BID.2G(eset=analysis.par$merge.ac.eset,G1=G1,G0=G0,G1_name='ALS',G0_name='Healthy')
# Save comparison result to list element in analysis.par, with comparison name
analysis.par$DE[[comp_name]] <- DE_gene_bid
analysis.par$DA[[comp_name]] <- DA_driver_bid

# Save Step 3 analysis.par as RData
NetBID.saveRData(analysis.par=analysis.par,step='act-DA')

# Visualize top drivers 以main_id为排序标准
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='ALS Vs Healthy', text_cex = 1, col_srt = 60)
draw.NetBID(DA_list=analysis.par$DA,DE_list=analysis.par$DE,main_id='ALS Vs Healthy', pdf_file=sprintf('%s/NetBID_Cluster_TOP.pdf',analysis.par$out.dir.PLOT),row_cex = 1, column_cex = 0.6,text_cex=0.8, col_srt = 45) # Save as PDF

############### Step 4: Generate a master table for drivers (ms-tab) ###############

# Reload analysis.par RData from Step 3
analysis.par <- list()
analysis.par$out.dir.DATA <- 'c9orf72/driver_Cerebellum/DATA'
NetBID.loadRData(analysis.par=analysis.par,step='act-DA')

# Reload data into R workspace, and saves it locally under db/ directory with specified species name and analysis level.
db.preload(use_level='gene',use_spe='human',update=FALSE)
# Get all comparison names
all_comp <- names(analysis.par$DE) # Users can use index or name to get target ones
# Prepare the conversion table (OPTIONAL)
use_genes <- unique(c(analysis.par$merge.network$network_dat$source.symbol,analysis.par$merge.network$network_dat$target.symbol))
transfer_tab <- get_IDtransfer2symbol2type(from_type = 'external_gene_name',use_genes=use_genes)
analysis.par$transfer_tab <- transfer_tab
# Create the final master table
analysis.par$final_ms_tab <- generate.masterTable(use_comp=all_comp,DE=analysis.par$DE,DA=analysis.par$DA,
                                                  target_list=analysis.par$merge.network$target_list,
                                                  tf_sigs=tf_sigs,z_col='Z-statistics',display_col=c('logFC','P.Value'),
                                                  main_id_type='external_gene_name')

# Path and file name of the output EXCEL file
out_file <- sprintf('%s/%s_ms_tab.xlsx',analysis.par$out.dir.DATA,analysis.par$project.name)
# Highlight marker genes
mark_gene <- list(ALS=c('SOD1', 'TARDBP', 'FUS', 'C9ORF72','HNRNPA1','OPTN','TBK1','KIF5A'),
                  Healthy=c('SOD1', 'TARDBP', 'FUS', 'C9ORF72','HNRNPA1','OPTN','TBK1','KIF5A'))
# Customize highlight color codes
#mark_col <- get.class.color(names(mark_gene)) # this randomly assign color codes
mark_col <- list('ALS'='red','Healthy'='green')
#因为前面markgene里的基因都是一样的，所以不同分组设置不同颜色最终只显示一种颜色。
# Save the final master table as EXCEL file
out2excel(analysis.par$final_ms_tab,out.xlsx = out_file,mark_gene,mark_col)

# Save Step 4 analysis.par as RData, ESSENTIAL
NetBID.saveRData(analysis.par=analysis.par,step='ms-tab')

