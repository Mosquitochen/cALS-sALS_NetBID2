## Complete Pipeline for Advanced Analysis and Visualization ##

############### Preparation ###############

# Load NetBID2 package
library(NetBID2)

# Give file path to reload `ms-tab` RData from driver estimation step
analysis.par <- list()
# Here we use demo data from NetBID2 package
analysis.par$out.dir.DATA <- 'c9orf72/driver_Cerebellum/DATA'
NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
# To see where the plots created later will be saved
print(analysis.par$out.dir.PLOT)
#analysis.par$out.dir.PLOT <- 'test/driver/PLOT' # Users can modify this path
analysis.par$out.dir.PLOT <- "c9orf72/driver_Cerebellum/PLOT/Visualization/" # Users can modify this path

############### Part I: More details about the top drivers  ###############

### QI.1: How to get the top drivers with significant differential activity (DA) in the comparison between G4 vs. other subtypes ?

ms_tab <- analysis.par$final_ms_tab ## get the master table data frame
ms_tab <- ms_tab[which(ms_tab$Size>=30 & ms_tab$Size <=1000),]
comp_name <- 'ALS Vs Healthy' ## get the comparison name

# sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name), #显示label名字
#                                Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.4,Pv_thre=1e-8,legend_cex=0.6,
#                                main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,label_cex = 1) #默认pdf_file为Null，生成PDF
# sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
#                                Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.4,Pv_thre=1e-8,
#                                main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=FALSE) #默认pdf_file为Null，生成PDF  不显示label
# sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
#                              Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1.5,Pv_thre=1e-4,
#                              main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,label_cex = 1)   #默认pdf_file为Null，生成PDF
# sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
#                              Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1.5,Pv_thre=1e-4,
#                              main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=FALSE) #默认pdf_file为Null，生成PDF

sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
                               Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.1,Pv_thre=1e-2,
                               main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=TRUE,
                               pdf_file=sprintf('%s/vocalno_label_DA.pdf',analysis.par$out.dir.PLOT),label_cex = 1) #默认pdf_file为Null，生成PDF，这里在指定文件加下生成PDF，默认参数为pdf(pdf_file, width = 8, height = 8) 
sig_driver <- draw.volcanoPlot(dat=ms_tab,label_col='gene_label',logFC_col=sprintf('logFC.%s_DA',comp_name),
                               Pv_col=sprintf('P.Value.%s_DA',comp_name),logFC_thre=0.1,Pv_thre=1e-2,
                               main=sprintf('Volcano Plot for %s_DA',comp_name),show_label=FALSE,
                               pdf_file=sprintf('%s/vocalno_nolabel_DA.pdf',analysis.par$out.dir.PLOT)) #默认pdf_file为Null，生成PDF
sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                             Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1,Pv_thre=1e-2,
                             main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=TRUE,
                             pdf_file=sprintf('%s/vocalno_label_DE.pdf',analysis.par$out.dir.PLOT),label_cex = 1)   #默认pdf_file为Null，生成PDF
sig_gene <- draw.volcanoPlot(dat=ms_tab,label_col='geneSymbol',logFC_col=sprintf('logFC.%s_DE',comp_name),
                             Pv_col=sprintf('P.Value.%s_DE',comp_name),logFC_thre=1,Pv_thre=1e-2,
                             main=sprintf('Volcano Plot for %s_DE',comp_name),show_label=FALSE,
                             pdf_file=sprintf('%s/vocalno_nolabel_DE.pdf',analysis.par$out.dir.PLOT)) #默认pdf_file为Null，生成PDF

### QI.2: How to interpret the significance of top DA drivers ?
# Get the DE data frame of target genes
DE <- analysis.par$DE[[comp_name]]
driver_list <- rownames(sig_driver) # The rownames is the originalID_label

draw.GSEA.NetBID(DE=DE,profile_col='logFC',profile_trend='neg2pos',name_col='ID',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in Healthy',right_annotation = 'high in ALS',
                 main= comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_logFC.pdf',analysis.par$out.dir.PLOT))

draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
                 driver_list = driver_list,
                 show_label=ms_tab[driver_list,'gene_label'],
                 driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                 driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
                 target_list=analysis.par$merge.network$target_list,
                 top_driver_number=30,target_nrow=2,target_col='RdBu',
                 left_annotation = 'high in ALS',right_annotation = 'high in Healthy',
                 main= comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
                 pdf_file=sprintf('%s/NetBID_GSEA_t.pdf',analysis.par$out.dir.PLOT))

# draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
#                  driver_list = driver_list,
#                  show_label=ms_tab[driver_list,'gene_label'],
#                  driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
#                  driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
#                  target_list=analysis.par$merge.network$target_list,
#                  top_driver_number=30,target_nrow=2,target_col='black',
#                  left_annotation = 'high in G4',right_annotation = 'high in others',
#                  main= comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
#                  pdf_file=sprintf('%s/NetBID_GSEA_demo3.pdf',analysis.par$out.dir.PLOT))
# 
# draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
#                  driver_list = driver_list,
#                  show_label=ms_tab[driver_list,'gene_label'],
#                  driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
#                  driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
#                  target_list=analysis.par$merge.network$target_list,
#                  top_driver_number=30,target_nrow=1,target_col='RdBu',
#                  left_annotation = 'high in G4',right_annotation = 'high in others',
#                  main=comp_name,target_col_type='DE',Z_sig_thre=1.64,profile_sig_thre = 1.64,
#                  pdf_file=sprintf('%s/NetBID_GSEA_demo4.pdf',analysis.par$out.dir.PLOT))
# 
# draw.GSEA.NetBID(DE=DE,profile_col='t',profile_trend='pos2neg',
#                  driver_list = driver_list,
#                  show_label=ms_tab[driver_list,'gene_label'],
#                  driver_DA_Z=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
#                  driver_DE_Z=ms_tab[driver_list,sprintf('Z.%s_DE',comp_name)],
#                  target_list=analysis.par$merge.network$target_list,
#                  top_driver_number=30,target_nrow=1,target_col='black',
#                  left_annotation = 'high in G4',right_annotation = 'high in others',
#                  main= comp_name,target_col_type='PN',Z_sig_thre=1.64,profile_sig_thre = 1.64,
#                  pdf_file=sprintf('%s/NetBID_GSEA_demo5.pdf',analysis.par$out.dir.PLOT))

### QI.3: What is the expression/activity pattern of these top DA drivers across sample subtypes?

exp_mat <- exprs(analysis.par$cal.eset) # expression matrix, the rownames must be the originalID
ac_mat <- exprs(analysis.par$merge.ac.eset) # activity matrix, the rownames must be the originalID_label
phe_info <- pData(analysis.par$cal.eset) # phenotype data frame

# Users can use the other paramters in Heatmap()
draw.heatmap(mat=exp_mat,use_genes=ms_tab[driver_list,'originalID'],use_gene_label=ms_tab[driver_list,'geneSymbol'],
             use_samples=colnames(exp_mat),use_sample_label=rownames(phe_info[colnames(exp_mat),]),
             phenotype_info=phe_info,use_phe=c('genotype','tissue','group'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,show_column_names=FALSE, clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12), pdf_file=sprintf('%s/heatmap_Cereblum_Exp.pdf',analysis.par$out.dir.PLOT))

draw.heatmap(mat=ac_mat,use_genes=driver_list,use_gene_label=ms_tab[driver_list,'gene_label'],
             use_samples=colnames(exp_mat),use_sample_label=rownames(phe_info[colnames(exp_mat),]),
             phenotype_info=phe_info,use_phe=c('genotype','tissue','group'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,show_column_names=FALSE, clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 12),pdf_file=sprintf('%s/heatmap_Cereblum_AC.pdf',analysis.par$out.dir.PLOT))


# Draw heatmaps using top DE genes
gene_list <- rownames(sig_gene)
draw.heatmap(mat=exp_mat,use_genes=ms_tab[gene_list,'originalID'],use_gene_label=ms_tab[gene_list,'geneSymbol'],
             use_samples=colnames(exp_mat),use_sample_label=rownames(phe_info[colnames(exp_mat),]),
             phenotype_info=phe_info,use_phe=c('genotype','tissue','group'),main='Expression for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 10),show_column_names=FALSE)

draw.heatmap(mat=ac_mat,use_genes= gene_list,use_gene_label=ms_tab[gene_list,'gene_label'],
             use_samples=colnames(exp_mat),use_sample_label=rownames(phe_info[colnames(exp_mat),]),
             phenotype_info=phe_info,use_phe=c('genotype','tissue','group'),main='Activity for Top drivers',scale='row',
             cluster_rows=TRUE,cluster_columns=TRUE,clustering_distance_rows='pearson',clustering_distance_columns='pearson',
             row_names_gp = gpar(fontsize = 10),show_column_names=FALSE)

### QI.4: What are the biological functions of these top DA drivers ?

# Download gene sets from MSigDB and save as RData, creat a global variable all_gs2gene
gs.preload(use_spe='Homo sapiens',update=FALSE)
print(all_gs2gene_info)

# Gene Set Enrichment Analysis
driver_list_up <- rownames(sig_driver)[which(sig_driver[,2]>0)] # up
driver_list_down <- rownames(sig_driver)[which(sig_driver[,2]<0)] # down
res_up <- funcEnrich.Fisher(input_list=ms_tab[driver_list_up,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs=c('H','CP:REACTOME','BP','CGP','CP:KEGG'),
                            Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
res_down <- funcEnrich.Fisher(input_list=ms_tab[driver_list_down,'geneSymbol'],bg_list=unique(ms_tab[,'geneSymbol']),use_gs=c('H','CP:REACTOME','BP','CGP','CP:KEGG'),
                              Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
# Save gene set enrichment analysis results as EXCEl
out2excel(list(up=res_up,down=res_down),out.xlsx=sprintf('%s/fisher_res.xlsx',analysis.par$out.dir.PLOT))

# Gene set enrichment analysis Barplot
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers',pdf_file=sprintf('%s/funcEnrich_bar_nogene.pdf',analysis.par$out.dir.PLOT))
draw.funcEnrich.bar(funcEnrich_res= res_up,top_number=30,main='Function Enrichment for Top drivers',display_genes = TRUE,gs_cex=0.6,
                    pdf_file=sprintf('%s/funcEnrich_bar_withgene.pdf',analysis.par$out.dir.PLOT))

# Gene set enrichment analysis Function Cluster Plot
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 1.4,gene_cex=1.5,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_clusterBOTH.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95) #对富集结果和基因都聚类
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_clusterGS.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=TRUE,cluster_gene = FALSE)  #对富集结果聚类
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 0.8,gene_cex=0.9,pv_cex=0.8,pdf_file = sprintf('%s/funcEnrich_clusterGENE.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=FALSE,cluster_gene = TRUE)  #对基因聚类
draw.funcEnrich.cluster(funcEnrich_res= res_up,top_number=30,gs_cex = 1.5,gene_cex=1.4,pv_cex=1.2,pdf_file = sprintf('%s/funcEnrich_clusterNO.pdf',analysis.par$out.dir.PLOT),
                        cluster_gs=FALSE,cluster_gene = FALSE) #对富集结果和基因都不聚类

# Get ID conversion table
transfer_tab <- analysis.par$transfer_tab

# Bubble Plot to show target genes enriched biological functions
draw.bubblePlot(driver_list= driver_list,show_label=ms_tab[driver_list,'gene_label'],
                Z_val=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                driver_type=ms_tab[driver_list,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=5,max_gs_size=500,use_gs=c('CP:KEGG','H'),
                top_geneset_number=10,top_driver_number=10,
                pdf_file = sprintf('%s/bubblePlot.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')

draw.bubblePlot(driver_list= driver_list,show_label=ms_tab[driver_list,'gene_label'],
                Z_val=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                driver_type=ms_tab[driver_list,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=10,max_gs_size=300,use_gs=c('CP:KEGG'),
                top_geneset_number=30,top_driver_number=30,
                pdf_file = sprintf('%s/bubblePlot_KEGG.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')
# Add marker genes

mark_gene <- intersect(rownames(sig_driver),rownames(sig_gene))
draw.bubblePlot(driver_list= driver_list,show_label=ms_tab[rownames(sig_driver),'gene_label'],
                Z_val=ms_tab[driver_list,'Z.ALS Vs Healthy_DA'],
                driver_type=ms_tab[driver_list,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab,
                bg_list=ms_tab[,'geneSymbol'],min_gs_size=10,max_gs_size=300,use_gs=c('CP:KEGG','CP:BIOCARTA','H'),
                top_geneset_number=30,top_driver_number=80,
                pdf_file = sprintf('%s/bubblePlot_combine.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets',
                mark_gene=ms_tab[which(rownames(ms_tab) %in% mark_gene),'originalID_label'],gs_cex = 1,driver_cex=1.2)

############### Part II: More details about the selected driver  ###############

### QII.1: How to interpret the significance of the selected driver ?

# Get the DE file
DE <- analysis.par$DE[[comp_name]]
DE_profile <- DE$`Z-statistics`;
names(DE_profile) <- rownames(DE)

# Use the second driver in the driver list as an example
use_driver <- driver_list[1]
use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
use_target_direction <- sign(analysis.par$merge.network$target_list[[use_driver]]$spearman) ## 1/-1
annot <- sprintf('P-value: %s',signif(ms_tab[use_driver,sprintf('P.Value.%s_DA',comp_name)],2))

# Draw classic GSEA plot for one driver
draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,use_direction=use_target_direction,
          main=sprintf('GSEA plot for driver %s',ms_tab[use_driver,'gene_label']),
          pdf_file = sprintf('%s/GSEA_each_direction.pdf',analysis.par$out.dir.PLOT),
          annotation=annot,annotation_cex=1.2,
          left_annotation='high in ALS',right_annotation='high in Healthy')

# GSEA plot without direction, without annotation
draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,use_direction=NULL,
          main=sprintf('GSEA plot for driver %s',ms_tab[use_driver,'gene_label']),
          pdf_file = sprintf('%s/GSEA_each.pdf',analysis.par$out.dir.PLOT),
          annotation=annot,annotation_cex=1.2,
          left_annotation='high in ALS',right_annotation='high in Healthy')


### QII.2: How to visualize the network structure of the selected driver ?

# Use the second driver in the driver list as an example
use_driver <- driver_list[1]
# Define edges of the network
edge_score <- analysis.par$merge.network$target_list[[use_driver]]$MI*sign(analysis.par$merge.network$target_list[[use_driver]]$spearman)
names(edge_score) <- analysis.par$merge.network$target_list[[use_driver]]$target

# Draw sub-network structure of selected driver
draw.targetNet(source_label=ms_tab[use_driver,'gene_label'],source_z=ms_tab[use_driver,sprintf('Z.%s_DA',comp_name)],
               edge_score = edge_score,pdf_file=sprintf('%s/targetNet_out.pdf',analysis.par$out.dir.PLOT),label_cex = 1,n_layer=4, alphabetical_order=TRUE)

draw.targetNet(source_label=ms_tab[use_driver,'gene_label'],source_z=ms_tab[use_driver,sprintf('Z.%s_DA',comp_name)],
               edge_score = edge_score,pdf_file=sprintf('%s/targetNet_in.pdf',analysis.par$out.dir.PLOT),label_cex = 1,arrow_direction = 'in',n_layer=4)

# Draw sub-network structure for two selected drivers
# Use the first two drivers in the driver list as an example
use_driver2 <- driver_list[2]
edge_score2 <- analysis.par$merge.network$target_list[[use_driver2]]$MI*sign(analysis.par$merge.network$target_list[[use_driver2]]$spearman)
names(edge_score2) <- analysis.par$merge.network$target_list[[use_driver2]]$target
use_genes <- unique(analysis.par$merge.network$network_dat$target.symbol)
draw.targetNet.TWO(source1_label=ms_tab[use_driver,'gene_label'],edge_score1 = edge_score,
                   source2_label=ms_tab[use_driver2,'gene_label'],edge_score2 = edge_score2,
                   source1_z=ms_tab[use_driver,sprintf('Z.%s_DA',comp_name)],source2_z=ms_tab[use_driver2,sprintf('Z.%s_DA',comp_name)],
                   pdf_file=sprintf('%s/targetNetTWO.pdf',analysis.par$out.dir.PLOT),total_possible_target=use_genes,show_test=TRUE,
                   label_cex = 1,n_layer=1,source_cex=1.2,alphabetical_order=FALSE)

# To show the overlapped target genes from two selected drivers
test.targetNet.overlap(source1_label=ms_tab[use_driver,'gene_label'],source2_label=ms_tab[use_driver2,'gene_label'],
                       target1 = names(edge_score),target2 = names(edge_score2),total_possible_target=use_genes)

### QII.3: What is the expression/activity of this selected driver across subtypes of sample?

# Creates a vector of each sample’s selected phenotype descriptive information 同时生成activity和表达值差异
use_driver <- driver_list[1]
use_obs_class <- get_obs_label(phe_info = phe_info,'group')
draw.categoryValue(ac_val=ac_mat[use_driver,],exp_val=exp_mat[ms_tab[use_driver,'originalID'],],use_obs_class=use_obs_class,
                   class_order=c('ALS','Healthy'),class_srt=30,main_ac = ms_tab[use_driver,'gene_label'],main_exp=ms_tab[use_driver,'geneSymbol'],
                   pdf_file=sprintf('%s/categoryValue_ADAMTSL1_AC&EXP.pdf',analysis.par$out.dir.PLOT))
# 只生成activity表达水平差异
draw.categoryValue(ac_val=ac_mat[use_driver,],exp_val=NULL,use_obs_class=use_obs_class,class_order=c('ALS','Healthy'),
                   pdf_file=sprintf('%s/categoryValue_ADAMTSL1_AC.pdf',analysis.par$out.dir.PLOT))

# group和genotype同时比较
use_obs_class <- get_obs_label(phe_info = phe_info,c('group','genotype'))
draw.categoryValue(ac_val=ac_mat[use_driver,],exp_val=exp_mat[ms_tab[use_driver,'originalID'],],use_obs_class=use_obs_class,
                   class_srt=50,main_ac = ms_tab[use_driver,'gene_label'],main_exp=ms_tab[use_driver,'geneSymbol'],
                   pdf_file=sprintf('%s/categoryValue_ADAMTSL1_group&genotype.pdf',analysis.par$out.dir.PLOT))

### QII.4: What are the functions of the target genes of this selected driver ?

# Use the first driver in the driver list as an example
use_driver <- driver_list[1]
use_target_genes <- analysis.par$merge.network$target_list[[use_driver]]$target
transfer_tab <- analysis.par$transfer_tab
#use_target_genes <- get_name_transfertab(use_genes= use_target_genes,transfer_tab=transfer_tab,ignore_version=TRUE)	# optional
# Converts the original gene IDs into target gene IDs, with conversion table provided
bg_list <- get_name_transfertab(use_genes= unique(analysis.par$merge.network$network_dat$target),transfer_tab=transfer_tab,ignore_version=TRUE)
res <- funcEnrich.Fisher(input_list= use_target_genes,bg_list= bg_list,use_gs=c('H','CP:REACTOME','BP','CGP','CP:KEGG'),
                         Pv_thre=0.1,Pv_adj = 'none',min_gs_size = 30, max_gs_size = 500)
draw.funcEnrich.cluster(funcEnrich_res= res,top_number=20,gs_cex = 1.2,gene_cex=1,pv_cex=1,Pv_thre=0.1,
                        pdf_file = sprintf('%s/funcEnrich_clusterBOTH_%s.pdf',analysis.par$out.dir.PLOT,use_driver),
                        cluster_gs=TRUE,cluster_gene = TRUE,h=0.95)


############### Part III: Other analyses NetBID2 can do  ###############

# ### QIII.1: What are the activities of the curated gene sets across all samples?
# 
# # Preload database files into R workspace
# db.preload(use_level='gene',use_spe='human',update=FALSE)
# # Get the expression matrix, row names must be gene symbols
# exp_mat <- exprs(analysis.par$cal.eset)
# # If original is gene-based expression matrix, just use the exp_mat
# exp_mat_gene <- exp_mat
# 
# # Calculate activity for all gene sets
# use_gs2gene <- merge_gs(all_gs2gene=all_gs2gene,use_gs=c('H','CP:BIOCARTA','CP:REACTOME','CP:KEGG','C5'))
# ac_gs <- cal.Activity.GS(use_gs2gene = use_gs2gene,cal_mat = exp_mat_gene)
# 
# # Calculate DA
# phe_info <- pData(analysis.par$cal.eset)
# G0  <- rownames(phe_info)[which(phe_info$`subgroup`!='G4')] # get sample list for G0
# G1  <- rownames(phe_info)[which(phe_info$`subgroup`=='G4')] # get sample list for G1
# DA_gs_bid <- getDE.BID.2G(eset=generate.eset(ac_gs),G1=G1,G0=G0,G1_name='G4',G0_name='others')
# 
# # Draw vocalno plot for top significant gene sets
# sig_gs <- draw.volcanoPlot(dat= DA_gs_bid,label_col='ID',logFC_col='logFC',
#                            Pv_col='P.Value',logFC_thre=0.2,Pv_thre=1e-3,
#                            main='Volcano Plot for gene sets',show_label=TRUE,label_type = 'distribute',label_cex = 0.5,
#                            pdf_file=sprintf('%s/vocalno_GS_DA.pdf',analysis.par$out.dir.PLOT))
# 
# # Draw heatmap for top significant gene sets
# draw.heatmap(mat=ac_gs[sig_gs$ID,],pdf_file=sprintf('%s/heatmap_GS.pdf',analysis.par$out.dir.PLOT),scale='row',
#              phenotype_info=phe_info,use_phe=c('gender','subgroup'),
#              pre_define=c('WNT'='blue','SHH'='red','G4'='green'))
# 
# # Draw GSEA plot for top significant gene sets
# DE <- analysis.par$DE[[comp_name]]
# 
# draw.GSEA.NetBID.GS(DE=DE,name_col='ID',profile_col='t',profile_trend='pos2neg',
#                     sig_gs_list = sig_gs$ID,
#                     gs_DA_Z= DA_gs_bid[sig_gs$ID,'Z-statistics'],
#                     use_gs2gene = use_gs2gene,
#                     top_gs_number=20,target_col='RdBu',
#                     left_annotation = 'high in others',right_annotation = 'high in G4',
#                     main= comp_name,Z_sig_thre=1.64,profile_sig_thre = 0,
#                     pdf_file=sprintf('%s/NetBID_GSEA_GS_demo1.pdf',analysis.par$out.dir.PLOT))
# 
# draw.GSEA.NetBID.GS(DE=DE,name_col='ID',profile_col='t',profile_trend='pos2neg',
#                     sig_gs_list = sig_gs$ID,
#                     gs_DA_Z= DA_gs_bid[sig_gs$ID,'Z-statistics'],
#                     use_gs2gene = use_gs2gene,
#                     top_gs_number=20,target_col='black',
#                     left_annotation = 'high in others',right_annotation = 'high in G4',
#                     main= comp_name,Z_sig_thre=1.64,profile_sig_thre = 0,
#                     pdf_file=sprintf('%s/NetBID_GSEA_GS_demo2.pdf',analysis.par$out.dir.PLOT))
# 
# # Draw GSEA plot for one gene set
# DE_profile <- DE$`Z-statistics`; names(DE_profile) <- rownames(DE)
# use_gs <- sig_gs$ID[1]
# use_target_genes <- rownames(DE)[which(DE$ID %in% use_gs2gene[[use_gs]])]
# draw.GSEA(rank_profile=DE_profile,use_genes=use_target_genes,
#           main=sprintf('GSEA plot for %s',use_gs),
#           pdf_file = sprintf('%s/GSEA_GS_each.pdf',analysis.par$out.dir.PLOT),
#           left_annotation='high in G4',right_annotation='high in others',
#           annotation=sprintf('P-value: %s',signif(sig_gs[use_gs,'P.Value'],2)))
# 
# # Draw category plot for one gene set
# use_obs_class <- get_obs_label(phe_info = phe_info,'subgroup')
# draw.categoryValue(ac_val=ac_gs[use_gs,],use_obs_class=use_obs_class,
#                    class_order=c('WNT','SHH','G4'),class_srt=30,pdf_file=sprintf('%s/categoryValue_GS_demo1.pdf',analysis.par$out.dir.PLOT),
#                    main_ac= use_gs,main_cex=0.8,
#                    pre_define=c('WNT'='blue','SHH'='red','G4'='green'))

### QIII.2: How to find drivers share significantly overlapped target genes ?

gs2gene_target <- analysis.par$merge.network$target_list[driver_list]
gs2gene_target <- lapply(gs2gene_target,function(x)x$target)
transfer_tab_fake <- data.frame(from=transfer_tab[,1],to=transfer_tab[,1],gene_biotype=transfer_tab[,3],stringsAsFactors=F)
draw.bubblePlot(driver_list= driver_list,show_label=ms_tab[driver_list,'gene_label'],
                Z_val=ms_tab[driver_list,sprintf('Z.%s_DA',comp_name)],
                driver_type=ms_tab[driver_list,'gene_biotype'],
                target_list=analysis.par$merge.network$target_list,transfer2symbol2type=transfer_tab_fake,
                bg_list=ms_tab[,'geneSymbol'],gs2gene=gs2gene_target,
                min_gs_size=10,max_gs_size=1000,
                use_gs='all',
                top_geneset_number=10,top_driver_number=10,
                pdf_file = sprintf('%s/bubblePlot_overlap.pdf',analysis.par$out.dir.PLOT),
                main='Bubbleplot for top driver targets')


# ### QIII.3: How to perform functional enrichment analysis to top drivers ?
# analysis.par <- list()
# analysis.par$out.dir.DATA <- system.file('demo1','driver/DATA/',package = "NetBID2")
# NetBID.loadRData(analysis.par=analysis.par,step='ms-tab')
# ms_tab <- analysis.par$final_ms_tab
# ## get DA profile for drivers, if one driver has both _TF and _SIG,
# # choose one with larger Z-statistics
# DA_profile <- processDriverProfile(Driver_name=ms_tab$gene_label,
#                                    Driver_profile=ms_tab$logFC.G4.Vs.others_DA,
#                                    choose_strategy='absmax',
#                                    return_type ='gene_statistics')
# ## perform GSEA for gene sets by DA profile
# res1 <- funcEnrich.GSEA(rank_profile=DA_profile,
#                         use_gs=c('H'),
#                         Pv_thre=0.1,Pv_adj = 'none')

# top_gs <- res1[1,'#Name'] ## draw for the top 1
# annot <- sprintf('NES: %s \nAdjusted P-value: %s',
#                  signif(res1[1,'NES'],2),
#                  signif(res1[1,'Adj_P'],2))
# draw.GSEA(rank_profile=DA_profile,
#           use_genes=all_gs2gene$H[[top_gs]],
#           main=sprintf('GSEA plot for gene set %s',
#                        top_gs),
#           annotation=annot,annotation_cex=1.2,
#           left_annotation='high in G4',
#           right_annotation='high in others')
# 
# 
# # QIII.3: How do I draw oncoPrint figures for samples with different mutation information ?
# #   The user can use draw.oncoprint() to draw oncoPrint figures.
# 
# all_sample <- sprintf('Sample%s',1:30)
# group1_sample <- sample(all_sample,18) ## demo sample for KRAS missense mutation
# group2_sample <- sample(all_sample,12) ## demo sample for MYC amplification
# group3_sample <- sample(all_sample,10) ## demo sample for MYC missense mutation
# group4_sample <- sample(all_sample,1) ## demo sample for MYC deletion
# phenotype_info_demo <-
#   data.frame(sample =sprintf('Sample%s',1:30),
#              KRAS_MIS=ifelse(all_sample %in% group1_sample,1,0),
#              MYC_AMP=ifelse(all_sample %in% group2_sample,1,0),
#              MYC_MIS=ifelse(all_sample %in% group3_sample,1,0),
#              MYC_DEL=ifelse(all_sample %in% group4_sample,1,0))
# draw.oncoprint(phenotype_info=phenotype_info_demo,
#                Missense_column=c('KRAS_MIS','MYC_MIS'),Missense_label=c('KRAS','MYC'),
#                Amplification_column=c('MYC_AMP'),Amplification_label=c('MYC'),
#                Deletion_column=c('MYC_DEL'),Deletion_label=c('MYC'),
#                main="OncoPrint for the demo dataset")
# 
# ##
write.table(sig_driver,'/Users/mosquito/Desktop/sig_driver_Cerebllum.csv')
write.table(sig_gene,'/Users/mosquito/Desktop/sig_gene_Cerebllum.csv')
