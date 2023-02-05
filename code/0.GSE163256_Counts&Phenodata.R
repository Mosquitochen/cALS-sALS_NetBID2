############################################################################===
#####  GSE67196 counts
############################################################################===

# 1 GSE67196_targets ----
library(GEOquery)
library(tidyverse)
library(limma)

GSE67196_sm <- getGEO(filename = "data/GSE67196_series_matrix.txt.gz",
                       getGPL = F) 

GSE67196_pd <- pData(GSE67196_sm)

count_raw <- read.table('data/GSE67196_Petrucelli2015_ALS_genes.rawcount.txt',header = T) 

#基因名重复
sum(duplicated(count_raw$GeneID)) #54个重复基因名

# 对重复基因名，挑选行平均值大的那一整行取最大值（通常最大值比较符合临床意义）
#计算行平均值，按降序排列
index=order(rowMeans(count_raw[,-(1:2)]),decreasing = T)
#调整表达谱的基因顺序
count_raw_ordered=count_raw[index,]
#对于有重复的基因，保留第一次出现的那个，即行平均值大的那个
keep=!duplicated(count_raw_ordered$GeneID)
#得到最后处理之后的表达谱矩阵
count_raw_max=count_raw_ordered[keep,]
count_raw_max

rownames(count_raw_max) <- count_raw_max$GeneID

count <- count_raw_max[,-(1:2)] 
count <- count[,sort(colnames(count))]


phenotype_info <- GSE67196_pd %>% separate(title, into = paste0("x", 1:3), sep = "_")

phenotype_info$x1 <- sprintf('%03d',as.numeric(phenotype_info$x1))

phenotype_info$x2 <- str_to_lower(phenotype_info$x2, locale = "en")

group <- str_replace_all(phenotype_info$`genotype:ch1`,c('c9ALS'='ALS','sALS'='ALS'))

phenotype_info <- mutate(phenotype_info, group=group, sample_name=paste0('ALS',str_c(phenotype_info$x1,phenotype_info$x2,sep = "_"))) %>%
  select(sample_name,geo_accession,`genotype:ch1`,`tissue:ch1`,group) %>% 
  setNames(c('sample_name','geo_accession','genotype','tissue','group')) %>% 
  arrange(.,sample_name)


identical(phenotype_info$sample_name,colnames(count))

colnames(count) <- rownames(phenotype_info)

sum(is.na(count))

counts <- rownames_to_column(count,'gene_name')

# 写出phenotype_info文件
write.table(phenotype_info,
            file = "data/phenotype_info.txt",
            quote = F,
            sep = "\t",
            row.names = F)

write.table(counts,
            file = "data/counts.txt",
            quote = F,
            sep = "\t",
            row.names = F)

# save
save(phenotype_info, count,counts,file = "data/phenotype_info.Rdata")



