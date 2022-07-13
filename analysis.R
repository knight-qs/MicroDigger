library(phyloseq)
library(qiime2R)
library(tidyverse)
library(ggplot2)
library(vegan)
library(ggpubr)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(ggtree)
phyloseq_jx_clo2 <- qza_to_phyloseq(
  features = "table.qza",
  taxonomy = "taxonomy.qza",
  metadata = "2021_jx_soil_clo2_metadata.txt",
  tree = "rooted-tree.qza"
)

sample_names(phyloseq_jx_clo2)
rank_names(phyloseq_jx_clo2)
sample_variables(phyloseq_jx_clo2)

#中位数测序深度归一化reads数目
total <- median(sample_sums(phyloseq_jx_clo2))
standf <- function(x, t=total) round(t * (x / sum(x)))
phyloseq_jx_clo2_rarefied <- transform_sample_counts(phyloseq_jx_clo2, standf)

sample_sums(phyloseq_jx_clo2)
sample_sums(phyloseq_jx_clo2_rarefied)

# 将reads数转换为相对丰度
# a <- transform_sample_counts(phyloseq_jx_clo2, function(x) x/sum(x)*100)
# a@otu_table

# #按最小reads数抽平
# set.seed(1)
# 
# phyloseq_jx_clo2_rarefied <- rarefy_even_depth(phyloseq_jx_clo2,sample.size = 1646)
# 
# sample_sums(phyloseq_jx_clo2)
# sample_sums(phyloseq_jx_clo2_rarefied)

########################>>>>>>>>>>>>>>>>>>>
# α多样性绘图
########################>>>>>>>>>>>>>>>>>>>
alpha_plot <- plot_richness(phyloseq_jx_clo2_rarefied,
                              x = "treatment1",
                              color = "additive_amount",
                              measures = c("Observed","Chao1","Shannon","Simpson","InvSimpson","Fisher")) +
  scale_color_gradient(low = "#74a9cf", high = "#0570b0") +
  geom_point(aes(size=additive_amount)) +
  theme_bw()

alpha_plot


########################>>>>>>>>>>>>>>>>>>>
# β多样性绘图
########################>>>>>>>>>>>>>>>>>>>

# pcoa_bray（效果差）
phyloseq_jx_clo2_rarefied_ord_pcoa_bray <- ordinate(phyloseq_jx_clo2_rarefied,
                                                method = "PCoA",
                                                distance = "bray")

sampledf  <-  data.frame(sample_data(phyloseq_jx_clo2_rarefied))

jx_clo2_bray <- distance(phyloseq_jx_clo2_rarefied,"bray")

jx_clo2_adonis_bray  <- adonis(jx_clo2_bray ~ treatment3, sampledf )

plot_ordination(phyloseq_jx_clo2_rarefied,
                phyloseq_jx_clo2_rarefied_ord_pcoa_bray,
                color = "treatment3",
                #shape = "treatment3",
                title = "PCoA bray
R2 = 0.66667
P = 1") +
  #stat_ellipse(aes(fill=treatment2),geom ="polygon",alpha=0.2,linetype = 'dashed',level = 0.90) +
  geom_point(aes(size=additive_amount)) +
  #geom_polygon(aes(fill=treatment))
  #geom_line() +
  xlab("PCoA1(33.3%)") +
  ylab("PCoA2(33.3%)") +
  theme_bw()

# NMDS_bray(效果差)
phyloseq_jx_clo2_rarefied_ord_NMDS_bray <- ordinate(phyloseq_jx_clo2_rarefied,
                                                method = "NMDS",
                                                distance = "bray")

jx_clo2_tre <- get_variable(phyloseq_jx_clo2_rarefied,"treatment3")
jx_clo2_anosim_bray  <- anosim(distance( phyloseq_jx_clo2_rarefied , "bray" ), jx_clo2_tre )

plot_ordination(phyloseq_jx_clo2_rarefied,
                phyloseq_jx_clo2_rarefied_ord_NMDS_bray,
                color = "treatment3" ,
                title = "NDMS bray
Stress = 0
R = -0
P = 1") +
  geom_point(aes(size=additive_amount)) +
  theme_bw()


# pcoa_wunifrac（效果相对较好）
phyloseq_jx_clo2_rarefied_ord_pcoa_wunifrac <- ordinate(phyloseq_jx_clo2_rarefied,
                                                    method = "PCoA",
                                                    distance = "wunifrac")

sampledf  <-  data.frame(sample_data(phyloseq_jx_clo2_rarefied))

jx_clo2_unifrac_W <- distance(phyloseq_jx_clo2_rarefied,"wunifrac")

jx_clo2_adonis_UFW  <- adonis(jx_clo2_unifrac_W ~ treatment3, sampledf )

jx_clo2_adonis_UFW

beta_plot_pcoa_wunifrac <- plot_ordination(phyloseq_jx_clo2_rarefied,
                phyloseq_jx_clo2_rarefied_ord_pcoa_wunifrac,
                color = "treatment3",
                #shape = "treatment3",
                title = "PCoA unifrac weighted") +
  #stat_ellipse(aes(fill=treatment2),geom ="polygon",alpha=0.2,linetype = 'dashed',level = 0.90) +
  geom_point(aes(size=additive_amount)) +
  #geom_polygon(aes(fill=treatment))
  #geom_line() +
  xlab("PCoA1(72.2%)") +
  ylab("PCoA2(19.8%)") +
  theme_bw()

beta_plot_pcoa_wunifrac

# NMDS_wunifrac（效果差，stress值为1，过大）
phyloseq_jx_clo2_rarefied_ord_NMDS_wunifrac <- ordinate(phyloseq_jx_clo2_rarefied,
                                                    method = "NMDS",
                                                    distance = "wunifrac")

jx_clo2_tre <- get_variable(phyloseq_jx_clo2_rarefied,"treatment3")
jx_clo2_anosim_UFW  <- anosim(distance( phyloseq_jx_clo2_rarefied , "wunifrac" ), jx_clo2_tre )

plot_ordination(phyloseq_jx_clo2_rarefied,
                phyloseq_jx_clo2_rarefied_ord_NMDS_wunifrac,
                color = "treatment3",
                #shape = "treatment3",
                title = "NDMS unifrac weighted
Stress = 0
R = 1
P = 0.16667") +
  #stat_ellipse(aes(fill=treatment2),geom ="polygon",alpha=0.2,linetype = 'dashed',level = 0.90) +
  #stat_ellipse(aes(fill=mcx_tre),geom ="polygon",alpha=0.2,linetype = 'dashed',level = 0.90) +
  geom_point(aes(size=additive_amount)) +
  #geom_polygon(aes(fill=treatment))
  #geom_line() +
  theme_bw()

########################>>>>>>>>>>>>>>>>>>>
# 树(画着玩，不用)
########################>>>>>>>>>>>>>>>>>>>
plot_tree(phyloseq_jx_clo2_rarefied, color="treatment3", shape="Phylum")
plot_tree(phyloseq_jx_clo2_rarefied, nodelabf=nodeplotboot(), ladderize="left", shape="Phylum",color="treatment3")


########################>>>>>>>>>>>>>>>>>>>
# 相对丰度热图
########################>>>>>>>>>>>>>>>>>>>
# Class纲水平
rank_names(phyloseq_jx_clo2_rarefied)

class_unique <- tax_glom(phyloseq_jx_clo2_rarefied,"Class")

dim(class_unique@tax_table)[1]

tax_class <- data.frame(class_unique@tax_table)

otu_class <- data.frame(class_unique@otu_table)

otu_tax_class <- otu_class

row.names(otu_tax_class) <- tax_class[,3]

colnames(otu_tax_class) <- c("0g","125g","250g","500g")

otu_tax_class <- as.matrix(otu_tax_class)

pheatmap(otu_tax_class,
         scale = "row",
         cluster_col = FALSE,
         cellwidth = 50,
         treeheight_row = 30,
         legend = FALSE,
         fontsize_row = 5,
         fontsize_col = 5,
         angle_col = 45,
         col = brewer.pal(9,"Blues")
         )

# Class纲水平前20

class_unique <- tax_glom(phyloseq_jx_clo2_rarefied,"Class")

class_unique_20 <- prune_taxa(names(sort(taxa_sums(class_unique),TRUE)[1:20]), class_unique)

dim(class_unique_20@tax_table)[1]

tax_class_20 <- data.frame(class_unique_20@tax_table)

otu_class_20 <- data.frame(class_unique_20@otu_table)

otu_tax_class_20 <- otu_class_20

row.names(otu_tax_class_20) <- tax_class_20[,3]

colnames(otu_tax_class_20) <- c("0g","125g","250g","500g")

otu_tax_class_20 <- as.matrix(otu_tax_class_20)

pheatmap(otu_tax_class_20,
         scale = "row",
         cluster_col = FALSE,
         cellwidth = 50,
         treeheight_row = 30,
         legend = FALSE,
         fontsize_row = 5,
         fontsize_col = 5,
         angle_col = 45,
         col = brewer.pal(9,"Blues")
)

# family(存在unculture重复，无法直接更改行名，因此需拼接字符串)

family_unique <- tax_glom(phyloseq_jx_clo2_rarefied,"Family")

dim(family_unique@tax_table)[1]

tax_family <- data.frame(family_unique@tax_table)

otu_family <- data.frame(family_unique@otu_table)

otu_tax_family <- otu_family

row.names(otu_tax_family) <- str_c(tax_family[,c("Phylum")],
                                  tax_family[,c("Class")],
                                  tax_family[,c("Order")],
                                  tax_family[,c("Family")],
                                  sep = " ; ")

colnames(otu_tax_family) <- c("0g","125g","250g","500g")

otu_tax_family <- as.matrix(otu_tax_family)

pheatmap(otu_tax_family,
         scale = "row",
         cluster_col = FALSE,
         cellwidth = 50,
         treeheight_row = 30,
         legend = FALSE,
         fontsize_row = 5,
         fontsize_col = 5,
         angle_col = 45,
         col = brewer.pal(9,"Blues")
)

########################>>>>>>>>>>>>>>>>>>>
# 样品聚类树+物种堆积柱状图
########################>>>>>>>>>>>>>>>>>>>

phylum_unique <- tax_glom(phyloseq_jx_clo2_rarefied,"Phylum")

dim(phylum_unique@tax_table)[1]

tax_phylum <- data.frame(phylum_unique@tax_table)

otu_phylum <- data.frame(phylum_unique@otu_table)

otu_tax_phylum <- otu_phylum

row.names(otu_tax_phylum) <- tax_phylum[,"Phylum"]

colnames(otu_tax_phylum) <- c("0g","125g","250g","500g")

otu_tax_phylum <- as.matrix(otu_tax_phylum)

# 按四个样品中的总丰度为phylum排序

otu_tax_phylum_sort <- as.data.frame(otu_tax_phylum)

otu_tax_phylum_sort$summ <- otu_tax_phylum_sort[,1]+otu_tax_phylum_sort[,2]+otu_tax_phylum_sort[,3]+otu_tax_phylum_sort[,4]

otu_tax_phylum_sort <- otu_tax_phylum_sort[order(otu_tax_phylum_sort$summ,decreasing = T),]

# 合并前10之外的分类为others

others <- vector(mode = "numeric", length = 5)
  
for (a in 11:dim(otu_tax_phylum_sort)[1]) {
  others <- others+otu_tax_phylum_sort[a,]
}

row.names(others) <- "others"

others

otu_tax_phylum_10 <- rbind(otu_tax_phylum_sort[1:10,],others)

otu_tax_phylum_10 <- otu_tax_phylum_10[,-5]

otu_tax_phylum_10_rela <- matrix(nrow = dim(otu_tax_phylum_10)[1],
                                 ncol = 0)

for (b in 1:dim(otu_tax_phylum_10)[2]) {
  c <- otu_tax_phylum_10[,b]/sum(otu_tax_phylum_10[,b])
  otu_tax_phylum_10_rela <- cbind(otu_tax_phylum_10_rela,c)
}

otu_tax_phylum_10_rela

row.names(otu_tax_phylum_10_rela) <- row.names(otu_tax_phylum_10)
colnames(otu_tax_phylum_10_rela) <- colnames(otu_tax_phylum_10)

otu_tax_phylum_10_rela

# 首先根据分类群丰度计算样本间距离
# 本示例以群落分析中常用的 Bray-curtis 距离为例，使用 vegan 包函数 vegdist() 计算
dis_bray_phylum <- vegdist(t(otu_tax_phylum_10_rela), method = 'canberra')

# 根据距离矩阵执行层次聚类，本示例以 UPGMA 为例
tree_phylum <- hclust(dis_bray_phylum, method = 'average')
tree_phylum

# 使用 ggtree 包提供的方法绘制一个简单的聚类树用作展示
ggtree(tree_phylum) +  #聚类树
  geom_tippoint(size = 2) +  #添加节点
  geom_tiplab(hjust = -0.5) + #添加名称
  theme_tree2()  #添加刻度线

# 转换树文件格式，便于添加分组
tree_phylum <- ape::as.phylo(tree_phylum)

# 读取样本分组并将分组信息添加在聚类树中
group <- read.delim('2021_jx_soil_clo2_metadata.txt', row.names = 1, sep = '\t')
group <- group[-1,]
group <- split(row.names(group), group$treatment1)

p <- ggtree(groupOTU(tree_phylum, group), aes(color = group), size = 1) + 
  # geom_tippoint(size = 2, show.legend = FALSE) +  #添加节点
  # geom_tiplab(hjust = -0.5, show.legend = FALSE) +  #添加名称
  theme_tree2() +  #添加刻度线
  scale_color_manual(values = c('#749A53','#E2AF4D','#D8854A','#B5342C')) +
  xlim_tree(0.05)  #可以控制聚类树长度范围

p

# 在聚类树的右侧添加分类群丰度的堆叠柱形图
taxo_phylum <- data.frame(otu_tax_phylum_10_rela)
colnames(taxo_phylum) <- c("0g","125g","250g","500g")
taxo_phylum$taxonomy <- factor(rownames(taxo_phylum), levels = rev(rownames(taxo_phylum)))
taxo_phylum_melt <- reshape2::melt(taxo_phylum, id = 'taxonomy')
taxo_phylum_melt <- taxo_phylum_melt[c(2, 3, 1)]
head(taxo_phylum_melt)  #注：作图数据的第一列必须为样本名称

p2 <- p + 
  geom_facet(panel = 'Relative abundance (%)', data = taxo_phylum_melt, geom = geom_bar, 
             mapping = aes(x = 100 * value, fill = taxonomy), color = 'gray30',
             orientation = 'y', width = 0.8, stat = 'identity') +  #绘制分类群丰度堆叠图
  scale_fill_manual(values = c('gray', '#CCEBC5', '#BC80BD', '#FCCDE5', '#B3DE69', '#FDB462', 
                               '#80B1D3', '#FB8072', '#BEBADA', '#FFFFB3', '#8DD3C7'))  #赋值分类群颜色

p2 <- facet_widths(p2, widths = c(1, 8)) 

p2

ggsave('taxo_bar_16s.pdf', p2, 'pdf', height = 2, width = 12, limitsize = F)

########################>>>>>>>>>>>>>>>>>>>
# lefse分析(做不了，可能是样品太少)
########################>>>>>>>>>>>>>>>>>>>

library(microbiomeMarker)

micbio_jx_clo2 <- import_qiime2(
  otu_qza = "table.qza", taxa_qza = "taxonomy.qza",
  sam_tab = "2021_jx_soil_clo2_metadata.txt", tree_qza = "rooted-tree.qza"
)

micbio_jx_clo2

micbio_jx_clo2_lefse <- run_lefse(
  micbio_jx_clo2,
  wilcoxon_cutoff = 0.01,
  norm = "CPM",
  group = "treatment1",
  kw_cutoff = 0.01,
  multigrp_strat = TRUE,
  lda_cutoff = 4
)
# Warning message:
# No marker was identified


########################>>>>>>>>>>>>>>>>>>>
# 微生物、代谢组相关性网络分析
########################>>>>>>>>>>>>>>>>>>>

jx_GCMS <- read.csv('jx_GCMS.csv',sep = ',',header = F)
#jx_GCMS <- read.table('jx_GCMS.txt',sep = '\t',,quote = "")#也可以，必须加quote = ""才能读取完整
#jx_GCMS <- read.table('jx_GCMS.csv',sep = ',')#最后会警告，串行且不完整，加quote还是报错
#看来read.table和read.csv还是不太一样
jx_GCMS_clo <- jx_GCMS[,2:24]
jx_GCMS_clo1 <- jx_GCMS_clo[-1:-5,]
jx_GCMS_clo2 <- jx_GCMS_clo1[,c(-3,-4,-8,-10,-11)]
jx_GCMS_clo3 <- jx_GCMS_clo2
colnames(jx_GCMS_clo3) <- c('Average_Rt(min)','Average_RI',
                         'English_name','molecular_formula','CAS',
                         'Chinese_name','0_1','0_2','0_3','125_1','125_2',
                         '125_3','250_1','250_2','250_3','500_1','500_2',
                         '500_3')
row.names(jx_GCMS_clo3) <- c(1:dim(jx_GCMS_clo3)[1])

mean_0g <- vector()
for (d in 1:dim(jx_GCMS_clo3)[1]) {
  mean_0g[d] <- mean(as.numeric(jx_GCMS_clo3[d,7:9]))
}
mean_125g <- vector()
for (e in 1:dim(jx_GCMS_clo3)[1]) {
  mean_125g[e] <- mean(as.numeric(jx_GCMS_clo3[e,10:12]))
}
mean_250g <- vector()
for (f in 1:dim(jx_GCMS_clo3)[1]) {
  mean_250g[f] <- mean(as.numeric(jx_GCMS_clo3[f,13:15]))
}
mean_500g <- vector()
for (g in 1:dim(jx_GCMS_clo3)[1]) {
  mean_500g[g] <- mean(as.numeric(jx_GCMS_clo3[g,16:18]))
}
jx_GCMS_clo4 <- cbind(jx_GCMS_clo3,mean_0g,mean_125g,mean_250g,mean_500g)

jx_GCMS_clo5 <- jx_GCMS_clo4[,c('English_name','mean_0g','mean_125g','mean_250g','mean_500g')]

jx_GCMS_clo6 <- t(jx_GCMS_clo5)

colnames(jx_GCMS_clo6) <- jx_GCMS_clo6[1,]

jx_GCMS_clo7 <- jx_GCMS_clo6[-1,]

# row.names(jx_GCMS_clo7) <- c('0g','125g','250g','500g')


otu_tax_family1 <- otu_tax_family
otu_tax_family1[otu_tax_family1>0] <- 1
freq <- vector()
for (h in 1:dim(otu_tax_family1)[1]) {
  freq[h] <- sum(otu_tax_family1[h,])
}
otu_tax_family2 <- cbind(otu_tax_family,freq)
otu_tax_family3 <- subset(otu_tax_family2,freq>1)
otu_tax_family4 <- otu_tax_family3[,-5]
otu_tax_family5 <- t(otu_tax_family4)

bact_gcms_family_clo <- cbind(otu_tax_family5,jx_GCMS_clo7)

# 相关性计算，准备边文件

speacor_bact_gcms_family_clo <- matrix(nrow=0,ncol=4)

for (i in 1:66) {
  for (j in 67:273) {
    test <- cor.test(as.numeric(bact_gcms_family_clo[,i]),
                     as.numeric(bact_gcms_family_clo[,j]),
                     method = "spearman")
    rho <- test$estimate
    p.value <- test$p.value
    new.row <- c(colnames(bact_gcms_family_clo)[i],
                 colnames(bact_gcms_family_clo)[j],
                 rho,p.value)
    speacor_bact_gcms_family_clo <- rbind(speacor_bact_gcms_family_clo,new.row)
  }
}

colnames(speacor_bact_gcms_family_clo) <- c('source','target','rho','p.value')

p_threshold <- 0.1

speacor_bact_gcms_family_clo1 <- subset(as.data.frame(speacor_bact_gcms_family_clo),
                                        p.value<=p_threshold & (abs(as.numeric(rho))>=0.8))
row.names(speacor_bact_gcms_family_clo1) <- c(1:dim(speacor_bact_gcms_family_clo1)[1])

rh_po <- speacor_bact_gcms_family_clo1$rho

rh_po[rh_po>0] <- 'pos'
rh_po[rh_po<0] <- 'neg'

weight <- abs(as.numeric(speacor_bact_gcms_family_clo1$rho))

speacor_bact_gcms_family_clo2 <- cbind(speacor_bact_gcms_family_clo1[,1:4],weight,rh_po)

write.csv(speacor_bact_gcms_family_clo2,
          file = 'network.edge_list.csv',
          sep = ',',col.names = T,row.names = F)

#准备点文件

network.node_list <- data.frame(rep(NA,length(unique(speacor_bact_gcms_family_clo2[,1]))),
                                rep(NA,length(unique(speacor_bact_gcms_family_clo2[,1]))),
                                rep(NA,length(unique(speacor_bact_gcms_family_clo2[,1]))),
                                rep(NA,length(unique(speacor_bact_gcms_family_clo2[,1]))),
                                rep(NA,length(unique(speacor_bact_gcms_family_clo2[,1]))))
colnames(network.node_list) <- c('ID','Phylum','Class','Order','Family')
network.node_list[,1] <- unique(speacor_bact_gcms_family_clo2[,1])

for (k in 1:dim(network.node_list)) {
  network.node_list[k,2:5] <- stringr::str_split(network.node_list[k,1],' ; ')[[1]]
}

# 修改'uncultured' family
network.node_list1 <- network.node_list
for (l in 1:dim(network.node_list1)[1]) {
  if (network.node_list1[l,5]=='uncultured')
    network.node_list1[l,5] <- stringr::str_c(network.node_list1[l,4],
                                              network.node_list1[l,5],
                                              sep='_')
}

metabolism <- data.frame(rep(NA,length(colnames(jx_GCMS_clo7))),
                         rep(NA,length(colnames(jx_GCMS_clo7))),
                         rep(NA,length(colnames(jx_GCMS_clo7))),
                         rep(NA,length(colnames(jx_GCMS_clo7))),
                         rep(NA,length(colnames(jx_GCMS_clo7))))

metabolism[,1] <- colnames(jx_GCMS_clo7)
row.names(metabolism) <- 1:dim(metabolism)[1]

metabolism[is.na(metabolism)] <- 'metabolite'

colnames(metabolism) <- c('ID','Phylum','Class','Order','Family')

network.node_list2 <- rbind(network.node_list1,metabolism)

phylum_lable <- c(network.node_list2[1:66,2],network.node_list2[,1][67:273])
class_lable <- c(network.node_list2[1:66,3],network.node_list2[,1][67:273])
order_lable <- c(network.node_list2[1:66,4],network.node_list2[,1][67:273])
family_lable <- c(network.node_list2[1:66,5],network.node_list2[,1][67:273])

network.node_list3 <- cbind(network.node_list2,
                            phylum_lable,class_lable,order_lable,family_lable)

write.csv(network.node_list3,
          file = 'network.node_list.csv',
          sep = ',',col.names = T,row.names = F)

# 还需手动修正物种注释里的一些错误


########################>>>>>>>>>>>>>>>>>>>
# 计算网络模块内、间连通度
########################>>>>>>>>>>>>>>>>>>>
# 获取邻接矩阵方法一：边列表 -> 邻接矩阵
# tidyr::pivot_wider
speacor_bact_gcms_family_clo3 <- speacor_bact_gcms_family_clo

for (m in 1:dim(speacor_bact_gcms_family_clo)[1]) {
  n <- speacor_bact_gcms_family_clo3[m,]
  if(n[4] >= p_threshold | (abs(as.numeric(n[3]))<0.8)){
    speacor_bact_gcms_family_clo3[m,3] <- 0
  }
}

speacor_bact_gcms_family_clo4 <- speacor_bact_gcms_family_clo3[,-4]

library(tidyr)

adjacent_metrix <- as.data.frame(speacor_bact_gcms_family_clo4) %>%
pivot_wider(id_cols = source,
            names_from = target,
            values_from = rho)

adjacent_metrix1 <- as.data.frame(adjacent_metrix)

colnames(adjacent_metrix1)[1] <- NA

rownames(adjacent_metrix1) <- adjacent_metrix1[,1]

adjacent_metrix2 <- adjacent_metrix1[,-1]

o <- as.numeric(unlist(adjacent_metrix2))
o[abs(o)>0] <- 1 #这一步很关键，否则无法直接将list内非0数字全部转化为1

adjacent_metrix2.1 <- data.frame(matrix(o, ncol = length(colnames(adjacent_metrix2))))
rownames(adjacent_metrix2.1) <- rownames(adjacent_metrix2)
colnames(adjacent_metrix2.1) <- colnames(adjacent_metrix2)

row <- length(rownames(adjacent_metrix2)) #66
col <- length(colnames(adjacent_metrix2)) #207

adjacent_metrix2.2.1 <- cbind(matrix(rep(0,row^2),nrow=row,ncol=row),
                              as.matrix(adjacent_metrix2)) #注意左右顺序

adjacent_metrix2.2.2 <- rbind(adjacent_metrix2.2.1,
                              matrix(rep(0,(row+col)*col),nrow=col,ncol=(row+col))) #注意上下顺序

adjacent_metrix2.2.3 <- matrix(as.numeric(adjacent_metrix2.2.2),
                               nrow=(row+col),ncol=(row+col))

adjacent_metrix2.2.4 <- adjacent_metrix2.2.3 + t(adjacent_metrix2.2.3) - 2*diag(diag(adjacent_metrix2.2.3)) #生成邻接矩阵，权有正负

rownames(adjacent_metrix2.2.4) <- c(rownames(adjacent_metrix2),colnames(adjacent_metrix2))
colnames(adjacent_metrix2.2.4) <- c(rownames(adjacent_metrix2),colnames(adjacent_metrix2))

adjacent_metrix2.2 <- abs(adjacent_metrix2.2.4) #含权，为正

write.table(adjacent_metrix2.2, 'adjacency_weight_16S_ClO2.txt', sep = '\t') # 用于与ITS比较网络稳定性
########
# 

# 只有精确到list内的对象才可以操作
# adjacent_metrix2.1[[1]][adjacent_metrix2.1[[1]]>0] <- -1


# 获取邻接矩阵方法二：边列表 -> igraph对象 -> 邻接矩阵
# graph_from_data_frame, as_adjacency_matrix
library(igraph)

igraph_clo_16s <- graph_from_data_frame(speacor_bact_gcms_family_clo2, directed =F)

adjacent_metrix3 <- as.matrix(as_adjacency_matrix(igraph_clo_16s, type = "both"))

adjacent_metrix4 <- adjacent_metrix3[which(rownames(adjacent_metrix3) %in% colnames(otu_tax_family5)), which(colnames(adjacent_metrix3) %in% colnames(jx_GCMS_clo7))]

adjacent_metrix5 <- adjacent_metrix4[rownames(adjacent_metrix2), colnames(adjacent_metrix2)]

# adjacent_metrix5和adjacent_metrix2.1一模一样，证明两种方法都可以获得邻接矩阵
all.equal(adjacent_metrix5,as.matrix(adjacent_metrix2.1))
all(adjacent_metrix5==as.matrix(adjacent_metrix2.1))
sum(adjacent_metrix5-adjacent_metrix2.1)


##igraph 包计算网络模块
library(igraph)

#邻接矩阵 -> igraph 的邻接列表，获得含权的无向网络（方法二的不含权）
igraph <- graph_from_adjacency_matrix(as.matrix(adjacent_metrix2.2), mode = 'undirected', weighted = T, diag = F)
igraph    #igraph 的邻接列表

#计算节点度
V(igraph)$degree <- degree(igraph)
#使用 E() 和 V() 函数，可以获取 igraph 对象的边和节点

#模块划分，详情 ?cluster_fast_greedy，有多种模型
#membership(): Functions to deal with the result of network community detection
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph)) #这里要求权重全为正

#输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)

rownames(nodes_list) <- nodes_list[,1]
nodes_list <- nodes_list[,-1]

head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块


##计算模块内连通度（Zi）和模块间连通度（Pi）;用的是白鱼龙哥校正版
source('zi_pi.r')

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacent_metrix2.2), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacent_metrix2.2, degree = 'degree', modularity_class = 'modularity')

classification <- vector()
for (p in 1:dim(zi_pi)) {
  if(any(colnames(otu_tax_family5) %in% zi_pi[p,1])) {
    classification[p] <- 'microbie'
  }else if(any(colnames(jx_GCMS_clo7) %in% zi_pi[p,1])) {
    classification[p] <- 'metabolite'
  }else {
    classification[p] <- NA
  }
}
zi_pi <- data.frame(zi_pi, classification)
head(zi_pi)

write.table(zi_pi, 'zi_pi_result.txt', sep = '\t', row.names = FALSE, quote = FALSE)



##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)

zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities, shape=classification)) +
  geom_point(aes(color = type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)

########################>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 计算网络节点度，加权度，接近中心性，介数中心性，特征向量中心性
########################>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# 取含正负权重的邻接矩阵
igraph1 <- graph_from_adjacency_matrix(as.matrix(adjacent_metrix2.2.4), mode = 'undirected', weighted = T, diag = F)
igraph1

E(igraph1)$corr <- E(igraph1)$weight
E(igraph1)$weight <- abs(E(igraph1)$weight)

length(V(igraph1)$name)

V(igraph1)$degree <- degree(igraph1)
V(igraph1)$degree

degree_dist <- degree.distribution(igraph1)[-1] #第一个是度为0的节点数
degree_num <- 1:max(V(igraph1)$degree)

# 查看度是否符合幂律分布

par(mfrow = c(1, 2))
hist(V(igraph1)$degree, xlab = 'Degree', ylab = 'Frequency',breaks = 8,
     main = 'Degree distribution')
plot(degree_num, degree_dist, log = 'xy', xlab = 'Log-degree',
     ylab = 'Log-intensity', main = 'Log-log degree distribution')

Degree_Frequency <- as.data.frame(V(igraph1)$degree)
colnames(Degree_Frequency)[1] <- 'Degree_Frequency'

deg_Fre <- ggplot(Degree_Frequency, aes(x = Degree_Frequency)) +
             geom_histogram(binwidth = 4, fill = '#D73027', alpha = 1) +
             theme_minimal()
deg_Fre
ggsave('deg_Fre_16S.pdf', deg_Fre, 'pdf', height = 5, width = 2.5, limitsize = F)
  

#查看节点度与其“邻居”的平均度的关系
#微生物网络中高度值的节点更倾向连接在一起，是普遍现象吗？
neighbor_degree <- graph.knn(igraph1, V(igraph1))$knn
plot(V(igraph1)$degree, neighbor_degree, log = 'xy',
     xlab = 'Log degree', ylab = 'Log average neighbor degree')

Deg_neighbor_deg_Fre <- data.frame(Degree_Frequency, neighbor_degree)
# Deg_neighbor_deg_Fre <- log10(Deg_neighbor_deg_Fre)
Deg_neighbor_deg <- ggplot(Deg_neighbor_deg_Fre, aes(Degree_Frequency, neighbor_degree)) +
                      geom_point(color = '#4393C3') +
                      stat_smooth(method = 'lm', formula = formula2) +
                      stat_poly_eq(
                        aes(label =paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = '~~')),
                        formula = formula2,  parse = TRUE,
                        family="serif",
                        size = 3.6,
                        #color="black", #在不需要分组的层中使用恒定值覆盖它
                        #label.x = 1,  #0-1之间的比例确定位置
                        #label.y = 1
                      ) +
                      scale_x_log10() +
                      scale_y_log10() +
                      xlab('Log degree') +
                      ylab('Log average neighbor degree') +
                      theme_minimal()
Deg_neighbor_deg
ggsave('Deg_neighbor_deg_16S.pdf', Deg_neighbor_deg, 'pdf', height = 5, width = 5, limitsize = F)

#加权度（Weighted degree）
V(igraph1)$weight_degree <- strength(igraph1)
V(igraph1)$weight_degree

#接近中心性（Closeness centrality）
V(igraph1)$closeness_centrality <- closeness(igraph1)
V(igraph1)$closeness_centrality

#介数中心性（Betweenness centrality）
V(igraph1)$betweenness_centrality <- betweenness(igraph1)
V(igraph1)$betweenness_centrality

#特征向量中心性（Eigenvector centrality）
V(igraph1)$eigenvector_centrality <- evcent(igraph1)$vector
V(igraph1)$eigenvector_centrality

library(car)
library(rgl)

scatter3d(V(igraph1)$closeness_centrality, V(igraph1)$betweenness_centrality, V(igraph1)$eigenvector_centrality,
          xlab =  'Closeness centrality', ylab = 'Betweenness centrality', zlab = 'Eigenvector centrality',
          surface = T)

#输出列表
nodes_list1 <- data.frame(
  node_id = V(igraph1)$name,
  degree = V(igraph1)$degree,
  weight_degree = V(igraph1)$weight_degree,
  closeness_centrality = V(igraph1)$closeness_centrality,
  betweenness_centrality = V(igraph1)$betweenness_centrality,
  eigenvector_centrality = V(igraph1)$eigenvector_centrality)

zi_pi1 <- zi_pi[,-1]
rownames(zi_pi1) <- zi_pi[,1]
nodes_list2 <- zi_pi1[nodes_list1[,1],]

nodes_list3 <- nodes_list[nodes_list1[,1],]

nodes_list4 <- data.frame(nodes_list1, nodes_list2, nodes_list3)
nodes_list4 <- nodes_list4[,c(-1,-9,-13,-14)]

scatter3d(nodes_list4$closeness_centrality, nodes_list4$betweenness_centrality,
          nodes_list4$eigenvector_centrality, groups = as.factor(nodes_list4$type),
          xlab =  'Closeness_centrality', ylab = 'Betweenness_centrality', zlab = 'Eigenvector centrality',
          surface = F, ellipsoid = F)

# PCOA

library(vegan)

nodes_list5 <- nodes_list4[,-8:-10]

groups <- nodes_list4[,8:10]
groups$vertex_name <- rownames(groups)

bray <- vegdist(nodes_list5)

b.pcoa <- cmdscale(bray, k=(nrow(nodes_list4)-1), eig=T)

pcoa_plot_data <-  data.frame(b.pcoa$points)[,1:2]

pcoa_plot_data$vertex_name <- rownames(pcoa_plot_data)

colnames(pcoa_plot_data)[1:2] <- c('PCoA1', 'PCoA2')

eig = b.pcoa$eig

pcoa_plot_data1 <- merge(pcoa_plot_data, groups[,c(-1,-3)], by = 'vertex_name', all.x = TRUE)

library(ggplot2)
library(ggpubr)

anosim_bray  <- anosim(nodes_list5, distance = 'bray', 
                       grouping = groups$classification )

anosim_bray

net_pcoa <- ggplot(data = pcoa_plot_data1, aes(x=PCoA1, y=PCoA2, color=classification)) +
          geom_point(alpha=.7, size=2) +
          ggtitle('Bray', 'R = 0.2978   P = 0.001') +
          #geom_text(aes(label=vertex_name),size=2) +
          stat_chull() +
          theme_minimal() +
          scale_color_manual(values=c('#4393c3','#d6604d')) +
          labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
               y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))
net_pcoa
ggsave('net_pcoa_16S.pdf', net_pcoa, 'pdf', height = 5, width = 6, limitsize = F)

#散点图
library(ggpmisc)
formula1 =  y ~ x
formula2 = y ~ poly(x, 3, raw = TRUE)

ggplot(nodes_list4, aes(degree, closeness_centrality, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()


ggplot(nodes_list4, aes(degree, betweenness_centrality, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(degree, eigenvector_centrality, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

##
within_deg <- ggplot(nodes_list4, aes(degree, within_module_connectivities, color=classification)) +
                stat_smooth(method = 'lm', formula = formula1) +
                stat_poly_eq(
                  aes(label =paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = '~~')),
                  formula = formula1,  parse = TRUE,
                  family="serif",
                  size = 3.6,
                  #color="black", #在不需要分组的层中使用恒定值覆盖它
                  #label.x = 1,  #0-1之间的比例确定位置
                  #label.y = 1
                ) +
                geom_point() +
                geom_rug() +
                scale_color_manual(values=c('#4393c3','#d6604d')) +
                theme_minimal()
within_deg
ggsave('within_deg_16S.pdf', within_deg, 'pdf', height = 5, width = 6, limitsize = F)

##
among_deg <- ggplot(nodes_list4, aes(degree, among_module_connectivities, color=classification)) +
               stat_smooth(method = 'lm', formula = formula1) +
               stat_poly_eq(
                 aes(label =paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = '~~')),
                 formula = formula1,  parse = TRUE,
                 family="serif",
                 size = 3.6,
                 #color="black", #在不需要分组的层中使用恒定值覆盖它
                 #label.x = 1,  #0-1之间的比例确定位置
                 #label.y = 1
               ) +
               geom_point() +
               geom_rug() +
               scale_color_manual(values=c('#4393c3','#d6604d')) +
               theme_minimal()
among_deg
ggsave('among_deg_16S.pdf', among_deg, 'pdf', height = 5, width = 6, limitsize = F)


ggplot(nodes_list4, aes(closeness_centrality, betweenness_centrality, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(closeness_centrality, eigenvector_centrality, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(closeness_centrality, within_module_connectivities, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(closeness_centrality, among_module_connectivities, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(betweenness_centrality, eigenvector_centrality, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(betweenness_centrality, within_module_connectivities, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(betweenness_centrality, among_module_connectivities, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(eigenvector_centrality, within_module_connectivities, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

ggplot(nodes_list4, aes(eigenvector_centrality, among_module_connectivities, color=classification)) +
  geom_point() +
  scale_color_manual(values=c('#4393c3','#d6604d')) +
  theme_minimal()

##
formula1 =  y ~ x
formula2 = y ~ poly(x, 3, raw = TRUE)
within_among_deg <- ggplot(nodes_list4, aes(within_module_connectivities, among_module_connectivities, color=classification)) +
                      stat_smooth(method = 'lm', formula = formula1) +
                      stat_poly_eq(
                        aes(label =paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = '~~')),
                        formula = formula1,  parse = TRUE,
                        family="serif",
                        size = 3.6,
                        #color="black", #在不需要分组的层中使用恒定值覆盖它
                        #label.x = 1,  #0-1之间的比例确定位置
                        #label.y = 1
                        ) +
                      geom_point() +
                      geom_rug() +
                      scale_color_manual(values=c('#4393c3','#d6604d')) +
                      geom_vline(xintercept = 2.5) +
                      geom_hline(yintercept = 0.62) +
                      theme_minimal()
within_among_deg
ggsave('within_among_deg_16S.pdf', within_among_deg, 'pdf', height = 5, width = 6, limitsize = F)



########################>>>>>>>>>>>>>>>>>>>
# iTOL注释文件准备
########################>>>>>>>>>>>>>>>>>>>

#准备物种注释文件
#注意：树文件FeatureID都是带单引号的，这里也需要加上
tax_species <- data.frame(phyloseq_jx_clo2_rarefied@tax_table)

tax_species1 <- tax_species
for (m in 1:dim(tax_species1)[1]) {
  tax_species1[m,c('Kingdom')]<-stringr::str_sub(tax_species1[m,c('Kingdom')],4)
}

tax_species2 <- tax_species1
tax_species2[is.na(tax_species2)] <- 'unassigned'

FeatureID_tax <- rownames(tax_species2)
tax_species3 <- data.frame(stringr::str_c("'",FeatureID_tax,"'"),tax_species)
colnames(tax_species3)[1] <- 'FeatureID'

write.table(tax_species3,file = "taxonomy.txt",sep = "\t",row.names = F,quote = FALSE)

#准备中位数抽平后的丰度文件，每格加1后取log10
#注意：树文件FeatureID都是带单引号的，这里也需要加上
otu_species <- log10(data.frame(phyloseq_jx_clo2_rarefied@otu_table)+1)

FeatureID_otu <- rownames(otu_species)
otu_species1 <- data.frame(stringr::str_c("'",FeatureID_otu,"'"),otu_species)
colnames(otu_species1) <- c('FeatureID','0g','125g','250g','500g')

write.table(otu_species1,file = "feature_table.txt",sep = "\t",row.names = F,quote =FALSE)

otu_species1_0g <- otu_species1[,c(1,2)]
otu_species1_125g <- otu_species1[,c(1,3)]
otu_species1_250g <- otu_species1[,c(1,4)]
otu_species1_500g <- otu_species1[,c(1,5)]

write.table(otu_species1_0g,file = "feature_table_0g.txt",sep = "\t",row.names = F,quote =FALSE)
write.table(otu_species1_125g,file = "feature_table_125g.txt",sep = "\t",row.names = F,quote =FALSE)
write.table(otu_species1_250g,file = "feature_table_250g.txt",sep = "\t",row.names = F,quote =FALSE)
write.table(otu_species1_500g,file = "feature_table_500g.txt",sep = "\t",row.names = F,quote =FALSE)


source("table2itol.R")

#终端执行代码为 ctrl+alt+enter
## 方案1. 外圈颜色、形状分类和丰度方案
# annotation.txt OTU对应物种注释和丰度，
#-a 找不到输入列将终止运行（默认不执行）-c 将整数列转换为factor或具有小数点的数字，-t 偏离提示标签时转换ID列，-w 颜色带，区域宽度等， -D输出目录，-i OTU列名，-l OTU显示名称如种/属/科名，
Rscript ./table2itol.R -a -c double -D plan1 -i FeatureID -l Genus -t %s -w 0.5 taxonomy.txt
# 生成注释文件中每列为单独一个文件

## 方案2. 生成丰度柱形图注释文件
Rscript ./table2itol.R -a -d -c none -D plan2 -b Phylum -i FeatureID -l Genus -t %s -w 0.5 taxonomy.txt

## 方案3. 生成热图注释文件
Rscript ./table2itol.R -c keep -D plan3 -i FeatureID -t %s feature_table.txt

## 方案4. 将整数转化成因子生成注释文件
Rscript ./table2itol.R -a -c factor -D plan4 -i FeatureID -l Genus -t %s -w 0 taxonomy.txt

## 方案5. 自定义颜色
Rscript ./table2itol.R -a -C ./colours_1.yml -c double -D plan5 -i FeatureID -l Genus -t %s -w 0.5 taxonomy.txt

## 方案6. 单独生成热图注释文件
Rscript ./table2itol.R -c keep -D plan6_0g -i FeatureID -t %s feature_table_0g.txt
Rscript ./table2itol.R -c keep -D plan6_125g -i FeatureID -t %s feature_table_125g.txt
Rscript ./table2itol.R -c keep -D plan6_250g -i FeatureID -t %s feature_table_250g.txt
Rscript ./table2itol.R -c keep -D plan6_500g -i FeatureID -t %s feature_table_500g.txt
