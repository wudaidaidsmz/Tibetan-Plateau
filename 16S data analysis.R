#library(DESeq2)
library(devtools)
#install.packages("rlang", version="1.1.0")
library(rlang)
library(phyloseq)
library(vegan)
library(ggplot2)
library(kableExtra)
#install.packages('kableExtra')
library(tibble)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggforce)
library(iNEXT)
library(ggpubr)
library(reshape2)
library(ggplot2)
library(ggstats)
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   #axis.line.x=element_line(size=1, colour="black"),
                   #axis.line.y=element_line(size=1, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=10),
                   axis.text.x=element_text(size=10,hjust=0,angle = 45,vjust=0),
                   axis.text.y = element_text(size=8),
                   legend.position="right",#right, none,left
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=10))


setwd('K:/16S_2025')
otu_table <- read.table("16S_ASV_filter_5.txt", header=T, row.names= 1, sep="\t")

otu_table_matrix<- as.matrix(otu_table)



OTU = otu_table(otu_table_matrix, taxa_are_rows = TRUE)


biom_table = phyloseq(OTU)
biom_table

class(biom_table)


metadata <- read.table("16S_metadata.txt", header=T, sep="\t")
table(metadata$Group)
head(metadata)
summary(metadata)
rownames(metadata) <- metadata$SampleID
mapfile <- sample_data(metadata)
class(mapfile)
biomtable_mapfile <- merge_phyloseq(biom_table, mapfile)
str(biomtable_mapfile)
biomtable_mapfile
biomtree <- read_tree("16S-rooted-tree.nwk")
summary(biomtree)
length(biomtree$tip.label)
biomtable_mapfile@otu_table
biomtable_mapfile@tax_table
physeqfull <- merge_phyloseq(biomtable_mapfile, biomtree)
physeqfull
class(physeqfull)
#is.rooted(physeqfull@phy_tree)
str(physeqfull)
physeqfull

# 计算Shannon指数
#shannon_index <- estimate_richness(physeqfull, measures = "Shannon")

# 查看计算结果
head(shannon_index)
# 获取OTU表数据
otu_table_data <- otu_table(physeqfull)























##perform alpha rarefaction cuve
# 获取测序深度范围
sample_depths <- seq(10, max(sample_sums(physeqfull)), by = 500)

# 创建空的数据框存储 rarefaction 结果
rarefaction_df <- data.frame()

# 计算不同深度的 alpha diversity
for (depth in sample_depths) {
  rarefied_physeq <- rarefy_even_depth(physeqfull, sample.size = depth, verbose = FALSE, rngseed = 123)
  
  # 计算 Shannon 和 物种丰富度（Richness）
  shannon_div <- estimate_richness(rarefied_physeq, measures = c("Shannon", "Observed"))
  
  # 合并样本信息
  shannon_div$Sample <- rownames(shannon_div)
  shannon_div$Depth <- depth
  
  rarefaction_df <- rbind(rarefaction_df, shannon_div)
}




# 绘制alpha稀释曲线
library(ggsci)
library(viridis)
library(scales)
library(RColorBrewer)
# 生成 96 种不同的颜色
my_colors <- colorRampPalette(brewer.pal(12, "Paired"))(96)



p=ggplot(rarefaction_df, aes(x = Depth, y = Observed, color = Sample)) +
  geom_smooth(method = "loess", span = 0.2, se = FALSE) +  # 使用 LOESS 平滑曲线
  #geom_point(size = 1) +
  scale_color_manual(values = my_colors) +
  #scale_color_jco() +  # 采用 JCO（Journal of Clinical Oncology）配色
  labs(title = "Alpha Rarefaction Curve",
       x = "Sequencing Depth",
       y = "Species Richness") +

  #theme_minimal() +
  theme(legend.title = element_blank())
p+main_theme


otu_table_nonrarefied <- as.data.frame(physeqfull@otu_table)
head(DataInfo(otu_table_nonrarefied))
##check sample coverage
DataInfo(otu_table_nonrarefied)
??DataInfo
write.table (DataInfo(otu_table_nonrarefied),file ="ASV_iNEXT_filter.xls", sep ="\t", row.names = T)
##all of the sample coverage are 1, so do not rarefy the ASV table
##Hill number
library(phytools)
library(hilldiv)
ASV_table_nonra <- as.data.frame(physeqfull@otu_table)
Hillnumber_diversity <- data.frame(q0=hilldiv::hill_div(ASV_table_nonra,qvalue=0),
                                   q1=hilldiv::hill_div(ASV_table_nonra,qvalue=1),
                                   q2=hilldiv::hill_div(ASV_table_nonra,qvalue=2))
Hillnumber_diversity
write.table(Hillnumber_diversity,file ='16S_Hillnumber_diversity.xls', sep ="\t", row.names = T)

#boxplot for hillnumber
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   #axis.line.x=element_line(size=1, colour="black"),#修改坐标轴线条的粗细
                   #axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),#修改坐标轴刻度线样式
                   axis.text=element_text(color="black", size=5),
                   axis.text.x=element_text(size=9, angle = 45, hjust=0,vjust=0),
                   axis.text.y = element_text(size=9),
                   legend.position="top",#right, none,left
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   axis.title.y = element_text(size=9),
                   legend.text= element_text(size=9))

data <- read.csv(file.choose())

data[1:3, 1:3]
d1 <- data[,c(2:7)]
d1 <- melt(d1[,], id.vars=c(1:2)) # reshape the data
head(d1)


p <- ggplot(d1, aes(x=Treatment, y=value, fill=Plant)) + 
  
  geom_boxplot(width =1)+
  #geom_jitter()+
  #facet_wrap(~variable, nrow=2, scale="free")+
  
  scale_color_manual(values=c('#C2C2C2','#FBE3C0','#6AAF2D','#FFE6FA',
                              '#FBE3C0','#C2C2C2','#6AAF2D','#FFE6FA',
                              '#FBE3C0','#C2C2C2','#6AAF2D','#FFE6FA'),
                     limits = c('July_nourea','July_urea','August_nourea','August_urea','October_nourea','October_urea'))+
  scale_fill_manual(values=c('#C2C2C2','#FBE3C0','#6AAF2D','#FFE6FA',
                             '#FBE3C0','#C2C2C2','#6AAF2D','#FFE6FA',
                             '#FBE3C0','#C2C2C2','#6AAF2D','#FFE6FA'))+ 
  scale_x_discrete(limits = c('July_nourea','July_urea','August_nourea','August_urea','October_nourea','October_urea'))+
  theme_bw()+
  theme(axis.title.x = element_blank(),strip.text = element_text(size = 10, face = "bold"))+
  main_theme
p + facet_wrap(vars(variable), scale="free_y", nrow = 2)


##linear mixed model for hillnumber
library(Matrix)
library(lme4)
library(car)
data <- read.csv(file.choose())
skewness(data$qpcr)
data$qpcr_log <- log(data$qpcr + 1)
skewness(data$qpcr_log)
model <- lmer(q2 ~ Plant*Treatment*Season + (1 |plot), data = data)
#model <- lmer(q2 ~ Plant*Treatment+Season + (1 |plot), data = data)
#model <- lmer(q2 ~ Plant*Treatment +(1|plot) + (1|Season),data = data)
summary(model)
fixed_effects <- as.data.frame(summary(model)$coefficients)
write.csv(fixed_effects, "fixed_effects_q2.csv")
Anova(model)
car::Anova(model,type=2)
AIC(model)
qqnorm(resid(model))
plot(model)
qqline(resid(model))
qqnorm(resid(model)); qqline(resid(model))
k.check(model)
appraise(model)


##beta diversity_NMDS
physeqfull_relativeabundance <- transform_sample_counts(physeqfull,function(x)100*x/sum(x) )
ASV_table_relativeabundance <- as.data.frame(physeqfull_relativeabundance@otu_table)
#write.table(ASV_table_relativeabundance,file ='ASV_table_relativeabundance.xls', sep ="\t", row.names = T)

ordu = ordinate(physeqfull_relativeabundance, "NMDS", 'unifrac', weighted=FALSE)

#提取PC1和PC2
nmds_scores <- scores(ordu, display = "sites")
head(nmds_scores)
write.table(nmds_scores,file ='NMDS1-NMDS2.xls', sep ="\t", row.names = T)
ordu = ordinate(physeqfull_relativeabundance, "NMDS", 'bray')
ordu$stress
p = plot_ordination(physeqfull_relativeabundance, ordu, color="Month_N",shape = 'Plant')

p = p + geom_point( size=4) +
  scale_shape_manual(values = c(15,16,17,18))+
  #stat_ellipse(aes(fill = Month_N), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)+
  #geom_polygon(aes(fill=Treatment),alpha=0.6)+
  #ggforce::geom_mark_hull(geom = "polygon", aes(group = Month_N, fill = Month_N), alpha = 0.2) +
  stat_chull(geom = "polygon", aes(group = Month_N, fill = Month_N), alpha = 0.2) +
  
  scale_color_manual(values=c('#FFE6FA','#FFC0F3','#C5E0AB','#6AAF2D','#FBE3C0','#FBC99A'))+
  scale_fill_manual(values=c('#FFE6FA','#FFC0F3','#C5E0AB','#6AAF2D','#FBE3C0','#FBC99A'))
p = p + theme_bw() +
  theme(panel.grid=element_blank())+main_theme
p

##adonis
metadata <- as(sample_data(physeqfull_relativeabundance), "data.frame")

adonis_result <- adonis2(distance(physeqfull_relativeabundance, method="unifrac") ~Plant*Treatment*Month,
                         data = metadata)

adonis_result


##NMDS for July
#extract  dataset
physeqfull_July <- subset_samples(physeqfull, Month=="July")
physeqfull_July
physeqfull_relativeabundance_July <- transform_sample_counts(physeqfull_July,function(x)100*x/sum(x) )
ordu = ordinate(physeqfull_relativeabundance_July, "NMDS", 'unifrac', weighted=TRUE)
#ordu = ordinate(physeqfull_July, "NMDS", 'bray')
ordu$stress
p = plot_ordination(physeqfull_relativeabundance_July, ordu, color="Plant",shape = 'Month_N')

p = p + geom_point( size=4) +
  scale_shape_manual(values = c(2,17,0,15,1,16))+
  #stat_ellipse(aes(fill = Month_N), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)+
  scale_color_manual(values=c("#ff717f",'#a1d5b9','#edd064','#6a73cf','red','blue'))+
  scale_fill_manual(values=c("#ff717f",'#a1d5b9','#edd064','#6a73cf','red','blue'))
#p = p + geom_point(size=4)+
# stat_ellipse(aes(fill = Plant_species), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)+
#scale_color_manual(values=c('#33a02c','#fb9a99'))+
#scale_fill_manual(values=c('#33a02c','#fb9a99'))
#scale_fill_manual(values=c("#a6cee3",'#b2df8a','#fb9a99','#cab2d6'))+
#scale_color_manual(values=c("#a6cee3",'#b2df8a','#fb9a99','#cab2d6'))
#scale_color_manual(values=c("#a6cee3","#1f78b4",'#b2df8a','#33a02c','#fb9a99','#e31a1c','#cab2d6','#fbdf6f'))+
#scale_fill_manual(values=c("#a6cee3","#1f78b4",'#b2df8a','#33a02c','#fb9a99','#e31a1c','#cab2d6','#fbdf6f'))
p = p + theme_bw() +
  theme(panel.grid=element_blank())
p



##adonis
metadata <- as(sample_data(physeqfull_relativeabundance_July), "data.frame")

adonis_result <- adonis2(distance(physeqfull_relativeabundance_July, method="unifrac") ~Plant*Treatment,
                         data = metadata)

adonis_result

##NMDS for August
#extract  dataset
physeqfull_August <- subset_samples(physeqfull, Month=="August")
physeqfull_August
physeqfull_relativeabundance_August <- transform_sample_counts(physeqfull_August,function(x)100*x/sum(x) )
ordu = ordinate(physeqfull_relativeabundance_August, "NMDS", 'unifrac', weighted=FALSE)
# 提取样本坐标
nmds_scores <- scores(ordu, display = "sites")
nmds_scores <- scores(ordu)
nmds_scores
nmds_stress <- ordu$stress
nmds_stress
print(nmds_stress)
# 查看NMDS1和NMDS2的得分
nmds1_scores <- nmds_scores[, 1]
nmds1_scores
nmds2_scores <- nmds_scores[, 2]
nmds2_scores

# 查看前几行
head(nmds_scores)
write.table(nmds_scores,file ='NMDS1-NMDS2_August.xls', sep ="\t", row.names = T)
#ordu = ordinate(physeqfull_July, "NMDS", 'bray')
ordu$stress
p = plot_ordination(physeqfull_relativeabundance_August, ordu, color="Plant",shape = 'Month_N')

p = p + geom_point( size=4) +
  scale_shape_manual(values = c(2,17,0,15,1,16))+
  #stat_ellipse(aes(fill = Month_N), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)+
  scale_color_manual(values=c("#ff717f",'#a1d5b9','#edd064','#6a73cf','red','blue'))+
  scale_fill_manual(values=c("#ff717f",'#a1d5b9','#edd064','#6a73cf','red','blue'))
#p = p + geom_point(size=4)+
# stat_ellipse(aes(fill = Plant_species), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)+
#scale_color_manual(values=c('#33a02c','#fb9a99'))+
#scale_fill_manual(values=c('#33a02c','#fb9a99'))
#scale_fill_manual(values=c("#a6cee3",'#b2df8a','#fb9a99','#cab2d6'))+
#scale_color_manual(values=c("#a6cee3",'#b2df8a','#fb9a99','#cab2d6'))
#scale_color_manual(values=c("#a6cee3","#1f78b4",'#b2df8a','#33a02c','#fb9a99','#e31a1c','#cab2d6','#fbdf6f'))+
#scale_fill_manual(values=c("#a6cee3","#1f78b4",'#b2df8a','#33a02c','#fb9a99','#e31a1c','#cab2d6','#fbdf6f'))
p = p + theme_bw() +
  theme(panel.grid=element_blank())
p



##adonis
metadata <- as(sample_data(physeqfull_relativeabundance_August), "data.frame")

adonis_result <- adonis2(distance(physeqfull_relativeabundance_August, method="unifrac") ~Plant*Treatment,
                         data = metadata)

adonis_result

##NMDS for October

#extract  dataset
physeqfull_October <- subset_samples(physeqfull, Month=="October")
physeqfull_October
physeqfull_relativeabundance_October <- transform_sample_counts(physeqfull_October,function(x)100*x/sum(x) )
ordu = ordinate(physeqfull_relativeabundance_October, "NMDS", 'unifrac', weighted=TRUE)
#ordu = ordinate(physeqfull_July, "NMDS", 'bray')
ordu$stress
p = plot_ordination(physeqfull_relativeabundance_October, ordu, color="Plant",shape = 'Month_N')

p = p + geom_point( size=4) +
  scale_shape_manual(values = c(2,17,0,15,1,16))+
  #stat_ellipse(aes(fill = Month_N), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)+
  scale_color_manual(values=c("#ff717f",'#a1d5b9','#edd064','#6a73cf','red','blue'))+
  scale_fill_manual(values=c("#ff717f",'#a1d5b9','#edd064','#6a73cf','red','blue'))
#p = p + geom_point(size=4)+
# stat_ellipse(aes(fill = Plant_species), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE)+
#scale_color_manual(values=c('#33a02c','#fb9a99'))+
#scale_fill_manual(values=c('#33a02c','#fb9a99'))
#scale_fill_manual(values=c("#a6cee3",'#b2df8a','#fb9a99','#cab2d6'))+
#scale_color_manual(values=c("#a6cee3",'#b2df8a','#fb9a99','#cab2d6'))
#scale_color_manual(values=c("#a6cee3","#1f78b4",'#b2df8a','#33a02c','#fb9a99','#e31a1c','#cab2d6','#fbdf6f'))+
#scale_fill_manual(values=c("#a6cee3","#1f78b4",'#b2df8a','#33a02c','#fb9a99','#e31a1c','#cab2d6','#fbdf6f'))
p = p + theme_bw() +
  theme(panel.grid=element_blank())
p



##adonis
metadata <- as(sample_data(physeqfull_relativeabundance_October), "data.frame")

adonis_result <- adonis2(distance(physeqfull_relativeabundance_October, method="unifrac") ~Plant*Treatment,
                         data = metadata)

adonis_result

##explained variation for vegetitave period

library(ggalluvial)
#install.packages('ggalluvial')
library(ggplot2)
library(reshape2)
library(vegan)
data <- read.csv(file.choose())


ggplot(data=data,aes(x=Month,y=value,fill=Condition),border=NA)+
  geom_bar(position = "stack",stat = 'identity',width = 0.4)+
  #geom_bar(position = "fill",stat = "identity",width = 0.7)+
  labs(x=NULL, y="Explained variation (%)")+#删掉x轴标签
  scale_fill_manual(values=c('#DEECF6','#AFC8E2','#E2F2CD','#B6DAA7','#F9D5D5','#E8E0EF','#FBC99A','#D2D2D2'))+#自定义颜色
  theme(legend.position = "right",
        panel.grid =element_blank())+
  scale_x_discrete(limits = c("July", "August", "October"))+
  guides(fill=guide_legend(title="",color="black",reverse=TRUE))+
  theme_classic()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=9))#字体设置


main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   #axis.line.x=element_line(size=1, colour="black"),
                   #axis.line.y=element_line(size=1, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=10),
                   axis.text.x=element_text(size=10,angle = 45,hjust=0,vjust=0),
                   axis.text.y = element_text(size=8),
                   legend.position="right",#right, none,left
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=10))

##
p <- ggplot(data, aes(x = Month, y = value, fill = Condition, 
                      stratum = Condition, alluvium = Condition)) +
  geom_stratum() +  #代替 geom_col() 绘制堆叠柱形图
  geom_flow(alpha = 0.5) +  #绘制同类别之间的连接线
  #facet_wrap(~group, scales = 'free_x', ncol = 2) +  #分面图
  scale_fill_manual(values=c('#FBE3C0','#6AAF2D','#FFE6FA','#C2C2C2'))+#自定义颜色
  theme(legend.position = "top",
        panel.grid =element_blank())+
  scale_x_discrete(limits = c("July", "August", "October"))+
  guides(fill=guide_legend(title="",reverse=TRUE))+
  
  theme_classic()+
  
  labs(x=NULL, y="Explained variation (%)")+
  theme(axis.text=element_text(colour='black',size=9))#字体设置
p+main_theme

##SEM analysis
setwd('K:/16S_2025')
library(lavaan)
library(plspm)
dat<-read.csv("SEM_16S_q1.csv",header = T,row.names = 1)
dat.scaled<-data.frame(scale(dat,center = F))

# model
model <- '
# regressions

pH~Removals+Treatment
NH4~Removals+Treatment
NO~Removals+Treatment
Plant~Removals+Treatment

q1~Plant+pH+NO+NH4+Removals+Treatment
'
Fit <- lavaan::sem(model, data=dat.scaled)
Fit <- plspm(model, data=dat.scaled)
summary(Fit, rsquare=T, standardized=T,fit.measures=TRUE)
summary(Fit, rsquare=T, standardized=T)

standardizedSolution(Fit)
residuals(Fit, type="cor")
modificationIndices(Fit,standardized=F)

#raw SEM plot
library(semPlot)

semPaths(Fit)
semPaths(Fit,"std",residuals=FALSE,nCharNodes=0,layout="tree2",edge.label.cex = 1,sizeMan =12)

##Effect size
data <-read.csv("effect_size_16S.csv",header = T,row.names = 1)
# 加载必要包
library(ggplot2)
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   #axis.line.x=element_line(size=1, colour="black"),#修改坐标轴线条的粗细
                   #axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),#修改坐标轴刻度线样式
                   axis.text=element_text(color="black", size=5),
                   axis.text.x=element_text(size=9,angle=45, hjust=0,vjust=0),
                   axis.text.y = element_text(size=9),
                   legend.position="bottom",#right, none,left
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   axis.title.y = element_text(size=9),
                   legend.text= element_text(size=9))

# 绘制图表
library(ggplot2)
ggplot(data, aes(x = Variable, y = Standardized_effect_size_based_on_SEM, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  facet_grid(. ~ Category, scales = "free_x") +
  scale_fill_manual(values = c("Richness" = '#6AAF2D', "Abundant_ASVs" = '#FBC99A',"Dominant_ASVs" = '#FFE6FA')) +
  theme_bw()+
  theme(axis.title.x = element_blank(),strip.text = element_text(size = 14, face = "bold"))+
  main_theme

##boxplot for 9 dominant class
setwd("K:/16S_2025")
library(ggalluvial)

#install.packages("readxl")
library(readxl)
# 替换以下路径为您的Excel文件路径
data <- read_excel("16S_ASV_tax.xlsx")

# 假设第二列的列名为 "Column2"
# 下面的代码将 "Column2" 分割成7个新的列
data_separated <- data %>%
  separate(Taxon, into = c("B", "phylum", "class", "oder", "family", "genus", "species"), sep = ";", remove = TRUE, extra = "merge")

library(dplyr)
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   #axis.line.x=element_line(size=1, colour="black"),#修改坐标轴线条的粗细
                   #axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),#修改坐标轴刻度线样式
                   axis.text=element_text(color="black", size=5),
                   axis.text.x=element_text(size=9, angle= 45,hjust=0,vjust=0),
                   axis.text.y = element_text(size=9),
                   legend.position="top",#right, none,left
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   axis.title.y = element_text(size=9),
                   legend.text= element_text(size=9))

# 计算每列的总和
column_sums <- data_separated %>%
  summarise(across(A1_0_1_J:C3_1_4_O, sum, na.rm = TRUE))

# 转换为丰度值
data_separated <- data_separated %>%
  mutate(across(A1_0_1_J:C3_1_4_O, ~ .x / column_sums[[cur_column()]]))


data_aggregated <- data_separated %>%
  group_by(phylum) %>%
  summarise(across(A1_0_1_J:C3_1_4_O, sum, na.rm = TRUE))

library(ggplot2)

# 1. 计算每个门的总丰度
total_abundance_by_phylum <- data_aggregated %>%
  rowwise() %>%
  mutate(TotalAbundance = sum(c_across(A1_0_1_J:C3_1_4_O), na.rm = TRUE)) %>%
  ungroup()

# 2. 选择除了NA之外的丰度最高的前9个门
top_8_phyla <- total_abundance_by_phylum %>%
  filter(!is.na(phylum)) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 8) %>%
  pull(phylum)

# 3. 标记其余的门和NA为 "Others"
data_aggregated <- data_aggregated %>%
  mutate(phylum = ifelse(phylum %in% top_8_phyla, phylum , "Others"))


data_aggregated <- data_aggregated %>%
  group_by(phylum) %>%
  summarise(across(A1_0_1_J:C3_1_4_O, sum, na.rm = TRUE))

write.table (data_aggregated,file ="top8_phylum.xls", sep ="\t", row.names = T)

data <- read.csv(file.choose())

data[1:3, 1:3]
d1 <- data[,c(2:11)]
d1 <- melt(d1[,], id.vars=c(1:2)) # reshape the data
head(d1)


p <- ggplot(d1, aes(x=Treatment, y=value, fill=Plant)) + 
  
  geom_boxplot(width =1)+
  #geom_jitter()+
  #facet_wrap(~variable, nrow=2, scale="free")+
  
  scale_color_manual(values=c('#C2C2C2','#FBE3C0','#6AAF2D','#FFE6FA',
                              '#FBE3C0','#C2C2C2','#6AAF2D','#FFE6FA',
                              '#FBE3C0','#C2C2C2','#6AAF2D','#FFE6FA'),
                     limits = c('July_nourea','July_urea','August_nourea','August_urea','October_nourea','October_urea'))+
  scale_fill_manual(values=c('#C2C2C2','#FBE3C0','#6AAF2D','#FFE6FA',
                             '#FBE3C0','#C2C2C2','#6AAF2D','#FFE6FA',
                             '#FBE3C0','#C2C2C2','#6AAF2D','#FFE6FA'))+ 
  scale_x_discrete(limits = c('July_nourea','July_urea','August_nourea','August_urea','October_nourea','October_urea'))+
  theme_bw()+
  theme(axis.title.x = element_blank(),strip.text = element_text(size = 9, face = "bold"))+
  main_theme
p + facet_wrap(vars(variable), scale="free_y", nrow = 4)+main_theme
##LMM for dominant class
library(Matrix)
library(lme4)
library(car)
data <- read.csv(file.choose())
model <- lmer(Pseudomonadota ~ Plant*Treatment*Season + (1 |plot), data = data)

summary(model)
fixed_effects <- as.data.frame(summary(model)$coefficients)
write.csv(fixed_effects, "fixed_effects_Pseudomonadota.csv")
car::Anova(model,type=2)
AIC(model)
## correlation analysis between env and ASVs
library(linkET)
library(vegan)
library(RColorBrewer)
library(tidyverse)
library(ggnewscale)
#install.packages('ggnewscale')
library(RColorBrewer)
#display.brewer.all()
setwd('K:/16S_2025')
varespec <-read.delim('16S_ASV.txt',row.names = 1,check.names = FALSE)
varechem <-read.delim('16S_env.txt',row.names = 1,check.names = FALSE)



mantel <- mantel_test(varespec, varechem) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.2, 0.3, Inf),
                  labels = c("< 0.2", "0.2 - 0.3", ">= 0.3")),
         pd = cut(p, breaks = c(-Inf, 0.005, 0.01, 0.05, Inf),
                  labels = c("< 0.005", "0.005 - 0.01", "0.01 - 0.05", ">= 0.05")))

p1 <- qcorrplot(correlate(varechem), 
                grid_col = "grey50",
                grid_size = 0.5,
                type = "upper", 
                diag = FALSE) +
  geom_square() +
  geom_mark(size = 4,
            only_mark = T,
            sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05,
            colour = 'white') +
  geom_couple(data = mantel,
              aes(color = pd, size = rd),  
              label.size = 3.88,
              label.family = "",
              label.fontface = 1,
              nudge_x = 0.2,
              curvature = nice_curvature(by = "from")) +    
  scale_fill_gradientn(limits = c(-0.8,0.8),
                       breaks = seq(-0.8,0.8,0.4),
                       colors = rev(brewer.pal(11, "RdBu"))) +
  scale_size_manual(values = c(0.5, 1.5, 3)) +
  scale_color_manual(values = color_pal(4, alpha = 0.6)) +  
  guides(size = guide_legend(title = "Mantel's r",                               
                             order = 2,                               
                             keyheight = unit(0.5, "cm")),           
         colour = guide_legend(title = "Mantel's p",                                  
                               order = 1,                                 
                               keyheight = unit(0.5, "cm")),           
         fill = guide_colorbar(title = "Pearson's r", 
                               keyheight = unit(2.2, "cm"),
                               keywidth = unit(0.5, "cm"),
                               order = 3)) + 
  theme(legend.box.spacing = unit(0, "pt"))
p1
##deseq2 analysis for Control_Kp
setwd('K:/16S_2025')
library(ggplot2)


# 读取第一个和第二个表格
ASV <- read.csv("16S_ASV_filter_5.csv")  # 包含 OTUID 和丰度
tax <- read.csv("tax.csv")  # 包含 OTUID 和物种注释

# 按照 OTUID 合并，保留 table1 中的 OTUID
merged_table <- merge(ASV, tax[, c("ASVID", "Taxon")], by = "ASVID", all.x = TRUE)

# 查看合并后的结果
head(merged_table)

# 将合并后的结果保存到 CSV 文件
write.csv(merged_table, "merged_table.csv", row.names = FALSE)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(SummarizedExperiment)
library(MatrixGenerics)
#install.packages('MatrixGenerics')

#install.packages('SummarizedExperiment')
#if (!require("BiocManager", quietly = TRUE))
#install.packages("MatrixGenerics")
# 分别读入OTU、metadata、tax表（数据关注公众号进群获取）
OTU <- read.table("16S_ASV_tax.txt", sep = "\t",  row.names = 1,stringsAsFactors =FALSE, check.names =FALSE,header=1)
metadata <- read.delim(file = "16S_metadata.txt", sep = '\t', stringsAsFactors = FALSE)
metadata$Plant <- as.factor(metadata$Plant)
metadata$Month <- as.factor(metadata$Month)
metadata$Treatment <- as.factor(metadata$Treatment)
# 获取OTU数据框的行名
row.names <- rownames(OTU)

# 获取OTU数据框的最后一列
last_column <- OTU[, ncol(OTU)]

# 将行名和最后一列合并成一个新的数据框
tax <- cbind(row.names, last_column)

# 给新数据框的列命名
colnames(tax) <- c("OTU", "tax")
tax <-data.frame(tax)

# 假设tax是数据框的名称，tax列包含完整的分类信息
# 首先，我们使用strsplit函数将tax列的每个条目分割成向量
tax_levels <- strsplit(as.character(tax$tax), "; ")

# 然后，我们使用lapply函数遍历每个向量，提取门级别的信息
tax_phyla <- sapply(tax_levels, function(x) {
  # 找到包含'p__'的元素
  phylum_info <- x[grep("p__", x)]
  # 返回门的名称，即'p__'后面的部分
  if (length(phylum_info) > 0) {
    return(strsplit(phylum_info, "p__")[[1]][2])
  } else {
    return(NA)  # 如果没有找到，返回NA
  }
})

# 将提取的门级别信息赋值回tax数据框的新列
tax$phylum <- tax_phyla
tax$phylum
# 删除最后一列的分类信息
OTU <- OTU[ , -ncol(OTU)]

# 过滤相对丰度低于阈值的OTU
otu_relative <- apply(OTU, 2, function(x){x/sum(x)})
threshold = 0.000# 设置阈值
idx <- rowSums(otu_relative > threshold) >= 1
otu <- as.data.frame(OTU[idx, ])
otu_relative <- as.data.frame(otu_relative[idx, ])
# DESeq2差异表达分析
# 构建DESeqDataSet对象 
dds <- DESeqDataSetFromMatrix(countData = otu, colData = metadata, design= ~ Month * Plant * Treatment ) 
#dds <- DESeqDataSetFromMatrix(countData = otu, colData = metadata, design= ~ Month + Plant + Treatment) 
metadata$Plant <- factor(metadata$Plant, levels = c("Control", "PsKpremoval"))
dds <- estimateSizeFactors(dds, type="poscounts")
# 对原始dds进行标准化
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, name="Plant_PsKpremoval_vs_Control")
# 显示dds信息
dds

# 对结果res利用order()函数按pvalue值进行排序
# order()函数先对数值排序，然后返回排序后各数值的索引，常用用法：V[order(V)]或者df[order(df$variable),]
# 对res进行排序
res = res[order(res$pvalue),]
res # 查看结果
summary(res)  # 简要统计结果

# 输出表格至本地
res <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
write.table(res, 'control_treat.DESeq_PsKp_Control.txt', col.names = NA, sep = '\t', quote = FALSE)

# 筛选差异OTU
# 首先对表格排个序，按 padj 值升序排序，相同 padj 值下继续按 log2FC 降序排序
res1 <- res[order(res$padj, res$log2FoldChange, decreasing = c(FALSE, TRUE)), ]

# log2FC≥1 & padj<0.01 标识 up，代表显著上调的OTU
# log2FC≤-1 & padj<0.01 标识 down，代表显著下调的OTU
# 其余标识 none，代表非差异的OTU
res1[which(res1$log2FoldChange >= 1 & res1$padj < 0.01),'sig'] <- 'up'
res1[which(res1$log2FoldChange <= -1 & res1$padj < 0.01),'sig'] <- 'down'
res1[which(abs(res1$log2FoldChange) <= 1 | res1$padj >= 0.01),'sig'] <- 'none'
write.table(res1, file = 'control_treat_PsKp_Control.txt', sep = '\t', col.names = NA, quote = FALSE)

# 输出选择的差异基因总表
res1_select <- subset(res1, sig %in% c('up', 'down'))
write.table(res1_select, file = 'control_treat.DESeq_PsKp_Control.select.txt', sep = '\t', col.names = NA, quote = FALSE)

# 根据 up 和 down 分开输出
res1_up <- subset(res1, sig == 'up')
res1_down <- subset(res1, sig == 'down')

write.table(res1_up, file = 'control_treat.DESeq2_PsKp_Control.up.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(res1_down, file = 'control_treat.DESeq2_PsKp_Control.down.txt', sep = '\t', col.names = NA, quote = FALSE)
# 查看fdr校正后的P<0.05的个数
table(res1$padj<0.05) 

# 提取差异OTU
diff_OTU_deseq2 <-subset(res1, padj < 0.01 & abs(log2FoldChange) > 1)
dim(diff_OTU_deseq2)
head(diff_OTU_deseq2)
write.csv(diff_OTU_deseq2, file= "DEOTU_Group1_vs_Group2_PsKp_Control.csv")

# 计算平均丰度值
abundance<-aggregate(t(otu_relative),by=list(metadata$Plant),FUN=mean)
abundance<-as.data.frame(t(abundance))
colnames(abundance)<-abundance[1,]
abundance<-abundance[-1,]
abundance<-as.data.frame(lapply(abundance, as.numeric))
res1$abundance <- apply(cbind(abundance$Control,abundance$PsKpremoval), 1, function(x){mean(x)})

# 合并数据
data <- merge(as.data.frame(res1), res1,by="row.names",sort=FALSE,all=F)

# 过滤相对丰度低于0.01%的
#data<-data[(data$abundance.x>0.0001),]

# 假设data是目标数据框，tax是包含分类信息的数据框
# 假设tax数据框中有一个名为phylum的列，其中包含门级别的分类信息
# 将data数据框的Row.names转换为字符型，如果它不是字符型的话
rownames(data) <- as.character(data$Row.names)

# 确保tax数据框的OTU列也是字符型
tax$OTU <- as.character(tax$OTU)

# 使用match函数找到data数据框的Row.names在tax的OTU列中的位置
indexes <- match(data$Row.names, tax$OTU)
indexes
# 使用这些位置索引来从tax数据框中提取phylum列的信息
data$phylum <- tax$phylum[indexes]
data$phylum
# 现在data数据框中应该有一个新的phylum列，包含了匹配的门级别分类信息

# 转换Pvalue为负对数
data$neglogp = -log(data$pvalue.x)

# 提取所需的列
data<-as.data.frame(cbind(data$Row.names, data$log2FoldChange.x, data$pvalue.x, data$phylum,data$abundance.x, data$neglogp))
colnames(data)<-c("otu","log2FoldChange","pvalue","phylum","Abundance","neglogp")

# 修改数据类型为数值型,numeric
# 查看数据类型sapply(data, mode)
data<-transform(data, Abundance = as.numeric(Abundance), neglogp = as.numeric(neglogp), pvalue= as.numeric(pvalue))
# 或者手动修改fix(data)

# 标记差异OTU类型
data$level = as.factor(ifelse(data$pvalue>=0.05, "nosig", ifelse(data$pvalue<0.05&data$log2FoldChange<0, "enriched","depleted")))

# 输出结果至本地
write.csv(data, file= 'OTU_Group1_vs_Group2_PsKp_Control.csv')
data <- read.table("deseq_PsKp_Control.txt", sep = "\t",  row.names = 1,stringsAsFactors =FALSE, check.names =FALSE,header=1)
# 调整标签顺序
label = unique(data$phylum)
label
label = label[!(label %in% "Low Abundance")] # Delete low abundance

# 设置为因子确保顺序
data$phylum = factor(data$phylum, levels = c(label, "Low Abundance"))
data$phylum
data$level = factor(data$level, levels = c("enriched","depleted","nosig"))

# 绘制曼哈顿图
# 调整y轴 使更美观 阈值自定
#data[data$neglogp>30,]$neglogp  = 30

library(ggplot2)
main_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),#修改坐标轴线条的粗细
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),#修改坐标轴刻度线样式
                   axis.text=element_text(color="black", size=5),
                   axis.text.x=element_text(size=15,angle = 45, hjust=0,vjust=0),
                   axis.text.y = element_text(size=15),
                   legend.position="none",#right, none,left
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   axis.title.y = element_text(size=15),
                   legend.text= element_text(size=15))
Title=paste('Differential ASV in' )

# OTU为x轴
p <- ggplot(data, aes(x=otu, y=neglogp, color=phylum, shape=level, size=Abundance)) +
  geom_hline(yintercept=-log(0.05), linetype=2, color="lightgrey") +
  geom_point(alpha=.6,position=position_jitter(0.5),stroke=2) +
  scale_shape_manual(values=c(17,25, 20))+
  scale_size(breaks=c(5,10,15))+
  labs(x="OTU", y="-log10(P)",title=Title)+
  theme_classic()+ 
  theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
        legend.position="right", 
        panel.grid = element_blank()) 
p

p<-ggplot(data, aes(x=phylum, y=neglogp, color=phylum, shape=level, size=Abundance)) +
  geom_hline(yintercept=-log(0.05), linetype=2, color="lightgrey") +
  geom_point(alpha=.6,position=position_jitter(0.5),stroke=2) +
  scale_shape_manual(values=c(17, 25, 20))+
  scale_color_manual(values=c('#BCD896','#F7D1E3','#F7BF6C','#C6C2DF','#2279B0','#D2DC38',
                              '#CE9FCA','#A9CCDE','#FAE93B','#C3B0D2','#EE9094','#90D0C2','#C0845C',
                              '#734E9C','#61B673',"#A6CEE3","#81B6D6","#5D9FC9","#3988BD","#297FB0","#519BA5",
                              "#78B69A","#A0D28F","#9FD57C","#7DC463","#5BB349","#39A330","#5E9E43","#949D61","#CA9B7E","#FA9695","#F37474","#ED5252",
                              "#E72F31","#E52A24","#EC563A","#F38250","#FAAE66","#FDB45C","#FDA23E"),
                     limits = c('Acidobacteriota','Actinobacteriota',"Bacteroidota","Chloroflexi","Firmicutes",'Gemmatimonadota','Planctomycetota','Proteobacteria'))+
  scale_fill_manual(values=c('#BCD896','#F7D1E3','#F7BF6C','#C6C2DF','#2279B0','#D2DC38',
                             '#CE9FCA','#A9CCDE','#FAE93B','#C3B0D2','#EE9094','#90D0C2','#C0845C',
                             '#734E9C','#61B673',"#A6CEE3","#81B6D6","#5D9FC9","#3988BD","#297FB0","#519BA5",
                             "#78B69A","#A0D28F","#9FD57C","#7DC463","#5BB349","#39A330","#5E9E43","#949D61","#CA9B7E","#FA9695","#F37474","#ED5252",
                             "#E72F31","#E52A24","#EC563A","#F38250","#FAAE66","#FDB45C","#FDA23E"))+ 
  scale_x_discrete(limits = c('Acidobacteriota','Actinobacteriota',"Bacteroidota","Chloroflexi","Firmicutes",'Gemmatimonadota','Planctomycetota','Proteobacteria'))+
  
  scale_size(breaks=c(5,10,20))+
  labs(x=NULL, y="-log10(P)",title=Title)+
  main_theme 
p 


