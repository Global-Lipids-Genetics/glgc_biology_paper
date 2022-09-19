rm(list = ls())
library(tidyverse)
library(data.table)
library(ggpubr)
library(gridExtra)
library(ggsci)
theme_set(theme_pubr())

####### Gold Standard Gene Set############
load("glgc/revision/glgc_mendelian_ind_97.Rdata")
df <- new_table[,c(1:4,7:18)]
ldl <- df[df$trait == "LDL",]
nonhdl <- df[df$trait == "nonHDL",]
tc <- df[df$trait == "TC",]
logTG <- df[df$trait == "logTG",]
hdl <- df[df$trait == "HDL",]

####### Fig 1B ###################
## precision, the fraction of prioritized genes in the gold standard set
df_precision <- apply(df[,c(6:10,13,15:16)], 2, function(x) sum(str_detect(x,df$Mendelian_gene),na.rm = T))/nrow(df)

# recall, the fraction of gold standard genes in the prioritized set.
all_df <- apply(df[,c(6:10,13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
df_recall <- apply(df[,c(6:10,13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(df[,5]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_df

methods <- c("eQTL","TWAS","PoPS","Closest","DEPICT","Coding Variants","PoPS+","DEPICT+")
group <- c("Locus-based","Locus-based","Similarity-based","Locus-based","Similarity-based","Locus-based","Combined","Combined")

df_plot_data <- data.frame(df_precision,df_recall,methods,group)
df_plot_data$methods <- factor(df_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","eQTL"))
df_plot_data$group <- factor(df_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p1b <- ggplot(df_plot_data,mapping = aes(x=df_recall,y=df_precision,colour=methods,shape=methods)) + 
  geom_point(size=10,stroke = 3)+
  labs(x="Proportion of correctly identified gold standard genes among all prioritized genes",y="Proportion of prioritized genes in the gold standard set")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,2))+
  theme_classic()+
  theme(text = element_text(size=40,family="Arial"),legend.position="top",
        axis.title=element_text(size=30))+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+ 
  scale_colour_d3(name = "Methods")

ggsave(plot = p1b, filename = "fig_1b_v2.png",width = 18, height =16,dpi = 300)
round(df_plot_data[,1:2] * 100,0)

####### Fig S1A ###################

## precision, the fraction of prioritized genes in the gold standard set
df_precision <- apply(df[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,df$Mendelian_gene),na.rm = T))/nrow(df)
ldl_precision <- apply(ldl[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,ldl$Mendelian_gene),na.rm = T))/nrow(ldl)
nonhdl_precision <- apply(nonhdl[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,nonhdl$Mendelian_gene),na.rm = T))/nrow(nonhdl)
tc_precision <- apply(tc[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,tc$Mendelian_gene),na.rm = T))/nrow(tc)
logTG_precision <- apply(logTG[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,logTG$Mendelian_gene),na.rm = T))/nrow(logTG)
hdl_precision <- apply(hdl[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,hdl$Mendelian_gene),na.rm = T))/nrow(hdl)

# recall, the fraction of gold standard genes in the prioritized set.
all_df <- apply(df[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_ldl <- apply(ldl[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_nonhdl <- apply(nonhdl[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_tc <- apply(tc[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_logTG <- apply(logTG[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_hdl <- apply(hdl[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))

# df_recall <- apply(df[,c(6:13,15:16)], 2, function(x) length(unique(str_subset(x,df$Mendelian_gene))))/all_df
# ldl_recall <- apply(ldl[,c(6:13,15:16)], 2, function(x) length(unique(str_subset(x,ldl$Mendelian_gene))))/all_ldl
# nonhdl_recall <- apply(nonhdl[,c(6:13,15:16)], 2, function(x) length(unique(str_subset(x,nonhdl$Mendelian_gene))))/all_nonhdl
# tc_recall <- apply(tc[,c(6:13,15:16)], 2, function(x) length(unique(str_subset(x,tc$Mendelian_gene))))/all_tc
# logTG_recall <- apply(logTG[,c(6:13,15:16)], 2, function(x) length(unique(str_subset(x,logTG$Mendelian_gene))))/all_logTG
# hdl_recall <-  apply(hdl[,c(6:13,15:16)], 2, function(x) length(unique(str_subset(x,hdl$Mendelian_gene))))/all_hdl

df_recall <- apply(df[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(df[,5]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_df
ldl_recall <- apply(ldl[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(ldl[,5]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_ldl
nonhdl_recall <- apply(nonhdl[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(nonhdl[,5]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_nonhdl
tc_recall <- apply(tc[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(tc[,5]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_tc
logTG_recall <- apply(logTG[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(logTG[,5]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_logTG
hdl_recall <-  apply(hdl[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(hdl[,5]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_hdl

methods <- c("eQTL","TWAS","PoPS","Closest","DEPICT","eQTL_Lipid","TWAS_Lipid","Coding Variants","PoPS+","DEPICT+")
group <- c("Locus-based","Locus-based","Similarity-based","Locus-based","Similarity-based","Locus-based","Locus-based","Locus-based","Combined","Combined")

df_plot_data <- data.frame(df_precision,df_recall,methods,group)
df_plot_data$methods <- factor(df_plot_data$methods, levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
df_plot_data$group <- factor(df_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

round(df_plot_data[,1:2] * 100,0)

p1 <- ggplot(df_plot_data,mapping = aes(x=df_recall,y=df_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="All Lipid Traits")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+ 
  scale_colour_d3(name = "Methods")


ldl_plot_data <- data.frame(ldl_precision,ldl_recall,methods,group)
ldl_plot_data$methods <- factor(ldl_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
ldl_plot_data$group <- factor(ldl_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p2 <- ggplot(ldl_plot_data,mapping = aes(x=ldl_recall,y=ldl_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="LDL")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+ 
  scale_colour_d3(name = "Methods")

nonhdl_plot_data <- data.frame(nonhdl_precision,nonhdl_recall,methods,group)
nonhdl_plot_data$methods <- factor(nonhdl_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
nonhdl_plot_data$group <- factor(nonhdl_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p3 <- ggplot(nonhdl_plot_data,mapping = aes(x=nonhdl_recall,y=nonhdl_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="nonHDL")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+ 
  scale_colour_d3(name = "Methods")

logTG_plot_data <- data.frame(logTG_precision,logTG_recall,methods,group)
logTG_plot_data$methods <- factor(logTG_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
logTG_plot_data$group <- factor(logTG_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p4 <- ggplot(logTG_plot_data,mapping = aes(x=logTG_recall,y=logTG_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="logTG")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+ 
  scale_colour_d3(name = "Methods")

tc_plot_data <- data.frame(tc_precision,tc_recall,methods,group)
tc_plot_data$methods <- factor(tc_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
tc_plot_data$group <- factor(tc_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p5 <- ggplot(tc_plot_data,mapping = aes(x=tc_recall,y=tc_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="TC")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+ 
  scale_colour_d3(name = "Methods")

hdl_plot_data <- data.frame(hdl_precision,hdl_recall,methods,group)
hdl_plot_data$methods <- factor(hdl_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
hdl_plot_data$group <- factor(hdl_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p6 <- ggplot(hdl_plot_data,mapping = aes(x=hdl_recall,y=hdl_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="HDL")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,1))+
  scale_y_continuous(limits = c(0,1))+ 
  scale_colour_d3(name = "Methods")

p <- ggarrange(p1, p2, p3, p5, p6, p4, ncol=3, nrow=2,common.legend = T)
p <- annotate_figure(p,
                bottom = text_grob("Proportion of correctly identified gold standard genes \n among all prioritized genes \n ",
                                   size = 50),
                left = text_grob(" \n Proportion of prioritized genes \n in the gold standard set", rot = 90,size = 50),
                fig.lab = "A. Gold standard gene set", fig.lab.face = "bold",fig.lab.size = 56)

ggsave(plot = p, filename = "fig_s1a_v2.png",width = 48, height =32,dpi = 300)



######## Silver standard gene set #####
rm(list = ls())

load("glgc/trans_ancestry_assign_ind.Rdata")
new_table <- trans_ancestry_assign_ind[!is.na(trans_ancestry_assign_ind$mouse_model_genes),]
df <- new_table[,c(1:4,7:18)]
ldl <- df[df$trait == "LDL",]
nonhdl <- df[df$trait == "nonHDL",]
tc <- df[df$trait == "TC",]
logTG <- df[df$trait == "logTG",]
hdl <- df[df$trait == "HDL",]

## precision, the fraction of prioritized genes in the silver standard set
df_precision <- apply(df[,c(6:10,13,15:16)], 2, function(x) sum(str_detect(x,str_replace(df$mouse_model_gene,":","|")),na.rm = T))/nrow(df)


# recall, the fraction of silver standard genes in the prioritized set.
all_df <- apply(df[,c(6:10,13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
df_recall <- apply(df[,c(6:10,13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(df[,14]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_df

methods <- c("eQTL","TWAS","PoPS","Closest","DEPICT","Coding Variants","PoPS+","DEPICT+")
group <- c("Locus-based","Locus-based","Similarity-based","Locus-based","Similarity-based","Locus-based","Combined","Combined")

df_plot_data <- data.frame(df_precision,df_recall,methods,group)
df_plot_data$methods <- factor(df_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","eQTL"))
df_plot_data$group <- factor(df_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

ps1b <- ggplot(df_plot_data,mapping = aes(x=df_recall,y=df_precision,colour=methods,shape=methods)) + 
  geom_point(size=12)+
  labs(x="Proportion of correctly identified gold standard genes among all prioritized genes",y="Proportion of prioritized genes in the gold standard set")+
  scale_shape_manual(name = "Methods",values = c(16,15,16,15,17,17,17,17,17))+
  theme_classic()+
  theme(text = element_text(size=40,family="Arial"),legend.position="top")+
  scale_x_continuous(limits = c(0,0.5))+
  scale_y_continuous(limits = c(0,0.5))+ 
  scale_colour_d3(name = "Methods")

ggsave(plot = ps1b, filename = "fig_s1b.png",width = 16, height =16,dpi = 300)
round(df_plot_data[,1:2] * 100,0)

## precision, the fraction of prioritized genes in the silver standard set
df_precision <- apply(df[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,str_replace(df$mouse_model_gene,":","|")),na.rm = T))/nrow(df)
ldl_precision <- apply(ldl[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,str_replace(ldl$mouse_model_gene,":","|")),na.rm = T))/nrow(ldl)
nonhdl_precision <- apply(nonhdl[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,str_replace(nonhdl$mouse_model_gene,":","|")),na.rm = T))/nrow(nonhdl)
tc_precision <- apply(tc[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,str_replace(tc$mouse_model_gene,":","|")),na.rm = T))/nrow(tc)
logTG_precision <- apply(logTG[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,str_replace(logTG$mouse_model_gene,":","|")),na.rm = T))/nrow(logTG)
hdl_precision <- apply(hdl[,c(6:13,15:16)], 2, function(x) sum(str_detect(x,str_replace(hdl$mouse_model_gene,":","|")),na.rm = T))/nrow(hdl)

# recall, the fraction of silver standard genes in the prioritized set.
all_df <- apply(df[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_ldl <- apply(ldl[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_nonhdl <- apply(nonhdl[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_tc <- apply(tc[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_logTG <- apply(logTG[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))
all_hdl <- apply(hdl[,c(6:13,15:16)], 2, function(x) length(unique(unlist(strsplit(na.omit(x),";")))))

df_recall <- apply(df[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(df[,14]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_df
ldl_recall <- apply(ldl[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(ldl[,14]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_ldl
nonhdl_recall <- apply(nonhdl[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(nonhdl[,14]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_nonhdl
tc_recall <- apply(tc[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(tc[,14]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_tc
logTG_recall <- apply(logTG[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(logTG[,14]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_logTG
hdl_recall <-  apply(hdl[,c(6:13,15:16)], 2, function(x) sum(unique(unlist(strsplit(na.omit(hdl[,14]),";"))) %in% unique(unlist(strsplit(na.omit(x),";")))))/all_hdl

methods <- c("eQTL","TWAS","PoPS","Closest","DEPICT","eQTL_Lipid","TWAS_Lipid","Coding Variants","PoPS+","DEPICT+")
group <- c("Locus-based","Locus-based","Similarity-based","Locus-based","Similarity-based","Locus-based","Locus-based","Locus-based","Combined","Combined")

df_plot_data <- data.frame(df_precision,df_recall,methods,group)
df_plot_data$methods <- factor(df_plot_data$methods, levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
df_plot_data$group <- factor(df_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

round(df_plot_data[,1:2] * 100,0)

p1 <- ggplot(df_plot_data,mapping = aes(x=df_recall,y=df_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="All Lipid Traits")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,0.5))+
  scale_y_continuous(limits = c(0,0.5))+ 
  scale_colour_d3(name = "Methods")


ldl_plot_data <- data.frame(ldl_precision,ldl_recall,methods,group)
ldl_plot_data$methods <- factor(ldl_plot_data$methods ,levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
ldl_plot_data$group <- factor(ldl_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p2 <- ggplot(ldl_plot_data,mapping = aes(x=ldl_recall,y=ldl_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="LDL")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,0.5))+
  scale_y_continuous(limits = c(0,0.5))+ 
  scale_colour_d3(name = "Methods")

nonhdl_plot_data <- data.frame(nonhdl_precision,nonhdl_recall,methods,group)
nonhdl_plot_data$methods <- factor(nonhdl_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
nonhdl_plot_data$group <- factor(nonhdl_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p3 <- ggplot(nonhdl_plot_data,mapping = aes(x=nonhdl_recall,y=nonhdl_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="nonHDL")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,0.5))+
  scale_y_continuous(limits = c(0,0.5))+ 
  scale_colour_d3(name = "Methods")

logTG_plot_data <- data.frame(logTG_precision,logTG_recall,methods,group)
logTG_plot_data$methods <- factor(logTG_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
logTG_plot_data$group <- factor(logTG_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p4 <- ggplot(logTG_plot_data,mapping = aes(x=logTG_recall,y=logTG_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="logTG")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,0.5))+
  scale_y_continuous(limits = c(0,0.5))+ 
  scale_colour_d3(name = "Methods")

tc_plot_data <- data.frame(tc_precision,tc_recall,methods,group)
tc_plot_data$methods <- factor(tc_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
tc_plot_data$group <- factor(tc_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p5 <- ggplot(tc_plot_data,mapping = aes(x=tc_recall,y=tc_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="TC")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,0.5))+
  scale_y_continuous(limits = c(0,0.5))+ 
  scale_colour_d3(name = "Methods")

hdl_plot_data <- data.frame(hdl_precision,hdl_recall,methods,group)
hdl_plot_data$methods <- factor(hdl_plot_data$methods , levels = c("PoPS","PoPS+","DEPICT","DEPICT+","Closest","Coding Variants","TWAS","TWAS_Lipid","eQTL","eQTL_Lipid"))
hdl_plot_data$group <- factor(hdl_plot_data$group , levels = c("Locus-based","Similarity-based","Combined"))

p6 <- ggplot(hdl_plot_data,mapping = aes(x=hdl_recall,y=hdl_precision,colour=methods,shape=methods)) + 
  geom_point(size=14,stroke = 3)+
  labs(title="HDL")+
  scale_shape_manual(name = "Methods",values = c(0,15,1,16,8,5,3,4,2,6))+
  theme_classic()+
  theme(text = element_text(size=44,family="Arial"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  scale_x_continuous(limits = c(0,0.5))+
  scale_y_continuous(limits = c(0,0.5))+ 
  scale_colour_d3(name = "Methods")

p <- ggarrange(p1, p2, p3, p5, p6, p4, ncol=3, nrow=2,common.legend = T)
p <- annotate_figure(p,
                     bottom = text_grob("Proportion of correctly identified silver standard genes \n among all prioritized genes \n ",
                                        size = 50),
                     left = text_grob(" \n Proportion of prioritized genes \n in the silver standard set", rot = 90,size = 50),
                     fig.lab = "B. Silver standard gene set", fig.lab.face = "bold",fig.lab.size = 56)

ggsave(plot = p, filename = "fig_s1b_v2.png",width = 48, height =32,dpi = 300)

######text mining
textmine_pops_baseline_count <- read_csv("pubmed/textmine/pops_baseline/textmine_pops_baseline_count.csv")
textmine_pops_baseline_count$Group <- "Reference_Gene_Set"
textmine_pops_plus_count <- read_csv("pubmed/textmine/pops_plus/textmine_pops_plus_count.csv")
textmine_pops_plus_count$Group <- "Lipid_Gene_PoPS+"

df <- rbind(textmine_pops_baseline_count,textmine_pops_plus_count)
df[df$Freq == 0,]$Freq <- df[df$Freq == 0,]$Freq  + 1

grob <- grobTree(textGrob("Mann-Whitney Test P-value < 2.2e-16", x=0.3,  y=0.9, hjust = 0, 
                          gp=gpar(col="black", fontsize=13, fontface="italic")))
p <- ggplot(df, aes(x = Freq)) 
p <- p + annotation_custom(grob)
p + geom_histogram(aes(color = Group, fill = Group),
                   alpha = 0.4, position = "identity",bins = 15) +
  scale_x_continuous(trans = log10_trans())+
  xlab("Number of Lipid Related Publications") +
  ylab("Number of Genes") + 
  theme(legend.position = "bottom") + 
  theme_classic() + 
  scale_colour_npg() + 
  scale_fill_nejm()


wilcox.test(textmine_pops_baseline_count$Freq,textmine_pops_plus_count$Freq,paired=FALSE)


# Create a text
# A Mann-Whitney U test showed that there was a significant difference (W= 7964, p = 0.017) between the leg ulcer free weeks for the Clinic group compared to the group receiving the standard treatment. The median ulcer free weeks was 20 weeks for the Clinic group compared to 3.1 weeks for those receiving the standard treatment at home.
# Plot


ggsave(filename = "glgc/revision/fig_s2.png",width = 8,height = 6)

