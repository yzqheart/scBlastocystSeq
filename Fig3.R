### 
rm(list = ls())
library(data.table)
library(plyr)
library(ggplot2)
library("scales")
library(ggsci)
library(tidyverse)
library(magrittr)

show_col(pal_nejm("default")(9))
io <- list()
#io$meta <- "/mnt/data/kongsiming/zhaifan/plot/fig3/"
#meta <- read.table(paste0(io$meta,"CNV.meta.day8.10.2021.11.1.txt"),header = T,sep="\t")
#head(meta)
#dim(meta)
io$meta <- "/mnt/data/kongsiming/zhaifan/plot/clinical.part/final_v1/fig4/"
meta <- read.csv(paste0(io$meta,"Supplementary_Table_2_Sample_Information.csv"),header = T,sep=",")
head(meta)
class(meta)
meta <- as.data.table(meta)
###################################################################################
data <- meta[meta$Day !="D6",] %>% 
  .[.$Day !="D12",] %>% 
  .[.$Day !="D14",] %>% 
  setDT(data) 
Aneu <- data[,.N,by=c("CNV")] 
Aneu
Aneu[,all := "embryo"]
Aneu[,ratio := c(round(Aneu[1,2]/(Aneu[1,2]+Aneu[2,2])*100,1),round(Aneu[2,2]/(Aneu[1,2]+Aneu[2,2])*100,1))]

pdf("/mnt/data/kongsiming/zhaifan/plot/final_v1/fig4/fig4a.RNA.aneu.ratio.in.all.embryos.pdf")
ggplot(data=Aneu,mapping = aes(x='CNV',y=ratio,fill=CNV))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values=c("#FEA842","#059D87"))+
  coord_polar(theta = 'y')+ labs(x = '', y = '', title = 'percentage of loss/gain in post implantation')+ 
  theme_bw() +
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(axis.text = element_blank())+
  theme(axis.ticks = element_blank()) 

dev.off()


## Aneuratio in split Day
data[,Day_Embryo := paste0(Day,"_",Embryo)]
DAneu <- data[,.N,by=c("Day_Embryo","CNV")] 
DAneu
long <- dcast(DAneu, Day_Embryo ~ CNV)  
long[is.na(long)] <- 0
long <- long %>%
  as.data.table() %>%
  .[,Aneuratio := Abnormal/(Abnormal+Normal)] %>%
  .[,Euratio := Normal/(Abnormal+Normal)]
long[,day := do.call(rbind, strsplit(as.character(long$Day_Embryo),"_"))[,1]]
long <-  long[long$Day_Embryo !="D10_E10",]
write.table(long,"/mnt/data/kongsiming/zhaifan/plot/final_v1/fig4/fig4a.RNA.embryo.in.boxplot.11.21.csv",append=FALSE,quote=FALSE, col.names=TRUE, row.names=FALSE)

## add pvalue
my_comparisons <- list( c("D8", "D10"))
short <- melt(long[,c(1,4,5,6)])
gastrulation <- data.frame(Day_Embryo=c("Gastrulation","Gastrulation"),day=c("D19","D19"),variable=c("Aneuratio","Euratio"),value = c(0.062,0.938))
short <- rbind(short,gastrulation)
name.order <- c("D8", "D10","D19")
pdf("/mnt/data/kongsiming/zhaifan/plot/final_v1/fig4/fig4a.RNA.embryo.in.split.day.11.21.pdf",width = 4,height = 4)
ggplot(short[short$variable=="Aneuratio",], aes(x=day, y = value, color = day)) + 
  geom_boxplot(outlier.shape = NA)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white", colour = "black")) +
  geom_dotplot(binaxis='y',dotsize=0.5,fill=NA)+
  scale_x_discrete(limits=name.order)+
  scale_color_manual(values=c("#EA5F5C","#37649E"))+
  labs(x = "Stage", y = "Aneuploid rate (%)") +
  ggtitle('gain ratio in different stages')+
  stat_compare_means(method="t.test",comparisons = my_comparisons,label.y = c(1, 0.9, 0.85))+ 
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),axis.text.x = element_text(size = 14,angle = 45, hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size = 14))
dev.off()

name.order <- c("D8_E1","D8_E2","D8_E3","D8_E4","D8_E5","D10_E1","D10_E2","D10_E3","D10_E4",
                "D10_E5","D10_E6","D10_E7","D10_E8","D10_E9","Gastrulation")
pdf("/mnt/data/kongsiming/zhaifan/plot/clinical.part/final_v1/fig4/fig4d.per.embryo.D8.D19.pdf",width = 6,height = 6)
ggplot(data = short, mapping = aes(x = Day_Embryo, y = value, fill = variable)) + 
  geom_bar(stat= 'identity', position = 'stack',width = 0.5)+
  scale_fill_manual(values=c("#BB0021FF","#008280FF"))+
  geom_col(width = 0.8)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  ggtitle('Pmat rate in different stages')+
  scale_x_discrete(limits=name.order)+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),axis.text.x = element_text(size = 14,angle = 45, hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size = 14))
dev.off()

### Day barplot

ratio <- data.frame(day=c("D8","D10","Gastrulation"),aneu.ratio=c(0.35294724,0.358362411,0.062))

library(Rmisc)
df2_count <- summarySE(df, measurevar = "len",
      
                                        groupvars = c("supp","dose"))
short.aneu <- short[short$variable == "Aneuratio",]
sd.E8 <-  sd(short.aneu[short.aneu$day=="D8",]$value, na.rm=TRUE)
sd.E10 <-  sd(short.aneu[short.aneu$day=="D10",]$value, na.rm=TRUE)

ratio$sd <- c(sd.E8,sd.E10,NA)
name.order <- c("D8","D10","Gastrulation")
pdf("/mnt/data/kongsiming/zhaifan/plot/clinical.part/final_v1/fig4/fig4d.per.embryo.D8.D19.2.pdf",width = 6,height = 6)
ggplot(data = ratio, mapping = aes(x = day, y = aneu.ratio, fill = day)) + 
  geom_bar(stat= 'identity',width = 0.5)+
  scale_fill_manual(values=c("#0072B5FF","#BC3C29FF","#E18727FF"))+
  geom_errorbar(aes(ymin = aneu.ratio - sd, ymax = aneu.ratio + sd), 
                width = 0.2, position = position_dodge(0.9),color="black")+
  
  geom_col(width = 0.8)+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  ggtitle('Pmat rate in different stages')+
  scale_x_discrete(limits=name.order)+
  scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.25))+
  theme(axis.title.x = element_text(size = 14),axis.title.y = element_text(size = 14),axis.text.x = element_text(size = 14,angle = 45, hjust = 0.5, vjust = 0.5),axis.text.y = element_text(size = 14))
dev.off()
