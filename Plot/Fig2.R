rm(list = ls())
'%!in%' <- function(x,y)!('%in%'(x,y))

df <- read.table("blastcyst.stat.data.txt",header = T,stringsAsFactors = F)
head(df)

####### aneuploid percent of each embryo
plot.data <- df[,c(1,3,10)]
No.mitotic.aneuploid <- plot.data$Embryo.cell.num - plot.data$Mitotic.aneuploidy.num
plot.data <- data.frame(plot.data,No.mitotic.aneuploid)
plot.data <- plot.data[,c(1,2,4)]

head(plot.data)
plot.data <- melt(plot.data)

p <- ggplot(plot.data, aes(x=Embryo,y=value,fill=variable))
p <- p + geom_bar(position="fill", stat="identity")
p <- p + scale_fill_manual(values=c("#BB0021FF","#008280FF"))
p <- p + theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                            axis.text.x = element_text(angle = 90))
p

pdf(file = "1.euploid.aneuploid.cell.number.in.each.embryo.pdf",height = 3,width = 6)
p
dev.off()

###### aneuploid rate of each embryo
head(df,n=2)
plot.data <- df[,c(1,11)]
head(plot.data)
plot.data$Mitotic.aneuploidy.per.embryo <- plot.data$Mitotic.aneuploidy.per.embryo * 100
head(plot.data)

p <- ggplot(plot.data, aes(y=Mitotic.aneuploidy.per.embryo))
p <- p + geom_boxplot(color="#BB0021FF",outlier.shape = NA)
p <- p + ylim(c(0,100))
p <- p + theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

pdf(file = "2.mitotic.aneuploidy.rate.of.each.embryo.pdf",height = 3,width = 1.5)
p
dev.off()

###### aneuploid rate of each embryo with and without meitic error
head(df,n=2)
plot.data <- df[,c(11,19)]
head(plot.data)
plot.data$Mitotic.aneuploidy.per.embryo <- plot.data$Mitotic.aneuploidy.per.embryo * 100
head(plot.data)
plot.data[plot.data$Special.event %!in% c("No"), 2] <- "withMeioticAneu"
plot.data$Special.event <- factor(plot.data$Special.event, levels = c("withMeioticAneu","No"))
plot.data

p <- ggplot(plot.data, aes(x=Special.event, y=Mitotic.aneuploidy.per.embryo, color=Special.event))
p <- p + geom_boxplot(outlier.shape = NA)
p <- p + ylim(c(0,100))
p <- p + scale_color_nejm()
p <- p + theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

pdf(file = "3.mitotic.aneuploidy.rate.of.each.embryo.withMeioticAneu.pdf",height = 3,width = 3)
p
dev.off()



##### harbor complementary embryo number

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

head(df,n=2)
table(df[,c(18)])

plot.data <- data.frame(harborComplement=14,noComplement=6)
plot.data <- melt(plot.data)
head(plot.data)

p <- ggplot(plot.data, aes(x="", y=value, fill=variable))
p <- p + geom_bar(width = 1, stat = "identity")
p <- p + coord_polar("y", start=0)
p <- p + scale_fill_manual(values=c("#6F99ADFF","#20854EFF"))
p <- p + blank_theme + theme(axis.text.x=element_blank())
p

pdf(file = "4.harbor.complementary.embryos.percent.pdf",height = 3,width = 3)
p
dev.off()
