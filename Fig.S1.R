rm(list = ls())
'%!in%' <- function(x,y)!('%in%'(x,y))

df <- read.table("blastcyst.pass.cell.vs.txt",header = T,stringsAsFactors = F)
head(df)

###### passed cells' vs distribution
head(df,n=2)
plot.data <- df
head(plot.data)
summary(plot.data$VS)

p <- ggplot(plot.data, aes(y=VS))
p <- p + geom_boxplot(color="#00008B")
p <- p + theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

pdf(file = "1.passed.cells.vs.distribution.pdf",height = 3,width = 1.5)
p
dev.off()


####### gain, loss rate of each embryo
rm(list = ls())
'%!in%' <- function(x,y)!('%in%'(x,y))

df <- read.table("../Fig2.blast/blastcyst.stat.data.txt",header = T,stringsAsFactors = F)
head(df)

head(df,n=2)
plot.data <- df[,c(15,16,17)]
head(plot.data)
plot.data <- plot.data * 100
head(plot.data)

plot.data <- melt(plot.data)
plot.data

p <- ggplot(plot.data, aes(y=value, color=variable))
p <- p + geom_boxplot(outlier.shape = NA)
p <- p + ylim(c(0,100))
p <- p + scale_color_nejm()
p <- p + theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

pdf(file = "3.gain.loss.cell.percent.of.each.embryo.pdf",height = 3,width = 4)
p
dev.off()



####### whole, seg rate of each embryo
head(df,n=2)
plot.data <- df[,c(12,13,14)]
head(plot.data)
plot.data <- plot.data * 100
head(plot.data)

plot.data <- melt(plot.data)
plot.data

p <- ggplot(plot.data, aes(y=value, color=variable))
p <- p + geom_boxplot(outlier.shape = NA)
p <- p + ylim(c(0,100))
p <- p + scale_color_nejm()
p <- p + theme_bw() + theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank())
p

pdf(file = "4.whole.seg.cell.percent.of.each.embryo.pdf",height = 3,width = 4)
p
dev.off()

