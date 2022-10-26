rm(list = ls())
'%!in%' <- function(x,y)!('%in%'(x,y))

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


plot.data <- data.frame(Abnormal=12,Normal=104)
plot.data <- melt(plot.data)
head(plot.data)

p <- ggplot(plot.data, aes(x="", y=value, fill=variable))
p <- p + geom_bar(width = 1, stat = "identity")
p <- p + coord_polar("y", start=0)
p <- p + scale_fill_manual(values=c("#A20056FF","#0072B5FF"))
p <- p + blank_theme + theme(axis.text.x=element_blank())
p

pdf(file = "1.normal.and.abnormal.pdf",height = 3,width = 3)
p
dev.off()


plot.data <- data.frame(Aneuploid=8,Mosaic=4)
plot.data <- melt(plot.data)
head(plot.data)

p <- ggplot(plot.data, aes(x="", y=value, fill=variable))
p <- p + geom_bar(width = 1, stat = "identity")
p <- p + coord_polar("y", start=0)
p <- p + scale_fill_manual(values=c("#DC0000FF","#EE4C97FF"))
p <- p + blank_theme + theme(axis.text.x=element_blank())
p

pdf(file = "2.aneuploid.and.mosaic.pdf",height = 3,width = 3)
p
dev.off()
