library(VennDiagram)
library(data.table)
library(eulerr)

d1 <- fread('f3_eur_genemapping.csv',data.table=F)

posmapped <- subset(d1,posMapSNPs>0)$ensg
eqtlmapped <- subset(d1,eqtlMapSNPs>0)$ensg
cimapped <- subset(d1,ciMap=="Yes")$ensg

# Generate 3 sets of 200 words

caro_orange <- rgb(249,164,24,maxColorValue=255)
beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)

# Chart
venn.diagram(
  x = list(posmapped,eqtlmapped,cimapped),
  category.names = c("Position" , "eQTL" , "Chromatin Interaction"),
  filename = 'genemapping_strategies.png',
  fill=c(adam_blue,beauty_red,caro_orange),
  output=TRUE,
          cex = 1.2,
        fontface = "bold",
        fontfamily = "sans",
        )

pdf('f3_genemapping.pdf',7,7)
png('f3_genemapping.png')

plot(euler(list(Position = posmapped, eQTL = eqtlmapped,CI = cimapped),shape='circle'),col=c(adam_blue,beauty_red,caro_orange),quantities=TRUE)
dev.off()
