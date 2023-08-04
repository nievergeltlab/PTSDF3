
library(data.table)

 blue=rgb(153,191,254,max=255)
  blue=rgb(153,191,254,max=255)
 
 
  blue2=rgb(0,0,255,max=255)
   blue3=rgb(0,0,205,max=255)
 
 
 red=rgb(230,185,184,max=255)
 green=rgb(147,205, 221,max=255)
 purple=rgb(179,162,199,max=255)

adam_blue <- rgb(41,153,208,maxColorValue=255)
adam_blue2 <-  rgb(60,160,255,maxColorValue=255)
caro_orange <- rgb(249,164,24,maxColorValue=255)

beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)

#read file


 esmr <- fread('eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.esmr.mh', header=T,data.table=F)
 
 
 pqtl <- fread('eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.pqtl.mh', header=T,data.table=F)
 
 twas <-  fread('eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.TWAS.mh', header=T,data.table=F)

 library(data.table)


 #esmr <- subset(esmr,P < 1e-4)
 
#Load manhattanplot script
 print("Plotting data")
 source('ManhattanPlotterFunction_colorfixed_max10ylim2_pgc_v2.R')


colorchoice='adam_blue'
 highlight_p <- 5e-15
 highlightboundary <- 20000

 
jpeg("eur_twas.jpg",width = 13.333 , height = 3.5, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(twas$CHR, twas$BP, twas$P, twas$SNP, genomesig = 0.05/14935, pchF=16,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = "ALL", pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,12))

 dev.off()
#topsnplabels = "ALL"
colorchoice='beauty_red'

jpeg("eur_smr.jpg",width = 13.333 , height = 3.5, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(esmr$CHR, esmr$BP, esmr$P, esmr$SNP, genomesig = 0.05/9903, pchF=16,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = "ALL", pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,10))

 dev.off()

colorchoice='caro_orange'

jpeg("eur_pqtl.jpg",width = 13.333 , height = 3.5, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(pqtl$CHR, pqtl$BP, pqtl$P, pqtl$SNP, genomesig = 0.05/1210, pchF=16,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = "ALL", pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,14))

 dev.off()



