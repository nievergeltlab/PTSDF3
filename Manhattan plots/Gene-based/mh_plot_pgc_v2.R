
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

#read file


 library(data.table)

 genebased <- fread('eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.GENEBASED.mh',data.table=F)

 #esmr <- subset(esmr,P < 1e-4)
 
#Load manhattanplot script
 print("Plotting data")
 source('ManhattanPlotterFunction_colorfixed_max10ylim2_pgc_v2.R')


colorchoice='adam_blue2'
 highlight_p <- 5e-10
 highlightboundary <- 20000

 
jpeg("eur_gene.jpg",width = 13.333 , height = 7, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(genebased$CHR, genebased$BP, genebased$P, genebased$SNP, genomesig = 0.05/19106, pchF=16,genomesug = 0.05/19106,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = "NO", pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,15))

 dev.off()


 
jpeg("eur_twas.jpg",width = 13.333 , height = 3.5, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(twas$CHR, twas$BP, twas$P, twas$SNP, genomesig = 0.05/14935, pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = adam_blue2, sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,10))

 dev.off()


jpeg("eur_esmr.jpg",width = 13.333 , height = 3.5, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(esmr$CHR, esmr$BP, esmr$P, esmr$SNP, genomesig = 0.05/9903, pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = adam_blue2, sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,10))

 dev.off()

jpeg("eur_pqtl.jpg",width = 13.333 , height = 3.5, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(pqtl$CHR, pqtl$BP, pqtl$P, pqtl$SNP, genomesig = 0.05/9903, pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = adam_blue2, sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,10))

 dev.off()



  jpeg("eur_pli.jpg",width = 13.333 , height = 3.5, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(pli$CHR, pli$BP, pli$P, pli$SNP, genomesig = 0.9, pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = adam_blue2, sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,1))

 dev.off()


    jpeg("eur_cadd.jpg",width = 13.333 , height = 7, units = "in", res = 500,quality = 100)

     ManhattanPlot_AJS_cut(cadd$CHR, cadd$BP, cadd$P, cadd$SNP, genomesig = 4.265795e-13, pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = adam_blue2, sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary,ylims=c(0,30))

 dev.off()




    ManhattanPlot_AJS_cut(cadd$CHR, cadd$BP, cadd$P, cadd$SNP, genomesig = 0.05/19106, pchF=1,genomesug = 0, replot=1,photoname = '', 
    outname = outfile, colors =c('red','red'), sigcolor = 'red', sugcolor = 'red', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary)
    
# c(rgb(25,25,25,max=255),eval(parse(text = colorchoice)))
    
 
