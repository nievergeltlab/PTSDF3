

 zcat ../eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz | awk '{if (NR==1) {print "CHR","BP","P","SNP"}; if (NR>1 && $9<=0.05) print $1,$2,$9,$3}' | sed 's/X/23/g' > eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.mh
 zcat ../trans_ptsd_pcs_v4_aug3_2021_1.tbl.allchr.filtered.gz | awk '{if (NR==1) {print "CHR","BP","P","SNP"}; if (NR>1 && $12<=0.05) print $1,$2,$12,$3}' | sed 's/X/23/g' > trans_ptsd_pcs_v4_aug3_2021_1.tbl.allchr.filtered.gz.mh
 zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/aam_ptsd_pcs_v5_jan4_2022.allchr.fuma.gz | awk '{if (NR==1) {print "CHR","BP","P","SNP"}; if (NR>1 && $9<=0.05) print $1,$2,$9,$3}' | sed 's/X/23/g' > aam_ptsd_pcs_v5_jan4_2022.allchr.fuma.gz.mh
 zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/hna_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz | awk '{if (NR==1) {print "CHR","BP","P","SNP"}; if (NR>1 && $9<=0.05) print $1,$2,$9,$3}' | sed 's/X/23/g' > hna_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.mh
 
 #Het analysis: merge with main results..
 zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz | awk '{print $3,$1,$2}' | LC_ALL=C sort -k1b,1 > /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.snpset
 zcat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/Reur_ptsd_pcs_v4_aug3_2021.allchr.het.fuma.gz | LC_ALL=C sort -k1b,1 > /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/Reur_ptsd_pcs_v4_aug3_2021.allchr.het.fuma.gz.temp
 
 LC_ALL=C join /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.snpset /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/Reur_ptsd_pcs_v4_aug3_2021.allchr.het.fuma.gz.temp  | sort -g -k 17 > /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/Reur_ptsd_pcs_v4_aug3_2021.allchr.het.fuma.gz.temp2
 
 cat /mnt/ukbb/adam/ptsd/freeze3/freeze3_sumstats/Reur_ptsd_pcs_v4_aug3_2021.allchr.het.fuma.gz.temp2 | awk '{if (NR==1) {print "CHR","BP","P","SNP"}; if (NR>1 && $9<=0.05) print $2,$3,$17,$1}' | sed 's/X/23/g' > Reur_ptsd_pcs_v4_aug3_2021.allchr.het.fuma.gz.mh 
 
 
 blue=rgb(153,191,254,max=255)
  blue=rgb(153,191,254,max=255)
 
 
  blue2=rgb(0,0,255,max=255)
   blue3=rgb(0,0,205,max=255)
 
 
 red=rgb(230,185,184,max=255)
 green=rgb(147,205, 221,max=255)
 purple=rgb(179,162,199,max=255)

beauty_red <- rgb(237,33,36,maxColorValue=255)
adam_blue2 <-  rgb(60,160,255,maxColorValue=255)
adam_blue <- rgb(41,153,208,maxColorValue=255)
#read file


 library(data.table)

 trans <- fread('trans_ptsd_pcs_v4_aug3_2021_1.tbl.allchr.filtered.gz.mh', header=T,data.table=F)
 eur <- fread('eur_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.mh',data.table=F)
 eurhet <- fread('Reur_ptsd_pcs_v4_aug3_2021.allchr.het.fuma.gz.mh',data.table=F)
 aam  <-  fread('aam_ptsd_pcs_v5_jan4_2022.allchr.fuma.gz.mh',data.table=F)
 lat  <-  fread('hna_ptsd_pcs_v4_aug3_2021.allchr.fuma.gz.mh',data.table=F)
 
#Load manhattanplot script
 print("Plotting data")
 source('ManhattanPlotterFunction_colorfixed_max10ylim2_pgc_v2.R')
# source('ManhattanPlotterFunction_colorfixed_max10ylim2_pgc_v2_highlight.R') #highlighted, beauty red built in..


colorchoice="adam_blue2"
 highlight_p <- 5e-8
 highlightboundary <- 20000
 #replot parametert puts manhanttan over manhattan
jpeg("transplot2.jpg",width = 16 , height = 9, units = "in", res = 500,quality = 100)

    ManhattanPlot_AJS_cut(trans$CHR, trans$BP, trans$P, trans$SNP, genomesig = 5*(10^-8), pchF=1,genomesug = 0, replot=0,photoname = '', 
  outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary)
    
     ManhattanPlot_AJS_cut(eur$CHR, eur$BP, eur$P, eur$SNP, genomesig = 5*(10^-8), pchF=16,genomesug = 0,replot=1,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = rgb(127,253,166,maxColorValue=255),highlight_p=highlight_p, highlightboundary=highlightboundary)
     
    dev.off()


jpeg("eurplot2.jpg",width = 16 , height = 9, units = "in", res = 500,quality = 100)


     ManhattanPlot_AJS_cut(eur$CHR, eur$BP, eur$P, eur$SNP, genomesig = 5*(10^-8), pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = rgb(127,253,166,maxColorValue=255),highlight_p=highlight_p, highlightboundary=highlightboundary)
     
    dev.off()


colorchoice="adam_blue2"
jpeg("eurhet.jpg",width = 16 , height = 9, units = "in", res = 500,quality = 100)


     ManhattanPlot_AJS_cut(eurhet$CHR, eurhet$BP, eurhet$P, eurhet$SNP, genomesig = 5*(10^-8), pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = rgb(127,253,166,maxColorValue=255),highlight_p=highlight_p, highlightboundary=highlightboundary)
     
    dev.off()


jpeg("aamplot.jpg",width = 16 , height = 9, units = "in", res = 500,quality = 100)

colorchoice="beauty_red"
     ManhattanPlot_AJS_cut(aam$CHR, aam$BP, aam$P, aam$SNP, genomesig = 5*(10^-8), pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = rgb(127,253,166,maxColorValue=255),highlight_p=highlight_p, highlightboundary=highlightboundary)
     
    dev.off()
    
    


jpeg("latplot.jpg",width = 16 , height = 9, units = "in", res = 500,quality = 100)

colorchoice="purple"
     ManhattanPlot_AJS_cut(lat$CHR, lat$BP, lat$P, lat$SNP, genomesig = 5*(10^-8), pchF=16,genomesug = 0,replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),eval(parse(text = colorchoice))), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = rgb(127,253,166,maxColorValue=255),highlight_p=highlight_p, highlightboundary=highlightboundary)
     
    dev.off()
    
    

 source('ManhattanPlotterFunction_colorfixed_max10ylim2_pgc_v2.R')

colorchoice="blue2"
 highlight_p <- 5e-8
 highlightboundary <- 20000
 #replot parametert puts manhanttan over manhattan
jpeg("doubleplot3.jpg",width = 16 , height = 9, units = "in", res = 500,quality = 100)

    ManhattanPlot_AJS_cut(trans$CHR, trans$BP, trans$P, trans$SNP, genomesig = 5*(10^-8), pchF=15,genomesug = 0, replot=0,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,130,max=255), rgb(0,0,255,125,max=255)) , sigcolor = 'darkred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = "red",highlight_p=highlight_p, highlightboundary=highlightboundary)
    
     ManhattanPlot_AJS_cut(eur$CHR, eur$BP, eur$P, eur$SNP, genomesig = 5*(10^-8), pchF=16,genomesug = 0,replot=1,photoname = '', 
    outname = outfile, colors = c(rgb(25,25,25,max=255),rgb(60,160,255,maxColorValue=255)), sigcolor = 'indianred', sugcolor = 'indianred', ncex = 1, 
    sugcex = 1, sigcex = 1, nonautosome = c(23,24,25,26),xlabel = 'Chromosomal Positions',ylabel = '-log10(p-value)', 
    pdf = 'FALSE',topsnplabels = 'FALSE', pvalue_miss = 'NA', sigsnpcolor = rgb(127,253,166,maxColorValue=255),highlight_p=highlight_p, highlightboundary=highlightboundary)
     
    dev.off()
