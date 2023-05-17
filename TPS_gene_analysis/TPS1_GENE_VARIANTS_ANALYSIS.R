
library(dplyr)
library(ggplot2)
library(trackViewer)
library(plotrix)

TPS1_GENE_VARIANTS_ANALYSIS <- function(direc,offset,tableFrame_snps,tableFrame_indel) 
{

  #=====================TPS1============================
  TPS_1_exon_start<-c(17431851,17431433,17431046,17430539,17430222,17429948,17429496)
  TPS_1_exon_end<-c(17431531,17431160,17430668,17430321,17430084,17429703,17428944)

  TPS_1_intron_start<-c(17431531,17431160,17430668,17430321,17430084,17429703)
  TPS_1_intron_end<-c(17431433,17431046,17430539,17430222,17429948,17429496)
  #=====================TPS1============================

  #=====================TPS2============================
  TPS_2_exon_start<-c(38840385,38847596,38848004,38856189,38856493)
  TPS_2_exon_end<-c(38840561,38847839,38848155,38856403,38856637)

  TPS_2_intron_start<-c(38840561,38847839,38848155,38856403)
  TPS_2_intron_end<-c(38847596,38848004,38856189,38856493)
  #=====================TPS2============================

  #=====================TPS3============================
  TPS_3_exon_start<-c(40013548,40013137,40012757,40012264,40011866,40011636,40010851)
  TPS_3_exon_end<-c(40013234,40012870,40012376,40012046,40011728,40011382,40010187)

  TPS_3_intron_start<-c(40013234,40012870,40012376,40012046,40011728,40011382)
  TPS_3_intron_end<-c(40013137,40012757,40012264,40011866,40011636,40010851)
  #=====================TPS3============================
  #=====================TPS1============================
  #locate and extract the variants
  direction=direc   #0--negative  1--positive
  chmoname="NC_039901.1"
  offsetValue=offset
  startPos=17431851
  endPos=17428944
  if(direction==0)
  {
    temp=startPos
    startPos=endPos
    endPos=temp+offsetValue

  }else{
    startPos=startPos-offsetValue
  }

  snps_tps1 <-filter(tableFrame_snps, tableFrame_snps$CHROM==chmoname &tableFrame_snps$POS <= endPos & tableFrame_snps$POS >=startPos)
  indel_tps1<-filter(tableFrame_indel, tableFrame_indel$CHROM==chmoname & tableFrame_indel$POS <= endPos & tableFrame_indel$POS >=startPos)
  
  #the position of each snps
  snps_Position_tps1<-c(snps_tps1$POS[])

  #the position of each indels
  indel_Position_tps1<-c(indel_tps1$POS[])
  # 
  # #draw graph SNP
  # SNP <- snps_Position_tps1
  # sample.gr <- GRanges("chr1", IRanges(SNP, width=1))
  # features <- GRanges("chr1", IRanges(TPS_1_exon_end,
  #                                     width=c(TPS_1_exon_start[1]-TPS_1_exon_end[1],TPS_1_exon_start[2]-TPS_1_exon_end[2],TPS_1_exon_start[3]-TPS_1_exon_end[3],TPS_1_exon_start[4]-TPS_1_exon_end[4],TPS_1_exon_start[5]-TPS_1_exon_end[5],TPS_1_exon_start[6]-TPS_1_exon_end[6],TPS_1_exon_start[7]-TPS_1_exon_end[7]),
  #                                     names=c("Exon1", "Exon2", "Exon3", "Exon4", "Exon5", "Exon6", "Exon7")))
  # features$fill <- c("#CCFFCC","#CCFF00","#00FF66","#0066FF","#CC00FF", "#FF0000","#666666")
  # sample.gr$color <- sample.int(85, length(SNP),replace=TRUE)
  # sample.gr$label <- as.character(1:length(sample.gr))
  # sample.gr$label.col <- "white"
  # lolliplot(sample.gr, features)
  
  #draw graph inDel
  # inDel <- indel_Position_tps1
  # sample.gr <- GRanges("chr1", IRanges(inDel, width=1, names=paste0("indel", inDel)))
  # features <- GRanges("chr1", IRanges(TPS_1_exon_end,
  #                                     width=c(TPS_1_exon_start[1]-TPS_1_exon_end[1],TPS_1_exon_start[2]-TPS_1_exon_end[2],TPS_1_exon_start[3]-TPS_1_exon_end[3],TPS_1_exon_start[4]-TPS_1_exon_end[4],TPS_1_exon_start[5]-TPS_1_exon_end[5],TPS_1_exon_start[6]-TPS_1_exon_end[6],TPS_1_exon_start[7]-TPS_1_exon_end[7]),
  #                                     names=c("Exon1", "Exon2", "Exon3", "Exon4", "Exon5", "Exon6", "Exon7")))
  # features$fill <- c("#CCFFCC","#CCFF00","#00FF66","#0066FF","#CC00FF", "#FF0000","#666666")
  # sample.gr$color <- sample.int(85, length(inDel),replace=TRUE)
  # sample.gr$label <- as.character(1:length(sample.gr))
  # sample.gr$label.col <- "white"
  # lolliplot(sample.gr, features)
  
  #draw Pie Charts 2D
  # x <- c(length(snps_Position_tps1), length(indel_Position_tps1))
  # labels <- c("SNPs", "Indel")
  # pie(x,labels,main="Variants distribution in TPS1",col = rainbow(length(x)))
  
  # #draw Pie Charts 3D
  # data <-  c(length(snps_Position_tps1), length(indel_Position_tps1))
  # text_lab <- c("SNPs: ", "Indel: ")
  # lab <-paste0(text_lab,round(data/sum(data) * 100, 2), "%")
  # pie3D(data,labels=lab,
  #       col = hcl.colors(length(data), "Spectral"),explode = 0.4)

  #=======================================snps============================================
  counter<-0

  index_exon<-0
  index_intron<-0

  index_array_exon<-c()
  pos_array_exon<-c()
  exon_id_array_exon<-c()
  ref_array_exon<-c()
  alt_array_exon<-c()

  index_array_intron<-c()
  pos_array_intron<-c()
  intron_id_array_intron<-c()
  ref_array_intron<-c()
  alt_array_intron<-c()

  df_exon <- data.frame()
  df_intron <- data.frame()

  variants_statisic_exon<-c()
  variants_statisic_intron<-c()

  basepair_statisic_ATCG<-c(0,0,0,0)

  for (Snp_Position in snps_Position_tps1) {

    counter<-counter+1

    indicator<-0

    i<-1
    while (i <=length(TPS_1_exon_start)) {
      if(Snp_Position<=TPS_1_exon_start[i]&&Snp_Position>=TPS_1_exon_end[i])
      {
        indicator<-1

        index_exon<-index_exon+1

        index_array_exon <- append(index_array_exon,index_exon)
        pos_array_exon <- append(pos_array_exon, Snp_Position)
        exon_id_array_exon<- append(exon_id_array_exon,i)
        ref_array_exon<-append(ref_array_exon, snps_tps1$REF[counter])
        alt_array_exon<-append(alt_array_exon, snps_tps1$ALT[counter])

        df_exon <- data.frame(index_array_exon,exon_id_array_exon,pos_array_exon,ref_array_exon,alt_array_exon)

        break

      }
      i<-i+1
    }

    if(indicator==0)
    {
      i<-1
      while (i <=length(TPS_1_intron_start)) {
        if(Snp_Position<=TPS_1_intron_start[i]&&Snp_Position>=TPS_1_intron_end[i])
        {
          index_intron<-index_intron+1

          index_array_intron <- append(index_array_intron, index_intron)
          pos_array_intron <- append(pos_array_intron, Snp_Position)
          intron_id_array_intron<- append(intron_id_array_intron,i)

          ref_array_intron<-append(ref_array_intron, snps_tps1$REF[counter])
          alt_array_intron<-append(alt_array_intron, snps_tps1$ALT[counter])
          df_intron <- data.frame(index_array_intron,intron_id_array_intron,pos_array_intron,ref_array_intron,alt_array_intron)

          break

        }
        i<-i+1
      }
    }
  }

  colnames(df_exon)<-c("index","exon_id","pos","ref","alt")
  colnames(df_intron)<-c("index","intron_id","pos","ref","alt")

  i<-1
  while(i<=7)
  {
    counter<-0
    for (pos in df_exon$exon_id)
    {
      if(i==pos)
      {
        counter<-counter+1
      }

    }

    variants_statisic_exon[i]<-counter

    i<-i+1
  }
# 
  # barplot(variants_statisic_exon,
  #         main = "SNPs in exon statistic activity for TPS1",
  #         xlab = "Exon id",
  #         ylab = "The number of SNPs",
  #         names.arg = c("Exon1", "Exon2", "Exon3", "Exon4", "Exon5", "Exon6", "Exon7"),
  #         col = "darkred",
  #         horiz = FALSE)
  
  
  #=================================draw Pie Charts===================
  #draw Pie Charts
  # labels <- c("Exon1", "Exon2", "Exon3", "Exon4", "Exon5", "Exon6", "Exon7")
  # pie(variants_statisic_exon,labels, main="SNPs distribution in exon of TPS1",col = rainbow(length(labels)))
  #draw Pie Charts 3D
  # data <-variants_statisic_exon
  # text_lab <- c("Exon1", "Exon2", "Exon3", "Exon4", "Exon5", "Exon6", "Exon7")
  # #lab <-paste0(text_lab,round(data/sum(data) * 100, 2), "%")
  # pie3D(data,labels=text_lab,
  #       col = hcl.colors(length(data), "Spectral"),explode = 0.1)
 
   #=================================draw Pie Charts===================
  i<-1
  while(i<=6)
  {
    counter<-0
    for (pos in df_intron$intron_id)
    {
      if(i==pos)
      {
        counter<-counter+1
      }

    }

    variants_statisic_intron[i]<-counter

    i<-i+1
  }

  # barplot(variants_statisic_intron,
  #         main = "SNPs in intron statistic activity for TPS1",
  #         xlab = "Intron id",
  #         ylab = "The number of SNPs",
  #         names.arg = c("Intron1", "Intron2", "Intron3", "Intron4", "Intron5", "Intron6"),
  #         col = "darkred",
  #         horiz = FALSE)
  #=================================draw Pie Charts===================
  #draw Pie Charts 2d
  # labels <- c("Intron1", "Intron2", "Intron3", "Intron4", "Intron5", "Intron6")
  # pie(variants_statisic_intron,labels, main="SNPs distribution in intron of TPS1",col = rainbow(length(labels)))
  #draw Pie Charts 3D
  # data <-variants_statisic_intron
  # labels <- c("Intron1", "Intron2", "Intron3", "Intron4", "Intron5", "Intron6")
  # pie3D(data,labels,col = hcl.colors(length(data), "Spectral"))
  #=================================draw Pie Charts===================

  for (Basepair in snps_tps1$ALT) {

    if(Basepair=="A")
    {
      basepair_statisic_ATCG[1]<-basepair_statisic_ATCG[1]+1
    }
    if(Basepair=="T")
    {
      basepair_statisic_ATCG[2]<-basepair_statisic_ATCG[2]+1
    }
    if(Basepair=="C")
    {
      basepair_statisic_ATCG[3]<-basepair_statisic_ATCG[3]+1
    }
    if(Basepair=="G")
    {
      basepair_statisic_ATCG[4]<-basepair_statisic_ATCG[4]+1
    }

  }

  # barplot(basepair_statisic_ATCG,
  #         main = "Alt base statistic activity in TPS1",
  #         xlab = "Alt base",
  #         ylab = "The number of Alt base",
  #         names.arg = c("A", "T", "C", "G"),
  #         col = "darkred",
  #         horiz = FALSE)
  #=======================================snps============================================
  #=======================================indel============================================
  counter<-0

  index_exon<-0
  index_intron<-0

  index_array_exon<-c()
  pos_array_exon<-c()
  exon_id_array_exon<-c()
  ref_array_exon<-c()
  alt_array_exon<-c()

  index_array_intron<-c()
  pos_array_intron<-c()
  intron_id_array_intron<-c()
  ref_array_intron<-c()
  alt_array_intron<-c()

  df_exon <- data.frame()
  df_intron <- data.frame()

  variants_statisic_exon<-c()
  variants_statisic_intron<-c()

  for (indel_Position in indel_Position_tps1) {

    counter<-counter+1
    indicator<-0

    i<-1
    while (i <=length(TPS_1_exon_start)) {
      if(indel_Position<=TPS_1_exon_start[i]&&indel_Position>=TPS_1_exon_end[i])
      {
        indicator<-1
        index_exon<-index_exon+1

        index_array_exon <- append(index_array_exon,index_exon)
        pos_array_exon <- append(pos_array_exon, indel_Position)
        exon_id_array_exon<- append(exon_id_array_exon,i)
        ref_array_exon<-append(ref_array_exon, indel_tps1$REF[counter])
        alt_array_exon<-append(alt_array_exon, indel_tps1$ALT[counter])

        df_exon <- data.frame(index_array_exon,exon_id_array_exon,pos_array_exon,ref_array_exon,alt_array_exon)

        break

      }
      i<-i+1
    }

    if(indicator==0)
    {
      i<-1
      while (i <=length(TPS_1_intron_start)) {
        if(indel_Position<=TPS_1_intron_start[i]&&indel_Position>=TPS_1_intron_end[i])
        {
          index_intron<-index_intron+1

          index_array_intron <- append(index_array_intron, index_intron)
          pos_array_intron <- append(pos_array_intron, indel_Position)
          intron_id_array_intron<- append(intron_id_array_intron, i)
          ref_array_intron<-append(ref_array_intron, indel_tps1$REF[counter])
          alt_array_intron<-append(alt_array_intron, indel_tps1$ALT[counter])

          df_intron <- data.frame(index_array_intron,intron_id_array_intron,pos_array_intron,ref_array_intron,alt_array_intron)

          break

        }
        i<-i+1
      }
    }
  }
  colnames(df_exon)<-c("index","exon_id","pos","ref","alt")
  colnames(df_intron)<-c("index","intron_id","pos","ref","alt")

  i<-1
  while(i<=7)
  {
    counter<-0
    for (pos in df_exon$exon_id)
    {
      if(i==pos)
      {
        counter<-counter+1
      }

    }

    variants_statisic_exon[i]<-counter

    i<-i+1
  }

    # barplot(variants_statisic_exon,
    #         main = "Indel in exon statistic activity for TPS1",
    #         xlab = "Exon id",
    #         ylab = "The number of Indel",
    #         names.arg = c("Exon1", "Exon2", "Exon3", "Exon4", "Exon5", "Exon6", "Exon7"),
    #         col = "darkred",
    #         horiz = FALSE)

  i<-1
  while(i<=6)
  {
    counter<-0
    for (pos in df_intron$intron_id)
    {
      if(i==pos)
      {
        counter<-counter+1
      }

    }

    variants_statisic_intron[i]<-counter

    i<-i+1
  }
 
  # barplot(variants_statisic_intron,
  #         main = "Indel in intron statistic activity for TPS1",
  #         xlab = "Intron id",
  #         ylab = "The number of Indel",
  #         names.arg = c("Intron1", "Intron2", "Intron3", "Intron4", "Intron5", "Intron6"),
  #         col = "darkred",
  #         horiz = FALSE)
#=======================================indel============================================
  return(0)


#=====================TPS1============================
}