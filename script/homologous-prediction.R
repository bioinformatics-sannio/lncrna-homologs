source("script/metrics.R")

load.hg38.promSeq <- function() {
  e <- load(url("https://www.dropbox.com/s/4loicpp90g4lrk2/promSeqs-hg38.Rdata?raw=1"))
  return(eval(parse(text=e)))
}

# set database genome nomenclature

ens.org = c("hsapiens_gene_ensembl","drerio_gene_ensembl", "mmusculus_gene_ensembl")
names(ens.org) = c("Human","Zebrafish","Mouse")
ucsc.org = c("hg38","danRer10","mm10")
names(ucsc.org) = names(ens.org)

species = list(HM=c("Human","Mouse"),MZ=c("Mouse","Zebrafish"),HZ=c("Human","Zebrafish"))

# define tx classes
longclass = c("lincRNA","sense_intronic","sense_overlapping","antisense")
small = c("rRNA","snRNA","miRNA","snoRNA")

datatab = NULL


for (sp in species) {
  i1=sp[1]
  i2=sp[2]
  cat(i1,i2,"\n")
  
  orthocol1 = gsub("_gene_ensembl","",ens.org[i2])
  
  # Load gene sequence data
  o1 = ucsc.org[i1]
  o2 = ucsc.org[i2]
  
  load(paste0("data/geneSeqs-",o1,".Rdata"))
  gSeqs1 = geneSeqs
  load(paste0("data/geneSeqs-",o2,".Rdata"))
  gSeqs2 = geneSeqs
  
  load(paste0("data/txSeqs-",o1,".Rdata"))
  txSeqs1 = txSeqs
  load(paste0("data/txSeqs-",o2,".Rdata"))
  txSeqs2 = txSeqs
  
  load(paste0("data/up1000Seqs-",o1,".Rdata"))
  upSeqs1 = upSeqs
  load(paste0("data/up1000Seqs-",o2,".Rdata"))
  upSeqs2 = upSeqs
  
  if(o1 == "hg38") {
    promSeqs1 = load.hg38.promSeq()
  } else {
    load(paste0("data/promSeqs-",o1,".Rdata"))
    promSeqs1 = promSeqs
  }
  load(paste0("data/promSeqs-",o2,".Rdata"))
  
  load(paste0("promSeqs-",o2,".Rdata"))
  promSeqs2 = promSeqs
  
  # Load transcript classes data
  load(paste0("data/txClass-",o1,".Rdata"))
  txClass1 = txClass
  load(paste0("data/txClass-",o2,".Rdata"))
  txClass2 = txClass
  
  # Only longs and existing seq
  txClass1 = txClass1[txClass1$GENEID %in% names(gSeqs1),]
  txClass2 = txClass2[txClass2$GENEID %in% names(gSeqs2),]
  txClass1 = txClass1[txClass1$TXNAME %in% names(txSeqs1),]
  txClass2 = txClass2[txClass2$TXNAME %in% names(txSeqs2),]
  
  # Load S1 - S2 Long coding ortho
  oracle = read.csv2(paste0("data/",i1,"-",i2,"-long-ortho.csv"),stringsAsFactors = FALSE)
  oracle = oracle[oracle$O2ENSID %in% names(gSeqs2),]
  oracle = oracle[oracle$O1ENSID %in% names(gSeqs1),]
  
  # Substitution Matrix 
  file = "data/NUC.4.4"
  mat = as.matrix(read.table(file, check.names=FALSE))
  tipo = "local"
  gapo = 10
  gape = 0.5
  
  for (seqt in c("transcript","promoter")) {
    for (gi in 1:nrow(oracle)) {
      for (gtype in c("long","pc")) {
        if (gtype == "long") {
          o1g = oracle[gi,"O1ENSID"]
          o2g = oracle[gi,"O2ENSID"]
          o1t = unlist(strsplit(oracle[gi,"O1ENSTID"],","))
          o1t = o1t[o1t %in% names(txSeqs1)]
          o2t = unlist(strsplit(oracle[gi,"O2ENSTID"],","))
          o2t = o2t[o2t %in% names(txSeqs2)]
          othergenes2 = sample(unique(txClass2[txClass2$GENEID!=o2g & txClass2$Class %in% longclass,c("GENEID")]),20)
          othergenes1 = sample(unique(txClass1[txClass1$GENEID!=o1g & txClass1$Class %in% longclass,c("GENEID")]),20)
          othertran2 = sample(unique(txClass2[txClass2$GENEID!=o2g & txClass2$Class %in% longclass,c("TXNAME")]),20)
          othertran1 = sample(unique(txClass1[txClass1$GENEID!=o1g & txClass1$Class %in% longclass,c("TXNAME")]),20)
        } else {
          pcgenes = txClass1$Class=="PCT" & !is.na(txClass1[,orthocol1]) & txClass1[,orthocol1]!="" & txClass1[,orthocol1] %in% txClass2$GENEID
          orthos = txClass1[pcgenes,c("GENEID",orthocol1)]
          srow = sample(1:nrow(orthos),1)
          o1g = orthos[srow,"GENEID"]
          o2g = orthos[srow,orthocol1]
          o1t=txClass1[txClass1$GENEID==o1g,"TXNAME"]
          o2t=txClass2[txClass2$GENEID==o2g,"TXNAME"]
          o1t = o1t[o1t %in% names(txSeqs1)]
          o2t = o2t[o2t %in% names(txSeqs2)]
          othergenes2 = sample(unique(txClass2[txClass2$GENEID!=o2g & txClass2$Class=="PCT",c("GENEID")]),20)
          othergenes1 = sample(unique(txClass1[txClass1$GENEID!=o1g & txClass1$Class=="PCT",c("GENEID")]),20)
          othertran2 = sample(unique(txClass2[txClass2$GENEID!=o2g & txClass2$Class=="PCT",c("TXNAME")]),20)
          othertran1 = sample(unique(txClass1[txClass1$GENEID!=o1g & txClass1$Class=="PCT",c("TXNAME")]),20)
        }
        cat(" ",gtype,seqt,o1g,"-",o2g,"\n")
        if (seqt == "promoter") {
          o1SeqL = subseq(promSeqs1[o1t],start=1, end=3000)
          o2SeqL = subseq(promSeqs2[o2t],start=1, end=3000)
          o1OtherSeq = subseq(promSeqs1[othertran1],start=1, end=3000)
          o2OtherSeq = subseq(promSeqs2[othertran2],start=1, end=3000)
        } else {
          o1SeqL = txSeqs1[o1t]
          o2SeqL = txSeqs2[o2t]
          o1OtherSeq = txSeqs1[othertran1]
          o2OtherSeq = txSeqs2[othertran2]
        }
        
        for(k in 1:length(o1SeqL)) {
          o1Seq = o1SeqL[[k]]
          sw12 = pairwiseAlignment(o2OtherSeq,o1Seq,substitutionMatrix=mat,type=tipo,gapOpening=gapo,gapExtension=gape,scoreOnly=TRUE)
          
          sdistm = matrix(0,nrow=length(o2OtherSeq),ncol=length(smetds))
          colnames(sdistm) = names(smetds)
          for(smet in names(smetds)) {
            sd = stringdistmatrix(as.character(o1Seq),as.character(o2OtherSeq),method=smetds[[smet]]["mth"],q=smetds[[smet]]["q"])
            sdistm[,smet]=sd
          }
          
          bbc12 = BBC(S1=o2OtherSeq,S2=o1Seq)
          jsd12 = JSD(S1=o2OtherSeq,S2=o1Seq)
          lz12= LZ(S1=o2OtherSeq,S2=o1Seq)
          acs12 = ACS(S1=o2OtherSeq,S2=o1Seq)
          
          datatab=rbind(datatab,data.frame(O1ENSID=o1g,O2ENSID=othergenes2,GENE=F,EXPID=gi,SEQ=seqt,GTYPE=gtype,
                                           O1Length=length(o1Seq),
                                           O2Length=width(o2OtherSeq),
                                           SW=sw12,BBC=bbc12, JSD=jsd12, LZ=lz12, ACS=acs12, sdistm, 
                                           Ortho=F,O1=i1,O2=i2))
        }
        for(k in 1:length(o2SeqL)) {
          o2Seq=o2SeqL[[k]]
          sw21 = pairwiseAlignment(o1OtherSeq,o2Seq,substitutionMatrix=mat,type=tipo,gapOpening=gapo,gapExtension=gape,scoreOnly=TRUE)
          bbc21 = BBC(S1=o1OtherSeq,S2=o2Seq)
          jsd21 = JSD(S1=o1OtherSeq,S2=o2Seq)
          lz21= LZ(S1=o1OtherSeq,S2=o2Seq)
          acs21 = ACS(S1=o1OtherSeq,S2=o2Seq)
          
          sw = pairwiseAlignment(o1SeqL,o2Seq,substitutionMatrix=mat,type=tipo,gapOpening=gapo,gapExtension=gape,scoreOnly=T)
          bbc = BBC(S1=o1SeqL,S2=o2Seq)
          jsd = JSD(S1=o1SeqL,S2=o2Seq)
          lz= LZ(S1=o1SeqL,S2=o2Seq)
          acs = ACS(S1=o1SeqL,S2=o2Seq)
          
          sdistm = matrix(0,nrow=length(o1OtherSeq),ncol=length(smetds))
          colnames(sdistm) = names(smetds)
          sdist = matrix(0,nrow=length(o1SeqL),ncol=length(smetds))
          colnames(sdist) = names(smetds)
          for(smet in names(smetds)) {
            sd = stringdistmatrix(as.character(o2Seq),as.character(o1OtherSeq),method=smetds[[smet]]["mth"],q=smetds[[smet]]["q"])
            sdistm[,smet]=sd
            sd = stringdistmatrix(as.character(o2Seq),as.character(o1SeqL),method=smetds[[smet]]["mth"],q=smetds[[smet]]["q"])
            sdist[,smet]=sd
          }
          
          datatab=rbind(datatab,data.frame(O1ENSID=othergenes1,O2ENSID=o2g,GENE=F,EXPID=gi,SEQ=seqt,GTYPE=gtype,
                                           O1Length=width(o1OtherSeq),
                                           O2Length=length(o2Seq),
                                           SW=sw21, BBC=bbc21, JSD=jsd21, LZ=lz21, ACS=acs21,sdistm, 
                                           Ortho=F,O1=i1,O2=i2))
          
          datatab=rbind(datatab,data.frame(O1ENSID=o1g,O2ENSID=o2g,GENE=F,EXPID=gi,SEQ=seqt,GTYPE=gtype,
                                           O1Length=width(o1SeqL),
                                           O2Length=length(o2Seq),
                                           SW=sw, BBC=bbc, JSD=jsd, LZ=lz, ACS=acs, sdist, 
                                           Ortho=T,O1=i1,O2=i2))
        }
      }
    }
  }
}