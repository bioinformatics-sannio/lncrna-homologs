
BBC = function(S1,S2){
  W.S2 = dinucleotideFrequency(S2)/(length(S2)-2 +1)
  H = NULL
  for(ss in 1:length(S1) ){
    W.S1 = dinucleotideFrequency(S1[[ss]])/(length(S1[[ss]])-2 +1)
    H = c(H,sum(abs(W.S1 - W.S2)))
  }
  return(H)
}

JSD = function(S1,S2){
  if(class(S2) != "DNAString") S2 = DNAString(S2)
  W.S2 = dinucleotideFrequency(S2,as.prob = T)
  if(class(S1) != "DNAStringSet"){
    S1 = DNAString(S1)
    S1 = DNAStringSet(S1)
  } 
  jsd = NULL
  for(ss in 1:length(S1) ){
    W.S1 = dinucleotideFrequency(S1[[ss]],as.prob = T)
    W.SM = (W.S1 + W.S2)/2 
    
    tmp1 = log(W.S1/W.SM,base = 2)
    tmp1 = tmp1[!tmp1 %in% c(NaN, -Inf,Inf)]
    KL1 = sum(W.S1[!W.S1 == 0]*tmp1)
    
    tmp2 = log(W.S2/W.SM,base = 2)
    tmp2 = tmp2[!tmp2 %in% c(NaN, -Inf,Inf)]
    KL2 = sum(W.S2[!W.S2 == 0]*tmp2)
    
    jsd = c(jsd,0.5*(KL1 +KL2))
  }
  return(jsd)
  
}

LZ = function(S1,S2){
  lz = NULL
  for(ss in 1:length(S1) ){
    tmp = DNAStringSet(c(as.character(S1[[ss]]),as.character(S2)),use.names=T)
    tmp = as(tmp, "XStringSet")
    writeXStringSet(tmp,filepath="lz.fasta",format="fasta")
    lz = c(lz,as.numeric(gsub(".* ","",system("decaf+py/compute-lz-distance.py --fasta lz.fasta",intern=T, ignore.stderr = T)[3])))
  }
  return(lz)
}

ACS = function(S1,S2){
  acs = NULL
  for(ss in 1:length(S1) ){
    tmp = DNAStringSet(c(as.character(S1[[ss]]),as.character(S2)),use.names=T)
    tmp = as(tmp, "XStringSet")
    writeXStringSet(tmp,filepath="acs.fasta",format="fasta")
    system("kmacs/kmacs acs.fasta k=0",ignore.stdout=T)
    acs = c(acs,read.table("DMat",skip=1)[2,1])
  }
  return(acs)
}

require(stringdist)
smetds = list(DLevDist=c(mth="dl",q=3))
qlist=c(3,4,5,6,7,8,9,10,11,12,13,14,15,16)
for(qq in qlist) {
  for(mm in c("qgram","cosine","jaccard")) {
    smetds[[paste0(mm,qq)]] = c(mth=mm,q=qq)
  }
}