#install.packages("stringdist")
library(stringdist)


#1#gene sequence
seq1="CCACGCGTCCGCCCACGCGTCCGGGAGGCACTAGGGATGGTCCGCAGGATTGGACTGATACAGAGGCCGCCACGGAGCCCGCCGGAGCCACCGTTCCTGCTGCTGCCGCCGCTGCCCGAATCGGAACCGTCGGGCCGCAGCCGCCGGCAATGCCGCGAAGGAAGAGGAATGCAGGCAGTAGTTCAGATGGAACCGAAGATTCCGATTTTCTACAGA"

#2#compliment of the seq
complement_DNA<-function(DNA_seq){
m=strsplit(DNA_seq,"")
m=m[[1]]
c_DNA_seq=unname(sapply(m, switch, "A"="T", "T"="A","G"="C","C"="G"))
c_DNA_seq=paste(c_DNA_seq, collapse="")
c_DNA_seq
}
c_seq1=complement_DNA(seq1)


#3#site enzimes
site_mut=list(HindIII = "AAGCTT",
EcoRI = "GAATTC",
Asp718I = "GGATACC",
KpnI = "GGTACC",
BamHI ="GGATCC" ,
EcoRV = "GATATC",
NotI = "CGGCCG" ,
XhoI = "CTCGAG" ,
XbaI ="TCTAGA")


#reverse enzimes
reverse_DNA<- function(DNA_seq){
  bp_split=strsplit(DNA_seq,"")
  rev_site_split=rev(bp_split[[1]])
  rev_site<- paste(rev_site_split, collapse="")
  rev_site
}

rev_site_mut=lapply(site_mut,reverse_DNA)



#DNA seq before start and stop
start= "ATG"
start_pos=regexpr("ATG", seq1)[1]
start_coding=substr(seq1, start_pos, nchar(seq1))
begining_seq=substr(seq1, 1, start_pos)
stop_pos=regexpr("TAG", start_coding)[1]
final_seq=substr(start_coding, stop_pos, nchar(start_coding))


#check for exact matches
#1st case: exact match
lapply(site_mut, regexpr,begining_seq)
lapply(site_mut, regexpr,final_seq)
#2nd case: exact reverse match
lapply(lapply(site_mut,reverse_DNA), regexpr,begining_seq)
lapply(lapply(site_mut,reverse_DNA), regexpr,final_seq)
#3rd case: complementary chain exact match
lapply(site_mut, regexpr,complement_DNA(begining_seq))
lapply(site_mut, regexpr,complement_DNA(final_seq))
#4th case: complementary chain  reverse match
lapply(lapply(site_mut,reverse_DNA), regexpr,complement_DNA(begining_seq))
lapply(lapply(site_mut,reverse_DNA), regexpr,complement_DNA(final_seq))


#levenshtain distance search FUNCTION START

test<- function(type){
  if(type=="yes"){
    
    "success"
  }else{
    "no success"
  }
  
}


find_DNA_match<- function (enzime, DNA_seq, distance, type){

if (type=="start"){
  combs=list()
  m=nchar(DNA_seq)-nchar(enzime)
for (i in 1:m){
  #first for amatch, create 1:6 combinations from begining to end-5
  sequence<-substr(DNA_seq,i, i+5)
  #same direction enzime, DNA strand
  if (!is.na(amatch(enzime, sequence, method='osa',maxDist=distance))){
    chk_mult=i-nchar(DNA_seq)
    if(chk_mult%%3==0){
  combs[[paste0("E_S", i)]]<-list("position"= seq(i,i+5), "sequence"=sequence, "enzime"=enzime, "type"= "same direction enzime, DNA strand")
    }
  }
  
  }
  
  #reverse direction enzime, DNA strand
  if (!is.na(amatch(reverse_DNA(enzime), sequence, method='osa',maxDist=distance))){
    chk_mult=i-nchar(DNA_seq)
    if(chk_mult%%3==0){
    combs[[paste0("RE_S", i)]]<-list("position"= seq(i,i+5), "sequence"=sequence, "enzime"=reverse_DNA(enzime), "type"= "reverse direction enzime, DNA strand")
    }
    }
  
  #same direction enzime, complementary DNA strand
  if (!is.na(amatch(enzime, complement_DNA(sequence), method='osa',maxDist=distance))){
    chk_mult=i-nchar(DNA_seq)
    if(chk_mult%%3==0){
    combs[[paste0("E_CS", i)]]<-list("position"= seq(i,i+5), "sequence"= complement_DNA(sequence), "enzime"=enzime, "type"= "same direction enzime, complement DNA strand")
    }
    }
  # reverse direction enzime, complement DNA strand
  if (!is.na(amatch(reverse_DNA(enzime), complement_DNA(sequence), method='osa',maxDist=distance))){
    chk_mult=i-nchar(DNA_seq)
    if(chk_mult%%3==0){
    
    combs[[paste0("RE_CS", i)]]<-list("position"= seq(i,i+5), "sequence"= complement_DNA(sequence), "enzime"=reverse_DNA(enzime), "type"= "reverse direction enzime, complement DNA strand")
    }
    
}

}



if (type=="final"){
#dna search function for the post-stop codon
#it HAS to cut at the 1 position
  combs=list()
  m=nchar(DNA_seq)-nchar(enzime)
  for (i in 1){
    #first for amatch, create 1:6 combinations from begining to end-5
    sequence<-substr(DNA_seq,i, i+5)
    #same direction enzime, DNA strand
    if (!is.na(amatch(enzime, sequence, method='osa',maxDist=distance))){
      #chk_mult=i-nchar(DNA_seq)
     # if(chk_mult%%3==0){
        combs[[paste0("E_S", i)]]<-list("position"= seq(i,i+5), "sequence"=sequence, "enzime"=enzime, "type"= "same direction enzime, DNA strand")
      #}
    }
    
  }
  
  #reverse direction enzime, DNA strand
  if (!is.na(amatch(reverse_DNA(enzime), sequence, method='osa',maxDist=distance))){
   # chk_mult=i-nchar(DNA_seq)
    #if(chk_mult%%3==0){
      combs[[paste0("RE_S", i)]]<-list("position"= seq(i,i+5), "sequence"=sequence, "enzime"=reverse_DNA(enzime), "type"= "reverse direction enzime, DNA strand")
    #}
  }
  
  #same direction enzime, complementary DNA strand
  if (!is.na(amatch(enzime, complement_DNA(sequence), method='osa',maxDist=distance))){
    #chk_mult=i-nchar(DNA_seq)
    #if(chk_mult%%3==0){
      combs[[paste0("E_CS", i)]]<-list("position"= seq(i,i+5), "sequence"= complement_DNA(sequence), "enzime"=enzime, "type"= "same direction enzime, complement DNA strand")
    #}
  }
  # reverse direction enzime, complement DNA strand
  if (!is.na(amatch(reverse_DNA(enzime), complement_DNA(sequence), method='osa',maxDist=distance))){
    #chk_mult=i-nchar(DNA_seq)
    #if(chk_mult%%3==0){
      
      combs[[paste0("RE_CS", i)]]<-list("position"= seq(i,i+5), "sequence"= complement_DNA(sequence), "enzime"=reverse_DNA(enzime), "type"= "reverse direction enzime, complement DNA strand")
    #}
    
  }
  
}
  combs

}

#Now, search in the begining and final sequences: begining_seq,final_seq
#find_DNA_match (enzime, DNA_seq, distance)
#for the begining seq, find the closest to the end
nchar(begining_seq)
str(lapply(site_mut, find_DNA_match, begining_seq, 1, "start"))
str(lapply(site_mut, find_DNA_match, begining_seq, 2, "start"))
#str(lapply(site_mut, find_DNA_match, begining_seq, 3))

#for the final seq, fin the closest to the begining
#find the closes to 1

str(lapply(site_mut, find_DNA_match, final_seq, 1, "final"))
str(lapply(site_mut, find_DNA_match, final_seq, 2, "final"))
#str(lapply(site_mut, find_DNA_match, final_seq, 3))






#in the begining seq, which is the closest to the end?
#in the final seq, which is the closest to the begining?


