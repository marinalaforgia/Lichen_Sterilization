write.csv(readRDS("taxa.rds"), 'taxa.csv')
write.csv(readRDS("seqtab_final.rds"), 'seqtab_final.csv')


library(dplyr)
library(stringr)

Taxa<-read.csv("taxa.csv")
Seqtab<-read.csv("seqtab_final.csv")

Counts<-Seqtab
rownames(Counts)<-Counts[,1]
Counts<-Counts[,-1]
Counts<-as.data.frame(t(Counts))
Counts$SPECIES_ID<-paste(substring(Taxa$Genus, 4), substring(Taxa$Species, 4))
Counts<-na.omit(Counts)
DedupCounts<-aggregate(.~SPECIES_ID, data=Counts, sum)
DedupCounts$GENUS_ID<-str_split_fixed(DedupCounts$SPECIES_ID, " ", 2)

FinalCounts<-DedupCounts

FinalCounts<-t(DedupCounts)

write.csv(FinalCounts, "Final_Counts_2020_Species.csv")
