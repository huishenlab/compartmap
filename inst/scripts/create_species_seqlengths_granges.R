##Create the GRanges objects for each supported species
#The goal is to create a lightweight GRanges object for seqlengths
#rather than having to keep larger packages (e.g. BSgenome) as
#dependencies. It also speeds up the build system.

#hg19
library(Homo.sapiens)

hg19.gr <- as(seqinfo(Homo.sapiens), "GRanges")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
hg19.tx.gr <- sort(genes(TxDb.Hsapiens.UCSC.hg19.knownGene))

#hg38
library(BSgenome.Hsapiens.UCSC.hg38)
hg38.gr <- as(seqinfo(BSgenome.Hsapiens.UCSC.hg38), "GRanges")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
hg38.tx.gr <- sort(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))

#mm9
library(BSgenome.Mmusculus.UCSC.mm9)
mm9.gr <- as(seqinfo(BSgenome.Mmusculus.UCSC.mm9), "GRanges")

library(TxDb.Mmusculus.UCSC.mm9.knownGene)
mm9.tx.gr <- sort(genes(TxDb.Mmusculus.UCSC.mm9.knownGene))

#mm10
library(Mus.musculus)

mm10.gr <- as(seqinfo(Mus.musculus), "GRanges")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
mm10.tx.gr <- sort(genes(TxDb.Mmusculus.UCSC.mm10.knownGene))

#save
save(hg19.gr, file = "../../data/hg19_gr.rda")
save(hg38.gr, file = "../../data/hg38_gr.rda")
save(mm9.gr, file = "../../data/mm9_gr.rda")
save(mm10.gr, file = "../../data/mm10_gr.rda")

save(hg19.tx.gr, file = "../../data/hg19_tx_gr.rda")
save(hg38.tx.gr, file = "../../data/hg38_tx_gr.rda")
save(mm9.tx.gr, file = "../../data/mm9_tx_gr.rda")
save(mm10.tx.gr, file = "../../data/mm10_tx_gr.rda")
