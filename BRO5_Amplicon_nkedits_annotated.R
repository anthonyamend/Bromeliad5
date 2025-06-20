#installing/loading packages:
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2") 
library(dada2); packageVersion("dada2")
#Version 1.32.0
library(ShortRead)
packageVersion("ShortRead")
#1.62.0
library(Biostrings)
packageVersion("Biostrings")
#2.72.2
library(here)


here()

path <- "/Users/Samira/Desktop/BRO5/R_Analysis/amend_bot662_miseq_demuxed/fq_gz" # CHANGE ME to the directory containing the raw fastq files


allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
        RevComp = Biostrings::reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}


list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME...

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))


list.files(path)
#46 - with 1 PCR negative pair, 1 sample/"true" negative pair, and 1 extraction negative pair
#40 samples total
#ignoring filtN folder in the count

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
#23 samples


#First, lets look at quality profile of R1 reads
plotQualityProfile(fnFs[c(1,2,3,4)])
plotQualityProfile(fnFs[c(20,21,22,23)])
#looks great up to 200


#Then look at quality profile of R2 reads
plotQualityProfile(fnRs[c(1,2,3,4)])
plotQualityProfile(fnRs[c(20,21,22,23)])
#keeping up to 200, though quality is worse than fwd


# Make directory and filenames for the filtered fastqs
filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(200,180), #leaves overlap
                     maxN=0, #DADA does not allow Ns
                     maxEE=c(2,2), #allow 2 errors
                     truncQ=2, 
                     #trimLeft=c(18,20), #N nucleotides to remove from the start of each read
                     rm.phix=TRUE, #remove reads matching phiX genome
                     matchIDs=TRUE, #enforce matching between id-line sequence identifiers of F and R reads
                     compress=TRUE, 
                     verbose=TRUE,
                     multithread=TRUE) # On Windows set multithread=FALSE

head(out)
tail(out)

#No samples had zero reads; this section is skipped
#out.zeroes <- out[out[,2]==0,]
#out.zeroes
##zeroes: 
#going back to the beginning now & moving these fastq files away manually, I'm sure there's an easier way

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#sanity check: visualize estimated error rates
#error rates should decline with increasing qual score
#red line is based on definition of quality score alone
#black line is estimated error rate after convergence
#dots are observed error rate for each quality score


plotErrors(errF, nominalQ=TRUE) 
plotErrors(errR, nominalQ=TRUE) 

#Errors look a little weird but generally decline as quality score goes up; pattern holds

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 


dadaFs[[1]]
dadaRs[[1]]

### Merge paired reads

#To further cull spurious sequence variants
#Merge the denoised forward and reverse reads
#Paired reads that do not exactly overlap are removed
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

summary((mergers[[1]]))

#We now have a data.frame for each sample with the merged $sequence, its $abundance, and the indices of the merged $forward and $reverse denoised sequences. Paired reads that did not exactly overlap were removed by mergePairs.

seqtab <- makeSequenceTable(mergers)
dim(seqtab)
#23 samples, 9973 ASVs

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

##save raw output, I doubt I'll change it but just in case
#saveRDS(seqtab,file="bro5.rawdada2.rds")

##visual
plot(table(nchar(getSequences(seqtab)))) 

#The core dada method removes substitution and indel errors, but chimeras remain. 
#Fortunately, the accuracy of the sequences after denoising makes identifying chimeras easier 
#than it is when dealing with fuzzy OTUs: all sequences which can be exactly reconstructed as 
#a bimera (two-parent chimera) from more abundant sequences.

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#Identified 4655 bimeras out of 9973 input sequences.
dim(seqtab.nochim)
#23 samples, 5318 sequences

sum(seqtab.nochim)/sum(seqtab)
#0.864187% reads

# Track Read Stats #
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
tail(track)

#write.csv(track,file="bro5.20240325.readstats.csv",row.names=TRUE,quote=FALSE)
#write.csv(seqtab.nochim,file="bro5.20240325.seqtab.nochim.csv")
# Save the R object to disk for later use
saveRDS(seqtab.nochim,file="bro5.seqtab.nochim.rds")


#read this table back in if it's not in your environment
#seqtab.nochim <- readRDS("bro5.20240325.seqtab.nochim.rds")

library(DECIPHER); packageVersion("DECIPHER")
#version 3.0.0
library("dada2"); packageVersion("dada2")
#version 1.32.0

dna <- DNAStringSet(getSequences(seqtab.nochim))                                 # Create a DNAStringSet from the ASVs
load("/Users/Samira/GenusListsBRO20/SILVA_SSU_r138_2019.RData")                  # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=TRUE)     # use all processors; Time difference of 1666.18 secs
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")   # ranks of interest

# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
        m <- match(ranks, x$rank)
        taxa <- x$taxon[m]
        taxa[startsWith(taxa, "unclassified_")] <- NA
        taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)

##making a column for sequences
taxid2 <- cbind(taxid,rownames(taxid))
colnames(taxid2)[8] <- "sequence"

##checking that sequences line up still between counts file & taxa file
taxid2[40,8] == colnames(seqtab.nochim)[40]
#should be true
taxid2[50,8] == colnames(seqtab.nochim)[51]
#should be false

##renaming with ASV IDs which will be easier to work with
asv.ids <- paste0("ASV",sprintf('%0.4d', 1:length(colnames(seqtab.nochim))))
#making output fasta file - did once then skipping
#uniquesToFasta(seqtab.nochim, fout="bro5.seqtab.nochim.ids.fasta", ids = asv.ids, mode = "w", width = 20000)

##making a copy before I overwrite with new ids
seqtab.nochim.ids <- seqtab.nochim
colnames(seqtab.nochim.ids) <- asv.ids
rownames(taxid2) <- asv.ids

##also the sample names still have extra info
rownames(seqtab.nochim.ids) <- sub("_F_filt.fastq.gz","",rownames(seqtab.nochim.ids))

##save files
#saveRDS(taxid2,file="bro5.seqtab.nochim.taxids.rds")
#saveRDS(seqtab.nochim.ids,file="bro5.seqtab.nochim.ids.rds")

#write.csv(taxid2,file="bro5.20240327.seqtab.nochim.taxids.csv")
#write.csv(seqtab.nochim.ids,file="bro5.20240327.seqtab.nochim.ids.csv")


library(stringr)

getwd()

#read in SeqTab.NoChim.TaxIDs as a dataframe
#tax <- data.frame(readRDS("bro5.20240327.seqtab.nochim.taxids.rds"))

# if you're just powering through you can use the code below instead of reading an R object back in
#as.data.frame(taxid2) -> tax

tax$asv.id <- row.names(tax)

tax.clean <- data.frame(row.names = row.names(tax),
                        Kingdom = str_replace(tax[,1], "D_0__",""),
                        Phylum = str_replace(tax[,2], "D_1__",""),
                        Class = str_replace(tax[,3], "D_2__",""),
                        Order = str_replace(tax[,4], "D_3__",""),
                        Family = str_replace(tax[,5], "D_4__",""),
                        Genus = str_replace(tax[,6], "D_5__",""),
                        #Species = str_replace(tax[,7], "D_6__",""),
                        Sequence = c(tax[,8]),
                        ASV_id = c(tax[,9]),
                        stringsAsFactors = FALSE)

tax.clean[is.na(tax.clean)] <- ""

for (i in 1:7){ tax.clean[,i] <- as.character(tax.clean[,i])}
####### Fill holes in the tax table
tax.clean[is.na(tax.clean)] <- ""
for (i in 1:nrow(tax.clean)){
  if (tax.clean[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
    tax.clean[i, 2:6] <- kingdom
  } else if (tax.clean[i,3] == ""){
    phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
    tax.clean[i, 3:6] <- phylum
  } else if (tax.clean[i,4] == ""){
    class <- paste("Class_", tax.clean[i,3], sep = "")
    tax.clean[i, 4:6] <- class
  } else if (tax.clean[i,5] == ""){
    order <- paste("Order_", tax.clean[i,4], sep = "")
    tax.clean[i, 5:6] <- order
  } else if (tax.clean[i,6] == ""){
    family <- paste("Family_", tax.clean[i,5], sep = "")
    tax.clean[i, 6:6] <- family
  }
}

# Save taxa list to a file
#write.csv(tax.clean,file="bro5.taxafull.csv")



library(readxl)
library(phyloseq)

#read back in the taxa full list
taxa.info <- read.csv("bro5.taxafull.csv",row.names=1)

#import metadata sheet
samdf_intermediate <- read_xlsx("/Users/Samira/Desktop/BRO5/ROL_MicrobiomeBootcamp_Metadata.xlsx", sheet = "sample_metadata")
samdf_intermediate$sampleID <- gsub('_','',samdf_intermediate$sampleID)
##nk: replaced sub with 'gsub' for this

#remove the last rows because they contain notes that are not relevant for this analysis
samdf_intermediate <- samdf_intermediate[-c(26, 27, 28, 29), ]

#there were extra filtrates collected that were not sequenced for metagenomics.
#these rows were removed as well from the metadata sheet with the command below.
#the rows removed are as follows: LY-BR1-W2; LY-BR2-W1; LY-BR2-W3; LY-BR3-W2; LY-BR5-W2
samdf_intermediate <- samdf_intermediate[-c(17, 18, 20, 22, 25), ]

#add rows for the controls
samdf_intermediate <- rbind(samdf_intermediate, list('EXTNEGLY', 'N', 0, '-', '-', '-', '-', '-', '-', '-', '-'))
samdf_intermediate <- rbind(samdf_intermediate, list('NTCLY', 'N', 0, '-', '-', '-', '-', '-', '-', '-', '-'))
samdf_intermediate <- rbind(samdf_intermediate, list('PCRposLY', 'P', 0, '-', '-', '-', '-', '-', '-', '-', '-'))
##nk: changed from 'N' to 'P'

samdf <- samdf_intermediate[order(samdf_intermediate$sampleID), ]

samdf <- as.data.frame(samdf)

rownames(samdf) <- samdf$sampleID
##nk: rearranged this so warning is gone now

# Remove extraneous characters in rownames.
gsub('-', '', rownames(seqtab.nochim.ids)) -> rownames(seqtab.nochim.ids)


#phyloseq object with new taxa ids
ps <- phyloseq(otu_table(seqtab.nochim.ids, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(as.matrix(taxa.info), taxa_names()))

ps 
##23 samples, 5318 taxa
samdf.ps <- data.frame(ps@sam_data)
##which ones are missing
#seqtab.check <- seqtab.nochim.ids[!rownames(seqtab.nochim.ids) %in% rownames(samdf.ps),]
#seqtab.check


library(decontam)
library(ggplot2)


df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
# Begin creating a plot with ggplot2
ggplot(data=df, aes(x=Index, y=LibrarySize)) + geom_point()

sample_data(ps)$lib_size <- sample_sums(ps)
sample_data(ps)$is.neg <- sample_data(ps)$sample_type == "N"
contamdf.prev <- isContaminant(ps, neg="is.neg",threshold=0.5)
table(contamdf.prev$contaminant)
# FALSE  TRUE 
#  5313     5 

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_type == "N", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_type != "N", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
# Begin creating a plot with ggplot2
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#remove from ps
ps.decontam1 <- prune_taxa(!contamdf.prev$contaminant,ps)
#below removes neg controls; I'm keeping them in
#ps.decontam2 <- subset_samples(ps.decontam1,(Mesocosm_type!="Negative"))
ps.decontam1

ps.decontam2 <- prune_samples(sample_sums(ps.decontam1)!=0,ps.decontam1)
ps.decontam2

ps.decontam <- prune_taxa(taxa_sums(ps.decontam2)!=0,ps.decontam2)
ps.decontam
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 5313 taxa and 22 samples ]
# sample_data() Sample Data:       [ 22 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 5313 taxa by 8 taxonomic ranks ]

tail(taxa_sums(ps.decontam))
tail(sort(sample_sums(ps.decontam),decreasing=T))

##save save
#saveRDS(ps.decontam,file="bro5.20240701.ps.decontam.rds")
otu.table <- data.frame(ps.decontam@otu_table)
##write.csv(otu.table, "nanoseq.dec23.otutab.decontam.csv")



#first, read in ASV table
#install.packages("remotes")
library("remotes")
#install_github("https://github.com/tobiasgf/lulu.git")
library("lulu")

#And match list
matchList <- read.table("match_list.txt")
head(matchList)

#Reformat ASV table to desired LULU format
ASVs <- data.frame(t(otu.table),check.names=FALSE)

#Now, run the LULU curation
##attempting to match CMAIKI steps
curated_result <- lulu(ASVs,matchList, minimum_ratio_type = "min", minimum_ratio = 1, minimum_match = 97, minimum_relative_cooccurence = 0.95)

summary(curated_result)
##nk: 568 discarded units

#Pull out the curated OTU list, re-transpose
seqtab.decontam.lulu <- data.frame(t(curated_result$curated_table),check.names=FALSE)

##fix the lulu order
#seqtab.nochim.lulu <- select(seqtab.nochim.lulu1, num_range("ASV", 0:2284))

#Continue on to your favorite analysis
#write.csv(seqtab.decontam.lulu,file="bro5.otutab.decontam.lulu.csv")

##back into phyloseq object
ps.decontam.lulu <- ps.decontam
ps.decontam.lulu@otu_table <- otu_table(seqtab.decontam.lulu,taxa_are_rows=FALSE)
ps.decontam.lulu
#saveRDS(ps.decontam.lulu,file="bro5.ps.decontam.lulu.rds")
#down to 4745 taxa, from 5313 taxa after decontam
##nk: 4745 taxa



ps.decontam.lulu.trim <- filter_taxa(ps.decontam.lulu, function (x) {sum(x > 0) > 1}, prune=TRUE)

##any 0 read samples?
sample_sums(ps.decontam.lulu.trim)==0
ps.decontam.lulu.trim.no0 <- prune_taxa(taxa_sums(ps.decontam.lulu.trim)!=0,ps.decontam.lulu.trim)
ps.decontam.lulu.trim.no0 #1995 ASVs, 22 samples
sample_sums(ps.decontam.lulu.trim.no0)==0

##save
#saveRDS(ps.decontam.lulu.trim.no0,file="bro5.ps.decontam.lulu.trim.rds")


tax.pre <- data.frame(ps.decontam.lulu.trim.no0@tax_table)
##to examine the stuff below manually

ps.mito <- subset_taxa(ps.decontam.lulu.trim.no0, Family=="Mitochondria")
ps.mito #6 taxa to remove
ps.chlor <- subset_taxa(ps.decontam.lulu.trim.no0, Order=="Chloroplast")
ps.chlor #6 taxa to remove
ps.notbact <- subset_taxa(ps.decontam.lulu.trim.no0, Kingdom!="Bacteria")
ps.notbact #341 taxa to remove

ps.nomito <- subset_taxa(ps.decontam.lulu.trim.no0, Family!="Mitochondria")
ps.nomito #1989 taxa
ps.nochlor <- subset_taxa(ps.nomito, Order!="Chloroplast")
ps.nochlor #1983 taxa
ps.clean.int <- subset_taxa(ps.nochlor, Kingdom=="Bacteria")
ps.clean.int #1642 taxa

#just archaea
ps.arch <- subset_taxa(ps.decontam.lulu.trim.no0, Kingdom=="Archaea")
ps.arch #21 taxa

#bye negatives & positive
ps.clean <- subset_samples(ps.clean.int,sample_type!="P"&sample_type!="N")
ps.clean
ps.clean.no0 <- prune_taxa(taxa_sums(ps.clean)!=0,ps.clean)
ps.clean.no0
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1642 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 1642 taxa by 8 taxonomic ranks ]

##final save
#saveRDS(ps.clean,"bro5.ps.decontam.lulu.trim.clean.rds")



ps.rel <- transform_sample_counts(ps.clean.no0, function(x) x / sum(x))

library(trekcolors)

plot_bar(ps.rel,fill="Genus")+
  theme(legend.position="none") +
  scale_fill_trek(palette = "enara")
  


library(DECIPHER)
library(dada2)
library(tidyverse)


as.data.frame(taxid) -> taxtab # the taxonomy table
t(seqtab.nochim) -> asvtab # the transposed asv count table



# Host sequences get in the way of alignment. Remove them.

nohost = (taxtab
          %>% mutate(Seq = rownames(.))
	%>% dplyr::filter(!is.na(phylum), 
	           domain %in% c('Bacteria', 'Archaea')))
# Get the sequences on their own for clustering
seqs = as.character(nohost$Seq)



# Start by making a multisequence alignment of the ASVs
aln = AlignSeqs(DNAStringSet(seqs))

# It's not a bad idea to look at the alignment and make sure it's not too bad
BrowseSeqs(aln)

# Create a distance matrix 
dmat = DistanceMatrix(aln, type = 'dist')

# I changed the line of code below to be able to use the latest version of DECIPHER.
# I want to cluster at 97%, so I set the cutoff 0.03 (1 - 0.97).
# Don't run this below; use TreeLine.
# clsts = Clusterize(DNAStringSet(seqs), cutoff = 0.03)
# head(clsts)

# I'm using the TreeLine output because that's what the authors recommend. The 'Clusterize' function is included above.
TreeLine(aln, dmat, cutoff = 0.03, showPlot = FALSE, type = "clusters", method = "complete", processors = 3) -> clsts
head(clsts)

# Each sequence has its own row, in the order they were originally listed in.
# The cluster number that each sequence belongs to is listed in its row in the "cluster" column

# Count the clusters
n_distinct(clsts$cluster)
#2467 distinct clusters identified


# Count the number of clusters of each size
OTUcounts <- clsts %>% dplyr::count(cluster, name = 'size') 
OTUcounts %>% dplyr::count(size, name = 'count') -> binsizes

# Add the sequences to the data frame so we can make consensus sequences
clsts$seqs = seqs
head(clsts)



# A function to convert a set of sequences into a consensus sequence
get_conseq = function(x){
	if (length(x) > 1){
		aln = AlignSeqs(DNAStringSet(x))
	} else {
		aln = DNAStringSet(x)
	}
	con = ConsensusSequence(aln)
	return(as.character(con))
}
conseq = (clsts
		  %>% group_by(cluster)
		  %>% summarize(consensus = get_conseq(seqs)))
head(conseq)



# For the refFasta, make sure that you are using the right file path for your own re-running.  
otu97.tax = assignTaxonomy(conseq$consensus, 
                        refFasta = "silva_nr99_v138.1_wSpecies_train_set.fa.gz",
                        tryRC = TRUE, multithread = FALSE)
head(otu97.tax)

#saveRDS(otu97.tax, file = "otu97_taxonomy.rds")



# Turn the ASV table into a data frame and add the sequences as a column

# Create data frame first.
asvdf <- as.data.frame((asvtab))
# Reassign the column names
colnames(asvdf) <- rownames(seqtab.nochim.ids)

# Add the sequences as a column
asvdf <- asvdf %>% mutate(seqs = rownames(.)) %>% select(seqs, everything())

# Reassign the row names to avoid confusion
rownames(asvdf) <- NULL

# Join with the cluster table
asvdf = (asvdf
		  %>% left_join(clsts, by = 'seqs')
		  %>% dplyr::select(cluster, everything(), -seqs))

# Create the OTU counts table
clstdf = (asvdf
		%>% reframe(across(where(is.numeric), sum), .by = cluster)
		%>% filter(!is.na(cluster))
		%>% data.frame())
dim(clstdf)

# Check to see that the OTU cluster assignments are in the same order in the Consensus Sequence data and the Cluster Dataframe data.
if(all(conseq$cluster == clstdf$cluster) == FALSE) {
  clstdf[order(clstdf$cluster),] -> clstdf
} else {
  print("everything looks good")
}

# Add the cluster consensus sequences as the rownames of the new table
rownames(clstdf) = conseq$consensus

# Make it a matrix just like dada2 produces
otu97.counts = as.matrix(select(clstdf, -cluster))

#saveRDS(otu97.counts, file = "otu97_counts.rds")




library(stringr)

#read in SeqTab.NoChim.TaxIDs as a dataframe
#tax.otu <- data.frame(readRDS("otu97_taxonomy.rds"))


# if you're just powering through you can use the code below instead of reading an R object back in
#as.data.frame(otu97.tax) -> tax.otu

# Count the number of times a null value appears in a column by first making a data frame.
nullcounts <- c(sum(is.na(tax.otu$Genus)), sum(is.na(tax.otu$Family)),sum(is.na(tax.otu$Order)), sum(is.na(tax.otu$Phylum)))
nullnames <- c("Genus", "Family", "Order", "Phylum")
nulldf <- data.frame(nullnames, nullcounts)
# Adds a column that gives you the percent of null cells in the dataframe
nulldf$pct <- (nulldf$nullcounts / nrow(tax.otu)) * 100

# Assign the sequences to a different column
tax.otu$seq <- row.names(tax.otu)

# Remove the species column; my OTUs are clustered at 97%, a threshold better suited for genus-level identification. Species names are not informative.
tax.otu <- tax.otu[, -which(colnames(tax.otu) == "Species")]

tax.otu[is.na(tax.otu)] <- ""

for (i in 1:7){ tax.otu[,i] <- as.character(tax.otu[,i])}
####### Fill holes in the tax table
tax.otu[is.na(tax.otu)] <- ""
for (i in 1:nrow(tax.otu)){
  if (tax.otu[i,2] == ""){
    kingdom <- paste("Kingdom_", tax.otu[i,1], sep = "")
    tax.otu[i, 2:6] <- kingdom
  } else if (tax.otu[i,3] == ""){
    phylum <- paste("Phylum_", tax.otu[i,2], sep = "")
    tax.otu[i, 3:6] <- phylum
  } else if (tax.otu[i,4] == ""){
    class <- paste("Class_", tax.otu[i,3], sep = "")
    tax.otu[i, 4:6] <- class
  } else if (tax.otu[i,5] == ""){
    order <- paste("Order_", tax.otu[i,4], sep = "")
    tax.otu[i, 5:6] <- order
  } else if (tax.otu[i,6] == ""){
    family <- paste("Family_", tax.otu[i,5], sep = "")
    tax.otu[i, 6:6] <- family
  }
}
# remove the i in the loop
rm(i)

rownames(tax.otu) <- paste0("OTU",sprintf('%0.4d', 1:length(rownames(tax.otu))))


# Save taxa list to a file
# write.csv(tax.otu,file="bro5.otu_taxfull.csv")




#read back in the taxa full list
#otu.info <- read.csv("bro5.otu_taxfull.csv",row.names=1)

#otu97.counts <- data.frame(readRDS("otu97_counts.rds"))

#import metadata sheet
samdf_intermediate <- read_xlsx("/Users/Samira/Desktop/BRO5/ROL_MicrobiomeBootcamp_Metadata.xlsx", sheet = "sample_metadata")
samdf_intermediate$sampleID <- gsub('_','',samdf_intermediate$sampleID)
##nk: replaced sub with 'gsub' for this

#remove the last rows because they contain notes that are not relevant for this analysis
samdf_intermediate <- samdf_intermediate[-c(26, 27, 28, 29), ]

#there were extra filtrates collected that were not sequenced for metagenomics.
#these rows were removed as well from the metadata sheet with the command below.
#the rows removed are as follows: LY-BR1-W2; LY-BR2-W1; LY-BR2-W3; LY-BR3-W2; LY-BR5-W2
samdf_intermediate <- samdf_intermediate[-c(17, 18, 20, 22, 25), ]

#add rows for the controls
samdf_intermediate <- rbind(samdf_intermediate, list('EXTNEGLY', 'N', 0, '-', '-', '-', '-', '-', '-', '-', '-'))
samdf_intermediate <- rbind(samdf_intermediate, list('NTCLY', 'N', 0, '-', '-', '-', '-', '-', '-', '-', '-'))
samdf_intermediate <- rbind(samdf_intermediate, list('PCRposLY', 'P', 0, '-', '-', '-', '-', '-', '-', '-', '-'))
##nk: changed from 'N' to 'P'

samdf <- samdf_intermediate[order(samdf_intermediate$sampleID), ]

samdf <- as.data.frame(samdf)

rownames(samdf) <- samdf$sampleID
##nk: rearranged this so warning is gone now

rownames(otu97.counts) <- rownames(otu.info)

as.data.frame(t(otu97.counts)) -> otu97.counts


#phyloseq object with new taxa ids
ps <- phyloseq(otu_table(otu97.counts, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(as.matrix(otu.info), taxa_names()))

ps 
##23 samples, 2467 OTUs
samdf.ps <- data.frame(ps@sam_data)
##which ones are missing
seqtab.check <- otu97.counts[!colnames(otu97.counts) %in% rownames(samdf.ps),]; seqtab.check




df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
# Begin creating a plot with ggplot2
ggplot(data=df, aes(x=Index, y=LibrarySize)) + geom_point()

sample_data(ps)$lib_size <- sample_sums(ps)
sample_data(ps)$is.neg <- sample_data(ps)$sample_type == "N"
contamdf.prev <- isContaminant(ps, neg="is.neg",threshold=0.5)
table(contamdf.prev$contaminant)
# FALSE  TRUE 
#  2461     6 

# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$sample_type == "N", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$sample_type != "N", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)
# Begin creating a plot with ggplot2
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

#remove from ps
ps.decontam1 <- prune_taxa(!contamdf.prev$contaminant,ps)
#below removes neg controls; I'm keeping them in
#ps.decontam2 <- subset_samples(ps.decontam1,(Mesocosm_type!="Negative"))
ps.decontam1

ps.decontam2 <- prune_samples(sample_sums(ps.decontam1)!=0,ps.decontam1)
ps.decontam2

ps.decontam <- prune_taxa(taxa_sums(ps.decontam2)!=0,ps.decontam2)
ps.decontam
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 2461 taxa and 22 samples ]
# sample_data() Sample Data:       [ 22 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 2461 taxa by 7 taxonomic ranks ]

tail(taxa_sums(ps.decontam))
tail(sort(sample_sums(ps.decontam),decreasing=T))

##save save
#saveRDS(ps.decontam,file="bro5.otus.ps.decontam.rds")
otu.tab <- data.frame(ps.decontam@otu_table)
#write.csv(otu.tab, "nanoseq.clust97.otutab.decontam.csv")




ps.decontam.trim <- filter_taxa(ps.decontam, function (x) {sum(x > 0) > 1}, prune=TRUE)

##any 0 read samples?
sample_sums(ps.decontam.trim)==0
ps.decontam.trim.no0 <- prune_taxa(taxa_sums(ps.decontam.trim)!=0,ps.decontam.trim)
ps.decontam.trim.no0 #1393 OTUs, 22 samples
sample_sums(ps.decontam.trim.no0)==0

##save
#saveRDS(ps.decontam.trim.no0,file="bro5otus.ps.decontam.trim.rds")




tax.pre <- data.frame(ps.decontam.trim.no0@tax_table)
##to examine the stuff below manually

ps.mito <- subset_taxa(ps.decontam.trim.no0, Family=="Mitochondria")
ps.mito #5 taxa to remove
ps.chlor <- subset_taxa(ps.decontam.trim.no0, Order=="Chloroplast")
ps.chlor #3 taxa to remove
ps.notbact <- subset_taxa(ps.decontam.trim.no0, Kingdom!="Bacteria")
ps.notbact #15 taxa to remove

ps.nomito <- subset_taxa(ps.decontam.trim.no0, Family!="Mitochondria")
ps.nomito #1388 taxa
ps.nochlor <- subset_taxa(ps.nomito, Order!="Chloroplast")
ps.nochlor #1385 taxa
ps.clean.int <- subset_taxa(ps.nochlor, Kingdom=="Bacteria")
ps.clean.int #1370 taxa

#just archaea
ps.arch <- subset_taxa(ps.decontam.trim.no0, Kingdom=="Archaea")
ps.arch #15 taxa

#bye negatives & positive
ps.clean <- subset_samples(ps.clean.int,sample_type!="P"&sample_type!="N")
ps.clean
ps.clean.no0 <- prune_taxa(taxa_sums(ps.clean)!=0,ps.clean)
ps.clean.no0
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 1370 taxa and 20 samples ]
# sample_data() Sample Data:       [ 20 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 1370 taxa by 7 taxonomic ranks ]

##final save
#saveRDS(ps.clean,"bro5.otus.rds")
