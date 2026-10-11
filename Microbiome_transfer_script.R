##########################################
## MICROBIOME TRANSFER IN A POISON FROG ##
##           ANGIOLANI ET AL            ##
##               SCRIPT                 ##
##########################################

# LOAD PACKAGES ----
library(mia)
library(tidyverse)
library(ggplot2)
library(vegan)
library(veganEx)
library(phyloseq)
library(decontam)
library(microViz)
library(metagMisc)
library(microbiome) 
library(MiscMetabar)
library(rstatix)
library(ggVennDiagram)
library(cowplot)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(dplyr)


# 1.- Rarefaction and coverage for sample cleaning ----
otufileraw <- read.table("full_otu_table_transposed.txt",header = T)

#Removing duplicate samples (samples with lower seq number) 
otufile1 <- otufileraw %>% 
  filter_out(SampleID %in% c("A011a","A045a","A047a","A050a","A055b",
                                    "A058a","A060a","A061a","A062a","A063a",
                                    "A065a","A066a","A067a","A068a","A077a",
                                    "A086b","S009a","S034a","S037a","S044b",
                                    "S064a","S067a","S070b","S083a","S085a",
                                    "S089a","S092a","S098a","S121a","S124a",
                                    "S128a","S130a","S139a","S177a","S180a",
                                    "S187a","S188a","S220a","S227a","S247a",
                                    "S252a","S253a","S265a","S271a","S289a",
                                    "S290a","S298a","S319a","S326b","S328a",
                                    "S330a","S337a","S343a","S345a","S354b",
                                    "S359a","S370a","S379a","S380a","S384a",
                                    "S388a","S402a","S405a","S414a","S415a",
                                    "S423a","S426a","S429a","S433a","S458a",
                                    "S459a","S461b","S472a","S473a","S477b",
                                    "S464a","S486a","S503a","S505a"))

#Analisis of sequence counts per sample
counts <- otufile1 %>% 
  pivot_longer(-SampleID, names_to = "otu", values_to= "values")

coverage <- counts%>%
  group_by(SampleID)%>% 
  summarise(n_seqs = sum(values))
write.table(coverage, "coverage.txt", sep="\t")

#Selecting a rarefaction threshold based on coverage percentage and Goods coverage
coverage%>%
  ggplot(aes(x = n_seqs))+
  geom_histogram(binwidth = 500)+
  coord_cartesian(xlim =c(0, 15000))

coverage%>%
  ggplot(aes(x=1, y= n_seqs))+
  geom_jitter()+
  scale_y_log10()

coverage%>%
  ggplot(aes(x=1, y= n_seqs))+
  geom_boxplot()

coverage%>%
  ggplot(aes(x=1, y= n_seqs))+
  geom_violin()+
  scale_y_log10()

coverage%>%
  arrange(n_seqs)

coverage%>%
  arrange(n_seqs)%>%
  ggplot(aes(x = 1:nrow(.), y = n_seqs))+
  geom_line()+
  coord_cartesian(xlim = c(20,30), ylim= c(200,15000))

coverage%>%
  arrange(n_seqs)%>%
  print(n=40)

#calculating Goods Coverage (metric to determine if threshold is valid)
coverage_stats <- 
  counts%>%
  group_by(SampleID)%>% 
  summarise(n_seqs = sum(values),
            n_sings = sum(values==1),
            goods = 100* (1- n_sings / n_seqs))%>% 
filter(n_seqs > 700)

#eliminate samples with less than 700 sequences
otufile1 <- otufile1 %>% 
  filter_out(SampleID %in% c("S155","S033","S502","S416","S427","S467"))

ggplot(coverage_stats, aes(x=n_seqs, y= goods))+
  geom_point()

coverage_stats%>%
  arrange(goods) #90% of coverage is reached at 700 sequences 

#rarefaction
rownames(otufile1) <- otufile1$SampleID
otufile1 <- otufile1[,-1]
OTUs <- (rarefy(otufile1, 700))
slopes<-as.data.frame(rareslope(otufile1, 322))
write.table(slopes, "slopes.txt", sep="\t")

map_dfr(rarecurve_data, bind_rows)%>%
  bind_cols(SampleID = rownames(subset_25),.)%>%
  pivot_longer(-SampleID)%>%
  drop_na()%>%
  mutate(n_seqs = as.numeric(str_replace(name,"N","")))%>%
  ggplot(aes(x=n_seqs, y=value, group=SampleID))+
  geom_vline(xintercept = 700, color= "red")+
  geom_line()+
  theme_classic()


#Removing  samples that don't reach plateau (less than 70% coverage)
otufile1 <- otufileraw %>% 
  filter_out(SampleID %in% c("S148","A068b","S070a","A058b","S127","S012",
                            "A057","S227b","A037","A011b","S319b","S044a","A055a",
                            "A051","A085","A054","S253b","S445","A052","S009b",
                            "S095","S034b","A077b","A050b","S124b","S099","S118",
                            "A056"))


## All trimmed out samples##
#A011a	A045a	A047a	A050a	A055b	A058a	A060a	A061a	A062a	A063a
#A065a	A066a	A067a	A068a	A077a	A086b	S009a	S034a	S037a	S044b
#S064a	S067a	S070b	S083a	S085a	S089a	S092a	S098a	S121a	S124a
#S128a	S130a	S139a	S177a	S180a	S187a	S188a	S220a	S227a	S247a
#S252a	S253a	S265a	S271a	S289a	S290a	S298a	S319a	S326b	S328a
#S330a	S337a	S343a	S345a	S354b	S359a	S370a	S379a	S380a	S384a
#S388a	S402a	S405a	S414a	S415a	S423a	S426a	S429a	S433a	S458a
#S459a	S461b	S472a	S473a	S477b	S464a	S486a	S503a	S505a	A056
#S148	A068b	S070a	A058b	S127	S012	S227b	A037	A011b	S319b S467
#A057	S044a	A055a	A085	A054	S253b	S445	A052	S009b	A051 S427
#S095	S034b	A077b	A050b	S124b	S099	S118 S033 S155 S416 S502



# 2.- Remove contamination from database ----
otufile     <- read.table("initial_otu_table.txt", sep = "\t", h=T, row.names = 1)
mapfile     <- read.table("initial_metadata.txt", sep = "\t", h=T, row.names = 1)
taxfile     <- read.table("tax_table.txt", h=T, sep = "\t", row.names = 1) 

OTU = otu_table(otufile, taxa_are_rows = TRUE)
taxfile <- as.matrix(taxfile)
TAX = tax_table(taxfile)
samples = sample_data(mapfile)

database.pso <- phyloseq(OTU, TAX, samples)

#rarefaction
set.seed(1000, kind = NULL, normal.kind = NULL, sample.kind = NULL)
database.pso_r <- rarefy_even_depth(database.pso)

#Identify Contaminants with witness with a treshhold of 0.1

sample_data(database.pso)$is.neg <- sample_data(database.pso_r)$Sample_Control == "Control"
contamdf.prev <- isContaminant(database.pso, method="prevalence", neg="is.neg", threshold=0.1)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

#with a treshhold of 0.5
contamdf.prev05 <- isContaminant(database.pso, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
head(which(contamdf.prev05$contaminant)) 

write.table(contamdf.prev05, "contaminant-0.5.txt", sep="\t")
write.table(contamdf.prev, "contaminant-0.1.txt", sep="\t")


# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa <- transform_sample_counts(database.pso, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.prev$contaminant)

## ELIMINATE CONTAMINANT zOTUs FROM DATA FILES WITH 0.1 THRESHOLD 

# 3.- Determine differences between dorsal and ventral adult samples----
#load data
mapfileAD    <- read.table("metadata_adults.txt", h=T, row.names = 1)
OTUads    <- read.table("otu_table-adults.txt", h=T, row.names = 1)
OTUads  <- as.data.frame(t(OTUads))

rownames(OTUads) <- OTUads$SampleID
OTUads <- OTUads[,-1]
OTUads <- (rarefy(OTUads, 700))

#compute bray-curtis distances
bray1 <- vegdist(OTUads, method="bray",binary=FALSE, diag=FALSE, upper=FALSE, na.rm=FALSE) 

#create bray-curtis dataframe
braydf1 = (as.data.frame(as.matrix(bray1))) 
write.table(braydf1, "bray-curtis_adults.txt", sep="\t")

#mMDS
set.seed(100)
nMDS<-metaMDS(bray1,k=2,stress=0.1)

# anosim to test the significance
anosim(bray1, mapfileAD$Side, permutations = 999, distance = "bray")

nMDS$points
scores(nMDS)

#NMDS Plot
scores(nMDS)%>%
  as_tibble(rownames="SampleID")%>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=mapfileAD$Side))+
  geom_point()

# 4.- Creating working-PSO, cleaning and rarefying ----

#load corresponding files: 1) otu_table_60.txt; 2) otu_table_75.txt
                          #1) metadata_60.txt; 2) metadata_75.txt

otufile     <- read.table("otu_table_60.txt", sep = "\t", h=T, row.names = 1)
mapfile     <- read.table("metadata_60.txt", sep = "\t", h=T, row.names = 1)
taxfile     <- read.table("tax_table_60.txt", h=T, sep = "\t", row.names = 1) 

OTU = otu_table(otufile, taxa_are_rows = TRUE)
taxfile <- as.matrix(taxfile)
TAX = tax_table(taxfile)
samples = sample_data(mapfile)

database.pso <- phyloseq(OTU, TAX, samples)

#clean#cleanfasta
database.pso <- tax_fix(database.pso, unknowns = 'NA' )
database.pso <- subset_taxa(database.pso, Kingdom!=" k__Eukaryota")
database.pso <- subset_taxa(database.pso, Kingdom!=" k__Archaea")
database.pso <- subset_taxa(database.pso, Family!=" f__Mitochondria")
database.pso <- subset_taxa(database.pso, Order!=" o__Chloroplast")

#rarefaction
set.seed(1000, kind = NULL, normal.kind = NULL, sample.kind = NULL)
database.pso_r <- rarefy_even_depth(database.pso)
database.pso <- database.pso_r 


# 5.- Differences in samples types ----
## Hellinger normalization
otu_hellinger <- t(otufile)
otu_hellinger <- as.matrix(otu_hellinger)
otu_hellinger <- decostand(otu_hellinger, method = "hellinger")

#mMDS
mMDS6 <- metaMDS(otu_hellinger, distance = "euclidean", k = 2, 
                 trymax = 100, stress = 0.1)

nmds_coords <- as.data.frame(scores(mMDS6, display = "sites"))
nmds_coords$Sample_type <- samples$Sample_type 
stress_valor <- round(mMDS6$stress, 3)

#NMDS Plot
ggplot(nmds_coords, aes(x = NMDS1, y = NMDS2, color = Sample_type)) +
  geom_point(size = 3, alpha = 0.8) +                  
  stat_ellipse(aes(fill = Sample_type),               
               geom = "polygon", type = "norm",alpha = 0.1, 
               show.legend = FALSE) + 
  labs(
    title = "NMDS Hellinger distances by Sample Type",
    subtitle = paste("Stress =", stress_valor), 
    x = "NMDS 1",
    y = "NMDS 2",
    color = "Sample type"
  ) +
  theme_bw(base_size = 14) +                  
  theme(
    panel.grid.minor = element_blank(),  
    legend.position = "right",       
    plot.title = element_text(face = "bold"))

# ANOSIM to test significance
anosim(otu_hellinger, samples$Sample_type, permutations = 999)

#pairwise post hoc
dist_hellinger <- vegdist(otu_hellinger, method = "euclidean")
anosim.pairwise(x = dist_hellinger, samples$Sample_type, perm = 999)


# 6.- Microbial composition ----

# Relative abundance
#same script should be use to extract Phylum and Genus data

#Males
Males <- subset_samples(database.pso, Sample_type == "Male")

Males <- tax_glom(Males, taxrank = "Genus")
Males <- transform_sample_counts(Males, function(x)100* x / sum(x))
OTUMa <- otu_table(Males)
TAXMa <- tax_table(Males)[,"Genus"]

AverageMa <- as.data.frame(rowMeans(OTUMa))
names(AverageMa) <- c("Mean")
TableMa <- merge(TAXMa, AverageMa, by=0, all=TRUE)
TableMa$Row.names = NULL
TableMa <- TableMa[order(desc(TableMa$Mean)),]

#Females
Females <- subset_samples(database.pso, Sample_type == "Female")

Females <- tax_glom(Females, taxrank = "Genus")
Females <- transform_sample_counts(Females, function(x)100* x / sum(x))
OTUF <- otu_table(Females)
TAXF <- tax_table(Females)[,"Genus"]

AverageF <- as.data.frame(rowMeans(OTUF))
names(AverageF) <- c("Mean")
TableF <- merge(TAXF, AverageF, by=0, all=TRUE)
TableF$Row.names = NULL
TableF <- TableF[order(desc(TableF$Mean)),]

#Substrate
Substrate <- subset_samples(database.pso, Sample_type == "ENVIRONMENT")

Substrate <- tax_glom(Substrate, taxrank = "Genus")
Substrate <- transform_sample_counts(Substrate, function(x)100* x / sum(x))
OTUS <- otu_table(Substrate)
TAXS <- tax_table(Substrate)[,"Genus"]

AverageS <- as.data.frame(rowMeans(OTUS))
names(AverageS) <- c("Mean")
TableS <- merge(TAXS, AverageS, by=0, all=TRUE)
TableS$Row.names = NULL
TableS <- TableS[order(desc(TableS$Mean)),]

#Water
Water<- subset_samples(database.pso, Sample_type == "WATER")

Water <- tax_glom(Water, taxrank = "Genus")
Water <- transform_sample_counts(Water, function(x)100* x / sum(x))
OTUW <- otu_table(Water)
TAXW <- tax_table(Water)[,"Genus"]

AverageW <- as.data.frame(rowMeans(OTUW))
names(AverageW) <- c("Mean")
TableW <- merge(TAXW, AverageW, by=0, all=TRUE)
TableW$Row.names = NULL
TableW <- TableW[order(desc(TableW$Mean)),]

#Clutch week 1
Clutch_W1 <- subset_samples(database.pso, Sample_type == "CLW1")

Clutch_W1 <- tax_glom(Clutch_W1, taxrank = "Genus")
Clutch_W1 <- transform_sample_counts(Clutch_W1, function(x)100* x / sum(x))
OTUW1 <- otu_table(Clutch_W1)
TAXW1 <- tax_table(Clutch_W1)[,"Genus"]

AverageW1 <- as.data.frame(rowMeans(OTUW1))
names(AverageW1) <- c("Mean")
TableW1 <- merge(TAXW1, AverageW1, by=0, all=TRUE)
TableW1$Row.names = NULL
TableW1 <- TableW1[order(desc(TableW1$Mean)),]

#Clutch week 2
Clutch_W2 <- subset_samples(database.pso, Sample_type == "CLW2")

Clutch_W2 <- tax_glom(Clutch_W2, taxrank = "Genus")
Clutch_W2 <- transform_sample_counts(Clutch_W2, function(x)100* x / sum(x))
OTUW2 <- otu_table(Clutch_W2)
TAXW2 <- tax_table(Clutch_W2)[,"Genus"]

AverageW2 <- as.data.frame(rowMeans(OTUW2))
names(AverageW2) <- c("Mean")
TableW2 <- merge(TAXW2, AverageW2, by=0, all=TRUE)
TableW2$Row.names = NULL
TableW2 <- TableW2[order(desc(TableW2$Mean)),]

#Tadpoles
Tadpoles <- subset_samples(database.pso, Sample_type == "TP")

Tadpoles <- tax_glom(Tadpoles, taxrank = "Genus")
Tadpoles <- transform_sample_counts(Tadpoles, function(x)100* x / sum(x))
OTUT <- otu_table(Tadpoles)
TAXT <- tax_table(Tadpoles)[,"Genus"]

AverageT <- as.data.frame(rowMeans(OTUT))
names(AverageT) <- c("Mean")
TableT <- merge(TAXT, AverageT, by=0, all=TRUE)
TableT$Row.names = NULL
TableT <- TableT[order(desc(TableT$Mean)),]

#Metamorphs
Metamorphs <- subset_samples(database.pso, Sample_type == "METAMORPH")

Metamorphs <- tax_glom(Metamorphs, taxrank = "Genus")
Metamorphs <- transform_sample_counts(Metamorphs, function(x)100* x / sum(x))
OTUM <- otu_table(Metamorphs)
TAXM <- tax_table(Metamorphs)[,"Genus"]

AverageM <- as.data.frame(rowMeans(OTUM))
names(AverageM) <- c("Mean")
TableM <- merge(TAXM, AverageM, by=0, all=TRUE)
TableM$Row.names = NULL
TableM <- TableM[order(desc(TableM$Mean)),]

# 7.- Intersection plots per treatment ---- 
NormalT <- subset_samples(database.pso, Posit_ctrl == "1")
ggvenn_pq(NormalT, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE)

NoT <- subset_samples(database.pso, Neg_ctrl == "1")
ggvenn_pq(NoT, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE)

TT <- subset_samples(database.pso, Test == "1")
ggvenn_pq(TT, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE)

# 8.- Intersection plots per trial ----
#Trial 1
T1 <- subset_samples(database.pso, S1 == "1")
ggvenn_pq(T1,fact = "Sample_type", 
          min_nb_seq = 0,print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 2
T2 <- subset_samples(database.pso, S2 == "1")
ggvenn_pq(T2, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 3
T3 <- subset_samples(database.pso, S3 == "1")
ggvenn_pq(T3, 
          fact = "FosterS3", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 5
T5 <- subset_samples(database.pso, S5 == "1")
ggvenn_pq(T5, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 6
T6 <- subset_samples(database.pso, S6 == "1")
ggvenn_pq(T6, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 9
T9 <- subset_samples(database.pso, S9 == "1")
ggvenn_pq(T9, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 10
T10 <- subset_samples(database.pso, S10 == "1")
ggvenn_pq(T10, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 12
T12 <- subset_samples(database.pso, S12 == "1")
ggvenn_pq(T12, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 15
T15 <- subset_samples(database.pso, S15 == "1")
ggvenn_pq(T15, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 16
T16 <- subset_samples(database.pso, S16 == "1")
ggvenn_pq(T16, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 17
T17 <- subset_samples(database.pso, S17 == "1")
ggvenn_pq(T17, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 18
T18 <- subset_samples(database.pso, S18 == "1")
ggvenn_pq(T18, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 23
T23 <- subset_samples(database.pso, S23 == "1")
ggvenn_pq(T23, 
          fact = "FosterS23", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 24
T24 <- subset_samples(database.pso, S24 == "1")
ggvenn_pq(T24, 
          fact = "FosterS24", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 26
T26 <- subset_samples(database.pso, S26 == "1")
ggvenn_pq(T26, 
          fact = "FosterS26", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 28
T28 <- subset_samples(database.pso, S28 == "1")
ggvenn_pq(T28, 
          fact = "FosterS28", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 29
T29 <- subset_samples(database.pso, S29 == "1")
ggvenn_pq(T29, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 34
T34 <- subset_samples(database.pso, S34 == "1")
ggvenn_pq(T34, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 35
T35 <- subset_samples(database.pso, S35 == "1")
ggvenn_pq(T35, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 37
T37 <- subset_samples(database.pso, S37 == "1")
ggvenn_pq(T37, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 39
T39 <- subset_samples(database.pso, S39 == "1")
ggvenn_pq(T39, 
          fact = "FosterS39", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 40
T40 <- subset_samples(database.pso, S40 == "1")
ggvenn_pq(T40, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 42
T42 <- subset_samples(database.pso, S42 == "1")
ggvenn_pq(T42, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 43
T43 <- subset_samples(database.pso, S43 == "1")
ggvenn_pq(T43, 
          fact = "Sample_type", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 48
T48 <- subset_samples(database.pso, S48 == "1")
ggvenn_pq(T48, 
          fact = "FosterS48", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 49
T49 <- subset_samples(database.pso, S49 == "1")
ggvenn_pq(T49, 
          fact = "FosterS49", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)

#Trial 52
T52 <- subset_samples(database.pso, S52 == "1")
ggvenn_pq(T52, 
          fact = "FosterS52", 
          min_nb_seq = 0, 
          print_values = TRUE,
          add_nb_samples = FALSE,force_upset = TRUE)


# 9.- Alpha-diversity analyses ----
richness <- estimate_richness(database.pso)

#### CHAO1
# Check normality
hist(richness$Chao1, main = "Chao1 index", xlab = "")
shapiro.test(richness$Chao1)
Sample_type = factor(sample_data(database.pso)$Sample_type)

# test
kruskal.test(richness$Chao1 ~ Sample_type)

# Perform pairwise Wilcoxon test with p-value adjustment
full_data_Chao1 <- na.omit(data.frame(Chao1 = richness$Chao1, Sample_type))

Chao1res <- full_data_Chao1 %>% 
  pairwise_wilcox_test(Chao1 ~ Sample_type, p.adjust.method = "bonferroni")

print(Chao1res)
write.table(Chao1res, "Chao1_stats_60.txt", sep="\t")

#### SHANNON

hist(richness$Shannon, main = "Shannon index", xlab = "")
shapiro.test(richness$Shannon)

kruskal.test(richness$Shannon ~ sample_data(database.pso)$Sample_type)

# Perform pairwise Wilcoxon test with p-value adjustment
full_data_Shan <- na.omit(data.frame(Shannon = richness$Shannon, Sample_type))

Shannres <- full_data_Shan %>% 
  pairwise_wilcox_test(Shannon ~ Sample_type, p.adjust.method = "bonferroni")

print(Shannres)
write.table(Shannres, "Shannon_stats_60.txt", sep="\t")

##### INVERSE SIMPSON

hist(richness$InvSimpson, main = "InvSimpson index", xlab = "")
shapiro.test(richness$InvSimpson)

kruskal.test(richness$InvSimpson ~ sample_data(database.pso)$Sample_type)

full_data_InvSim <- na.omit(data.frame(InvSimpson = richness$InvSimpson,Sample_type))

InvSimpres <- full_data_InvSim %>% 
  pairwise_wilcox_test(InvSimpson ~ Sample_type, p.adjust.method = "bonferroni")

print(InvSimpres)
write.table(InvSimpres, "InvSimpson_stats_60.txt", sep="\t")


## ESTIMATING ALPHA DIVERSITY PER OFFSPRING SAMPLE TYPE

offspring <- subset_samples(database.pso, GruopID == "offspring")

#Clutch week 1
CW1 <- subset_samples(offspring, Sample_type == "CLW1")
richness_CW1 <- estimate_richness(CW1)

hist(richness_CW1$Shannon, main="Shannon index", xlab="")
shapiro.test(richness_CW1$Shannon)
summary(aov(richness_CW1$Shannon ~ sample_data(CW1)$transport))

hist(richness_CW1$Chao1, main="Chao1 index", xlab="")
shapiro.test(richness_CW1$Chao1)
summary(aov(richness_CW1$Chao1 ~ sample_data(CW1)$transport))

hist(richness_CW1$InvSimpson, main="Inverse Simpson index", xlab="")
shapiro.test(richness_CW1$InvSimpson)
kruskal.test(richness_CW1$InvSimpson ~ sample_data(CW1)$transport)

#Clutch week 2
CW2 <- subset_samples(offspring, Sample_type == "CLW2")
richness_CW2 <- estimate_richness(CW2)

hist(richness_CW2$Shannon, main="Shannon index", xlab="")
shapiro.test(richness_CW2$Shannon)
summary(aov(richness_CW2$Shannon ~ sample_data(CW2)$transport))

hist(richness_CW2$Chao1, main="Chao1 index", xlab="")
shapiro.test(richness_CW2$Chao1)
summary(aov(richness_CW2$Chao1 ~ sample_data(CW2)$transport))

hist(richness_CW2$InvSimpson, main="Inverse Simpson index", xlab="")
shapiro.test(richness_CW2$InvSimpson)
summary(aov(richness_CW2$InvSimpson ~ sample_data(CW2)$transport))

#Tadpoles
TP <- subset_samples(offspring, Sample_type == "TP")
richness_TP <- estimate_richness(TP)

hist(richness_TP$Shannon, main="Shannon index", xlab="")
shapiro.test(richness_TP$Shannon)
kruskal.test(richness_TP$Shannon ~ sample_data(TP)$transport)

hist(richness_TP$Chao1, main="Chao1 index", xlab="")
shapiro.test(richness_TP$Chao1)
kruskal.test(richness_TP$Chao1 ~ sample_data(TP)$transport)

hist(richness_TP$InvSimpson, main="Inverse Simpson index", xlab="")
shapiro.test(richness_TP$InvSimpson)
kruskal.test(richness_TP$InvSimpson ~ sample_data(TP)$transport)

#Metamorphs
Meta <- subset_samples(offspring, Sample_type == "METAMORPH")
richness_Meta <- estimate_richness(Meta)

hist(richness_Meta$Shannon, main="Shannon index", xlab="")
shapiro.test(richness_Meta$Shannon)
kruskal.test(richness_Meta$Shannon ~ sample_data(Meta)$transport)

hist(richness_Meta$Chao1, main="Chao1 index", xlab="")
shapiro.test(richness_Meta$Chao1)
summary(aov(richness_Meta$Chao1 ~ sample_data(Meta)$transport))

hist(richness_Meta$InvSimpson, main="Inverse Simpson index", xlab="")
shapiro.test(richness_Meta$InvSimpson)
kruskal.test(richness_Meta$InvSimpson ~ sample_data(Meta)$transport)


# 10.- Contribution proportion of different sources ----
# Load data
contribs   <- read.table("Contributions-parent_offspring.txt", h=T)
contrib   <- subset(contribs, !source == "Unknown")

#Model 1
#Prepare data with proportions of contribution by all source types
contrib$source <- as.factor(contrib$source) 
contrib$treatment <- as.factor(contrib$treatment) 
contrib$Trial <- as.factor(contrib$Trial) 

#model for testing differences in contribution by type of source
model1 <- glmmTMB(contribution_median ~ treatment + source + Offspring + (1|Trial),
                  data = contrib, family=beta_family(), 
                  dispformula = ~ treatment + source + Offspring)

glmmTMB:::Anova.glmmTMB(model1)
emmeansM1 <- emmeans(model1, pairwise ~ source )

#diagnostics of model
testDispersion(simulateResiduals(model1, plot = TRUE))


#Model 2 
#Prepare data only for foster treatment with proportions of contribution by all source types
contribFOSTER   <- subset(contrib, treatment == "Foster_T")

#model for testing differences in contribution by type of source in Foster Transport treatment
model2 <- glmmTMB(contribution_median ~ source + Offspring + (1|Trial),
                  data = contribFOSTER, family=beta_family(),
                  dispformula = ~ source + Offspring)

glmmTMB:::Anova.glmmTMB(model2)
emmeansM2 <- emmeans(model2, pairwise ~ source)

#diagnostics of model
testDispersion(simulateResiduals(model2, plot = TRUE))

#Model 3 
#Prepare data only for no transfer treatment with proportions of contribution by all source types
contribNO   <- subset(contrib, treatment == "No_T")

#model for testing differences in contribution by type of source
model3 <- glmmTMB(contribution_median ~ source + (1|Trial),
                  data = contribNO, family=beta_family(),
                  dispformula = ~ source)

glmmTMB:::Anova.glmmTMB(model3)
emmeansM3 <- emmeans(model3, pairwise ~ source)

#diagnostics of model
testDispersion(simulateResiduals(model3, plot = TRUE))

# 11.- Contribution proportion of developmental steps ----
# Load data
devo   <- read.table("Contributions-dev_stage.txt", h=T)
devo$dev_step <- as.factor(devo$dev_step) 
devo$treatment <- as.factor(devo$treatment) 
devo$Trial <- as.factor(devo$Trial) 

#model 4 testing differences in contribution from developmental stages
model4 <- glmmTMB(contribution ~ treatment + dev_step + (1|Trial),
               data = devo, family=beta_family(), 
               dispformula = ~treatment+ dev_step) #to control for heterocedasticity/overdispersion of residuals

glmmTMB:::Anova.glmmTMB(model4)
emmeans(model4, pairwise ~  dev_step)

simulateResiduals(model4, plot = TRUE)

# 12.- Extraction of acquired taxa by transportation ----

grupo_A <- "Foster"
grupo_B <- "TP"
grupo_C <- "METAMORPH"
grupo_D <- "Biological"
grupo_E <- "Male"
taxa_table <- as.data.frame(tax_table(database.pso))

#S1
S1 <- subset_samples(database.pso, S1 == "1")
pS1 <- merge_samples(S1, "Sample_type")
binaryS1 <- as.data.frame(t(otu_table(pS1)))
binaryS1[binaryS1 > 0] <- 1

#Tadpoles - Biological
sharedST1 <- rownames(binaryS1[
  binaryS1[[grupo_E]] == 1 & 
    binaryS1[[grupo_B]] == 1 & 
    rowSums(binaryS1) == 2, ])

ST1 <- c("ZOTU4609", "ZOTU4195", "ZOTU4749", "ZOTU5908") 
taxaST1 <- taxa_table[ST1, ]

#Metamorphs - Biological
sharedSM1 <- rownames(binaryS1[
  binaryS1[[grupo_E]] == 1 & 
    binaryS1[[grupo_C]] == 1 & 
    rowSums(binaryS1) == 2, ])

SM1 <- c("ZOTU10303", "ZOTU532", "ZOTU18855", "ZOTU1382", "ZOTU12547") 
taxaSM1 <- taxa_table[SM1, ]

#S2
S2 <- subset_samples(database.pso, S2 == "1")
pS2 <- merge_samples(S2, "Sample_type")
binaryS2 <- as.data.frame(t(otu_table(pS2)))
binaryS2[binaryS2 > 0] <- 1

#Tadpoles - Biological
sharedST2 <- rownames(binaryS2[
  binaryS2[[grupo_E]] == 1 & 
    binaryS2[[grupo_B]] == 1 & 
    rowSums(binaryS2) == 2, ])

ST2 <- c("ZOTU5425","ZOTU348", "ZOTU3977","ZOTU827","ZOTU15696",
         "ZOTU852","ZOTU11051","ZOTU5335","ZOTU7934","ZOTU19823",
         "ZOTU7937","ZOTU20359","ZOTU263") 
taxaST2 <- taxa_table[ST2, ]

#Metamorphs - Biological
sharedSM2 <- rownames(binaryS2[
  binaryS2[[grupo_E]] == 1 & 
    binaryS2[[grupo_C]] == 1 & 
    rowSums(binaryS2) == 2, ])

SM2 <- c("ZOTU20629","ZOTU3268","ZOTU3366") 
taxaSM2 <- taxa_table[SM2, ]

#S3
S3 <- subset_samples(database.pso, S3 == "1")
pS3 <- merge_samples(S3, "FosterS3")
binaryS3 <- as.data.frame(t(otu_table(pS3)))
binaryS3[binaryS3 > 0] <- 1

#Tadpoles - Foster
sharedSTF3 <- rownames(binaryS3[
  binaryS3[[grupo_A]] == 1 & 
    binaryS3[[grupo_B]] == 1 & 
    rowSums(binaryS3) == 2, ])

STF3 <- c("ZOTU133") 
taxaSTF3 <- taxa_table[STF3, ]

#Metamorphs - Foster
sharedSMF3 <- rownames(binaryS3[
  binaryS3[[grupo_A]] == 1 & 
    binaryS3[[grupo_C]] == 1 & 
    rowSums(binaryS3) == 2, ])

SMF3 <- c("ZOTU1207","ZOTU435" ) 
taxaSMF3 <- taxa_table[SMF3, ]

#Tadpoles - Biological
sharedST3 <- rownames(binaryS3[
  binaryS3[[grupo_D]] == 1 & 
    binaryS3[[grupo_B]] == 1 & 
    rowSums(binaryS3) == 2, ])

ST3 <- c("ZOTU2098","ZOTU52","ZOTU506","ZOTU2047","ZOTU803",
         "ZOTU99","ZOTU20323","ZOTU937","ZOTU7742","ZOTU5558",
         "ZOTU5832","ZOTU51") 
taxaST3 <- taxa_table[ST3, ]

#Metamorphs - Biological
sharedSM3 <- rownames(binaryS3[
  binaryS3[[grupo_D]] == 1 & 
    binaryS3[[grupo_C]] == 1 & 
    rowSums(binaryS3) == 2, ])

SM3 <- c("ZOTU18693","ZOTU4238","ZOTU1220") 
taxaSM3 <- taxa_table[SM3, ]

#S6
S6 <- subset_samples(database.pso, S6 == "1")
pS6 <- merge_samples(S6, "Sample_type")
binaryS6 <- as.data.frame(t(otu_table(pS6)))
binaryS6[binaryS6 > 0] <- 1

#Tadpoles - Biological
sharedST6 <- rownames(binaryS6[
  binaryS6[[grupo_E]] == 1 & 
    binaryS6[[grupo_B]] == 1 & 
    rowSums(binaryS6) == 2, ])

ST6 <- c("ZOTU12021","ZOTU1961","ZOTU12232","ZOTU663","ZOTU6669",
         "ZOTU14083","ZOTU7175","ZOTU6411") 
taxaST6 <- taxa_table[ST6, ]

#Metamorphs - Biological
sharedSM6 <- rownames(binaryS6[
  binaryS6[[grupo_E]] == 1 & 
    binaryS6[[grupo_C]] == 1 & 
    rowSums(binaryS6) == 2, ])

SM6 <- c("ZOTU87","ZOTU2534","ZOTU2042","ZOTU4238") 
taxaSM6 <- taxa_table[SM6, ]

#S9
S9 <- subset_samples(database.pso, S9 == "1")
pS9 <- merge_samples(S9, "Sample_type")
binaryS9 <- as.data.frame(t(otu_table(pS9)))
binaryS9[binaryS9 > 0] <- 1

#Tadpoles - Biological
sharedST9 <- rownames(binaryS9[
  binaryS9[[grupo_E]] == 1 & 
    binaryS9[[grupo_B]] == 1 & 
    rowSums(binaryS9) == 2, ])

ST9 <- c("ZOTU6152","ZOTU3386") 
taxaST9 <- taxa_table[ST9, ]

#Metamorphs - Biological
sharedSM9 <- rownames(binaryS9[
  binaryS9[[grupo_E]] == 1 & 
    binaryS9[[grupo_C]] == 1 & 
    rowSums(binaryS9) == 2, ])

SM9 <- c("ZOTU17506","ZOTU6151") 
taxaSM9 <- taxa_table[SM9, ]

#S10
S10 <- subset_samples(database.pso, S10 == "1")
pS10 <- merge_samples(S10, "Sample_type")
binaryS10 <- as.data.frame(t(otu_table(pS10)))
binaryS10[binaryS10 > 0] <- 1

#Tadpoles - Biological
sharedST10 <- rownames(binaryS10[
  binaryS10[[grupo_E]] == 1 & 
    binaryS10[[grupo_B]] == 1 & 
    rowSums(binaryS10) == 2, ])

ST10 <- c("ZOTU380","ZOTU615","ZOTU12227","ZOTU4537","ZOTU11007",
          "ZOTU4904","ZOTU2571") 
taxaST10 <- taxa_table[ST10, ]

#Metamorphs - Biological
sharedSM10 <- rownames(binaryS10[
  binaryS10[[grupo_E]] == 1 & 
    binaryS10[[grupo_C]] == 1 & 
    rowSums(binaryS10) == 2, ])

SM10 <- c("ZOTU4387","ZOTU1853","ZOTU503","ZOTU9069") 
taxaSM10 <- taxa_table[SM10, ]

#S12
S12 <- subset_samples(database.pso, S12 == "1")
pS12 <- merge_samples(S12, "Sample_type")
binaryS12 <- as.data.frame(t(otu_table(pS12)))
binaryS12[binaryS12 > 0] <- 1

#Metamorphs - Biological
sharedSM12 <- rownames(binaryS12[
  binaryS12[[grupo_E]] == 1 & 
    binaryS12[[grupo_C]] == 1 & 
    rowSums(binaryS12) == 2, ])

SM12 <- c("ZOTU212","ZOTU18693","ZOTU300","ZOTU17506","ZOTU10902",
          "ZOTU5561","ZOTU7467") 
taxaSM12 <- taxa_table[SM12, ]

#S23
S23 <- subset_samples(database.pso, S23 == "1")
pS23 <- merge_samples(S23, "FosterS23")
binaryS23 <- as.data.frame(t(otu_table(pS23)))
binaryS23[binaryS23 > 0] <- 1

#Tadpoles - Foster
sharedSTF23 <- rownames(binaryS23[
  binaryS23[[grupo_A]] == 1 & 
    binaryS23[[grupo_B]] == 1 & 
    rowSums(binaryS23) == 2, ])

STF23 <- c("ZOTU21763","ZOTU5486","ZOTU21200") 
taxaSTF23 <- taxa_table[STF23, ]

#Metamorphs - Foster
sharedSMF23 <- rownames(binaryS23[
  binaryS23[[grupo_A]] == 1 & 
    binaryS23[[grupo_C]] == 1 & 
    rowSums(binaryS23) == 2, ])

SMF23 <- c("ZOTU687","ZOTU12021","ZOTU287","ZOTU827","ZOTU4238",
           "ZOTU19775","ZOTU3090","ZOTU5976") 
taxaSMF23 <- taxa_table[SMF23, ]

#Tadpoles - Biological
sharedST23 <- rownames(binaryS23[
  binaryS23[[grupo_D]] == 1 & 
    binaryS23[[grupo_B]] == 1 & 
    rowSums(binaryS23) == 2, ])

ST23 <- c("ZOTU1157","ZOTU2473","ZOTU212","ZOTU18459") 
taxaST23 <- taxa_table[ST23, ]

#Metamorphs - Biological
sharedSM23 <- rownames(binaryS23[
  binaryS23[[grupo_D]] == 1 & 
    binaryS23[[grupo_C]] == 1 & 
    rowSums(binaryS23) == 2, ])

SM23 <- c("ZOTU659","ZOTU18467","ZOTU6781","ZOTU12227","ZOTU267",
          "ZOTU6260","ZOTU691","ZOTU668","ZOTU11007") 
taxaSM23 <- taxa_table[SM23, ]

#S24
S24 <- subset_samples(database.pso, S24 == "1")
pS24 <- merge_samples(S24, "FosterS24")
binaryS24 <- as.data.frame(t(otu_table(pS24)))
binaryS24[binaryS24 > 0] <- 1

#Tadpoles - Foster
sharedSTF24 <- rownames(binaryS24[
  binaryS24[[grupo_A]] == 1 & 
    binaryS24[[grupo_B]] == 1 & 
    rowSums(binaryS24) == 2, ])

STF24 <- c("ZOTU1212","ZOTU659","ZOTU2473","ZOTU6260","ZOTU16416","ZOTU550") 
taxaSTF24 <- taxa_table[STF24, ]

#Metamorphs - Foster
sharedSMF24 <- rownames(binaryS24[
  binaryS24[[grupo_A]] == 1 & 
    binaryS24[[grupo_C]] == 1 & 
    rowSums(binaryS24) == 2, ])

#Tadpoles - Biological
sharedST24 <- rownames(binaryS24[
  binaryS24[[grupo_D]] == 1 & 
    binaryS24[[grupo_B]] == 1 & 
    rowSums(binaryS24) == 2, ])

ST24 <- c("ZOTU687","ZOTU4258","ZOTU4660","ZOTU15060","ZOTU11243","ZOTU19422") 
taxaST24 <- taxa_table[ST24, ]

#Metamorphs - Biological
sharedSM24 <- rownames(binaryS24[
  binaryS24[[grupo_A]] == 1 & 
    binaryS24[[grupo_D]] == 1 & 
    rowSums(binaryS24) == 2, ])

SM24 <- c("ZOTU579","ZOTU19002","ZOTU781","ZOTU20894","ZOTU725") 
taxaSM24 <- taxa_table[SM24, ]

#S26
S26 <- subset_samples(database.pso, S26 == "1")
pS26 <- merge_samples(S26, "FosterS26")
binaryS26 <- as.data.frame(t(otu_table(pS26)))
binaryS26[binaryS26 > 0] <- 1

#Tadpoles - Foster
sharedSTF26 <- rownames(binaryS26[
  binaryS26[[grupo_A]] == 1 & 
    binaryS26[[grupo_B]] == 1 & 
    rowSums(binaryS26) == 2, ])

STF26 <- c("ZOTU2135","ZOTU12535","ZOTU16279") 
taxaSTF26 <- taxa_table[STF26, ]

#Metamorphs - Foster
sharedSMF26 <- rownames(binaryS26[
  binaryS26[[grupo_A]] == 1 & 
    binaryS26[[grupo_C]] == 1 & 
    rowSums(binaryS26) == 2, ])

SMF26 <- c("ZOTU6970","ZOTU2098","ZOTU10469","ZOTU4100","ZOTU4379",
           "ZOTU4642","ZOTU5521","ZOTU20009","ZOTU19791","ZOTU2168",
           "ZOTU7175") 
taxaSMF26 <- taxa_table[SMF26, ]

#Tadpoles - Biological
sharedST26 <- rownames(binaryS26[
  binaryS26[[grupo_D]] == 1 & 
    binaryS26[[grupo_B]] == 1 & 
    rowSums(binaryS26) == 2, ])

ST26 <- c("ZOTU9339","ZOTU11346","ZOTU294") 
taxaST26 <- taxa_table[ST26, ]

#Metamorphs - Biological
sharedSM26 <- rownames(binaryS26[
  binaryS26[[grupo_D]] == 1 & 
    binaryS26[[grupo_C]] == 1 & 
    rowSums(binaryS26) == 2, ])

SM26 <- c("ZOTU4782","ZOTU7423","ZOTU8344","ZOTU14213","ZOTU16517",
          "ZOTU4181","ZOTU18855","ZOTU1571","ZOTU181","ZOTU7537",
          "ZOTU4626","ZOTU12261") 
taxaSM26 <- taxa_table[SM26, ]

#S28
S28 <- subset_samples(database.pso, S28 == "1")
pS28 <- merge_samples(S28, "FosterS28")
binaryS28 <- as.data.frame(t(otu_table(pS28)))
binaryS28[binaryS28 > 0] <- 1

#Tadpoles - Foster
sharedSTF28 <- rownames(binaryS28[
  binaryS28[[grupo_A]] == 1 & 
    binaryS28[[grupo_B]] == 1 & 
    rowSums(binaryS28) == 2, ])

STF28 <- c("ZOTU3790","ZOTU5980") 
taxaSTF28 <- taxa_table[STF28, ]

#Metamorphs - Foster
sharedSMF28 <- rownames(binaryS28[
  binaryS28[[grupo_A]] == 1 & 
    binaryS28[[grupo_C]] == 1 & 
    rowSums(binaryS28) == 2, ])

SMF28 <- c("ZOTU707","ZOTU4517","ZOTU4723","ZOTU4318","ZOTU11206") 
taxaSMF28 <- taxa_table[SMF28, ]

#Tadpoles - Biological
sharedST28 <- rownames(binaryS28[
  binaryS28[[grupo_D]] == 1 & 
    binaryS28[[grupo_B]] == 1 & 
    rowSums(binaryS28) == 2, ])

ST28 <- c("ZOTU2529") 
taxaST28 <- taxa_table[ST28, ]

#Metamorphs - Biological
sharedSM28 <- rownames(binaryS28[
  binaryS28[[grupo_D]] == 1 & 
    binaryS28[[grupo_C]] == 1 & 
    rowSums(binaryS28) == 2, ])

SM28 <- c("ZOTU18459","ZOTU3386") 
taxaSM28 <- taxa_table[SM28, ]

#S34
S34 <- subset_samples(database.pso, S34 == "1")
pS34 <- merge_samples(S34, "Sample_type")
binaryS34 <- as.data.frame(t(otu_table(pS34)))
binaryS34[binaryS34 > 0] <- 1

#Tadpoles - Biological
sharedST34 <- rownames(binaryS34[
  binaryS34[[grupo_E]] == 1 & 
    binaryS34[[grupo_B]] == 1 & 
    rowSums(binaryS34) == 2, ])

ST34 <- c("ZOTU18459","ZOTU14083","ZOTU5542","ZOTU1014") 
taxaST34 <- taxa_table[ST34, ]

#Metamorphs - Biological
sharedSM34 <- rownames(binaryS34[
  binaryS34[[grupo_E]] == 1 & 
    binaryS34[[grupo_C]] == 1 & 
    rowSums(binaryS34) == 2, ])

SM34 <- c("ZOTU493","ZOTU17455") 
taxaSM34 <- taxa_table[SM34, ]

#S42
S42 <- subset_samples(database.pso, S42 == "1")
pS42 <- merge_samples(S42, "Sample_type")
binaryS42 <- as.data.frame(t(otu_table(pS42)))
binaryS42[binaryS42 > 0] <- 1

#Tadpoles - Biological
sharedST42 <- rownames(binaryS42[
  binaryS42[[grupo_E]] == 1 & 
    binaryS42[[grupo_B]] == 1 & 
    rowSums(binaryS42) == 2, ])

ST42 <- c("ZOTU7681","ZOTU16715","ZOTU6117")
taxaST42 <- taxa_table[ST42, ]

#Metamorphs - Biological
sharedSM42 <- rownames(binaryS42[
  binaryS42[[grupo_E]] == 1 & 
    binaryS42[[grupo_C]] == 1 & 
    rowSums(binaryS42) == 2, ])

SM42 <- c("ZOTU11875") 
taxaSM42 <- taxa_table[SM42, ]


#S52
S52 <- subset_samples(database.pso, S52 == "1")
pS52 <- merge_samples(S52, "FosterS52")
binaryS52 <- as.data.frame(t(otu_table(pS52)))
binaryS52[binaryS52 > 0] <- 1

#Tadpoles - Foster
sharedSTF52 <- rownames(binaryS52[
  binaryS52[[grupo_A]] == 1 & 
    binaryS52[[grupo_B]] == 1 & 
    rowSums(binaryS52) == 2, ])

STF52 <- c("ZOTU827","ZOTU1382","ZOTU20359","ZOTU20061")
taxaSTF52 <- taxa_table[STF52, ]

#Metamorphs - Foster
sharedSMF52 <- rownames(binaryS52[
  binaryS52[[grupo_A]] == 1 & 
    binaryS52[[grupo_C]] == 1 & 
    rowSums(binaryS52) == 2, ])

SMF52 <- c("ZOTU21763","ZOTU6669","ZOTU20759")
taxaSMF52 <- taxa_table[SMF52, ]

#Tadpoles - Biological
sharedST52 <- rownames(binaryS52[
  binaryS52[[grupo_D]] == 1 & 
    binaryS52[[grupo_B]] == 1 & 
    rowSums(binaryS52) == 2, ])

#Metamorphs - Biological
sharedSM52 <- rownames(binaryS52[
  binaryS52[[grupo_D]] == 1 & 
    binaryS52[[grupo_C]] == 1 & 
    rowSums(binaryS52) == 2, ])

SMB52 <- c("ZOTU12334")
taxaSMB52 <- taxa_table[SMB52, ]

# 13.- Relative abundance of shared zOTUs ----

#Trials from normal transport and no transport
#Replace numbers in the code according to the trial number
pS42m <- subset_samples(S42, Sample_type == "Male")
pS42T <- subset_samples(S42, Sample_type == "TP")
pS42T <- merge_samples(pS42T, "Sample_type")
pS42M <- subset_samples(S42, Sample_type == "METAMORPH")
pS42M <- merge_samples(pS42M, "Sample_type")
pSO42t <- merge_phyloseq(pS42m,pS42T)
pSO42m <- merge_phyloseq(pS42m,pS42M)

#Tadpoles
Tr42t <- tax_glom(pSO42t, taxrank = "Genus")
Tr42t <- transform_sample_counts(Tr42t, function(x)100* x / sum(x))
OTUTr42t <- otu_table(Tr42t)
TAXTr42t <- tax_table(Tr42t)[,"Genus"]
AverageTr42t <- as.data.frame(rowMeans(OTUTr42t))
names(AverageTr42t) <- c("Mean")
TableTr42t <- merge(TAXTr42t, AverageTr42t, by=0, all=TRUE)
TableTr42t$Row.names = NULL
TableTr42t <- TableTr42t[order(desc(TableTr42t$Mean)),]
head(TableTr42t)

#Metamorphs
Tr42m <- tax_glom(pSO42m, taxrank = "Genus")
Tr42m <- transform_sample_counts(Tr42m, function(x)100* x / sum(x))
OTUTr42m <- otu_table(Tr42m)
TAXTr42m <- tax_table(Tr42m)[,"Genus"]
AverageTr42m <- as.data.frame(rowMeans(OTUTr42m))
names(AverageTr42m) <- c("Mean")
TableTr42m <- merge(TAXTr42m, AverageTr42m, by=0, all=TRUE)
TableTr42m$Row.names = NULL
TableTr42m <- TableTr42m[order(desc(TableTr42m$Mean)),]
head(TableTr42m)

#Trials for foster transport
#Replace numbers in the code according to the trial number
pS52mb <- subset_samples(S52, FosterS52 == "Biological")
pS52mf <- subset_samples(S52, FosterS52 == "Foster")
pS52T <- subset_samples(S52, Sample_type == "TP")
pS52T <- merge_samples(pS52T, "Sample_type")
pS52M <- subset_samples(S52, Sample_type == "METAMORPH")
pS52M <- merge_samples(pS52M, "Sample_type")
pSO52tb <- merge_phyloseq(pS52mb,pS52T)
pSO52mb <- merge_phyloseq(pS52mb,pS52M)
pSO52tf <- merge_phyloseq(pS52mf,pS52T)
pSO52mf <- merge_phyloseq(pS52mf,pS52M)

#Foster
#Tadpoles
Tr52tf <- tax_glom(pSO52tf, taxrank = "Genus")
Tr52tf <- transform_sample_counts(Tr52tf, function(x)100* x / sum(x))
OTUTr52tf <- otu_table(Tr52tf)
TAXTr52tf <- tax_table(Tr52tf)[,"Genus"]
AverageTr52tf <- as.data.frame(rowMeans(OTUTr52tf))
names(AverageTr52tf) <- c("Mean")
TableTr52tf <- merge(TAXTr52tf, AverageTr52tf, by=0, all=TRUE)
TableTr52tf$Row.names = NULL
TableTr52tf <- TableTr52tf[order(desc(TableTr52tf$Mean)),]
head(TableTr52tf)

#Metamorphs
Tr52mf <- tax_glom(pSO52mf, taxrank = "Genus")
Tr52mf <- transform_sample_counts(Tr52mf, function(x)100* x / sum(x))
OTUTr52mf <- otu_table(Tr52mf)
TAXTr52mf <- tax_table(Tr52mf)[,"Genus"]
AverageTr52mf <- as.data.frame(rowMeans(OTUTr52mf))
names(AverageTr52mf) <- c("Mean")
TableTr52mf <- merge(TAXTr52mf, AverageTr52mf, by=0, all=TRUE)
TableTr52mf$Row.names = NULL
TableTr52mf <- TableTr52mf[order(desc(TableTr52mf$Mean)),]
head(TableTr52mf)

#Biological
#Tadpoles
Tr52tb <- tax_glom(pSO52tb, taxrank = "Genus")
Tr52tb <- transform_sample_counts(Tr52tb, function(x)100* x / sum(x))
OTUTr52tb <- otu_table(Tr52tb)
TAXTr52tb <- tax_table(Tr52tb)[,"Genus"]
AverageTr52tb <- as.data.frame(rowMeans(OTUTr52tb))
names(AverageTr52tb) <- c("Mean")
TableTr52tb <- merge(TAXTr52tb, AverageTr52tb, by=0, all=TRUE)
TableTr52tb$Row.names = NULL
TableTr52tb <- TableTr52tb[order(desc(TableTr52tb$Mean)),]
head(TableTr52tb)

#Metamorphs
Tr52mb <- tax_glom(pSO52mb, taxrank = "Genus")
Tr52mb <- transform_sample_counts(Tr52mb, function(x)100* x / sum(x))
OTUTr52mb <- otu_table(Tr52mb)
TAXTr52mb <- tax_table(Tr52mb)[,"Genus"]
AverageTr52mb <- as.data.frame(rowMeans(OTUTr52mb))
names(AverageTr52mb) <- c("Mean")
TableTr52mb <- merge(TAXTr52mb, AverageTr52mb, by=0, all=TRUE)
TableTr52mb$Row.names = NULL
TableTr52mb <- TableTr52mb[order(desc(TableTr52mb$Mean)),]
head(TableTr52mb)


# 14.- Putative inhibitory zOTUs transferred - models----
#Load data with number of putative Bd inhibitory zOTUs acquired from fathers for all treatments
transf <- read.table("Transfer_by_treatment.txt", h=T)
transf$Treatment <- as.factor(transf$Treatment) 
transf$Offspring <- as.factor(transf$Offspring) 
transf$Trial <- as.factor(transf$Trial) 

#model for testing proportion of transfer by treatment
model5 <- glmmTMB(numer_transferred  ~ Treatment  + Offspring + (1|Trial ),  
                 data = transf, family = nbinom2())
glmmTMB:::Anova.glmmTMB(model5)
simulateResiduals(model5, plot = TRUE)

ggplot(transf, aes(x = Treatment, y = numer_transferred, fill =Offspring)) +
  geom_boxplot(position = position_dodge(width = 0.8), alpha = 0.7) +
  labs(title = "Shared zOTUs pooled by treatment",
       x = "Treatment",
       y = "Median zOTUs per treatment") +
  theme_minimal()

#Model transfer putative inhibitory by treatment
model6 <- glmmTMB(proportion_Bd_transf ~ Treatment + Offspring + (1|Trial), 
                 family = ordbeta(), data = transf)
glmmTMB:::Anova.glmmTMB(model6)
simulateResiduals(model6, plot = TRUE)

#Model transfer putative inhibitory of foster transport
fosterBd   <- subset(transf, Treatment == "Foster_T")

model7 <- glmmTMB(proportion_Bd_transf ~ Source + Offspring + (1|Trial), 
                  family = beta_family(), ziformula = ~ 1, data = fosterBd)
glmmTMB:::Anova.glmmTMB(model7)
simulateResiduals(model7, plot = TRUE)

