#Data analysis workflow
#Xiaoqing Tan
#xvt5028@psu.edu
#The Pennsylvania State University

#Microbiome sequence data analyses (carried out in Mothur v 1.39.5)
#From the raw sequencing files (.fastq)
#Set the inputdirectory and number of processors

#Making.files which include list of samples and the sequences associated with those samples (paired R1, R2)
mothur > make.file(inputdir=., type=fastq, prefix=16s)

#Combine the paired-end reads together
#Extract the sequence and quality score data from fastq files, create the contigs
mothur > make.contigs(file=16s.files)

#Summary.seqs provide descriptive statistics of all sequences in .fasta file ##
mothur > summary.seqs(fasta=16s.trim.contigs.fasta)

#Filter sequences
##Remove any sequences with ambiguous bases ("N"), shorter than 292, longer than 292
mothur > screen.seqs(fasta=16s.trim.contgis.fasta,
                     group=16s.contigs.groups, summary=16s.trim.contigs.summary,
                     minlength=292, maxlenth=292, maxambig=0)

#Remove identical (grouped) sequences; representative sequence will be picked and stored in .fasta and corresponding sequences will be saved as sequence names to reducw computational work
mothur > unique.seqs(fasta=16s.trim.contigs.good.fasta)

#Create a count table of current unique sequences
mothur > count.seqs(name=16s.trim.contigs.good.names,
                    group=16s.contigs.good.groups)

#Analyze the target 16S rRNA region (V4)
#Customize database to target V4 region
#Define start and end positions within a 16S rRNA sequence
#Set keepdots=F to false to remove output trailing dots from fragments
mothur > pcr.seqs(fasta=silva.nr_v132.align, start=11894, end=25319,
                  keepdots=F, processors=8)

#Rename reference files
mothur > rename.file(input = silva.bacteria.pcr.fasta, new=
                       silva.v4.fasta)

#Align the target region to silva reference database
#Ksize can be determined (ksize=); default is 8
#Allow reverse matching to the reference using flip=T
mothur > align.seqs(fasta=16s.trim.contigs.good.unique.fasta,
                    reference=silva.v4.fasta, flip=T)

#Provide descriptive statistics of all sequences in .fasta file using summary.seqs
mothur > summary.seqs(fasta=16s.trim.contigs.good.unique.align,
                      count=16s.trim.contigs.good.count_table)

#Remove sequences that are before or after the sites of alignment from the previous step
mothur > screen.seqs(fasta=16s.trim.contgis.good.unique.align,
                     group=16s.contigs.good.count_table,
                     summary=16s.trim.contigs.good.unique.summary, minlength=Undecided,
                     maxlenth=Undecided, maxhomop=8)

#Remove overhangs
#Remove the alignment characters that only consist of "-", using vertical=T
#Remove the sequences containing '.' by using trump=.
mothur > filter.seqs(fasta=16s.trim.contigs.good.unique.good.align,
                     vertical=T, trump=.)

#Rerun unique.seqs in case new redundant sequences were created by filtering
mothur > unique.seqs(fasta=16s.trim.contigs.good.unique.good.filter.fasta,
                     count = 16s.trim.contigs.good.good.count_table)

#De-noise sequences
#Diffs=2 means threshold of mismatches in the sequence
mothur > pre.cluster(fasta=16s.trim.contigs.good.unique.good.filter.unique.fast
                     a, count=16s.trim.contigs.good.unique.good.filter.count_table,
                     diffs=2)

#Read a fasta and count file to chimera sequences
mothur > chimera.vsearch(fasta=16s.trim.contigs.good.unique.good.filter.unique.p
                         recluster.fasta,
                         count=16s.trim.contigs.good.unique.good.filter.unique.precluster.count
                         _table, dereplicate=t)

#Remove chimera
mothur > remove.seqs(fasta=16s.trim.contigs.good.unique.good.filter.unique.prec
                     luster.fasta,
                     accnos=16s.trim.contigs.good.unique.good.filter.unique.precluster.deno
                     vo.vsearch.accnos)

#Assign taxonomy
mothur > classify.seqs(fasta=16s.trim.contigs.good.unique.good.filter.unique.pr
                       ecluster.pick.fasta,
                       count=16s.trim.contigs.good.unique.good.filter.unique.precluster.denov
                       o.uchime.pick.count_table, reference=silva.nr_v123.align,
                       taxonomy=silva.nr_v123.tax)

#Remove chloroplast and mitochondria
mothur > remove.lineage(fasta=16s.trim.contigs.good.unique.good.filter.unique.p
                        recluster.pick.fasta,
                        count=16s.trim.contigs.good.unique.good.filter.unique.precluster.denov
                        o.uchime.pick.count_table,
                        taxonomy=16s.trim.contigs.good.unique.good.filter.unique.precluster.pi
                        ck.nr_v132.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-
                          Eukaryota)

#Calculate uncorrected pairwise distances between aligned DNA sequences. By default, a gap is only penalized; cutoff value indicates that distances larger than 0.03 (>97% similarity) will not be saved 
mothur > dist.seqs(fasta=16s.trim.contigs.good.unique.good.filter.unique.preclu
                   ster.pick.pick.fasta, cutoff=0.03)

#Assign sequences to OTUs, Mothur provides three different methods of alignment. By default, opticlust method is used.
mothur > cluster(column=16s.trim.contigs.good.unique.good.filter.unique.preclus
                 ter.pick.pick.dist, count =
                   16s.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsea
                 rch.pick.pick.count_table)

#Determine how many sequences are in each OTU at the 0.03 cutoff level. Distribute OTUs into groups
#The output shared file is used as an OTU table
mothur > make.shared(list=16s.trim.contigs.good.unique.good.filter.unique.precl
                     uster.pick.pick.opti_mcc.list,
                     count=16s.trim.contigs.good.unique.good.filter.unique.precluster.denov
                     o.vsearch.pick.pick.count_table, label=0.03)

#Determine taxonomy for all OTUs. Outcome of this command is a taxonomy file
mothur > classify.otu(list=16s.trim.contigs.good.unique.good.filter.unique.prec
                      luster.pick.pick.opti_mcc.list,
                      count=16s.trim.contigs.good.unique.good.filter.unique.precluster.denov
                      o.vsearch.pick.pick.count_table,
                      taxonomy=16s.trim.contigs.good.unique.good.filter.unique.precluster.pi
                      ck.pds.wang.pick.taxonomy, label=0.03)

####Downstream analaysis###

#Plot L.monocytogenes occurrence in three facilities using a bar plot
#Load the file containing data about facility and L.monocytogenes occurrence
lm <- read.csv(file.choose(), header = T)

lmoccurrence <- ggplot(data = lm, aes(x=Facility, y=Number.of.samples,fill=L..monocytogenes))+
  geom_bar(stat = "identity",width = 0.5)+
  geom_text(aes(label=Number.of.samples, size=3), hjust=0.5, vjust=3) +
  theme_bw(base_size = 12)+
  theme(legend.text=element_text(size=15), legend.title= element_text(size=15)) +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=15)) +
  scale_x_discrete(name = "Facility") +
  scale_y_discrete(name = "Number of Samples")

ggsave("lmoccurrence.pdf", plot = lmoccurrence, device="pdf", width=10, height=7, units="in",dpi=600)

#Plot the beta diversity using PCoA
#Obtain required R packages
library(phyloseq)
library(ape)
library(vegan)
library(ggplot2)

#Import data
set.seed(336)
otus <- import_mothur(mothur_shared_file= file.choose())
otus2 <- as.data.frame(otus)
otus.t <- t(otus)
min(rowSums(otus.t))
otus.r <- rrarefy(otus.t,4501)
OTU <- otu_table(otus.r , taxa_are_rows=FALSE)

taxon <- import_mothur(mothur_constaxonomy_file = file.choose())
taxon <- as.data.frame(taxon)
colnames(taxon) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX = tax_table(as.matrix(taxon)) 
metadat <- read.table(file.choose(), sep=",", header=T, row.names=1)

META = sample_data(metadat)

#Plot PCoA for rarefied samples (16S rRNA data)
phyloseq = phyloseq(OTU, TAX, META)  
TREE = rtree(ntaxa(phyloseq), rooted=TRUE, tip.label = taxa_names(phyloseq))  
phyloseq = phyloseq(OTU,TAX,META,TREE)  
phyloseq

ord = ordinate (phyloseq, "PCoA", "unifrac", weighted = TRUE)
po = plot_ordination(phyloseq, ord, color="Facility",shape="L..monocytogenes")
PCOA16S <- po + 
  geom_point(size=4)+theme_classic() +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=15)) +
  scale_x_continuous(name = "PC1 (10.7%)") +
  scale_y_continuous(name = "PC2 (6.6%)")

PCOA16S

#Import data for ITS
set.seed(336)
otus_ITS <- import_mothur(mothur_shared_file= file.choose())
otus2_ITS <- as.data.frame(otus_ITS)
otus.t_ITS <- t(otus_ITS)
min(rowSums(otus.t_ITS))
otus.r_ITS <- rrarefy(otus.t_ITS,5323)
OTU_ITS <- otu_table(otus.r_ITS , taxa_are_rows=FALSE)

taxon_ITS <- import_mothur(mothur_constaxonomy_file = file.choose())
taxon_ITS <- as.data.frame(taxon_ITS)
colnames(taxon_ITS) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
TAX_ITS = tax_table(as.matrix(taxon_ITS)) 
metadat_ITS <- read.table(file.choose(), sep=",", header=T, row.names=1)

META_ITS = sample_data(metadat_ITS)

#Plot PCoA for rarefied samples (ITS data)
phyloseq_ITS = phyloseq(OTU, TAX, META)  
TREE_ITS = rtree(ntaxa(phyloseq_ITS), rooted=TRUE, tip.label = taxa_names(phyloseq_ITS))  
phyloseq_ITS = phyloseq(OTU,TAX,META,TREE_ITS)  
phyloseq_ITS

ord_ITS = ordinate (phyloseq_ITS, "PCoA", "unifrac", weighted = TRUE)
po_ITS = plot_ordination(phyloseq_ITS, ord_ITS, color="Facility",shape="L.monocytogenes")
PCOAITS <- po_ITS + 
  geom_point(size=4)+theme_classic() +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=15)) +
  scale_x_continuous(name = "PC1 (43.1%)") +
  scale_y_continuous(name = "PC2 (20.5%)")

PCOAITS

#Combine the 16S rRNA and ITS PCoA plots
a = plot_grid(PCOA16S + theme(legend.position= "none") , PCOAITS + theme(legend.position = "none") ,
              ncol=1, nrow=2, labels=c("A", "B"), label_size = 20)
b = get_legend(PCOAITS)
c = plot_grid(a, b, ncol=2, rel_widths = c(3,1))
c
ggsave("PCOAcombined.pdf", plot =c, device="pdf", width=10, height=10, units="in",dpi=600)
ggsave("PCOAcombined.png", plot =c, device="png", width=10, height=10, units="in",dpi=600)

#Stack barplot for bacterial and fungal communities at a family level
#Use phyoseq objects for 16S rRNA and ITS

#Melt to long format (for ggploting) 
family_16 <- phyloseq %>%
  tax_glom(taxrank = "Family") %>%                      #agglomerate at a family level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  #Transform to relative abundance
  psmelt() %>%                                          #Melt to long format
  arrange(Family)                                       #Sort data frame alphabetically by phylum

family_ITS <- phyloseq_ITS %>%
  tax_glom(taxrank = "Rank5") %>%                      
  transform_sample_counts(function(x) {x/sum(x)} ) %>%  
  psmelt() %>%                                          
  arrange(Family)     

write.csv(family_16, "combined_family_16.csv")
write.csv(family_ITS, "combined_family_ITS.csv")         #Write the filtered file into csv format

#Filter the 'Abundance' column to 'less than' 0.10 (in Excel)
#All of the outcome rows are representative families with abundance lower than 0.10 
#Change all column labels under 'Family' to "Other" 
#Save as .csv file
#Open file in R

#Figure 1: Facility vs. L. monocytogenes
combined_family_16s <- read.csv(file.choose(), sep=",", header=T, row.names=1)
Family_colors <- c("#FFF5F0","#525252","#CB181D","#99000D","#EF3B2C","#FB6A4A","#FC9272","#FEE0D2","#4292C6","#084594","#FFF5F0","#EF3B2C","#9ECAE1","#C6DBEF","#DEEBF7","#6BAED6","#2171B5","#737373")

#Plot Figure 1
fig1 <- ggplot(combined_family_16s, aes(x = SampleOrder, y = Abundance , fill = Family))  + 
  facet_grid(Facility~Lmono)+
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = Family_colors) +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(panel.background = element_rect(fill="transparent", color =NA),
                                     plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))

fig1
ggsave("figure1.pdf", plot =fig1, device="pdf", width=8, height=5, units="in",dpi=600)


#Figure 2: Facility vs. section
combined_family_16s <- read.csv(file.choose(), sep=",", header=T, row.names=1)
Family_colors <- c("#FFF5F0","#525252","#CB181D","#99000D","#EF3B2C","#FB6A4A","#FC9272","#FEE0D2","#4292C6","#084594","#FFF5F0","#EF3B2C","#9ECAE1","#C6DBEF","#DEEBF7","#6BAED6","#2171B5","#737373")

#Plot Figure 2 
fig2 <- ggplot(combined_family_16s, aes(x = SampleOrder, y = Abundance , fill = Family))  + 
  facet_grid(Facility~Section)+
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = Family_colors) +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(panel.background = element_rect(fill="transparent", color =NA),
                                     plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black",fill=NA))

fig2
ggsave("figure2.pdf", plot =fig1, device="pdf", width=8, height=5, units="in",dpi=600)
ggsave("figure2.png", plot =fig1, device="png", width=8, height=5, units="in",dpi=600)

library(cowplot)
a = plot_grid(fig1 + theme(legend.position= "none") , fig2 + theme(legend.position = "none") ,
              ncol=1, nrow=2, labels=c("A", "B"), label_size = 20)
b = get_legend(fig2)
c = plot_grid(a, b, ncol=3, rel_widths = c(10,1))
c

ggsave("figre_combine.pdf", plot=c, device="pdf", width=10, height=7, units="in", dpi=600)
ggsave("figre_combine.png", plot=c, device="png", width=10, height=7, units="in", dpi=600)


#Make stack barplot for ITS data
combined_family_ITS <- read.csv(file.choose(), sep=",", header=T, row.names=1)
Family_colors <- c("#084594", "#2171B5", "#4292C6","#9ECAE1","#FFF5F0","#C6DBEF", "#DEEBF7","#6BAED6" ,"#99000D","#EF3B2C","#FC9272","#F7FBFF","#FB6A4A", "#FCBBA1", "#FEE0D2","#525252","#737373","#CB181D","#FFDAB9","#E6E6FA")

#Plot facility vs. L.monocytogenes
ITSfacilitystack <- ggplot(combined_family_ITS, aes(x = SampleOrder, y = Abundance, fill = Family))  + facet_grid(Facility~L.monocytogenes)+
  geom_bar(stat = "identity") + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = Family_colors) +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(panel.background = element_rect(fill="transparent", color =NA),
                                     plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))

#Plot Facility vs. Sections
ITSsectionstack <- ggplot(combined_family_ITS, aes(x = SampleOrder, y = Abundance, fill = Family))  + facet_grid(Facility~Section) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = Family_colors) +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_text(size=15)) + 
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1, ncol=1)) +
  ylab("Relative Abundance") + theme(panel.background = element_rect(fill="transparent", color =NA),
                                     plot.background = element_rect(fill="transparent", color =NA)) +
  theme(strip.background= element_blank(), strip.text = element_text(size=15),
        panel.border = element_rect(color="black", fill=NA))  

a_ITS = plot_grid(ITSfacilitystack + theme(legend.position= "none") , ITSsectionstack + theme(legend.position = "none") ,
                  ncol=1, nrow=2, labels=c("A", "B"), label_size = 20)
b_ITS = get_legend(ITSfacilitystack)
c_ITS = plot_grid(a_ITS, b_ITS, ncol=3, rel_widths = c(5,1))
c_ITS

#Save and export the figure
ggsave("ITSstackcombined.pdf", plot=c, device="pdf", width=12, height=10, units="in", dpi=600)

#Making phyloseq object for rarefaction curve before normalization
phyloseq_rare_16s = phyloseq(otu_table(otus.t, taxa_are_rows=FALSE), TAX, META)
phyloseq_rare_ITS = phyloseq(otu_table(otus.t_ITS, taxa_are_rows=FALSE), TAX_ITS, META_ITS)

#Rarefaction curves
rare_16s_apple_plot <- ggrare(phyloseq_rare_16s, step = 100, se= TRUE, color="Facility")
rare_16s_byfacility_plot <- rare_16s_apple_plot + facet_grid(Facility~.) + 
  theme(strip.text.y=element_blank()) +xlab("Number of OTUs") + ylab("Number of unique OTUs") +
  scale_x_continuous(breaks= seq(0,180000, 10000)) + theme(axis.text.x = element_text(size=10, angle=90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf)

rare_ITS_apple_plot <- ggrare(phyloseq_rare_ITS, step = 100, se= TRUE, color="Facility")
rare_ITS_byfacility_plot <- rare_ITS_applot_plot + facet_grid(Facility ~ .) + 
  theme(strip.text.y=element_blank()) +xlab("Number of OTUs") + ylab("Number of unique OTUs") +
  scale_x_continuous(breaks= seq(0,400000, 20000)) + theme(axis.text.x = element_text(size=10, angle=90)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend= -Inf)

rarefig <- plot_grid(rare_16s_byfacility_plot, rare_ITS_byfacility_plot, nrow=1, ncol=2, labels=c("A", "B"), label_size = 20)

ggsave("rarefig.pdf", plot=rarefig, device="pdf", width=11, height=6, units="in", dpi=600)
ggsave("rarefig.png", plot=rarefig, device="png", width=11, height=6, units="in", dpi=600)

#Alpha diversity

alpha <-estimate_richness(phyloseq, measures=c("Shannon", "InvSimpson", "Chao1"))
estimate_richness(phyloseq, split= TRUE, measures=c("Chao1", "Shannon", "InvSimpson"))

#Import 16S rRNA data
alpha_16s <- read.csv(file.choose(), sep = ",", header = T, row.names = 1)
alpha_ITS <- read.csv(file.choose(), sep = ",", header = T, row.names = 1)

#Pairwise.t.test for alpha diversity using Shannon and Inverse Simpson indices
pairwise.t.test(alpha_16s$Shannon, alpha_16s$Facility, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_16s$InvSimpson, alpha_16s$Facility, p.adjust.method = "bonferroni")

#Import ITS data
alpha_ITS <- read.csv(file.choose(), sep = ",", header = T, row.names = 1)

#Pairwise.t.test for alpha diversity using Shannon and Inverse Simpson indices
pairwise.t.test(alpha_ITS$Shannon, alpha_ITS$Facility, p.adjust.method = "bonferroni")
pairwise.t.test(alpha_ITS$InvSimpson, alpha_ITS$Facility, p.adjust.method = "bonferroni")

#Violin plots for alpha diversity
library(ggpubr)

#16S alpha diversity violin plots
alpha_16s <- read.csv(file.choose(), sep = ",", header = T, row.names = 1)
alpha16s1 <- ggviolin(alpha_16s, x = "Facility", y = "Shannon", add = "boxplot", 
                      fill= "Facility" ) +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=15))
alpha16s2 <- ggviolin(alpha_16s, x = "Facility", y = "InvSimpson", add = "boxplot", 
                      fill  = "Facility" ) +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=15))


#ITS alpha diversity violinplots
alpha_ITS <- read.csv(file.choose(), sep = ",", header = T, row.names = 1)
alphaITS1 <- ggviolin(alpha_ITS, x = "Facility", y = "Shannon", add = "boxplot", 
                      fill= "Facility" ) +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=15))
alphaITS2 <- ggviolin(alpha_ITS, x = "Facility", y = "InvSimpson", add = "boxplot", 
                      fill  = "Facility" ) +
  theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13)) +
  theme(axis.title = element_text(size=15))

#Combined 16S and ITS alpha diversity plots

a = plot_grid(alpha16s1 + theme(legend.position= "none") , 
              alpha16s2 + theme(legend.position = "none") , 
              alphaITS1 + theme(legend.position = "none"), 
              alphaITS2 + theme(legend.position = "none"),
              ncol=2, nrow=2, labels=c("A", "B","C","D"), label_size = 20)

ggsave("alphadiversity.pdf", plot =a, device="pdf", width=12, height=10, units="in",dpi=600)
ggsave("alphadiversity.png", plot =a, device="png", width=12, height=10, units="in",dpi=600)

#Pairwise PERMANOVA
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

#Run pairwise PERMANOVA for 16S rRNA data
permanova_data_16s <- data.frame(sample_data(phyloseq))
pairwise_perm_16s_f <- pairwise.adonis(otu_table(phyloseq), permanova_data_16s$Facility)
pairwise_perm_16s_s <- pairwise.adonis(otu_table(phyloseq), permanova_data_16s$Section)

#Export the file in .csv
write.csv(pairwise_perm_16s_f, "microbiome_pairwise_Facility.csv")
write.csv(pairwise_perm_16s_s, "microbiome_pairwise_Section.csv")

#Run pairwise PERMANOVA for ITS data
permanova_data_ITS <- data.frame(sample_data(phyloseq))
pairwise_perm_ITS_f <- pairwise.adonis(otu_table(phyloseq), permanova_data_ITS$Facility)
pairwise_perm_ITS_s <- pairwise.adonis(otu_table(phyloseq), permanova_data_ITS$Section)

#Export the file in .csv
write.csv(pairwise_perm_ITS_f, "mycobiome_pairwise_Facility.csv")
write.csv(pairwise_perm_ITS_s, "mycobiome_pairwise_Section.csv")

#PICRUSt analysis plot
#Import csv file for picrust, actural abundance
picrustfuntion <- read.csv(file.choose(), sep = ",", header = T)

#Make boxplot based on picrust information
allfunctionabun <- ggplot(picrustfuntion, aes(x=Facility, y=Abundance, fill=Pathway.category)) + 
  theme_bw() +geom_boxplot()+
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.title.x = element_text(size=15), 
        axis.title.y=element_text(size=15))

#Make a plot based on function abundance
functionabundance <- ggplot(picrustfuntion, aes(x=Facility, y=Abundance, fill=Pathway.category)) +theme_bw() +geom_boxplot()

#Import .csv file for picurst, relative abundance, all combined
picrustrelabun <- read.csv(file.choose(), sep = ",", header = T)
refunctionabun <- ggplot(picrustrelabun, aes(x=Facility, y=relative.abundance, fill=Category)) +theme_bw() + 
  geom_boxplot() +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.title.x = element_text(size=15), 
        axis.title.y=element_text(size=15))

#Create a plot for functional categories based on relative abundance
re_all_plot <-ggplot(picrustrelabun, aes(x=Facility, y=relative.abundance, fill=Category)) +theme_bw() + geom_boxplot()

#Import csv file for picrust, relative abundance, by category
picrustfuntioncate <- read.csv(file.choose(), sep = ",", header = T)

#Pairwise.t.test for significant difference between categories
pairwise.t.test(picrustfuntioncate$Metabolism, picrustfuntioncate$Facility, p.adjust.method = "bonferroni")
ggplot(picrustfuntioncate, aes(x=Facility, y=Cellular.Processes)) + theme_bw() + geom_col()

#Plot metabolism and environment functional categories for each facility
metabolism <- ggplot(picrustfuntioncate, aes(x=Facility, y=Metabolism)) + theme_bw() + 
  geom_boxplot() +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.title.x = element_text(size=15), 
        axis.title.y=element_text(size=15))

Environment <- ggplot(picrustfuntioncate, aes(x=Facility, y=Environmental.Information.Processing)) + 
  theme_bw() + geom_boxplot() +
  theme(legend.text=element_text(size=13), legend.title= element_text(size=15)) +
  theme(axis.title.x = element_text(size=15), 
        axis.title.y=element_text(size=15))

a = plot_grid(allfunctionabun, refunctionabun, nrow=2, labels=c("A", "B"), label_size = 20)
b =plot_grid(metabolism, Environment,nrow = 2, labels = c("C","D"), label_size = 20)
c = plot_grid(a, b, ncol=2, rel_widths = c(6,3))
c

ggsave("picrust.pdf", plot=c, device="pdf", width=12, height=10, units="in", dpi=600)
ggsave("picrust.png", plot=c, device="png", width=12, height=10, units="in", dpi=600)


