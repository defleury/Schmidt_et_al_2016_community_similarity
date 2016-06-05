#!/usr/bin/Rscript
################################################################################
#Interaction-adjusted beta diversity (community similarity) estimation.
#
#Test phylogenetic signal in between (urogenital) body subsites
#=> is there a statistically significant correlation of phylogeny and body subsite?
#
#2016-06-01
#sebastian.schmidt@imls.uzh.ch
################################################################################

################################################################################
################################################################################
# Load Packages
library("foreach")
library("doMC")
library("doParallel")
library("parallel")
library("Matrix")
library("gtools")
library("gplots")
library("ggplot2")
library("picante")
library("ggtree")
library("RColorBrewer")
library("vegan")
library("phyloseq")
################################################################################
################################################################################

################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$file.input <- list();
PARAM$folder.data <- "PUT DATA FOLDER HERE";
PARAM$folder.output <- "PUT RESULTS/OUTPUT FOLDER HERE"
PARAM$file.functions <- "PUT FUNCTIONS TO SOURCE HERE (funtions.communtiy_similarity.R)"
PARAM$file.data <- paste0(PARAM$folder.data, "hmp_samples.otu_table.RData");
PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$use.cores <- 40;
################################################################################
################################################################################

################################################################################
################################################################################
#Load data
load(PARAM$file.data);
#Prune data to subselect only urogenital samples
curr.samples <- rownames(sample.data)[sample.data$Body_Site == "Urogenital_tract"];
curr.sample.data <- sample.data[curr.samples,];
curr.ot <- ot[rowSums(ot[,curr.samples]) > 50, curr.samples];
curr.otus <- rownames(curr.ot);
curr.tree <- drop.tip(phy_tree(my.ps), rownames(ot)[! rownames(ot) %in% curr.otus]);
###################################
#Calculate relative penetrance of each OTU per body subsites
#=> in how many samples per body subsite was a given OTU observed?
curr.otu_penetrance <- matrix(nrow=length(curr.otus), ncol=3, dimnames=list(curr.otus, c("Vaginal_introitus", "Mid_vagina", "Posterior_fornix")));
for (s in colnames(curr.otu_penetrance)) {
	s.smpl <- rownames(curr.sample.data)[curr.sample.data$Body_Subsite == s];
	curr.otu_penetrance[, s] <- apply(curr.ot[, s.smpl], 1, function(o) {length(which(o > 0)) / length(s.smpl)});
}
#Calculate distribution of (relative) abundance of OTUs per body subsite
curr.otu_rel_abd <- matrix(nrow=length(curr.otus), ncol=3, dimnames=list(curr.otus, c("Vaginal_introitus", "Mid_vagina", "Posterior_fornix")));
curr.ot.rel <- t(t(curr.ot) / colSums(curr.ot));
for (s in colnames(curr.otu_rel_abd)) {
	s.smpl <- rownames(curr.sample.data)[curr.sample.data$Body_Subsite == s];
	curr.otu_rel_abd[, s] <- rowSums(curr.ot.rel[,s.smpl]) / rowSums(curr.ot.rel);
}
#Calculate phylogenetic signal of body subsite on subtree
curr.test.penetrance <- multiPhylosignal(curr.otu_penetrance, curr.tree);
curr.test.rel_abd <- multiPhylosignal(curr.otu_rel_abd, curr.tree);
###################################

###################################
#Cluster urogenital samples and plot heatmaps next to tree
#=> within each subsite, to determine order in subsequent heatmap plotting
#=> by Jaccard index
###################################
#Plot tree
curr.phylo_plot <- ggtree(curr.tree, layout="rectangular");
ggsave(curr.phylo_plot, width=10, height=20, filename=paste(PARAM$folder.output, "phylogenetic_signal.urogenital_subsites.pdf", sep = ""), useDingbats=F);
###################################
#Vaginal introitus
s.samples <- curr.samples[curr.sample.data$Body_Subsite == "Vaginal_introitus"];
s.ot <- curr.ot[, s.samples];
#Calculate pairwise Jaccard similarities
s.jac <- community.similarity.par(s.ot, distance="jaccard", use.cores=4);
#Cluster
s.clust <- hclust(as.dist(s.jac), method="complete");
#Get clustered sample order and pass data
s.ot.ordered <- curr.ot.rel[curr.phylo_plot$data$label[1:nrow(s.ot)], s.samples[s.clust$order]];
#Plot heatmap
pdf(file=paste0(PARAM$folder.output, "phylogenetic_signal.urogenital_subsites.heatmap.Vaginal_introitus.pdf"), width=8, height=40, useDingbats=F);
heatmap.2(
	as.matrix(s.ot.ordered ^ (1/3)),
	Rowv=F, Colv=F, dendrogram="none",
	revC=T,
	col=colorRampPalette(c("white", brewer.pal(9, "Purples")))(1000),
	trace="none", key=T, density.info="none"
);
dev.off();
###################################
#Mid vagina
s.samples <- curr.samples[curr.sample.data$Body_Subsite == "Mid_vagina"];
s.ot <- curr.ot[, s.samples];
#Calculate pairwise Jaccard similarities
s.jac <- community.similarity.par(s.ot, distance="jaccard", use.cores=4);
#Cluster
s.clust <- hclust(as.dist(s.jac), method="complete");
#Get clustered sample order and pass data
s.ot.ordered <- curr.ot.rel[curr.phylo_plot$data$label[1:nrow(s.ot)], s.samples[s.clust$order]];
#Plot heatmap
pdf(file=paste0(PARAM$folder.output, "phylogenetic_signal.urogenital_subsites.heatmap.Mid_vagina.pdf"), width=8, height=40, useDingbats=F);
heatmap.2(
	as.matrix(s.ot.ordered ^ (1/3)),
	Rowv=F, Colv=F, dendrogram="none",
	revC=T,
	col=colorRampPalette(c("white", brewer.pal(9, "Greens")))(1000),
	trace="none", key=T, density.info="none"
);
dev.off();
###################################
#Posterior fornix
s.samples <- curr.samples[curr.sample.data$Body_Subsite == "Posterior_fornix"];
s.ot <- curr.ot[, s.samples];
#Calculate pairwise Jaccard similarities
s.jac <- community.similarity.par(s.ot, distance="jaccard", use.cores=4);
#Cluster
s.clust <- hclust(as.dist(s.jac), method="complete");
#Get clustered sample order and pass data
s.ot.ordered <- curr.ot.rel[curr.phylo_plot$data$label[1:nrow(s.ot)], s.samples[s.clust$order]];
#Plot heatmap
pdf(file=paste0(PARAM$folder.output, "phylogenetic_signal.urogenital_subsites.heatmap.Posterior_fornix.pdf"), width=8, height=40, useDingbats=F);
heatmap.2(
	as.matrix(s.ot.ordered ^ (1/3)),
	Rowv=F, Colv=F, dendrogram="none",
	revC=T,
	col=colorRampPalette(c("white", brewer.pal(9, "Oranges")))(1000),
	trace="none", key=T, density.info="none"
);
dev.off();
###################################
#Plot circular tree with OTU occurrence heatmap "halo"
curr.phylo_plot <- ggtree(curr.tree, layout="circular");
curr.plot <- gheatmap(p=curr.phylo_plot, data=as.data.frame(curr.otu_rel_abd), width=0.2, low="white", high="#3f007d");
ggsave(curr.plot, width=40, height=40, filename=paste(PARAM$folder.output, "phylogenetic_signal.urogenital_subsites.otu_rel_abd.pdf", sep = ""), useDingbats=F);
################################################################################
################################################################################






q()







