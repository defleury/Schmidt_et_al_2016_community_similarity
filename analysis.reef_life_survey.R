#!/usr/bin/Rscript
################################################################################
#Ecologically-informed beta diversity (community similarity) estimation.
#Test adapted beta diversity estimators on the Reef Life Survey dataset
#
#=>	load Reef Life Survey Data
#=> filter data down to relevant samples & taxa
#=>	calculate global cooccurrence network
#=>	for these datasets, estimate beta diversity (community similarity) using either "vanilla" or corrected formulas
#
#
#2016-02-15
#sebastian.schmidt@imls.uzh.ch
################################################################################


################################################################################
################################################################################
# Load Packages
library("foreach")
library("doMC")
library("iterators")
library("doParallel")
library("parallel")
library("Matrix")
library("bigmemory")
library("biganalytics")
library("gRbase")
library("gplots")
library("ggplot2")
library("maps")
library("grid")
library("gridExtra")
library("data.table")
library("plyr")
library("igraph")
library("ape")
library("phyloseq")
library("vegan")
library("RColorBrewer")
################################################################################
################################################################################


################################################################################
################################################################################
#Preallocate global data structures
PARAM <- list();
PARAM$file.input <- list();
PARAM$folder.input <- #PUT INPUT FOLDER HERE
PARAM$folder.data <- #PUT DATA FOLDER HERE
PARAM$folder.output <- #PUT RESULTS/OUTPUT FOLDER HERE
PARAM$file.functions <- #PUT FUNCTIONS TO SOURCE HERE ("funtions.communtiy_similarity.R")
PARAM$folder.input <- "/mnt/mnemo3/sebastian/projects/beta_div/1-data/Reef_Life_Survey/";
PARAM$folder.data <- "/mnt/mnemo3/sebastian/projects/beta_div/1-data/"
PARAM$folder.output <- "/mnt/mnemo3/sebastian/projects/beta_div/5-results/"
PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$use.cores <- 20;
###################################

###################################
#Set parameters for data processing
###################################
#Minimum OTU size required across all samples
PARAM$thresh.taxa_size <- 1;
#Minimum sample size required
PARAM$thresh.sample_size <- 100;
#Minimum SparCC threshold for network export
PARAM$thresh.SparCC <- 0.5;
###################################

###################################
#Import functions
source(PARAM$file.functions)

#Define functions for geodesic distance calculation
deg2rad <- function(deg) return(deg*pi/180)
geodesic <- function(long1, lat1, long2, lat2, R=6371) {acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R}
################################################################################
################################################################################


################################################################################
################################################################################
#Read and pre-process data
#=> survey (meta-)data
#=> fish count data
#=> cryptic fish count data
#=> invertebrate count data
#=> habitat quadrant data
#
#Data and description available online: http://reeflifesurvey.com
#Data descriptor in Nature Scientific Data: http://www.nature.com/articles/sdata20147
################################################################################
################################################################################
#Load survey (meta-)data
data.survey <- read.csv(file=paste0(PARAM$folder.input, "Reef_Life_Survey#_Surveys.csv"), row.names=1);

#Load fish count data (M1)
data.m1 <- read.csv(file=paste0(PARAM$folder.input, "Reef_Life_Survey#_Global_reef_fish_dataset.csv"), row.names=1);

#Load cryptic fish count data (M2)
data.m2 <- read.csv(file=paste0(PARAM$folder.input, "Reef_Life_Survey#_Cryptic_Fish.csv"), row.names=1);

#Load invertebrate count data (M2)
data.iv <- read.csv(file=paste0(PARAM$folder.input, "Reef_Life_Survey#_Invertebrates.csv"), row.names=1);

###################################
#Generate integrated taxa count table
#=> taxa (rows) by sites (columns)
###################################
#Get list of unique sites
unique.sites <- sort(unique(data.survey$SiteCode));
unique.surveys <- sort(unique(data.survey$SurveyID));
###################################
#Get table for fish count data (M1)
#=> by surveys, will later be integrated
unique.taxa.m1 <- sort(unique(data.m1$Taxon));
table.m1 <- matrix(data=0, nrow=length(unique.taxa.m1), ncol=length(unique.surveys));
dimnames(table.m1) <- list(unique.taxa.m1, as.character(unique.surveys));
#Populate table
for (t in as.character(unique.taxa.m1)) {
	curr.data <- data.m1[data.m1$Taxon == t, ];
	#Sum up counts across blocks
	curr.unique.surveys <- unique(curr.data$SurveyID);
	table.m1[t, as.character(curr.unique.surveys)] <- sapply(curr.unique.surveys, function(s) {sum(curr.data$Total[curr.data$SurveyID == s])});
}
###################################
#Get table for cryptic fish count data (M2)
#=> by surveys, will later be integrated
unique.taxa.m2 <- sort(unique(data.m2$Taxon));
table.m2 <- matrix(data=0, nrow=length(unique.taxa.m2), ncol=length(unique.surveys));
dimnames(table.m2) <- list(unique.taxa.m2, as.character(unique.surveys));
#Populate table
for (t in as.character(unique.taxa.m2)) {
	curr.data <- data.m2[data.m2$Taxon == t, ];
	#Sum up counts across blocks
	curr.unique.surveys <- unique(curr.data$SurveyID);
	table.m2[t, as.character(curr.unique.surveys)] <- sapply(curr.unique.surveys, function(s) {sum(curr.data$Total[curr.data$SurveyID == s])});
}
###################################
#Get table for invertebrate count data (M2)
#=> by surveys, will later be integrated
unique.taxa.iv <- sort(unique(data.iv$Taxon));
table.iv <- matrix(data=0, nrow=length(unique.taxa.iv), ncol=length(unique.surveys));
dimnames(table.iv) <- list(unique.taxa.iv, as.character(unique.surveys));
#Populate table
for (t in as.character(unique.taxa.iv)) {
	curr.data <- data.iv[data.iv$Taxon == t, ];
	#Sum up counts across blocks
	curr.unique.surveys <- unique(curr.data$SurveyID);
	table.iv[t, as.character(curr.unique.surveys)] <- sapply(curr.unique.surveys, function(s) {sum(curr.data$Total[curr.data$SurveyID == s])});
}
###################################
#Integrate taxa counts across tables
tmp.taxa <- sort(unique(c(rownames(table.m1), rownames(table.m2), rownames(table.iv))));
tmp.table <- matrix(data=0, nrow=length(tmp.taxa), ncol=ncol(table.m1));
dimnames(tmp.table) <- list(tmp.taxa, colnames(table.m1));
for (t in rownames(table.m1)) {tmp.table[t, ] <- tmp.table[t, ] + table.m1[t, ]}
for (t in rownames(table.m2)) {tmp.table[t, ] <- tmp.table[t, ] + table.m2[t, ]}
for (t in rownames(table.iv)) {tmp.table[t, ] <- tmp.table[t, ] + table.iv[t, ]}
#Integrate surveys per site
#=> sum up survey data per site
tmp.table.sites <- matrix(data=0, nrow=nrow(tmp.table), ncol=length(unique.sites));
dimnames(tmp.table.sites) <- list(rownames(tmp.table), unique.sites);
for (s in as.character(unique.sites)) {
	site.idx <- as.character(data.survey$SurveyID[data.survey$SiteCode == s]);
	if (length(site.idx) > 1) {tmp.table.sites[,s] <- rowSums(tmp.table[, site.idx])} else {tmp.table.sites[,s] <- tmp.table[, site.idx]}
}
###################################
#Filter for minimum taxa abundances and sample sizes
use.sites <- colnames(tmp.table.sites)[colSums(tmp.table.sites) > PARAM$thresh.sample_size];
use.taxa <- sort(unique(rownames(tmp.table.sites)[rowSums(tmp.table.sites[, use.sites]) > PARAM$thresh.taxa_size]));
taxa.table <- Matrix(tmp.table.sites[use.taxa, use.sites], sparse=T);
###################################
#Filter survey metadata and aggregate depths by site
data.sample <- data.frame(row.names=use.sites, SiteName=NA, Ecoregion=NA, Realm=NA, CountryRegion=NA, StateArea=NA, Location=NA, SiteLat=NA, SiteLon=NA, Depth_Mean=NA, Depth_SD=NA, No_of_Surveys=NA, Size=colSums(taxa.table));
#Populate metadata table
for (s in use.sites) {
	data.sample[s, "SiteName"] <- as.character(data.survey$SiteName[data.survey$SiteCode == s][1]);
	data.sample[s, "Ecoregion"] <- as.character(data.m1$Ecoregion[data.m1$SiteCode == s][1]);
	data.sample[s, "Realm"] <- as.character(data.m1$Realm[data.m1$SiteCode == s][1]);
	data.sample[s, "CountryRegion"] <- as.character(data.survey$CountryRegion[data.survey$SiteCode == s][1]);
	data.sample[s, "StateArea"] <- as.character(data.survey$StateArea[data.survey$SiteCode == s][1]);
	data.sample[s, "Location"] <- as.character(data.survey$Location[data.survey$SiteCode == s][1]);
	data.sample[s, "SiteLat"] <- data.survey$SiteLat[data.survey$SiteCode == s][1];
	data.sample[s, "SiteLon"] <- data.survey$SiteLon[data.survey$SiteCode == s][1];
	curr.depths <- data.survey$Depth[data.survey$SiteCode == s];
	data.sample[s, "Depth_Mean"] <- mean(curr.depths, na.rm=T);
	data.sample[s, "Depth_SD"] <- sd(curr.depths, na.rm=T);
	data.sample[s, "No_of_Surveys"] <- length(which(data.survey$SiteCode == s));
}
data.sample$SiteName <- as.factor(data.sample$SiteName);
data.sample$Ecoregion <- as.factor(data.sample$Ecoregion);
data.sample$Realm <- as.factor(data.sample$Realm);
data.sample$CountryRegion <- as.factor(data.sample$CountryRegion);
data.sample$StateArea <- as.factor(data.sample$StateArea);
data.sample$Location <- as.factor(data.sample$Location);
###################################
#Filter taxa metadata
data.taxa <- data.frame(row.names=use.taxa, Phylum=NA, Class=NA, Family=NA, Size=rowSums(taxa.table));
data.tax.raw <- rbind(data.m1, data.m2, data.iv);
for (t in use.taxa) {
	data.taxa[t, "Phylum"] <- as.character(data.tax.raw$Phylum[data.tax.raw$Taxon == t][1]);
	data.taxa[t, "Class"] <- as.character(data.tax.raw$Class[data.tax.raw$Taxon == t][1]);
	data.taxa[t, "Family"] <- as.character(data.tax.raw$Family[data.tax.raw$Taxon == t][1]);
}
data.taxa$Phylum <- as.factor(data.taxa$Phylum);
data.taxa$Class <- as.factor(data.taxa$Class);
data.taxa$Family <- as.factor(data.taxa$Family);
#Tidy up
rm(data.tax.raw); gc();
#Save data
save(taxa.table, data.sample, data.taxa, use.taxa, use.sites, file=paste0(PARAM$folder.data, "reef_life_survey.data.RData"));
################################################################################
################################################################################


################################################################################
################################################################################
#Process and analyze data
#=> calculate SparCC network
#=> calculate correlations of SparCC correlations
#=> export SparCC interaction network
#=> calculate pairwise community dissimilarities according to different indices
#=> test data structuring according to various factors
################################################################################
################################################################################
#SparCC
if ("SparCC" == "SparCC") {
	########################
	#SparCC correlation analysis
	#=> following Friedman & Alm, PLOS CB, 2012
	########################
	#Set parameters
	size.thresh <- 1;
	pseudocount <- 10^-6;
	nblocks <- 200;
	use.cores <- 40;
	########################
	
	########################
	#Load packages for parallel processing
	require("foreach");
	require("bigmemory");
	library("doMC", quietly=T);
	#Register cluster
	registerDoMC(cores=use.cores);
	########################
	
	########################
	#Filter OTU table by removing all OTUs that are observed less then <size.thresh> times across all samples
	#=> their correlations to all other OTUs will be (manually) set to 0 later on
	#Add pseudocount to all observed counts (to avoid issues with zeroes in log-space)
	cat(paste("SparCC preallocation steps =>", Sys.time()), sep="\n");
	remove.otus <- rowSums(taxa.table) < size.thresh; keep.otus <- ! remove.otus;
	o.t <- taxa.table[! remove.otus, ] + pseudocount;
	otus <- rownames(o.t);
	n.otu <- length(otus);
	########################
	#Preallocate blocks for parallel processing & Aitchinson's T matrix
	#=> based on https://gist.github.com/bobthecat/5024079
	size.split <- floor(n.otu / nblocks);
	if (size.split < 1) {size.split <- 1}
	my.split <- list(); length(my.split) <- nblocks;
	my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split));
	my.split[[nblocks]] <- (size.split*(nblocks-1)):n.otu;
	dat.split <- mclapply(my.split, function(g) {o.t[g,]}, mc.cores=use.cores);
	#Get combinations of splits
	my.combs <- expand.grid(1:length(my.split), 1:length(my.split));
	my.combs <- t(apply(my.combs, 1, sort));
	my.combs <- unique(my.combs);
	#Preallocate Aitchinson's T matrix as big.matrix ("shared" in memory, so accessible from w/in foreach loop)
	mat.T <- big.matrix(nrow=n.otu, ncol=n.otu, dimnames=list(otus, otus), shared=T);
	mat.T.desc <- describe(mat.T);
	cat(paste("Done with preallocations =>", Sys.time()), sep="\n");
	########################
	
	########################
	#Compute Aitchinson's T matrix
	#=> iterate through each block combination, calculate matrix
	#		between blocks and store them in the preallocated matrix on both
	#		symmetric sides of the diagonal
	cat(paste("Starting parallel T matrix calculation =>", Sys.time()), sep="\n");
	results <- foreach(i = 1:nrow(my.combs)) %dopar% {
		#Get current combination
		curr.comb <- my.combs[i, ];
		#Get current data
		g.1 <- my.split[[curr.comb[1]]];
		g.2 <- my.split[[curr.comb[2]]];
		dat.1 <- dat.split[[curr.comb[1]]];
		dat.2 <- dat.split[[curr.comb[2]]];
		#Get current part of Aitchinson's matrix
		curr.T <- apply(dat.1, 1, function(x) {apply(dat.2, 1, function(y) {var(log(x/y))})});
		#Store
		curr.mat.T <- attach.big.matrix(mat.T.desc);
		curr.mat.T[g.2, g.1] <- curr.T;
		curr.mat.T[g.1, g.2] <- t(curr.T);
		#Return
		TRUE
	}
	cat(paste("Done with parallel T matrix calculation =>", Sys.time()), sep="\n");
	########################
	
	########################
	#Compute component variations ("t_i")
	cat(paste("Computing component variations =>", Sys.time()), sep="\n");
	var.t <- colsum(mat.T);
	cat(paste("Done with component variation calculations =>", Sys.time()), sep="\n");
	#Estimate component variances ("omega_i") from t_i by solving a linear equation system
	cat(paste("Estimating component variances from linear equation system =>", Sys.time()), sep="\n");
	mat.a <- matrix(data=1, nrow=n.otu, ncol=n.otu); diag(mat.a) <- n.otu-1;
	omega <- sqrt(solve(a=mat.a, b=var.t));
	cat(paste("Done with component variation estimation =>", Sys.time()), sep="\n");
	#Estimate pairwise correlations based on these values
	cat(paste("Estimating correlations =>", Sys.time()), sep="\n");
	global.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
	rownames(global.sparcc) <- colnames(global.sparcc) <- otus;
	cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
	########################
	
	########################
	#Tidy up
	save(global.sparcc, file=paste0(PARAM$folder.data, "reef_life_survey.SparCC_global.RData"));
	rm(my.combs, dat.split, mat.T, mat.T.desc, var.t, mat.a, omega);
	gc();
	########################
}
###################################
#Correlation of SparCC correlations
cat(paste("Calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
tmp.S <- cor.par(global.sparcc, method="pearson", use=PARAM$cor.use, use.cores=40);
cat(paste("Done =>", Sys.time()), sep="\n");
S.sparcc <- 0.5 * (tmp.S + 1);
###################################
#Save correlations & tidy up
save(S.sparcc, file=paste(PARAM$folder.data, "reef_life_survey.S_SparCC_global.RData", sep=""));
rm(tmp.S); gc();
################################################################################
################################################################################
#Export SparCC network
###################################
#Filter SparCC correlations
#=> get top 5 interactions per taxon
keep.idx <- apply(global.sparcc, 1, function(t.vec) {order(t.vec, decreasing=T)[1:6]});
curr.sparcc <- matrix(data=0, nrow=nrow(global.sparcc), ncol=ncol(global.sparcc)); dimnames(curr.sparcc) <- dimnames(global.sparcc);
for (t in rownames(curr.sparcc)) {curr.sparcc[t, keep.idx[, t]] <- global.sparcc[t, keep.idx[, t]]}
curr.sparcc[curr.sparcc < PARAM$thresh.SparCC] <- 0;
#Coerce into igraph object
curr.graph <- graph.adjacency(curr.sparcc, mode="undirected", weighted=T, diag=F);
#Remove unconnected vertices
remove.vertices <- which(degree(curr.graph) == 0);
keep.vertices <- which(degree(curr.graph) > 0);
#keep.vertices <- seq(1, vcount(curr.graph));
curr.graph <- delete.vertices(curr.graph, remove.vertices);
############################
#Get some stats and metadata on current graph
curr.graph.dat <- list();
curr.graph.dat$size <- rowSums(taxa.table[use.taxa[keep.vertices],]);
############################
#Color by taxonomy
curr.graph.dat$color.tax <- rep.int("#bdbdbd", length(use.taxa));
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Actinopterygii"] <- "#80b1d3";
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Asteroidea"] <- "#fb8072";
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Crinoidea"] <- "#fdb462";
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Echinoidea"] <- "#ffffb3";
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Elasmobranchii"] <- "#8dd3c7";
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Gastropoda"] <- "#fccde5";
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Holothuroidea"] <- "#bc80bd";
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Malacostraca"] <- "#bebada";
curr.graph.dat$color.tax[data.taxa[use.taxa[keep.vertices], "Class"] == "Mammalia"] <- "#b3de69";
############################
#Compute a layout
curr.layout <- layout.kamada.kawai(curr.graph);
############################
#Export graph, coloured by taxonomy
pdf(file = paste(PARAM$folder.output, "reef_life_survey.full_network.pdf", sep = ""), width=20, height=20, useDingbats=F);
plot(curr.graph, vertex.size=log10(curr.graph.dat$size), vertex.label=NA, vertex.color=curr.graph.dat$color.tax, edge.arrow.size=0, edge.width=0.8, edge.color=adjustcolor("darkgrey", alpha.f=0.7), edge.curved=0.3, layout=curr.layout);
dev.off();
############################
#Rank all taxa by their overall interaction strength
data.taxa$summed_SparCC <- colSums(global.sparcc);
################################################################################
################################################################################
#Calculate pairwise community (dis)similarities
############################
#Preallocate
get.cs <- list(); t <- 1;
############################
#Classical indices
############################
#Jaccard index, classical, unweighted
get.cs[[t]] <- list(); get.cs[[t]]$name <- "jaccard.classical.unweighted"; get.cs[[t]]$call <- "jaccard"; t <- t + 1;
#Jaccard index, classical, weighted, fraction-based
get.cs[[t]] <- list(); get.cs[[t]]$name <- "jaccard.classical.weighted.frac"; get.cs[[t]]$call <- "jaccard.abd.frac"; t <- t + 1;
#Jaccard index, classical, Chao' abundance correction
get.cs[[t]] <- list(); get.cs[[t]]$name <- "jaccard.classical.weighted.chao"; get.cs[[t]]$call <- "jaccard.abd.chao"; t <- t + 1;
#Bray-Curtis, classical
get.cs[[t]] <- list(); get.cs[[t]]$name <- "bray_curtis.classical"; get.cs[[t]]$call <- "bray_curtis"; t <- t + 1;
#Morisita-Horn, classical
get.cs[[t]] <- list(); get.cs[[t]]$name <- "morisita_horn.classical"; get.cs[[t]]$call <- "morisita_horn"; t <- t + 1;
#Iterate through indices and calculate pairwise sample similarities
for (t in seq(1, length(get.cs))) {
	print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
	#Calculate all pairwise distances
	my.cs <- community.similarity.par(taxa.table, distance= get.cs[[t]]$call, use.cores=20);
	#Manually correct for rounding errors
	my.cs[my.cs < 0] <- 0;
	#Store current similarities
	get.cs[[t]]$cs <- my.cs;
	#Report
	print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
}
############################
#Corrected indices
############################
t <- start.t <- length(get.cs) + 1;
#TINA, unweighted
get.cs[[t]] <- list(); get.cs[[t]]$name <- "tina.uw"; get.cs[[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
#TINA, weighted
get.cs[[t]] <- list(); get.cs[[t]]$name <- "tina.w"; get.cs[[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
############################
#Iterate through (corrected) indices and calculate pairwise community similarities
for (t in seq(start.t, length(get.cs))) {
	print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
	#Calculate all pairwise similarities
	my.cs <- community.similarity.corr.par(taxa.table, S=S.sparcc, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
	#Manually correct for rounding errors
	my.cs[my.cs < 0] <- 0;
	#Store current similarities
	rownames(my.cs) <- colnames(my.cs) <- rownames(data.sample);
	get.cs[[t]]$cs <- my.cs;
	#Report
	print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
	rm(my.cs)
}
############################
#Save
save(get.cs, file=paste(PARAM$folder.data, "reef_life_survey.get_cs.RData", sep=""));
################################################################################
################################################################################
#Ordination & data structuring by factors
#=> PCoA plots & PERMANOVA analyses
#=> test effect of ecoregion, realm
############################
#Preallocate
all.methods <- unlist(lapply(get.cs, function(i) i$name));
collect.F <- collect.R2 <- collect.P <- matrix(nrow=length(all.methods), ncol=2, dimnames=list(all.methods, c("Realm", "Ecoregion")));
#Select samples in 8 largest realms for plotting only
use.realms <- names(sort(table(data.sample$Realm, useNA="no"), decreasing=T))[1:8];
curr.samples <- rownames(data.sample)[data.sample$Realm %in% use.realms];
curr.data.sample <- data.sample[curr.samples, ];
############################
#Iterate through community similarity metrics
for (t in seq(1, length(get.cs))) {
	curr.cs <- get.cs[[t]]$cs[curr.samples, curr.samples];
	#PCoA on water layer and region
	curr.pcoa <- pcoa(curr.cs, rn=rownames(curr.data.sample));
	#curr.pcoa <- ordinate(my.ps, method = "PCoA", distance=as.dist(curr.cs));
	############################
	#Ordinate by realm
	plot.data <- data.frame(curr.data.sample, Axis.1=curr.pcoa$vectors[, "Axis.1"], Axis.2=curr.pcoa$vectors[, "Axis.2"], Group=curr.data.sample$Realm);
	plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"]/sum(curr.pcoa$values[curr.pcoa$values[, "Relative_eig"] > 0, "Relative_eig"]), digits=1);
	plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
	plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
	curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
	#Make plot
	curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", color="Realm")) +
		geom_point(size=5) +
		geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
		geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.75) +
		geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group), shape=15, size=8) +
		scale_color_brewer(type="qual", palette="Accent") +
		scale_fill_brewer(type="qual", palette="Accent") +
		xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
		ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep=""));
	ggsave(plot=curr.plot,height=20,width=20,dpi=300, filename=paste0(PARAM$folder.output, "reef_life_survey.", get.cs[[t]]$name, ".PCoA.pdf"), useDingbats=FALSE);
	#Calculate PERMANOVA
	curr.adonis <- adonis(curr.cs ~ Realm + Ecoregion, data=curr.data.sample, parallel=PARAM$use.cores);
	#Store
	get.cs[[t]]$permanova <- curr.adonis;
	#F values
	collect.F[t, "Realm"] <- curr.adonis$aov.tab$F.Model[1];
	collect.F[t, "Ecoregion"] <- curr.adonis$aov.tab$F.Model[2];
	#R2 values
	collect.R2[t, "Realm"] <- curr.adonis$aov.tab$R2[1];
	collect.R2[t, "Ecoregion"] <- curr.adonis$aov.tab$R2[2];
	#P value
	collect.P[t, "Realm"] <- curr.adonis$aov.tab[1, "Pr(>F)"];
	collect.P[t, "Ecoregion"] <- curr.adonis$aov.tab[2, "Pr(>F)"];
}
write.table(collect.F, file=paste0(PARAM$folder.output, "reef_life_survey.adonis.F.tsv"), col.names=NA, sep="\t");
write.table(collect.R2, file=paste0(PARAM$folder.output, "reef_life_survey.adonis.R2.tsv", sep = ""), col.names=NA, sep="\t");
write.table(collect.P, file=paste0(PARAM$folder.output, "reef_life_survey.adonis.P.tsv", sep = ""), col.names=NA, sep="\t");
################################################################################
################################################################################
#Biogeographic analysis
#=> correlation of community dissimilarity to geographic distance and latitudinal distance
############################
#Plot all samples on a map of the world
world_map <- map_data("world");
#Empty base plot
p <- ggplot() + coord_fixed();
#Add map to base plot
base_world <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), fill="gray");
#Add data to world map plot
curr.plot <- base_world + geom_point(data=data.sample, aes(x=SiteLon, y=SiteLat, color=Realm), size=4);
ggsave(plot=curr.plot, height=20,width=40,dpi=300, filename=paste0(PARAM$folder.output, "reef_life_survey.all_samples.world_map.pdf"), useDingbats=FALSE);
############################
#Calculate pairwise geodesic distances between samples
geodesic.list <- mclapply(seq(1, nrow(data.sample)), function(r.1) {s.1 <- data.sample[r.1,]; lat.1 <- unname(deg2rad(as.numeric(s.1["SiteLat"]))); lon.1 <- unname(deg2rad(as.numeric(s.1["SiteLon"]))); unlist(lapply(seq(1, nrow(data.sample)), function(r.2) {s.2 <- data.sample[r.2,]; lat.2 <- unname(deg2rad(as.numeric(s.2["SiteLat"]))); lon.2 <- unname(deg2rad(as.numeric(s.2["SiteLon"]))); geodesic(lon.1, lat.1, lon.2, lat.2)}))}, mc.cores=PARAM$use.cores);
geodesic.dist <- do.call(rbind, geodesic.list);
geodesic.dist.num <- as.numeric(as.dist(geodesic.dist));
#Get pairwise delta(latitude)
latitude.list <- lapply(data.sample$SiteLat, function(l.1) {abs(l.1 - data.sample$SiteLat)});
latitude.dist <- do.call(rbind, latitude.list);
latitude.dist.num <- as.numeric(as.dist(latitude.dist));
dimnames(geodesic.dist) <- dimnames(latitude.dist) <- list(use.sites, use.sites);
#Preallocate correlations collector
cor.dist <- cor.lat <- matrix(data=NA, nrow=length(get.cs), ncol=10, dimnames=list(all.methods, c("All", "Western Indo-Pacific", "Central Indo-Pacific", "Eastern Indo-Pacific", "Temperate Australasia", "Temperate Northern Atlantic", "Temperate Northern Pacific", "Temperate South America", "Tropical Atlantic", "Tropical Eastern Pacific")));
sample.idx <- sample(length(geodesic.dist.num), 10000);
#Preallocate plot collection lists
plots.geodesic <- plots.latitude <- list();
jj <- 1;
#Iterate through community similarity metrics
for (t in seq(1, length(get.cs))) {
	curr.cs <- get.cs[[t]]$cs;
	curr.dist <- as.numeric(as.dist(curr.cs));
	#Get correlation across all samples
	cor.dist[t, "All"] <- cor(curr.dist, geodesic.dist.num, method="spearman");
	cor.lat[t, "All"] <- cor(curr.dist, latitude.dist.num, method="spearman");
	#Plot against geodesic distances
	plot.data <- data.frame(Geographic_Distance=geodesic.dist.num[sample.idx], Community_Distance=curr.dist[sample.idx]);
	plots.geodesic[[jj]] <- ggplot(plot.data, aes(y=Geographic_Distance, x=Community_Distance)) + geom_point(size=3, alpha=0.4) + geom_smooth(method = "lm", se=FALSE, formula = y ~ x) + ggtitle(paste(all.methods[t], "All,", "Cor", cor.dist[t, "All"])) + xlim(c(0,1)) + ylim(c(0, 20000)) + theme_bw();
	#Plot against delta(lat)
	plot.data <- data.frame(Latitude_Distance=latitude.dist.num[sample.idx], Community_Distance=curr.dist[sample.idx]);
	plots.latitude[[jj]] <- ggplot(plot.data, aes(y=Latitude_Distance, x=Community_Distance)) + geom_point(size=3, alpha=0.4) + geom_smooth(method = "lm", se=FALSE, formula = y ~ x) + ggtitle(paste(all.methods[t], "All,", "Cor", cor.lat[t, "All"])) + xlim(c(0,1)) + ylim(c(0, 150)) + theme_bw();
	#Increment plot counter
	jj <- jj + 1;
	#Iterate through relevant realms
	for (r in c("Western Indo-Pacific", "Central Indo-Pacific", "Eastern Indo-Pacific", "Temperate Australasia", "Temperate Northern Atlantic", "Temperate Northern Pacific", "Temperate South America", "Tropical Atlantic", "Tropical Eastern Pacific")) {
		curr.samples <- rownames(data.sample)[data.sample$Realm == r & !is.na(data.sample$Realm)];
		curr.dist.r <- as.numeric(as.dist(curr.cs[curr.samples, curr.samples]));
		curr.geodesic <- as.numeric(as.dist(geodesic.dist[curr.samples, curr.samples]));
		curr.latitude <- as.numeric(as.dist(latitude.dist[curr.samples, curr.samples]));
		#Get current correlations across samples
		cor.dist[t, r] <- cor(curr.dist.r, curr.geodesic, method="spearman");
		cor.lat[t, r] <- cor(curr.dist.r, curr.latitude, method="spearman");
		#Get data for plotting
		if (length(curr.dist.r) > 10000) {
			s.idx <- sample(length(curr.dist.r), 10000);
			plot.data.geo <- data.frame(Geographic_Distance=curr.geodesic[s.idx], Community_Distance=curr.dist.r[s.idx]);
			plot.data.lat <- data.frame(Latitude_Distance=curr.latitude[s.idx], Community_Distance=curr.dist.r[s.idx]);
		} else {
			plot.data.geo <- data.frame(Geographic_Distance=curr.geodesic, Community_Distance=curr.dist.r);
			plot.data.lat <- data.frame(Latitude_Distance=curr.latitude, Community_Distance=curr.dist.r);
		}
		#Generate plots
		plots.geodesic[[jj]] <- ggplot(plot.data.geo, aes(y=Geographic_Distance, x=Community_Distance)) + geom_point(size=3, alpha=0.4) + geom_smooth(method = "lm", se=FALSE, formula = y ~ x) + ggtitle(paste(all.methods[t], r, "Cor", cor.dist[t, r])) + xlim(c(0,1)) + ylim(c(0, 20000)) + theme_bw();
		plots.latitude[[jj]] <- ggplot(plot.data.lat, aes(y=Latitude_Distance, x=Community_Distance)) + geom_point(size=3, alpha=0.4) + geom_smooth(method = "lm", se=FALSE, formula = y ~ x) + ggtitle(paste(all.methods[t], r, "Cor", cor.dist[t, r])) + xlim(c(0,1)) + ylim(c(0, 60)) + theme_bw();
		#Increment plot counter
		jj <- jj + 1
	}
}
#Export plots
pdf(file = paste0(PARAM$folder.output, "reef_life_survey.geo.geodesic_distance_correlations.pdf"), width=100, height=70);
do.call(grid.arrange, c(plots.geodesic, nrow=7, ncol=10));
dev.off();
pdf(file = paste0(PARAM$folder.output, "reef_life_survey.geo.latitude_distance_correlations.pdf"), width=100, height=70);
do.call(grid.arrange, c(plots.latitude, nrow=7, ncol=10));
dev.off();
#ggsave(plot=curr.plot,height=20,width=20,dpi=300, filename=paste0(PARAM$folder.output, "reef_life_survey.", get.cs[[t]]$name, ".geography_correlations.pdf"), useDingbats=FALSE);





################################################################################
################################################################################






q()