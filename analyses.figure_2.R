#!/usr/bin/Rscript
################################################################################
#Interaction-adjusted beta diversity (community similarity) estimation.
#
#Prepare analyses shown in Figure 2
#
#=>	load HMP data (V35) for all samples, de novo 97% AL OTUs
#=> load community similarities (previously calculated in "prepare.community_similarity_data.R")
#=> perform PERMANOVA using vegan::adonis
#
#
#2015-08-31
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
library("grid")
library("gridExtra")
library("data.table")
library("plyr")
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
PARAM$folder.input <- #PUT INPUT FOLDER HERE
	PARAM$folder.data <- #PUT DATA FOLDER HERE
	PARAM$folder.output <- #PUT RESULTS/OUTPUT FOLDER HERE
	PARAM$file.functions <- #PUT FUNCTIONS TO SOURCE HERE ("funtions.communtiy_similarity.R")
	PARAM$file.sample_list <- #PUT RAW LIST OF SAMPLES FILE HERE
	PARAM$file.sample_metadata <- #PUT RAW SAMPLE METADATA FILE HERE
	PARAM$file.otu_table <- #PUT OTU TABLE FILE HERE
	PARAM$file.tree <- #PUT OTU PHYLO TREE FILE HERE
	PARAM$file.otu_data <- #PUT OTU DATA FILE HERE
	PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$sample.steps <- c(0.01, 0.02, 0.05, 0.1, 0.15, seq(0.2, 0.9, by=0.1));
PARAM$use.cores <- 40;
###################################

###################################
#Import functions
source(PARAM$file.functions)
###################################

###################################
#Set parameters for data processing
###################################
#Ordered factor of body subsite names
order.body_subsite <- factor(c(
	"Stool",
	"Anterior_nares",
	"Left_Antecubital_fossa", "Right_Antecubital_fossa", "Left_Retroauricular_crease", "Right_Retroauricular_crease",
	"Vaginal_introitus", "Mid_vagina", "Posterior_fornix",
	"Saliva", "Subgingival_plaque", "Supragingival_plaque", "Tongue_dorsum", "Buccal_mucosa", "Hard_palate", "Attached_Keratinized_gingiva", "Throat", "Palatine_Tonsils"
));
#Minimum OTU size required across all samples
PARAM$thresh.otu_size <- 3;
#Minimum sample size required
PARAM$thresh.sample_size <- 1000;
#Minimum relative number of sequences retained
PARAM$thresh.retained_seqs <- 0.5;
#Get list of colors to use for body sites
#=> based on the R Color Brewer qualitative palette "Accent"
PARAM$colors.body_site <- c("#33a02c", "#bf5b17", "#386cb0", "#fdc086", "#beaed4");
names(PARAM$colors.body_site) <- c("Airways", "Gastrointestinal_tract", "Oral", "Skin", "Urogenital_tract");
###################################

###################################
#Read list of samples to process
samples.raw <- as.vector(read.delim(PARAM$file.sample_list, header=F, sep="\n"))$V1;
n.sample.raw <- length(PARAM$use.samples);
###################################
#Read sample raw metadata
sample.data.raw <- read.table(PARAM$file.sample_metadata, header=T, sep="\t");
rownames(sample.data.raw) <- paste("SN_", as.character(sample.data.raw$SN), sep="");
################################################################################
################################################################################


################################################################################
################################################################################
#Load data
load(file=paste(PARAM$folder.data, "hmp_samples.otu_table.RData", sep=""));
load(file=paste(PARAM$folder.data, "hmp_samples.beta_div.get_cs.RData", sep = ""));

#Get lists of communtiy similarity methods
all.methods <- sapply(get.cs, function(i) {i$name});
all.calls <- sapply(get.cs, function(i) {i$call});
################################################################################
################################################################################


################################################################################
################################################################################
#Ordination & tests for separability by body (sub)site
#=> PCoA plots
#=> PERMANOVA (vegan::adnois) analysis per body sites and body subsites
################################################################################
#Preallocate results collector matrices
collect.F <- collect.R <- matrix(nrow=length(get.cs), ncol=4, dimnames=list(names.m, c("Body_Site", "Skin", "Urogenital_tract", "Oral")));
#Get combinations of samples
combn.samples <- t(combn(samples, 2));
###################################
#Process body sites
###################################
#Iterate through similarity / distance measures
for (t in seq(1, length(get.cs))) {
	curr.dist <- as.dist(get.cs[[t]]$cs);
	curr.sample_data <- as(sample_data(my.ps), "data.frame");
	
	############################
	#Ordinate on body sites
	############################
	#PCoA
	#=> using the Lingoes correction for negative eigenvalues (see ape::pcoa documentation for details)
	curr.pcoa <- pcoa(get.cs[[t]]$cs, correction="lingoes", rn=rownames(curr.sample_data));
	#Prepare data for plotting
	plot.data <- data.frame(curr.sample_data, Axis.1=curr.pcoa$vectors[, "Axis.1"], Axis.2=curr.pcoa$vectors[, "Axis.2"], Group=curr.sample_data$Body_Site);
	plot.var_explained <- round(100*curr.pcoa$values[1:2, "Rel_corr_eig"]/sum(curr.pcoa$values[, "Rel_corr_eig"]), digits=1);
	plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
	plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
	curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
	#Make plot
	curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", color="Body_Site")) +
		stat_ellipse(level=0.95, segments=101, alpha=0.8) + 
		geom_point(size=5, alpha=0.6) +
		#geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
		scale_color_manual(values = PARAM$colors.body_site) +
		scale_fill_manual(values = PARAM$colors.body_site) +
		xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
		ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
		theme_bw();
	ggsave(plot=curr.plot, height=20, width=20, dpi=300, filename=paste(PARAM$folder.output, "hmp_samples.", get.cs[[t]]$name, ".PCoA.pdf", sep = ""), useDingbats=FALSE);
	############################
	
	############################
	#Perform PERMANOVA (vegan::adonis) analysis on body sites
	tmp.aov <- adonis(curr.dist ~ Body_Site, data=curr.sample_data, parallel=PARAM$use.cores);
	capture.output(tmp.aov, file = paste(PARAM$folder.output, "hmp_samples.", get.cs[[t]]$name, ".PERMANOVA.Body_Site.txt", sep = ""));
	get.cs[[t]]$permanova.body_site <- tmp.aov;
	collect.F[t, "Body_Site"] <- tmp.aov$aov.tab$F.Model[1];
	collect.R[t, "Body_Site"] <- tmp.aov$aov.tab$R2[1];
	############################
	
	############################
	#Plot distributions of intra- vs inter-group distances
	############################
	#Extract distances for all pairs
	combn.dist <- get.cs[[t]]$cs[combn.samples];
	combn.ranks <- rank(combn.dist) / length(combn.dist);
	#Annotate pairs as "intra" or "inter"
	combn.bs <- cbind(sample.data[combn.samples[,1], "Body_Site"], sample.data[combn.samples[,2], "Body_Site"])
	combn.intra <- which(combn.bs[,1] == combn.bs[,2]);
	combn.inter <- which(combn.bs[,1] != combn.bs[,2]);
	combn.label <- character(length(combn.dist));
	combn.label[combn.intra] <- "Intra";
	combn.label[combn.inter] <- "Inter"
	#Plot density histograms of both intra- and inter-group distances
	plot.intra_inter.frame <- data.frame(Label=combn.label, Distance=combn.dist, Rank=combn.ranks);
	dens.plot <- ggplot(plot.intra_inter.frame, aes(x=Distance, fill=Label)) + geom_density(alpha=0.6) + xlim(c(0,1)) + theme_bw();
	ggsave(plot=dens.plot, height=10, width=20, dpi=300, filename=paste(PARAM$folder.output, "hmp_samples.", get.cs[[t]]$name, ".distance_histograms.pdf", sep = ""), useDingbats=FALSE);
	############################
	#Iterate though body sites and export intra- vs inter-group distance histograms
	for (bs in c("Airways", "Gastrointestinal_tract", "Oral", "Skin", "Urogenital_tract")) {
		bs.idx <- which(levels(sample.data$Body_Site) == bs);
		#Get intra- and inter-group distances, focussing on current body site only
		combn.curr_bs <- combn.bs[,1] == bs.idx | combn.bs[,2] == bs.idx;
		combn.intra <-  which(combn.curr_bs & combn.bs[,1] == combn.bs[,2]);
		combn.inter <- which(combn.curr_bs & combn.bs[,1] != combn.bs[,2]);
		combn.label <- rep.int(NA, length(combn.dist));
		combn.label[combn.intra] <- "Intra";
		combn.label[combn.inter] <- "Inter"
		#Plot density histograms of both intra- and inter-group distances
		plot.intra_inter.frame <- data.frame(Label=combn.label[c(combn.intra, combn.inter)], Distance=combn.dist[c(combn.intra, combn.inter)]);
		dens.plot <- ggplot(plot.intra_inter.frame, aes(x=Distance, fill=Label)) + geom_density(alpha=0.6, y = ..scaled..) + xlim(c(0,1)) + theme_bw();
		ggsave(plot=dens.plot, height=10, width=20, dpi=300, filename=paste(PARAM$folder.output, "hmp_samples.", get.cs[[t]]$name, ".distance_histograms.", bs, ".pdf", sep = ""), useDingbats=FALSE);
	}
	############################
}
###################################

###################################
#Process body subsites
###################################
for (bs in c("Skin", "Urogenital_tract", "Oral")) {
	############################
	#Sub-set OTU table
	bs.samples <- rownames(sample.data)[sample.data$Body_Site == bs];
	tmp.ot <- ot[, bs.samples]; bs.ot <- tmp.ot[rowSums(tmp.ot) > 0, ]; rm(tmp.ot);
	bs.sample.data <- sample.data[bs.samples, ];
	bs.otus <- rownames(bs.ot); bs.otu_sizes <- rowSums(bs.ot); bs.sample_sizes <- colSums(bs.ot);
	#Prune tree
	bs.tree <- drop.tip(phy_tree(my.ps), rownames(ot)[ !rownames(ot) %in% bs.otus]);
	#Report
	cat(paste("Body site", bs, "has", ncol(bs.ot), "samples and", nrow(bs.ot), "OTUs =>", Sys.time()), sep="\n");
	############################
	
	############################
	#Calculate S matrices
	############################
	#SparCC
	if ("SparCC" == "SparCC") {
		########################
		#Set parameters
		size.thresh <- 1;
		pseudocount <- 10^-6;
		nblocks <- 400;
		use.cores <- 20;
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
		remove.otus <- rowSums(bs.ot) < size.thresh; keep.otus <- ! remove.otus;
		#o.t <- otu.table[1:1000, ] + pseudocount;
		o.t <- bs.ot[! remove.otus, ] + pseudocount;
		otus <- rownames(bs.ot);
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
		curr.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
		rownames(curr.sparcc) <- colnames(curr.sparcc) <- otus;
		cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
		########################
		
		########################
		#Tidy up
		rm(my.combs, dat.split, mat.T, mat.T.desc, var.t, mat.a, omega);
		########################
	}
	#Calculate derived OTU similarity matrix S
	cat(paste(bs, "calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
	tmp.S <- cor.par(curr.sparcc, method="pearson", use=PARAM$cor.use, use.cores=40);
	cat(paste("Done =>", Sys.time()), sep="\n");
	bs.S[[bs]]$S.sparcc <- 0.5 * (tmp.S + 1);
	#Tidy up
	rm(tmp.S, curr.sparcc);
	############################
	#Phylogenetic distances
	#Calculate cophenetic phylogenetic distances
	curr.cophenetic_tree <- cophenetic(bs.tree);
	#Reorder to match order in OTU table
	curr.cophenetic_tree <- curr.cophenetic_tree[rownames(bs.ot), rownames(bs.ot)];
	#Calculate derived OTU similarity matrix S
	cat(paste(bs, "calculating correlations of cophenetic phylogenetic distances =>", Sys.time()), sep="\n");
	tmp.S <- cor.par(curr.cophenetic_tree, method="pearson", use=PARAM$cor.use, use.cores=40);
	cat(paste("Done =>", Sys.time()), sep="\n");
	bs.S[[bs]]$S.phylo <- 0.5 * (tmp.S + 1);
	#Tidy up
	rm(tmp.S, curr.cophenetic_tree);
	############################
	
	############################
	#Calculate community similarities
	bs.cs[[bs]] <- list(); t <- 1;
	############################
	#Classical indices
	############################
	#Jaccard index, classical, unweighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical"; bs.cs[[bs]][[t]]$call <- "jaccard"; t <- t + 1;
	#Jaccard index, classical, weighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, weighted, no fractions"; bs.cs[[bs]][[t]]$call <- "jaccard.abd"; t <- t + 1;
	#Jaccard index, classical, weighted, fraction-based
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, weighted"; bs.cs[[bs]][[t]]$call <- "jaccard.abd.frac"; t <- t + 1;
	#Jaccard index, classical, weighted, fraction-based, Chao's version
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, Chao fraction-based"; bs.cs[[bs]][[t]]$call <- "jaccard.abd.frac_chao"; t <- t + 1;
	#Jaccard index, classical, weighted, alternative formula
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, alternative formula"; bs.cs[[bs]][[t]]$call <- "jaccard.abd.alt"; t <- t + 1;
	#Jaccard index, classical, Chao' abundance correction
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Jaccard, classical, Chao"; bs.cs[[bs]][[t]]$call <- "jaccard.abd.chao"; t <- t + 1;
	#Bray-Curtis, classical
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Bray-Curtis"; bs.cs[[bs]][[t]]$call <- "bray_curtis"; t <- t + 1;
	#Morisita-Horn, classical
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "Morisita-Horn"; bs.cs[[bs]][[t]]$call <- "morisita_horn"; t <- t + 1;
	#Iterate through indices and calculate pairwise sample similarities
	for (t in seq(1, length(bs.cs[[bs]]))) {
		print(paste(Sys.time(), bs, "starting", bs.cs[[bs]][[t]]$name));
		#Calculate all pairwise distances
		my.cs <- community.similarity.par(bs.ot, distance= bs.cs[[bs]][[t]]$call, use.cores=PARAM$use.cores);
		#Manually correct for rounding errors
		my.cs[my.cs < 0] <- 0;
		#Store current similarities
		bs.cs[[bs]][[t]]$cs <- my.cs;
		#Report
		print(paste(Sys.time(), "Done with", bs.cs[[bs]][[t]]$name));
	}
	############################
	#UniFrac
	############################
	t <- length(bs.cs[[bs]]) + 1;
	#Unweighted UniFrac
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "UniFrac, unweighted"; bs.cs[[bs]][[t]]$call <- "UniFrac_unweighted";
	print(paste(Sys.time(), bs, "starting", bs.cs[[bs]][[t]]$name));
	bs.cs[[bs]][[t]]$cs <- fasterUniFrac(bs.ot, bs.tree, distance=bs.cs[[bs]][[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
	print(paste(Sys.time(), "Done with", bs.cs[[bs]][[t]]$name));
	t <- t + 1;
	#Weighted UniFrac
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "UniFrac, weighted"; bs.cs[[bs]][[t]]$call <- "UniFrac_weighted";
	print(paste(Sys.time(), bs, "starting", bs.cs[[bs]][[t]]$name));
	bs.cs[[bs]][[t]]$cs <- fasterUniFrac(bs.ot, bs.tree, distance=bs.cs[[bs]][[t]]$call, blocksize=1000, use.cores=PARAM$use.cores);
	print(paste(Sys.time(), "Done with", bs.cs[[bs]][[t]]$name));
	t <- t + 1;
	############################
	#Corrected indices
	############################
	t <- start.t <- length(bs.cs[[bs]]) + 1;
	#TINA, unweighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "TINA, unweighted"; bs.cs[[bs]][[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
	#TINA, weighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "TINA, weighted"; bs.cs[[bs]][[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
	#PINA, unweighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "PINA, unweighted"; bs.cs[[bs]][[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
	#PINA, weighted
	bs.cs[[bs]][[t]] <- list(); bs.cs[[bs]][[t]]$name <- "PINA, weighted"; bs.cs[[bs]][[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
	############################
	#Iterate through (corrected) indices and calculate pairwise community similarities
	for (t in seq(start.t, length(bs.cs[[bs]]))) {
		print(paste(Sys.time(), bs, "starting", bs.cs[[bs]][[t]]$name));
		#Calculate all pairwise similarities
		if (bs.cs[[bs]][[t]]$name %in% c("TINA, unweighted", "TINA, weighted")) {
			curr.cs <- community.similarity.corr.par(bs.ot, S=bs.S[[bs]]$S.sparcc, distance=bs.cs[[bs]][[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
		} else {
			curr.cs <- community.similarity.corr.par(bs.ot, S=bs.S[[bs]]$S.phylo, distance=bs.cs[[bs]][[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
		}
		#Manually correct for rounding errors
		curr.cs[curr.cs < 0] <- 0;
		#Store current similarities
		bs.cs[[bs]][[t]]$cs <- curr.cs;
		#Report
		print(paste(Sys.time(), "Done with", bs.cs[[bs]][[t]]$name));
		rm(curr.cs)
	}
	############################
	
	############################
	#Export PCoA ordination plots on subsites
	############################
	for (t in seq(1, length(bs.cs))) {
		#PCoA
		#=> using the Lingoes correction for negative eigenvalues (see ape::pcoa documentation for details)
		curr.pcoa <- pcoa(bs.cs[[bs]][[t]]$cs, rn=rownames(bs.sample.data));
		#Prepare data for plotting
		plot.data <- data.frame(bs.sample.data, Axis.1=curr.pcoa$vectors[, "Axis.1"], Axis.2=curr.pcoa$vectors[, "Axis.2"], Group=bs.sample.data$Body_Subsite);
		plot.var_explained <- round(100*curr.pcoa$values[1:2, "Rel_corr_eig"]/sum(curr.pcoa$values[, "Rel_corr_eig"]), digits=1);
		plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
		plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
		curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
		#Make plot
		curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", color="Body_Subsite")) +
			stat_ellipse(level=0.95, segments=101, alpha=0.8) + 
			geom_point(size=5, alpha=0.8) +
			#geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
			scale_color_brewer(palette="Paired") +
			xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
			ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
			theme_bw();
		ggsave(plot=curr.plot, height=20, width=20, dpi=300, filename=paste(PARAM$folder.output, "hmp_samples.", get.cs[[t]]$name, ".", bs, ".retrained.PCoA.pdf", sep = ""), useDingbats=FALSE);
	}
	############################
	
	############################
	#Calculate PERMANOVA scores, in parallel
	#=> as a baseline for subsequent downsampling experiments
	############################
	registerDoMC(cores=length(bs.cs[[bs]]));
	bs.adonis[[bs]] <- foreach(t = seq(1, length(bs.cs[[bs]]))) %dopar% {curr.dist <- as.dist(bs.cs[[bs]][[t]]$cs); adonis(curr.dist ~ Body_Subsite, data=bs.sample.data)}
	#Store current adonis statistics
	for (t in seq(1, length(bs.cs[[bs]]))) {
		collect.F[t, "bs"] <- bs.adonis[[bs]][[t]]$aov.tab$F.Model[1];
		collect.R[t, "Body_Site"] <- bs.adonis[[bs]][[t]]$aov.tab$R2[1];
	}
	############################
	
	############################
	#Save current progress
	save(bs.S, file=paste(PARAM$folder.data, "hmp_samples.data.body_subsites.S.RData", sep=""))
	save(bs.cs, bs.adonis, file=paste(PARAM$folder.data, "hmp_samples.data.body_subsites.RData", sep=""))
	############################
}
###################################

#Export F and R2 statistic results tables
write.table(collect.F, file=paste(PARAM$folder.output, "hmp_samples.adonis.F_statistic.tsv", sep = ""), sep="\t", col.names=NA);
write.table(collect.R, file=paste(PARAM$folder.output, "hmp_samples.adonis.R2_statistic.tsv", sep = ""), sep="\t", col.names=NA);
################################################################################
################################################################################


################################################################################
################################################################################
#Ordination & tests for separability for selected body (sub)site
#=> Figure 2C
#=> "harder" problems, pairwise separability
#=> PERMANOVA (vegan::adnois) analysis per pair of body (sub)sites
################################################################################
#Generate pairs for testing
tp <- list(); p <- 1;
############################
#All combinations of oral subsites
curr.comb <- combn(as.character(sort(unique(sample.data$Body_Subsite[sample.data$Body_Site == "Oral"]))), 2, simplify = F);
for (i.c in seq(1, length(curr.comb))) {
	comb <- curr.comb[[i.c]];
	tp[[p]] <- list();
	tp[[p]]$comb <- comb;
	tp[[p]]$name <- paste(comb[1], "vs", comb[2], sep="_");
	tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% comb];
	tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
	tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {as.character(i)});
	p <- p + 1;
}
############################
#All combinations of vaginal sites
tp[[p]] <- list();
tp[[p]]$comb <- c("Vaginal_introitus", "Mid_vagina");
tp[[p]]$name <- paste(tp[[p]]$comb[1], "vs", tp[[p]]$comb[2], sep="_");
tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% tp[[p]]$comb];
tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {as.character(i)});
p <- p + 1;
tp[[p]] <- list();
tp[[p]]$comb <- c("Vaginal_introitus", "Posterior_fornix");
tp[[p]]$name <- paste(tp[[p]]$comb[1], "vs", tp[[p]]$comb[2], sep="_");
tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% tp[[p]]$comb];
tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {as.character(i)});
p <- p + 1;
tp[[p]] <- list();
tp[[p]]$comb <- c("Mid_vagina", "Posterior_fornix");
tp[[p]]$name <- paste(tp[[p]]$comb[1], "vs", tp[[p]]$comb[2], sep="_");
tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% tp[[p]]$comb];
tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {as.character(i)});
p <- p + 1;
############################
#Combinations of skin sites
tp[[p]] <- list();
tp[[p]]$comb <- c("Left_Antecubital_fossa", "Right_Antecubital_fossa");
tp[[p]]$name <- paste(tp[[p]]$comb[1], "vs", tp[[p]]$comb[2], sep="_");
tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% tp[[p]]$comb];
tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {as.character(i)});
p <- p + 1;
tp[[p]] <- list();
tp[[p]]$comb <- c("Left_Retroauricular_crease", "Right_Retroauricular_crease");
tp[[p]]$name <- paste(tp[[p]]$comb[1], "vs", tp[[p]]$comb[2], sep="_");
tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% tp[[p]]$comb];
tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {as.character(i)});
p <- p + 1;
tp[[p]] <- list();
tp[[p]]$comb <- c("Left_Antecubital_fossa", "Right_Antecubital_fossa", "Left_Retroauricular_crease", "Right_Retroauricular_crease");
tp[[p]]$name <- "Antecubital_fossa_vs_Retroauricular_crease";
tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% tp[[p]]$comb];
tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {if (i %in% tp[[p]]$comb[1:2]) {"Antecubital_fossa"} else {"Retroauricular_crease"}});
p <- p + 1;
############################
#Airways vs skin subsites
tp[[p]] <- list();
tp[[p]]$comb <- c("Left_Antecubital_fossa", "Right_Antecubital_fossa", "Anterior_nares");
tp[[p]]$name <- "Antecubital_fossa_vs_Anterior_nares";
tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% tp[[p]]$comb];
tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {if (i %in% tp[[p]]$comb[1:2]) {"Antecubital_fossa"} else {"Anterior_nares"}});
p <- p + 1;
tp[[p]] <- list();
tp[[p]]$comb <- c("Left_Retroauricular_crease", "Right_Retroauricular_crease", "Anterior_nares");
tp[[p]]$name <- "Retroauricular_crease_vs_Anterior_nares";
tp[[p]]$samples <- rownames(sample.data)[sample.data$Body_Subsite %in% tp[[p]]$comb];
tp[[p]]$sample.data <- sample.data[tp[[p]]$samples, ];
tp[[p]]$sample.data$Group <- sapply(tp[[p]]$sample.data$Body_Subsite, function(i) {if (i %in% tp[[p]]$comb[1:2]) {"Retroauricular_crease"} else {"Anterior_nares"}});
p <- p + 1;
############################

############################
#Preallocate results collectors
pairwise.F <- pairwise.R2 <- pairwise.P <- matrix(nrow=length(tp), ncol=length(all.methods), dimnames=list(sapply(tp, function(i) {i$name}), all.methods));
#Iterate through pairwise tests
for (p in seq(1, length(tp))) {
	############################
	#Unpack current data
	curr.name <- tp[[p]]$name; curr.samples <- tp[[p]]$samples; curr.sample.data <- tp[[p]]$sample.data;
	cat(paste("Processing", curr.name , "=>", Sys.time()), sep="\n");
	#Subset OTU table
	curr.ot <- ot[rowSums(ot[, curr.samples]) > 0, curr.samples];
	curr.tree <- drop.tip(phy_tree(my.ps), rownames(ot)[ !rownames(ot) %in% rownames(curr.ot)]);
	############################
	
	############################
	#Calculate sub-setted SparCC correlations
	if ("SparCC" == "SparCC") {
		########################
		#Set parameters
		size.thresh <- 1;
		pseudocount <- 10^-6;
		nblocks <- 400;
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
		o.t <- curr.ot + pseudocount;
		otus <- rownames(curr.ot);
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
		curr.sparcc <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
		rownames(curr.sparcc) <- colnames(curr.sparcc) <- otus;
		cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
		########################
		
		########################
		rm(my.combs, dat.split, mat.T, mat.T.desc, var.t, mat.a, omega); gc()
		########################
	}
	#Calculate derived OTU similarity matrix S
	cat(paste("Calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
	tmp.S <- cor.par(curr.sparcc, method="pearson", use=PARAM$cor.use, use.cores=40);
	cat(paste("Done =>", Sys.time()), sep="\n");
	curr.S.sparcc <- 0.5 * (tmp.S + 1);
	#Tidy up
	rm(tmp.S);
	############################
	#Calculate cophenetic distances on current tree
	curr.cophenetic_tree <- cophenetic(curr.tree);
	#Reorder to match order in OTU table
	curr.cophenetic_tree <- curr.cophenetic_tree[rownames(curr.ot), rownames(curr.ot)];
	############################
	#Calculate derived OTU similarity matrix S
	cat(paste("Calculating correlations of cophenetic phylogenetic distances =>", Sys.time()), sep="\n");
	tmp.S <- cor.par(curr.cophenetic_tree, method="pearson", use=PARAM$cor.use, use.cores=40);
	cat(paste("Done =>", Sys.time()), sep="\n");
	curr.S.phylo <- 0.5 * (tmp.S + 1);
	#Tidy up
	rm(tmp.S);
	############################
	
	############################
	#Calculate community similarities
	curr.cs <- list();
	for (t in seq(1, length(all.methods))) {
		m <- all.methods[t];
		t.call <- all.calls[t];
		curr.cs[[t]] <- list();
		cat(paste("Calculating", m , "=>", Sys.time()), sep="\n");
		if (m %in% c(
			"Jaccard, classical",
			"Jaccard, classical, weighted, no fractions",
			"Jaccard, classical, weighted",
			"Jaccard, classical, Chao fraction-based",
			"Jaccard, classical, alternative formula",
			"Jaccard, classical, Chao",
			"Bray-Curtis",
			"Morisita-Horn"
		)) {
			iter.cs <- community.similarity.par(curr.ot, distance=t.call, use.cores=PARAM$use.cores)
		} else if (m %in% c("UniFrac, unweighted", "UniFrac, weighted")) {
			iter.cs <- fasterUniFrac(curr.ot, curr.tree, distance=t.call, blocksize=1000, use.cores=PARAM$use.cores);
		} else if (m %in% c("TINA, unweighted", "TINA, weighted")) {
			iter.cs <- community.similarity.corr.par(curr.ot, S=curr.S.sparcc, distance=t.call, blocksize=1000, use.cores=PARAM$use.cores)
		} else if (m %in% c("PINA, unweighted", "PINA, weighted")) {
			iter.cs <- community.similarity.corr.par(curr.ot, S=curr.S.phylo, distance=t.call, blocksize=1000, use.cores=PARAM$use.cores)
		}
		#Manually correct for rounding errors
		iter.cs[iter.cs < 0] <- 0;
		#Pass results
		curr.cs[[t]]$cs <- iter.cs;
	}
	############################
	
	############################
	#Calculate PERMANOVA statistics on all different community similarities
	curr.adonis <- mclapply(seq(1, length(all.methods)), function (t) {adonis(as.dist(curr.cs[[t]]$cs) ~ Group, data=curr.sample.data)}, mc.cores=length(all.methods));
	#Store results
	for (t in seq(1, length(all.methods))) {
		pairwise.F[curr.name, all.methods[t]] <- curr.adonis[[t]]$aov.tab$F.Model[1];
		pairwise.R2[curr.name, all.methods[t]] <- curr.adonis[[t]]$aov.tab$R2[1];
		pairwise.P[curr.name, all.methods[t]] <- curr.adonis[[t]]$aov.tab[1,"Pr(>F)"];
	}
	############################
	#Save
	write.table(pairwise.F, file=paste(PARAM$folder.output, "hmp_samples.pairwise_adonis.F_statistic.tsv", sep = ""), sep="\t", col.names=NA);
	write.table(pairwise.R2, file=paste(PARAM$folder.output, "hmp_samples.pairwise_adonis.R2_statistic.tsv", sep = ""), sep="\t", col.names=NA);
	write.table(pairwise.P, file=paste(PARAM$folder.output, "hmp_samples.pairwise_adonis.P_statistic.tsv", sep = ""), sep="\t", col.names=NA);
	save(pairwise.F, pairwise.R2, pairwise.P, file=paste(PARAM$folder.output, "hmp_samples.pairwise_adonis.RData", sep = ""));
}
############################
#Visualize data as heatmaps
############################
#F statistic
pairwise.F <- as.matrix(read.table("hmp_samples.pairwise_adonis.F_statistic.tsv", header=T, row.names=1))
curr.data <- pairwise.F;
curr.data <- curr.data[, colnames(curr.data) %in% c("Jaccard, classical", "Jaccard, classical, weighted", "Jaccard, classical, Chao", "Bray-Curtis", "Moristia-Horn", "UniFrac, unweighted", "UniFrac, weighted", "TINA, unweighted", "TINA, weighted", "PINA, unweighted", "PINA, weighted")];
colnames(curr.data) <- c("Jaccard, classical", "Jaccard, weighted", "Jaccard, Chao", "Bray-Curtis", "Morisita-Horn", "UniFrac, unweighted", "UniFrac, weighted", "TINA, unweighted", "TINA, weighted", "PINA, unweighted", "PINA, weighted");
pdf(file = paste(PARAM$folder.output, "hmp_samples.pairwise_adonis.F_statistic.heatmap.pdf", sep = ""), width=20, height=20);
heatmap.2(
	log2(curr.data),
	Rowv = F, Colv = F,
	na.rm = T,
	revC = F,
	col=colorRampPalette(brewer.pal(9, "BuGn")),
	trace="none",
	density.info="none",
	margins=c(20, 20)
)
dev.off();
#R2 statistic
pairwise.R2 <- as.matrix(read.table("hmp_samples.pairwise_adonis.R2_statistic.tsv", header=T, row.names=1))
curr.data <- pairwise.R2;
curr.data <- curr.data[, colnames(curr.data) %in% c("Jaccard, classical", "Jaccard, classical, weighted", "Jaccard, classical, Chao", "Bray-Curtis", "Moristia-Horn", "UniFrac, unweighted", "UniFrac, weighted", "TINA, unweighted", "TINA, weighted", "PINA, unweighted", "PINA, weighted")];
colnames(curr.data) <- c("Jaccard, classical", "Jaccard, weighted", "Jaccard, Chao", "Bray-Curtis", "Morisita-Horn", "UniFrac, unweighted", "UniFrac, weighted", "TINA, unweighted", "TINA, weighted", "PINA, unweighted", "PINA, weighted");
#Oral sample pairs
my.data <- curr.data[1:36,];
pdf(file = paste(PARAM$folder.output, "hmp_samples.pairwise_adonis.R2_statistic.heatmap.oral.pdf", sep = ""), width=20, height=20);
heatmap.2(
	my.data,
	Rowv = F, Colv = F,
	na.rm = T,
	revC = F,
	breaks = seq(0, 1, by=0.001),
	col=colorRampPalette(brewer.pal(9, "PuBu")),
	trace="none",
	density.info="none",
	margins=c(20, 20)
)
dev.off();
#Urogenital sample pairs
my.data <- curr.data[37:39,];
pdf(file = paste(PARAM$folder.output, "hmp_samples.pairwise_adonis.R2_statistic.heatmap.urogenital.pdf", sep = ""), width=20, height=20);
heatmap.2(
	my.data,
	Rowv = F, Colv = F,
	na.rm = T,
	revC = F,
	breaks = seq(0, 0.01, by=0.00001),
	col=colorRampPalette(brewer.pal(9, "PuRd")),
	trace="none",
	density.info="none",
	margins=c(20, 20)
)
dev.off();
#Skin sample pairs
my.data <- curr.data[40:44,];
pdf(file = paste(PARAM$folder.output, "hmp_samples.pairwise_adonis.R2_statistic.heatmap.skin.pdf", sep = ""), width=20, height=20);
heatmap.2(
	my.data,
	Rowv = F, Colv = F,
	na.rm = T,
	revC = F,
	breaks = seq(0, 0.4, by=0.0001),
	col=colorRampPalette(brewer.pal(9, "OrRd")),
	trace="none",
	density.info="none",
	margins=c(20, 20)
)
dev.off();
################################################################################
################################################################################


q()