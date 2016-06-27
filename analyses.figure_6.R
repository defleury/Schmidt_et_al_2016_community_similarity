#!/usr/bin/Rscript
################################################################################
#Interaction-adjusted beta diversity (community similarity) estimation.
#
#Simulation study to test TINA index performance
#=> Lotka-Volterra modeling of community dynamics
#=> with or without added interactions
#=> how is TINA performance influenced by different parameters?
#
#This script relies on code kindly provided by Stefanie Widder and David Berry, 
#of the university of Vienna, first used in Berry & Widder, Front Microbiol 2014
#(http://dx.doi.org/10.3389/fmicb.2014.00219). The respective code snippets are
#marked throughout the script.
#
#
#2016-06-27
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
library("reshape2")
library("plyr")
library("gtools")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("vegan")
library("ape")
library("deSolve")
library("igraph")
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
PARAM$bin.microbedmm <- "/home/sebastian/bin/MicrobeDMMv1.0/DirichletMixtureGHPFit";
PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$use.cores <- 40;
###################################

###################################
#Set parameters for data processing
###################################
#Number of samples in simulation
#=> use an even number
PARAM$n.samples <- 50;
#Number of taxa in simulations
#=> use an even number
PARAM$n.taxa <- 200;
#Number of simulation iterations to go through
PARAM$n.iter <- 20;
#Size of simulated samples
PARAM$size.sample <- 10000;
#Modeling time parameters
PARAM$time.end <- 500;
PARAM$time.steps <- 5000;
#Min and max F values for plotting
PARAM$F.lim <- c(0.1, 10^7);
###################################
#Source functions
source(PARAM$file.functions);
################################################################################
################################################################################


################################################################################
################################################################################
#Define functions for Lotka-Volterra community dynamics
#=> code modified from Berry & Widder, 2014, Front Microbiol
################################################################################
################################################################################
#Multispecies Lotka-Volterra model
lvm <- function(t,x,parms){
	with(as.list(parms,c(x)),{
		x[x<0.001]<-0 # If abundance is less than 0.001 round to zero
		dx <- ((r*x) * (1 - (a %*% x)/k - (sum(x) / kg)))
		list(dx)
	})
}

#Numerical integration of mLVM
n.integrate <- function(time=time,init.x=init.x,model=model, parms=parms){
	t.out <- seq(time$start,time$end,length=time$steps)
	as.data.frame(lsoda(init.x,t.out,lvm,parms, maxsteps=10000))
}
################################################################################
################################################################################


################################################################################
################################################################################
#Fit Dirichlet multinomial mixture model
#=> for a set of randomly selected samples from two habitats (oral vs gut)
#=> for the top 200 OTUs (by abundance) characterising these two habitats
#=> using the "microbedmm" package by Holmes et al., 2012 (http://dx.doi.org/10.1371/journal.pone.0030126)
###################################
#Load OTU data
load(PARAM$file.data);
#Randomly sub-select gut and oral samples, to a total of PARAM$n.samples
curr.samples <- c(sample(rownames(sample.data)[sample.data$Body_Site == "Oral"], PARAM$n.samples / 2), sample(rownames(sample.data)[sample.data$Body_Site == "Gastrointestinal_tract"], PARAM$n.samples / 2))
#Filter OTU table to current selected samples
tmp.ot <-  ot[rownames(ot)[rowSums(ot[,curr.samples]) > 0], curr.samples];
#Filter OTU table to PARAM$n.taxa OTUs, sampled from top-ranking OTUs by abundance from current selected samples
curr.ot <- tmp.ot[sample(order(rowSums(tmp.ot), decreasing=T)[1:(PARAM$n.taxa * 5)], PARAM$n.taxa), ];
#Export csv-table for Dirichlet multinomial fitting
write.table(as.matrix(curr.ot), file=paste0(PARAM$folder.data, "hmp_samples.otu_table.filtered.csv"), sep=",", col.names=NA, quote=F);
###################################
#Perform Dirichlet mixture model fit
#=> from the command line
#=> choose two classes (k=2) a priori, because samples are drawn from two environments (oral vs gut)
system(paste(PARAM$bin.microbedmm, "-in", paste0(PARAM$folder.data, "hmp_samples.otu_table.filtered.csv"), "-out", paste0(PARAM$folder.data, "hmp_samples.otu_table.dirichlet"), "-k 2"));
#"./DirichletMixtureGHPFit -in [...]/hmp_samples.otu_table.filtered.csv -out [...]/hmp_samples.otu_table.dirichlet -k 2"
###################################
#Read clustering results (which sample belongs to which group?)
dirichlet.groups <- read.csv(file=paste0(PARAM$folder.data, "hmp_samples.otu_table.dirichlet.z"), header=F, skip=1, row.names=1);
#Read Dirichlet mixture pseudo-counts
dirichlet.mixture <- read.csv(file=paste0(PARAM$folder.data, "hmp_samples.otu_table.dirichlet.mixture"), header=F, skip=2, row.names=1);
#Define factor of sample-to-habitat mappings
curr.habitat <- factor(c(rep("Oral", PARAM$n.samples / 2), rep("Gastrointestinal_tract", PARAM$n.samples / 2)));
#Cluster OTUs by abundances across all selected samples to define preferences for oral and gut habitats
#=> assign each OTU to the habitat in which the majority of its abundance (counts) are observed
otu.preferences <- apply(as.matrix(curr.ot), 1, function(v) {
	N.tot <- sum(v);
	N.oral <- sum(v[1:(PARAM$n.samples / 2)]);
	N.gut <- sum(v[((PARAM$n.samples / 2)+1) : PARAM$n.samples]);
	if (N.oral > N.gut) {return("oral")} else if (N.oral == N.gut) {return(NA)} else {return("gut")}
});
oral.otus <- names(otu.preferences)[otu.preferences == "oral"]; n.otu.oral <- length(oral.otus);
gut.otus <- names(otu.preferences)[otu.preferences == "gut"]; n.otu.gut <- length(gut.otus);
################################################################################
################################################################################


################################################################################
################################################################################
#Simulate data and test separation by community similarity metrics
#=> setting 1: Sampling from a Dirichlet Multinomial Mixture model (with implicit habitat preferences and taxa co-abundances)
#=> setting 2: Generalized Lotka-Volterra dynamics, no habitat preferences and neutral interactions (negative control setup)
#=> setting 3: Generalized Lotka-Volterra dynamics, with habitat preferences and neutral interactions
#=> setting 4: Generalized Lotka-Volterra dynamics, with habitat preferences and non-neutral interactions
#=> setting 5: Generalized Lotka-Volterra dynamics, no habitat preferences and non-neutral interactions
################################################################################
################################################################################
#Setup 1: Dirichlet sampling
#=> generate count tables by sampling from previously fitted Dirichlet models
#=> test separation between groups by different metrics
###################################
#Preallocate
collect.F.sim_1 <- collect.R.sim_1 <- collect.F.null_oral.sim_1 <- collect.R.null_oral.sim_1 <- collect.F.null_gut.sim_1 <- collect.R.null_gut.sim_1 <- matrix(nrow=6, ncol=PARAM$n.iter, dimnames=list(c("jaccard_classical", "jaccard_weighted", "bray_curtis", "morisita_horn", "TINA_uw", "TINA_w"), paste0("Iteration_", sprintf("%03i", 1:PARAM$n.iter))));
#Iterate simulation runs
for (i in 1:PARAM$n.iter) {
	sim.ot <- matrix(data=NA, nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), curr.samples));
	#Simulate oral and gut samples
	#=> per group, first get probabilities by sampling from Dirichlet alpha vector
	#=> get counts by multinomial sampling given probabilities
	curr.dirichlet <- t(rdirichlet(PARAM$n.samples / 2, dirichlet.mixture[,1]));
	sim.ot[, 1:(PARAM$n.samples / 2)] <- apply(curr.dirichlet, 2, function(prob) {rmultinom(1, PARAM$size.sample, prob)});
	curr.dirichlet <- t(rdirichlet(PARAM$n.samples / 2, dirichlet.mixture[,2]));
	sim.ot[, ((PARAM$n.samples / 2)+1) : PARAM$n.samples] <- apply(curr.dirichlet, 2, function(prob) {rmultinom(1, PARAM$size.sample, prob)});
	#Filter simulated OTU table
	sim.ot <- sim.ot[rowSums(sim.ot) > 0,];
	###################################
	#Calculate SparCC correlations and derived correlation similarities
	sim.sparcc <- sparcc(sim.ot, size.thresh=1, pseudocount=10^-6, nblocks=10, use.cores=40);
	sim.S_sparcc <- 0.5 * (cor(sim.sparcc, method="pearson", use=PARAM$cor.use) + 1);
	###################################
	#Calculate community similarities
	###################################
	t <- 1;
	sim.cs <- list();
	#Jaccard, classical
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, classical";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard", use.cores=20);
	t <- t + 1;
	#Jaccard, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, weighted";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard.abd.frac", use.cores=20)
	t <- t + 1;
	#Bray-Curtis
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Bray-Curtis";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="bray_curtis", use.cores=20)
	t <- t + 1;
	#Morisita-Horn
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Morisita-Horn";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="morisita_horn", use.cores=20)
	t <- t + 1;
	#TINA, unweighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, unweighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.uw.norm", blocksize=100, use.cores=20)
	t <- t + 1;
	#TINA, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, weighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.w.norm", blocksize=100, use.cores=20)
	###################################
	#Test group separation
	for (t in 1:length(sim.cs)) {
		#Habitat separation
		curr.aov <- adonis(as.dist(sim.cs[[t]]$cs) ~ curr.habitat, permutations=2);
		collect.F.sim_1[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.sim_1[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for oral samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[1:(PARAM$n.sample / 2), 1:(PARAM$n.sample / 2)])
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), PARAM$n.sample / 2, replace=T)), permutations=2);
		collect.F.null_oral.sim_1[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_oral.sim_1[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for gastrointestinal samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[((PARAM$n.sample/2)+1):PARAM$n.sample, ((PARAM$n.sample/2)+1):PARAM$n.sample]);
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), PARAM$n.sample / 2, replace=T)), permutations=2);
		collect.F.null_gut.sim_1[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_gut.sim_1[t, i] <- curr.aov$aov.tab$R2[1];
	}
	###################################
	#Export PCoA plot (based on Euclidean distances => i.e., a PCA plot)
	curr.dist <- dist(t(sim.ot^(1/2)));
	#Generate PCoA
	curr.pcoa <- pcoa(as.matrix(curr.dist), rn=colnames(sim.ot));
	#Collect plotting data
	plot.data <- data.frame(
		Group=curr.habitat,
		Axis.1=curr.pcoa$vectors[, "Axis.1"],
		Axis.2=curr.pcoa$vectors[, "Axis.2"]
	);
	#Get percent variation explained
	if (curr.pcoa$correction[1] == "none") {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"]/sum(curr.pcoa$values[, "Relative_eig"]), digits=1)} else {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Rel_corr_eig"]/sum(curr.pcoa$values[, "Rel_corr_eig"]), digits=1)}
	#Get convex hulls around points
	plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
	#Get centroids
	plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
	curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
	#Plot ordination
	curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", colour="Group")) +
		#Add group centroids
		geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group), shape=15, size=8) +
		#Add segments of points to centroids
		geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.75) +
		#Add points (per sample)
		geom_point(size=5, alpha=0.6) +
		#Add ellipses (at 95%)
		#stat_ellipse(level=0.95, segments=101, alpha=0.8) +
		#Add convex hull around all points per group
		geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
		#Add x & y axis labels
		xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
		ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
		#Set overall plot theme
		theme_bw();
	#Save plot
	ggsave(plot=curr.plot, height=10, width=10, dpi=300, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_1.PCA.iteration_", i, ".pdf"), useDingbats=FALSE);
	###################################
	#Report
	if ((i %% 10) == 0) {print(paste(Sys.time(), "=> Done with iteration", i));}
}
###################################
#Export simulation results
write.table(collect.F.sim_1, file=paste0(PARAM$folder.output, "hmp_samples.simulation_1.F.tsv"), sep="\t", col.names=NA);
write.table(collect.R.sim_1, file=paste0(PARAM$folder.output, "hmp_samples.simulation_1.R.tsv"), sep="\t", col.names=NA);
write.table(collect.F.null_oral.sim_1, file=paste0(PARAM$folder.output, "hmp_samples.simulation_1.F.null_oral.tsv"), sep="\t", col.names=NA);
write.table(collect.R.null_oral.sim_1, file=paste0(PARAM$folder.output, "hmp_samples.simulation_1.R.null_oral.tsv"), sep="\t", col.names=NA);
write.table(collect.F.null_gut.sim_1, file=paste0(PARAM$folder.output, "hmp_samples.simulation_1.F.null_gut.tsv"), sep="\t", col.names=NA);
write.table(collect.R.null_gut.sim_1, file=paste0(PARAM$folder.output, "hmp_samples.simulation_1.R.null_gut.tsv"), sep="\t", col.names=NA);

collect.F.sim_1 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_1.F.tsv"), header=T, row.names=1);
collect.R.sim_1 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_1.R.tsv"), header=T, row.names=1);
collect.F.null_oral.sim_1 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_1.F.null_oral.tsv"), header=T, row.names=1);
collect.R.null_oral.sim_1 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_1.R.null_oral.tsv"), header=T, row.names=1);
collect.F.null_gut.sim_1 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_1.F.null_gut.tsv"), header=T, row.names=1);
collect.R.null_gut.sim_1 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_1.R.null_gut.tsv"), header=T, row.names=1);
###################################
#Prepare data for plotting
df.F.sim_1 <- data.frame(t(collect.F.sim_1)); df.F.sim_1$val_type <- "F"; df.F.sim_1$comparison <- "oral_vs_gut";
df.F.null_oral.sim_1 <- data.frame(t(collect.F.null_oral.sim_1)); df.F.null_oral.sim_1$val_type <- "F"; df.F.null_oral.sim_1$comparison <- "oral_null";
df.F.null_gut.sim_1 <- data.frame(t(collect.F.null_gut.sim_1)); df.F.null_gut.sim_1$val_type <- "F"; df.F.null_gut.sim_1$comparison <- "gut_null";
df.R.sim_1 <- data.frame(t(collect.R.sim_1)); df.R.sim_1$val_type <- "R"; df.R.sim_1$comparison <- "oral_vs_gut";
df.R.null_oral.sim_1 <- data.frame(t(collect.R.null_oral.sim_1)); df.R.null_oral.sim_1$val_type <- "R"; df.R.null_oral.sim_1$comparison <- "oral_null";
df.R.null_gut.sim_1 <- data.frame(t(collect.R.null_gut.sim_1)); df.R.null_gut.sim_1$val_type <- "R"; df.R.null_gut.sim_1$comparison <- "gut_null";
#Melt into long format
plot.df <- rbind(
	melt(df.F.sim_1, id.vars=c("val_type", "comparison")),
	melt(df.F.null_oral.sim_1, id.vars=c("val_type", "comparison")),
	melt(df.F.null_gut.sim_1, id.vars=c("val_type", "comparison")),
	melt(df.R.sim_1, id.vars=c("val_type", "comparison")),
	melt(df.R.null_oral.sim_1, id.vars=c("val_type", "comparison")),
	melt(df.R.null_gut.sim_1, id.vars=c("val_type", "comparison"))
);
###################################
#Plot simulation results (F statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "F", ], aes(x=variable, y=value, fill=comparison)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	#geom_point(position=position_jitterdodge(), size=1, alpha=0.4) +
	#Scale log10
	scale_y_log10(breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), limits=PARAM$F.lim) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=10, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_1.F_statistic.pdf"), useDingbats=F);
#Plot simulation results (R2 statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "R", ], aes(x=variable, y=value, fill=comparison)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitterdodge(), size=1, alpha=0.4) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=10, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_1.R2_statistic.pdf"), useDingbats=F);
################################################################################
################################################################################


################################################################################
################################################################################
#Setup 2: Generalized Lotka-Volterra dynamics with neutral (i.e., no) taxa interactions
#=> generate starting counts per taxon by sampling abundances from combined oral and gut samples (sampling without habitat preference)
#=> simulate community dynamics according to an interaction-neutral Lotka-Volterra model
#=> test separation between groups by different metrics
###################################
#Preallocate
collect.F.sim_2 <- collect.R.sim_2 <- collect.F.null_oral.sim_2 <- collect.R.null_oral.sim_2 <- collect.F.null_gut.sim_2 <- collect.R.null_gut.sim_2 <- matrix(nrow=6, ncol=PARAM$n.iter, dimnames=list(c("jaccard_classical", "jaccard_weighted", "bray_curtis", "morisita_horn", "TINA_uw", "TINA_w"), paste0("Iteration_", sprintf("%03i", 1:PARAM$n.iter))));
###################################
#Set general model starting parameters
#=> parameter values modified from Berry & Widder, 2014
###################################
#Modeling time
time <- list(start=0, end=PARAM$time.end, steps=PARAM$time.steps);
#Number of sites to model
sites <- PARAM$n.samples;
#Average number of taxa per site
ot.occ <- apply(curr.ot, 2, function(x) {x[x>0] <- 1; x});
siteRich <- round(mean(colSums(ot.occ)));
#Average proportion of species per site shared with at least one other site
#=> calculated from the number of OTUs and the average richness per site
sharedSp <- floor(nrow(curr.ot) / siteRich);
#Average number of interactions per species
NumInter=0;
###################################
#Iterate through simulation runs
for (i in 1:PARAM$n.iter) {
	start.ot <- sim.ot <- matrix(data=NA, nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), curr.samples));
	###################################
	#Simulate oral and gut samples
	#=> per group, by shuffling OTU counts per row
	#=> normalize and re-assign initial population sizes at 10% of total sample size
	curr.shuffled <- t(apply(as.matrix(curr.ot), 1, sample, size=PARAM$n.samples));
	curr.shuffled.rel <- apply(curr.shuffled, 2, function(x) {x/sum(x)});
	start.ot <- apply(curr.shuffled.rel, 2, function(prob) {rmultinom(1, PARAM$size.sample / 10, prob)});
	dimnames(start.ot) <- dimnames(curr.ot);
	###################################
	#Plot (clustered) count heatmap of start.ot
	log.start.ot <- log10(start.ot[c(oral.otus, gut.otus),]); log.start.ot[is.infinite(log.start.ot)] <- 0;
	log.start.ot[, curr.habitat == "Gastrointestinal_tract"] <- 0 - log.start.ot[, curr.habitat == "Gastrointestinal_tract"];
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_2.start_ot_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(log.start.ot, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(-3, 3, by=0.01), col=colorRampPalette(brewer.pal(9, "PuOr")), trace="none", density.info="none");
	dev.off();
	###################################
	#Get (randomized) parameters for Lotka-Volterra model
	#=> modified from Berry & Widder, 2014
	###################################
	#Assign random growth rates to OTUs
	#=> draw from a uniform distribuion on [0,1]
	r <- runif(PARAM$n.taxa, min=0, max=1);
	#Set all carrying capacities for all OTUs in all samples to total carrying capacity (sample size)
	K <- matrix(data=PARAM$size.sample, nrow=PARAM$n.taxa, ncol=PARAM$n.samples);
	#Generate neutral interaction matrix
	#=> all interactions are 0, except self-competition (set diagonal to 1)
	alpha <- matrix(data=0, nrow=PARAM$n.taxa, ncol=PARAM$n.taxa);
	diag(alpha) <- 1;
	###################################
	#Preallocate results collectors
	sim.ot <- time.end <- rsd <- matrix(nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), paste0("Site_", 1:PARAM$n.samples)));
	#Iterate through sites
	for (s in 1:PARAM$n.samples) {
		k <- K[,s];
		#Simulation code adapted from Berry & Widder, 2014
		parms <- list(r=r, a=alpha, k=k, kg=PARAM$size.sample);
		#Run simulation
		out <- n.integrate(time=time, init.x=start.ot[,s], model=lvm, parms=parms)
		#Store simulation results
		sim.ot[, s] <- as.numeric(out[nrow(out),2:ncol(out)]);
	}
	#Filter results to remove all OTUs that do not occur in any sample
	sim.ot <- sim.ot[rowSums(sim.ot) > 0, ];
	###################################
	#Calculate SparCC correlations and derived correlation similarities
	sim.sparcc <- sparcc(sim.ot, size.thresh=0, pseudocount=10^-6, nblocks=10, use.cores=20);
	sim.S_sparcc <- 0.5 * (cor(sim.sparcc, method="pearson", use=PARAM$cor.use) + 1);
	###################################
	#Calculate community similarities
	###################################
	t <- 1;
	sim.cs <- list();
	#Jaccard, classical
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, classical";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard", use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#Jaccard, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, weighted";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard.abd.frac", use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#Bray-Curtis
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Bray-Curtis";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="bray_curtis", use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#Morisita-Horn
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Morisita-Horn";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="morisita_horn", use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, unweighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, unweighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.uw.norm", blocksize=10, use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, weighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.w.norm", blocksize=10, use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	###################################
	#Test group separation
	for (t in 1:length(sim.cs)) {
		curr.aov <- adonis(as.dist(sim.cs[[t]]$cs) ~ curr.habitat, permutations=2);
		collect.F.sim_2[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.sim_2[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for oral samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[1:(PARAM$n.sample / 2), 1:(PARAM$n.sample / 2)])
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), PARAM$n.sample / 2, replace=T)), permutations=2);
		collect.F.null_oral.sim_2[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_oral.sim_2[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for gastrointestinal samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[((PARAM$n.sample/2)+1):PARAM$n.sample, ((PARAM$n.sample/2)+1):PARAM$n.sample]);
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), PARAM$n.sample / 2, replace=T)), permutations=2);
		collect.F.null_gut.sim_2[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_gut.sim_2[t, i] <- curr.aov$aov.tab$R2[1];
	}
	###################################
	#Export heatmap of current interaction matrix
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_2.alpha_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(alpha, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(-1,1,by=0.01), col=colorRampPalette(rev(brewer.pal(9, "PRGn"))), trace="none", margins=c(20,20), density.info="none");
	dev.off()
	###################################
	#Export PCoA plot (based on Euclidean distances => i.e., a PCA plot)
	curr.dist <- dist(t(sim.ot^(1/2)));
	#Generate PCoA
	curr.pcoa <- pcoa(as.matrix(curr.dist), rn=colnames(sim.ot));
	#Collect plotting data
	plot.data <- data.frame(
		Group=curr.habitat,
		Axis.1=curr.pcoa$vectors[, "Axis.1"],
		Axis.2=curr.pcoa$vectors[, "Axis.2"]
	);
	#Get percent variation explained
	if (curr.pcoa$correction[1] == "none") {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"]/sum(curr.pcoa$values[, "Relative_eig"]), digits=1)} else {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Rel_corr_eig"]/sum(curr.pcoa$values[, "Rel_corr_eig"]), digits=1)}
	#Get convex hulls around points
	plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
	#Get centroids
	plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
	curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
	#Plot ordination
	curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", colour="Group")) +
		#Add group centroids
		geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group), shape=15, size=8) +
		#Add segments of points to centroids
		geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.75) +
		#Add points (per sample)
		geom_point(size=5, alpha=0.6) +
		#Add ellipses (at 95%)
		#stat_ellipse(level=0.95, segments=101, alpha=0.8) +
		#Add convex hull around all points per group
		geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
		#Add x & y axis labels
		xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
		ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
		#Set overall plot theme
		theme_bw();
	#Save plot
	ggsave(plot=curr.plot, height=10, width=10, dpi=300, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_2.PCA.iteration_", i, ".pdf"), useDingbats=FALSE);
	###################################
	#Report
	print(paste(Sys.time(), "=> Done with iteration", i));
}
###################################
#Export simulation results
write.table(collect.F.sim_2, file=paste0(PARAM$folder.output, "hmp_samples.simulation_2.F.tsv"), sep="\t", col.names=NA);
write.table(collect.R.sim_2, file=paste0(PARAM$folder.output, "hmp_samples.simulation_2.R.tsv"), sep="\t", col.names=NA);
write.table(collect.F.null_oral.sim_2, file=paste0(PARAM$folder.output, "hmp_samples.simulation_2.F.null_oral.tsv"), sep="\t", col.names=NA);
write.table(collect.R.null_oral.sim_2, file=paste0(PARAM$folder.output, "hmp_samples.simulation_2.R.null_oral.tsv"), sep="\t", col.names=NA);
write.table(collect.F.null_gut.sim_2, file=paste0(PARAM$folder.output, "hmp_samples.simulation_2.F.null_gut.tsv"), sep="\t", col.names=NA);
write.table(collect.R.null_gut.sim_2, file=paste0(PARAM$folder.output, "hmp_samples.simulation_2.R.null_gut.tsv"), sep="\t", col.names=NA);

collect.F.sim_2 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_2.F.tsv"), header=T, row.names=1);
collect.R.sim_2 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_2.R.tsv"), header=T, row.names=1);
collect.F.null_oral.sim_2 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_2.F.null_oral.tsv"), header=T, row.names=1);
collect.R.null_oral.sim_2 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_2.R.null_oral.tsv"), header=T, row.names=1);
collect.F.null_gut.sim_2 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_2.F.null_gut.tsv"), header=T, row.names=1);
collect.R.null_gut.sim_2 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_2.R.null_gut.tsv"), header=T, row.names=1);
###################################
#Prepare data for plotting
df.F.sim_2 <- data.frame(t(collect.F.sim_2)); df.F.sim_2$val_type <- "F"; df.F.sim_2$comparison <- "oral_vs_gut";
df.F.null_oral.sim_2 <- data.frame(t(collect.F.null_oral.sim_2)); df.F.null_oral.sim_2$val_type <- "F"; df.F.null_oral.sim_2$comparison <- "oral_null";
df.F.null_gut.sim_2 <- data.frame(t(collect.F.null_gut.sim_2)); df.F.null_gut.sim_2$val_type <- "F"; df.F.null_gut.sim_2$comparison <- "gut_null";
df.R.sim_2 <- data.frame(t(collect.R.sim_2)); df.R.sim_2$val_type <- "R"; df.R.sim_2$comparison <- "oral_vs_gut";
df.R.null_oral.sim_2 <- data.frame(t(collect.R.null_oral.sim_2)); df.R.null_oral.sim_2$val_type <- "R"; df.R.null_oral.sim_2$comparison <- "oral_null";
df.R.null_gut.sim_2 <- data.frame(t(collect.R.null_gut.sim_2)); df.R.null_gut.sim_2$val_type <- "R"; df.R.null_gut.sim_2$comparison <- "gut_null";
#Melt into long format
plot.df <- rbind(
	melt(df.F.sim_2, id.vars=c("val_type", "comparison")),
	melt(df.F.null_oral.sim_2, id.vars=c("val_type", "comparison")),
	melt(df.F.null_gut.sim_2, id.vars=c("val_type", "comparison")),
	melt(df.R.sim_2, id.vars=c("val_type", "comparison")),
	melt(df.R.null_oral.sim_2, id.vars=c("val_type", "comparison")),
	melt(df.R.null_gut.sim_2, id.vars=c("val_type", "comparison"))
);
###################################
#Plot simulation results (F statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "F" & plot.df$comparison == "oral_vs_gut", ], aes(x=variable, y=value, fill=variable, colour=variable)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitter(), size=5, alpha=0.4) +
	#Scale log10
	scale_y_log10(breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), limits=PARAM$F.lim) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=8, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_2.F_statistic.pdf"), useDingbats=F);
#Plot simulation results (R2 statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "R", ], aes(x=variable, y=value, fill=comparison)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitterdodge(), size=1, alpha=0.4) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=10, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_2.R2_statistic.pdf"), useDingbats=F);
################################################################################
################################################################################


################################################################################
################################################################################
#Setup 3: Generalized Lotka-Volterra dynamics with habitat preferences and neutral (i.e., no) taxa interactions
#=> generate starting counts per taxon by sampling abundances from oral and gut samples (two groups, with habitat preference)
#=> simulate community dynamics according to an interaction-neutral Lotka-Volterra model
#=> test separation between groups by different metrics
###################################
#Preallocate
collect.F.sim_3 <- collect.R.sim_3 <- collect.F.null_oral.sim_3 <- collect.R.null_oral.sim_3 <- collect.F.null_gut.sim_3 <- collect.R.null_gut.sim_3 <- matrix(nrow=6, ncol=PARAM$n.iter, dimnames=list(c("jaccard_classical", "jaccard_weighted", "bray_curtis", "morisita_horn", "TINA_uw", "TINA_w"), paste0("Iteration_", sprintf("%03i", 1:PARAM$n.iter))));
###################################
#Set general model starting parameters
#=> adapted from code by Berry & Widder, 2014
###################################
#Modeling time
time <- list(start=0, end=PARAM$time.end, steps=PARAM$time.steps);
#Number of sites to model
sites <- PARAM$n.samples;
#Average number of taxa per site
ot.occ <- apply(curr.ot, 2, function(x) {x[x>0] <- 1; x});
siteRich <- round(mean(colSums(ot.occ)));
#Average proportion of species per site shared with at least one other site
#=> calculated from the number of OTUs and the average richness per site
sharedSp <- floor(nrow(curr.ot) / siteRich);
#Average number of interactions per species
NumInter=0;
###################################
#Iterate through simulation runs
for (i in 1:PARAM$n.iter) {
	start.ot <- sim.ot <- matrix(data=NA, nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), curr.samples));
	###################################
	#Simulate oral and gut samples
	#=> per group, by shuffling OTU counts per row
	#=> normalize and re-assign initial population sizes at 10% of total sample size
	curr.shuffled <- cbind(
		t(apply(as.matrix(curr.ot[, 1:(PARAM$n.samples/2)]), 1, sample, size=PARAM$n.samples/2)),
		t(apply(as.matrix(curr.ot[, (PARAM$n.samples/2 + 1):PARAM$n.samples]), 1, sample, size=PARAM$n.samples/2))
	);
	curr.shuffled.rel <- apply(curr.shuffled, 2, function(x) {x/sum(x)});
	start.ot <- apply(curr.shuffled.rel, 2, function(prob) {rmultinom(1, PARAM$size.sample / 10, prob)});
	dimnames(start.ot) <- dimnames(curr.ot);
	###################################
	#Plot (clustered) count heatmap of start.ot
	log.start.ot <- log10(start.ot[c(oral.otus, gut.otus),]); log.start.ot[is.infinite(log.start.ot)] <- 0;
	log.start.ot[, curr.habitat == "Gastrointestinal_tract"] <- 0 - log.start.ot[, curr.habitat == "Gastrointestinal_tract"];
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_3.start_ot_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(log.start.ot, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(-3, 3, by=0.01), col=colorRampPalette(brewer.pal(9, "PuOr")), trace="none", density.info="none");
	dev.off();
	###################################
	#Get (randomized) parameters for Lotka-Volterra model
	#=> adapted from Berry & Widder, 2014
	###################################
	#Assign random growth rates to OTUs
	#=> draw from a uniform distribuion on [0,1]
	r <- runif(PARAM$n.taxa, min=0, max=1);
	#Set all carrying capacities for all OTUs in all samples to total carrying capacity (sample size)
	K <- matrix(data=PARAM$size.sample, nrow=PARAM$n.taxa, ncol=PARAM$n.samples);
	#Generate neutral interaction matrix
	#=> all interactions are 0, except self-competition (set diagonal to 1)
	alpha <- matrix(data=0, nrow=PARAM$n.taxa, ncol=PARAM$n.taxa);
	diag(alpha) <- 1;
	###################################
	#Preallocate results collectors
	sim.ot <- time.end <- rsd <- matrix(nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), paste0("Site_", 1:PARAM$n.samples)));
	#Iterate through sites
	for (s in 1:PARAM$n.samples) {
		k <- K[,s];
		#Simulation code, adapted from Berry & Widder, 2014
		parms <- list(r=r, a=alpha, k=k, kg=PARAM$size.sample);
		#Run simulation
		out <- n.integrate(time=time, init.x=start.ot[,s], model=lvm, parms=parms)
		#Store simulation results
		sim.ot[, s] <- as.numeric(out[nrow(out),2:ncol(out)]);
	}
	#Filter results to remove all OTUs that do not occur in any sample
	sim.ot <- sim.ot[rowSums(sim.ot) > 0, ];
	###################################
	#Calculate SparCC correlations and derived correlation similarities
	sim.sparcc <- sparcc(sim.ot, size.thresh=0, pseudocount=10^-6, nblocks=10, use.cores=20);
	sim.S_sparcc <- 0.5 * (cor(sim.sparcc, method="pearson", use=PARAM$cor.use) + 1);
	###################################
	#Calculate community similarities
	###################################
	t <- 1;
	sim.cs <- list();
	#Jaccard, classical
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, classical";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard", use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#Jaccard, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, weighted";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard.abd.frac", use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#Bray-Curtis
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Bray-Curtis";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="bray_curtis", use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#Morisita-Horn
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Morisita-Horn";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="morisita_horn", use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, unweighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, unweighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.uw.norm", blocksize=10, use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, weighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.w.norm", blocksize=10, use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	###################################
	#Test group separation
	for (t in 1:length(sim.cs)) {
		curr.aov <- adonis(as.dist(sim.cs[[t]]$cs) ~ curr.habitat, permutations=2);
		collect.F.sim_3[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.sim_3[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for oral samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[1:(PARAM$n.sample / 2), 1:(PARAM$n.sample / 2)])
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), PARAM$n.sample / 2, replace=T)), permutations=2);
		collect.F.null_oral.sim_3[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_oral.sim_3[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for gastrointestinal samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[((PARAM$n.sample/2)+1):PARAM$n.sample, ((PARAM$n.sample/2)+1):PARAM$n.sample]);
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), PARAM$n.sample / 2, replace=T)), permutations=2);
		collect.F.null_gut.sim_3[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_gut.sim_3[t, i] <- curr.aov$aov.tab$R2[1];
	}
	###################################
	#Export heatmap of current interaction matrix
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_3.alpha_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(alpha, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(-1,1,by=0.01), col=colorRampPalette(rev(brewer.pal(9, "PRGn"))), trace="none", margins=c(20,20), density.info="none");
	dev.off()
	###################################
	#Export PCoA plot (based on Euclidean distances => i.e., a PCA plot)
	curr.dist <- dist(t(sim.ot^(1/2)));
	#Generate PCoA
	curr.pcoa <- pcoa(as.matrix(curr.dist), rn=colnames(sim.ot));
	#Collect plotting data
	plot.data <- data.frame(
		Group=curr.habitat,
		Axis.1=curr.pcoa$vectors[, "Axis.1"],
		Axis.2=curr.pcoa$vectors[, "Axis.2"]
	);
	#Get percent variation explained
	if (curr.pcoa$correction[1] == "none") {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"]/sum(curr.pcoa$values[, "Relative_eig"]), digits=1)} else {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Rel_corr_eig"]/sum(curr.pcoa$values[, "Rel_corr_eig"]), digits=1)}
	#Get convex hulls around points
	plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
	#Get centroids
	plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
	curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
	#Plot ordination
	curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", colour="Group")) +
		#Add group centroids
		geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group), shape=15, size=8) +
		#Add segments of points to centroids
		geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.75) +
		#Add points (per sample)
		geom_point(size=5, alpha=0.6) +
		#Add ellipses (at 95%)
		#stat_ellipse(level=0.95, segments=101, alpha=0.8) +
		#Add convex hull around all points per group
		geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
		#Add x & y axis labels
		xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
		ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
		#Set overall plot theme
		theme_bw();
	#Save plot
	ggsave(plot=curr.plot, height=10, width=10, dpi=300, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_3.PCA.iteration_", i, ".pdf"), useDingbats=FALSE);
	###################################
	#Report
	print(paste(Sys.time(), "=> Done with iteration", i));
}
###################################
#Export simulation results
write.table(collect.F.sim_3, file=paste0(PARAM$folder.output, "hmp_samples.simulation_3.F.tsv"), sep="\t", col.names=NA);
write.table(collect.R.sim_3, file=paste0(PARAM$folder.output, "hmp_samples.simulation_3.R.tsv"), sep="\t", col.names=NA);
write.table(collect.F.null_oral.sim_3, file=paste0(PARAM$folder.output, "hmp_samples.simulation_3.F.null_oral.tsv"), sep="\t", col.names=NA);
write.table(collect.R.null_oral.sim_3, file=paste0(PARAM$folder.output, "hmp_samples.simulation_3.R.null_oral.tsv"), sep="\t", col.names=NA);
write.table(collect.F.null_gut.sim_3, file=paste0(PARAM$folder.output, "hmp_samples.simulation_3.F.null_gut.tsv"), sep="\t", col.names=NA);
write.table(collect.R.null_gut.sim_3, file=paste0(PARAM$folder.output, "hmp_samples.simulation_3.R.null_gut.tsv"), sep="\t", col.names=NA);

collect.F.sim_3 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_3.F.tsv"), header=T, row.names=1);
collect.R.sim_3 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_3.R.tsv"), header=T, row.names=1);
collect.F.null_oral.sim_3 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_3.F.null_oral.tsv"), header=T, row.names=1);
collect.R.null_oral.sim_3 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_3.R.null_oral.tsv"), header=T, row.names=1);
collect.F.null_gut.sim_3 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_3.F.null_gut.tsv"), header=T, row.names=1);
collect.R.null_gut.sim_3 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_3.R.null_gut.tsv"), header=T, row.names=1);
###################################
#Prepare data for plotting
df.F.sim_3 <- data.frame(t(collect.F.sim_3)); df.F.sim_3$val_type <- "F"; df.F.sim_3$comparison <- "oral_vs_gut";
df.F.null_oral.sim_3 <- data.frame(t(collect.F.null_oral.sim_3)); df.F.null_oral.sim_3$val_type <- "F"; df.F.null_oral.sim_3$comparison <- "oral_null";
df.F.null_gut.sim_3 <- data.frame(t(collect.F.null_gut.sim_3)); df.F.null_gut.sim_3$val_type <- "F"; df.F.null_gut.sim_3$comparison <- "gut_null";
df.R.sim_3 <- data.frame(t(collect.R.sim_3)); df.R.sim_3$val_type <- "R"; df.R.sim_3$comparison <- "oral_vs_gut";
df.R.null_oral.sim_3 <- data.frame(t(collect.R.null_oral.sim_3)); df.R.null_oral.sim_3$val_type <- "R"; df.R.null_oral.sim_3$comparison <- "oral_null";
df.R.null_gut.sim_3 <- data.frame(t(collect.R.null_gut.sim_3)); df.R.null_gut.sim_3$val_type <- "R"; df.R.null_gut.sim_3$comparison <- "gut_null";
#Melt into long format
plot.df <- rbind(
	melt(df.F.sim_3, id.vars=c("val_type", "comparison")),
	melt(df.F.null_oral.sim_3, id.vars=c("val_type", "comparison")),
	melt(df.F.null_gut.sim_3, id.vars=c("val_type", "comparison")),
	melt(df.R.sim_3, id.vars=c("val_type", "comparison")),
	melt(df.R.null_oral.sim_3, id.vars=c("val_type", "comparison")),
	melt(df.R.null_gut.sim_3, id.vars=c("val_type", "comparison"))
);
###################################
#Plot simulation results (F statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "F" & plot.df$comparison == "oral_vs_gut", ], aes(x=variable, y=value, fill=variable, colour=variable)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitter(), size=5, alpha=0.4) +
	#Scale log10
	scale_y_log10(breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), limits=PARAM$F.lim) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=8, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_3.F_statistic.pdf"), useDingbats=F);
#Plot simulation results (R2 statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "R", ], aes(x=variable, y=value, fill=comparison)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitterdodge(), size=1, alpha=0.4) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=10, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_3.R2_statistic.pdf"), useDingbats=F);
################################################################################
################################################################################


################################################################################
################################################################################
#Setup 4: Generalized Lotka-Volterra dynamics with habitat preferences and taxa interactions
#=> generate starting counts per taxon by sampling abundances from oral and gut samples (two groups, with habitat preference)
#=> simulate community dynamics according to a Lotka-Volterra model
#=> interactions are all negative between OTUs between groups, and randomized (simulated) for OTUs within group
#=> test separation between groups by different metrics
###################################
#Preallocate
collect.F.sim_4 <- collect.R.sim_4 <- collect.F.null_oral.sim_4 <- collect.R.null_oral.sim_4 <- collect.F.null_gut.sim_4 <- collect.R.null_gut.sim_4 <- matrix(nrow=8, ncol=PARAM$n.iter, dimnames=list(c("jaccard_classical", "jaccard_weighted", "bray_curtis", "morisita_horn", "TINA_uw", "TINA_w", "TINA_uw.int", "TINA_w.int"), paste0("Iteration_", sprintf("%03i", 1:PARAM$n.iter))));
###################################
#Set general model starting parameters
#=> adapted from Berry & Widder, 2014
###################################
#Modeling time
time <- list(start=0, end=PARAM$time.end, steps=PARAM$time.steps);
#Number of sites to model
sites <- PARAM$n.samples;
#Average number of taxa per site
ot.occ <- apply(curr.ot, 2, function(x) {x[x>0] <- 1; x});
siteRich <- round(mean(colSums(ot.occ)));
#Average proportion of species per site shared with at least one other site
#=> calculated from the number of OTUs and the average richness per site
sharedSp <- floor(nrow(curr.ot) / siteRich);
#Average number of interactions per species
NumInter=10;
###################################
#Iterate through simulation runs
for (i in 1:PARAM$n.iter) {
	start.ot <- sim.ot <- matrix(data=NA, nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), curr.samples));
	###################################
	#Simulate oral and gut samples
	#=> per group, by shuffling OTU counts per row
	#=> normalize and re-assign initial population sizes at 10% of total sample size
	curr.shuffled <- cbind(
		t(apply(as.matrix(curr.ot[, 1:(PARAM$n.samples/2)]), 1, sample, size=PARAM$n.samples/2)),
		t(apply(as.matrix(curr.ot[, (PARAM$n.samples/2 + 1):PARAM$n.samples]), 1, sample, size=PARAM$n.samples/2))
	);
	curr.shuffled.rel <- apply(curr.shuffled, 2, function(x) {x/sum(x)});
	start.ot <- apply(curr.shuffled.rel, 2, function(prob) {rmultinom(1, PARAM$size.sample / 10, prob)});
	dimnames(start.ot) <- dimnames(curr.ot);
	###################################
	#Plot (clustered) count heatmap of start.ot
	log.start.ot <- log10(start.ot[c(oral.otus, gut.otus),]); log.start.ot[is.infinite(log.start.ot)] <- 0;
	log.start.ot[, curr.habitat == "Gastrointestinal_tract"] <- 0 - log.start.ot[, curr.habitat == "Gastrointestinal_tract"];
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.start_ot_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(log.start.ot, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(-3, 3, by=0.01), col=colorRampPalette(brewer.pal(9, "PuOr")), trace="none", density.info="none");
	dev.off();
	###################################
	#Get (randomized) parameters for Lotka-Volterra model
	#=> code modified from Berry & Widder, 2014
	###################################
	#Assign random growth rates to OTUs
	#=> draw from a uniform distribuion on [0,1]
	r <- runif(PARAM$n.taxa, min=0, max=1);
	#Set all carrying capacities for all OTUs in all samples to total carrying capacity (sample size)
	K <- matrix(data=PARAM$size.sample, nrow=PARAM$n.taxa, ncol=PARAM$n.samples);
	###################################
	#Generate interaction matrix
	#=> code modified from Berry & Widder, 2014
	#=> set all interactions of OTUs between groups (habitats) as completely negative
	#=> randomize interactions between OTUs within groups
	#=> generate network with scale-free (Barabasi) properties (see Berry & Widder, 2014)
	#=> set diagonal to 1
	###################################
	#Preallocate
	alpha <- matrix(data=0, nrow=PARAM$n.taxa, ncol=PARAM$n.taxa, dimnames=list(rownames(curr.ot), rownames(curr.ot)));
	###################################
	#Oral samples
	curr.otus <- which(otu.preferences == "oral");
	#Generate uniformly distributed random interaction coefficients to replace rows of interaction matrix with corresponding indices
	m1 <- matrix(runif(length(curr.otus)^2, min=-1, max=1), ncol=length(curr.otus));
	#Generate mean probability of interaction based on number of species and mean number of interactions.
	prob <- NumInter / length(curr.otus);
	#Cap max number of interactions at number of species.
	if (prob>1) prob <- 1; 
	#Generate matrix of 0's and 1's based on probability of interaction
	m2 <- matrix(rbinom(length(curr.otus)^2, 1, prob), ncol=length(curr.otus));
	#Play the Barabasi game
	TMP0 <- barabasi.game(n=length(curr.otus), power=1, out.seq=rep(NumInter, length(curr.otus)), directed=F);
	TMP0 <- simplify(TMP0);
	#Turn into an adjacency matrix
	m2 <- as.matrix(get.adjacency(TMP0));
	#Replace the corresponding rows in the interaction matrix
	alpha[curr.otus, curr.otus] <- m1*m2;
	###################################
	#Gut samples
	curr.otus <- which(otu.preferences == "gut");
	#Generate uniformly distributed random interaction coefficients to replace rows of interaction matrix with corresponding indices
	m1 <- matrix(runif(length(curr.otus)^2, min=-1, max=1), ncol=length(curr.otus));
	#Generate mean probability of interaction based on number of species and mean number of interactions.
	prob <- NumInter / length(curr.otus);
	#Cap max number of interactions at number of species.
	if (prob>1) prob <- 1; 
	#Generate matrix of 0's and 1's based on probability of interaction
	m2 <- matrix(rbinom(length(curr.otus)^2, 1, prob), ncol=length(curr.otus));
	#Play the Barabasi game
	TMP0 <- barabasi.game(n=length(curr.otus), power=1, out.seq=rep(NumInter, length(curr.otus)), directed=F);
	TMP0 <- simplify(TMP0);
	#Turn into an adjacency matrix
	m2 <- as.matrix(get.adjacency(TMP0));
	#Replace the corresponding rows in the interaction matrix
	alpha[curr.otus, curr.otus] <- m1*m2;
	###################################
	#Set interactions between groups to "competitive"
	alpha[otu.preferences == "gut", otu.preferences == "oral"] <- matrix(runif(n.otu.oral*n.otu.gut, min=0.25, max=1), nrow=n.otu.gut, ncol=n.otu.oral);
	alpha[otu.preferences == "oral", otu.preferences == "gut"] <- matrix(runif(n.otu.oral*n.otu.gut, min=0.25, max=1), nrow=n.otu.oral, ncol=n.otu.gut);
	#Make sure that diagonal elements are 1
	diag(alpha) <- 1;
	###################################
	#Export heatmap of current interaction matrix
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.alpha_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(alpha[c(oral.otus, gut.otus), c(oral.otus, gut.otus)], Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(-1,1,by=0.01), col=colorRampPalette(rev(brewer.pal(9, "PRGn"))), trace="none", margins=c(20,20), density.info="none");
	dev.off()
	###################################
	#Preallocate results collectors
	sim.ot <- time.end <- rsd <- matrix(nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), paste0("Site_", 1:PARAM$n.samples)));
	#Iterate through sites
	for (s in 1:PARAM$n.samples) {
		k <- K[,s];
		#Simulation code modified from Berry & Widder, 2014
		parms <- list(r=r, a=alpha, k=k, kg=PARAM$size.sample);
		#Run simulation
		out <- n.integrate(time=time, init.x=start.ot[,s], model=lvm, parms=parms)
		#Store simulation results
		sim.ot[, s] <- as.numeric(out[nrow(out),2:ncol(out)]);
		#Report
		print(paste(Sys.time(), "=> Done with iteration", i, "sample", s));
	}
	#Filter results to remove all OTUs that do not occur in any sample
	sim.ot[sim.ot < 0.001] <- 0;
	raw.sim.ot <- sim.ot;
	sim.ot <- sim.ot[rowSums(sim.ot) > 0, ];
	keep.otus <- rownames(sim.ot);
	###################################
	#Plot count heatmap of sim.ot
	log.sim.ot <- log10(raw.sim.ot[c(oral.otus, gut.otus),]); log.sim.ot[is.infinite(log.sim.ot)] <- 0;
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.sim_ot_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(log.sim.ot, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(0, 4, by=0.01), col=colorRampPalette(brewer.pal(9, "Blues")), trace="none", density.info="none");
	dev.off();
	###################################
	#Calculate SparCC correlations and derived correlation similarities
	sim.sparcc <- sparcc(sim.ot, size.thresh=0, pseudocount=10^-6, nblocks=10, use.cores=20);
	#Filter out OTUs for which no correlation could be computed
	keep.otus.sparcc <- rownames(sim.sparcc)[apply(sim.sparcc, 1, function(x) {!all(is.na(x))})];
	sim.sparcc <- sim.sparcc[keep.otus.sparcc, keep.otus.sparcc];
	sim.ot <- sim.ot[keep.otus.sparcc, ];
	sim.S_sparcc <- 0.5 * (cor(sim.sparcc, method="pearson", use=PARAM$cor.use) + 1);
	#Calculate correlations of interactions in interaction matrix
	sim.S_int <- 0.5 * (cor(alpha[keep.otus.sparcc, keep.otus.sparcc]) + 1);
	sim.S_int[is.na(sim.S_int)] <- 0.5;
	rownames(sim.S_int) <- colnames(sim.S_int) <- rownames(sim.ot);
	###################################
	#Calculate community similarities
	###################################
	t <- 1;
	sim.cs <- list();
	#Jaccard, classical
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, classical";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard", use.cores=20);
	t <- t + 1;
	#Jaccard, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, weighted";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard.abd.frac", use.cores=20);
	t <- t + 1;
	#Bray-Curtis
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Bray-Curtis";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="bray_curtis", use.cores=20);
	t <- t + 1;
	#Morisita-Horn
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Morisita-Horn";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="morisita_horn", use.cores=20);
	t <- t + 1;
	#TINA, unweighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, unweighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.uw.norm", blocksize=10, use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, weighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.w.norm", blocksize=10, use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, unweighted, interaction correlations
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, unweighted, interaction correlations";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_int, distance="jaccard.corr.uw.norm", blocksize=10, use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, weighted, interaction correlations
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, weighted, interaction correlations";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_int, distance="jaccard.corr.w.norm", blocksize=10, use.cores=20);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	###################################
	#Test group separation
	for (t in 1:length(sim.cs)) {
		curr.aov <- adonis(as.dist(sim.cs[[t]]$cs) ~ curr.habitat, permutations=2);
		collect.F.sim_4[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.sim_4[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for oral samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[1:(PARAM$n.sample / 2), 1:(PARAM$n.sample / 2)])
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), PARAM$n.sample / 2, replace=T)), permutations=2);
		collect.F.null_oral.sim_4[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_oral.sim_4[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for gastrointestinal samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[((PARAM$n.sample/2)+1):PARAM$n.sample, ((PARAM$n.sample/2)+1):PARAM$n.sample]);
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), PARAM$n.sample / 2, replace=T)), permutations=2);
		collect.F.null_gut.sim_4[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_gut.sim_4[t, i] <- curr.aov$aov.tab$R2[1];
	}
	###################################
	#Export PCoA plot (based on Euclidean distances => i.e., a PCA plot)
	curr.dist <- dist(t(sim.ot^(1/2)));
	#Generate PCoA
	curr.pcoa <- pcoa(as.matrix(curr.dist), rn=colnames(sim.ot));
	#Collect plotting data
	plot.data <- data.frame(
		Group=curr.habitat,
		Axis.1=curr.pcoa$vectors[, "Axis.1"],
		Axis.2=curr.pcoa$vectors[, "Axis.2"]
	);
	#Get percent variation explained
	if (curr.pcoa$correction[1] == "none") {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"]/sum(curr.pcoa$values[, "Relative_eig"]), digits=1)} else {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Rel_corr_eig"]/sum(curr.pcoa$values[, "Rel_corr_eig"]), digits=1)}
	#Get convex hulls around points
	plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
	#Get centroids
	plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
	curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
	#Plot ordination
	curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", colour="Group")) +
		#Add group centroids
		geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group), shape=15, size=8) +
		#Add segments of points to centroids
		geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.75) +
		#Add points (per sample)
		geom_point(size=5, alpha=0.6) +
		#Add ellipses (at 95%)
		#stat_ellipse(level=0.95, segments=101, alpha=0.8) +
		#Add convex hull around all points per group
		geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
		#Add x & y axis labels
		xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
		ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
		#Set overall plot theme
		theme_bw();
	#Save plot
	ggsave(plot=curr.plot, height=10, width=10, dpi=300, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_4.PCA.iteration_", i, ".pdf"), useDingbats=FALSE);
	###################################
	#Export simulation results
	write.table(collect.F.sim_4, file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.F.tsv"), sep="\t", col.names=NA);
	write.table(collect.R.sim_4, file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.R.tsv"), sep="\t", col.names=NA);
	write.table(collect.F.null_oral.sim_4, file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.F.null_oral.tsv"), sep="\t", col.names=NA);
	write.table(collect.R.null_oral.sim_4, file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.R.null_oral.tsv"), sep="\t", col.names=NA);
	write.table(collect.F.null_gut.sim_4, file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.F.null_gut.tsv"), sep="\t", col.names=NA);
	write.table(collect.R.null_gut.sim_4, file=paste0(PARAM$folder.output, "hmp_samples.simulation_4.R.null_gut.tsv"), sep="\t", col.names=NA);
	###################################
	#Report
	print(paste(Sys.time(), "=> Done with iteration", i));
}
###################################
collect.F.sim_4 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_4.F.tsv"), header=T, row.names=1);
collect.R.sim_4 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_4.R.tsv"), header=T, row.names=1);
collect.F.null_oral.sim_4 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_4.F.null_oral.tsv"), header=T, row.names=1);
collect.R.null_oral.sim_4 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_4.R.null_oral.tsv"), header=T, row.names=1);
collect.F.null_gut.sim_4 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_4.F.null_gut.tsv"), header=T, row.names=1);
collect.R.null_gut.sim_4 <- read.table(paste0(PARAM$folder.output, "hmp_samples.simulation_4.R.null_gut.tsv"), header=T, row.names=1);
#Prepare data for plotting
df.F.sim_4 <- data.frame(t(collect.F.sim_4)); df.F.sim_4$val_type <- "F"; df.F.sim_4$comparison <- "oral_vs_gut";
df.F.null_oral.sim_4 <- data.frame(t(collect.F.null_oral.sim_4)); df.F.null_oral.sim_4$val_type <- "F"; df.F.null_oral.sim_4$comparison <- "oral_null";
df.F.null_gut.sim_4 <- data.frame(t(collect.F.null_gut.sim_4)); df.F.null_gut.sim_4$val_type <- "F"; df.F.null_gut.sim_4$comparison <- "gut_null";
df.R.sim_4 <- data.frame(t(collect.R.sim_4)); df.R.sim_4$val_type <- "R"; df.R.sim_4$comparison <- "oral_vs_gut";
df.R.null_oral.sim_4 <- data.frame(t(collect.R.null_oral.sim_4)); df.R.null_oral.sim_4$val_type <- "R"; df.R.null_oral.sim_4$comparison <- "oral_null";
df.R.null_gut.sim_4 <- data.frame(t(collect.R.null_gut.sim_4)); df.R.null_gut.sim_4$val_type <- "R"; df.R.null_gut.sim_4$comparison <- "gut_null";
#Melt into long format
plot.df <- rbind(
	melt(df.F.sim_4, id.vars=c("val_type", "comparison")),
	melt(df.F.null_oral.sim_4, id.vars=c("val_type", "comparison")),
	melt(df.F.null_gut.sim_4, id.vars=c("val_type", "comparison")),
	melt(df.R.sim_4, id.vars=c("val_type", "comparison")),
	melt(df.R.null_oral.sim_4, id.vars=c("val_type", "comparison")),
	melt(df.R.null_gut.sim_4, id.vars=c("val_type", "comparison"))
);
###################################
#Plot simulation results (F statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "F" & plot.df$comparison == "oral_vs_gut", ], aes(x=variable, y=value, fill=variable, colour=variable)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitter(), size=5, alpha=0.4) +
	#Scale log10
	scale_y_log10(breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), limits=PARAM$F.lim) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=10, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_4.F_statistic.pdf"), useDingbats=F);
#Plot simulation results (R2 statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "R", ], aes(x=variable, y=value, fill=comparison)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitterdodge(), size=1, alpha=0.4) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=10, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_4.R2_statistic.pdf"), useDingbats=F);
################################################################################
################################################################################


################################################################################
################################################################################
#Setup 5: Generalized Lotka-Volterra dynamics without habitat preferences, but with taxa interactions
#=> generate uniform (but noisy) starting counts for all taxa
#=> simulate community dynamics according to a Lotka-Volterra model
#=> interactions are all negative (but noisy) between OTUs between groups, and randomized (simulated) for OTUs within group
#=> test separation between groups by different metrics
###################################
#Preallocate
collect.F.sim_5 <- collect.R.sim_5 <- collect.F.null_oral.sim_5 <- collect.R.null_oral.sim_5 <- collect.F.null_gut.sim_5 <- collect.R.null_gut.sim_5 <- matrix(nrow=8, ncol=PARAM$n.iter, dimnames=list(c("jaccard_classical", "jaccard_weighted", "bray_curtis", "morisita_horn", "TINA_uw", "TINA_w", "TINA_uw.int", "TINA_w.int"), paste0("Iteration_", sprintf("%03i", 1:PARAM$n.iter))));
###################################
#Set general model starting parameters
#=> modified from Berry & Widder, 2014
###################################
#Modeling time
time <- list(start=0, end=PARAM$time.end, steps=PARAM$time.steps);
#Number of sites to model
sites <- PARAM$n.samples;
#Average number of taxa per site
ot.occ <- apply(curr.ot, 2, function(x) {x[x>0] <- 1; x});
siteRich <- round(mean(colSums(ot.occ)));
#Average proportion of species per site shared with at least one other site
#=> calculated from the number of OTUs and the average richness per site
sharedSp <- floor(nrow(curr.ot) / siteRich);
#Average number of interactions per species
NumInter=10;
#Generate (synthetic and symmetric) OTU habitat preferences
otu.pref.sim <- c(rep("oral", PARAM$n.taxa/2), rep("gut", PARAM$n.taxa/2));
names(otu.pref.sim) <- rownames(curr.ot);
###################################
#Iterate through simulation runs
for (i in 1:PARAM$n.iter) {
	start.ot <- sim.ot <- matrix(data=NA, nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), curr.samples));
	#Simulate oral and gut samples
	#=> per sample, generate random counts for OTUs from Poisson distribution
	for (s in 1:PARAM$n.samples) {start.ot[,s] <- rpois(PARAM$n.taxa, round((PARAM$size.sample/10)/PARAM$n.taxa))}
	###################################
	#Plot (clustered) count heatmap of start.ot
	log.start.ot <- log10(start.ot); log.start.ot[is.infinite(log.start.ot)] <- 0;
	log.start.ot[, ((PARAM$n.samples/2)+1):PARAM$n.samples] <- 0 - log.start.ot[, ((PARAM$n.samples/2)+1):PARAM$n.samples];
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.start_ot_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(log.start.ot, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(-3, 3, by=0.01), col=colorRampPalette(brewer.pal(9, "PuOr")), trace="none", density.info="none");
	dev.off();
	###################################
	#Get (randomized) parameters for Lotka-Volterra model
	#=> modified from Berry & Widder, 2014
	###################################
	#Set fixed growth rate for all taxa
	r <- rep(0.5, PARAM$n.taxa)
	#Set all carrying capacities for all OTUs in all samples to total carrying capacity (sample size)
	K <- matrix(data=PARAM$size.sample, nrow=PARAM$n.taxa, ncol=PARAM$n.samples);
	###################################
	#Generate interaction matrix
	#=> set all interactions of OTUs between groups (habitats) as completely negative
	#=> randomize interactions between OTUs within (synthetic and symmetric) groups
	#=> generate networks with scale-free (Barabasi) properties within groups (see Berry & Widder, 2014)
	#=> set diagonal to 1
	###################################
	#Preallocate
	alpha <- matrix(data=0, nrow=PARAM$n.taxa, ncol=PARAM$n.taxa, dimnames=list(rownames(curr.ot), rownames(curr.ot)));
	###################################
	#"Oral" samples
	curr.otus <- 1:(PARAM$n.taxa/2);
	#Generate uniformly distributed random interaction coefficients to replace rows of interaction matrix with corresponding indices
	m1 <- matrix(runif(length(curr.otus)^2, min=-1, max=1), ncol=length(curr.otus));
	#Generate mean probability of interaction based on number of species and mean number of interactions.
	prob <- NumInter / length(curr.otus);
	#Cap max number of interactions at number of species.
	if (prob>1) prob <- 1; 
	#Generate matrix of 0's and 1's based on probability of interaction
	m2 <- matrix(rbinom(length(curr.otus)^2, 1, prob), ncol=length(curr.otus));
	#Play the Barabasi game
	TMP0 <- barabasi.game(n=length(curr.otus), power=1, out.seq=rep(NumInter, length(curr.otus)), directed=F);
	TMP0 <- simplify(TMP0);
	#Turn into an adjacency matrix
	m2 <- as.matrix(get.adjacency(TMP0));
	#Replace the corresponding rows in the interaction matrix
	alpha[curr.otus, curr.otus] <- m1*m2;
	###################################
	#"Gut" samples
	curr.otus <- (PARAM$n.taxa/2 + 1):PARAM$n.taxa;
	#Generate uniformly distributed random interaction coefficients to replace rows of interaction matrix with corresponding indices
	m1 <- matrix(runif(length(curr.otus)^2, min=-1, max=1), ncol=length(curr.otus));
	#Generate mean probability of interaction based on number of species and mean number of interactions.
	prob <- NumInter / length(curr.otus);
	#Cap max number of interactions at number of species.
	if (prob>1) prob <- 1; 
	#Generate matrix of 0's and 1's based on probability of interaction
	m2 <- matrix(rbinom(length(curr.otus)^2, 1, prob), ncol=length(curr.otus));
	#Play the Barabasi game
	TMP0 <- barabasi.game(n=length(curr.otus), power=1, out.seq=rep(NumInter, length(curr.otus)), directed=F);
	TMP0 <- simplify(TMP0);
	#Turn into an adjacency matrix
	m2 <- as.matrix(get.adjacency(TMP0));
	#Replace the corresponding rows in the interaction matrix
	alpha[curr.otus, curr.otus] <- m1*m2;
	###################################
	#Set interactions between groups to "competitive" (with noise)
	alpha[otu.pref.sim == "gut", otu.pref.sim == "oral"] <- matrix(runif(length(curr.otus)^2, min=0.25, max=1), ncol=length(curr.otus));
	alpha[otu.pref.sim == "oral", otu.pref.sim == "gut"] <- matrix(runif(length(curr.otus)^2, min=0.25, max=1), ncol=length(curr.otus));
	#Make sure that diagonal elements are 1
	diag(alpha) <- 1;
	###################################
	#Preallocate results collectors
	sim.ot <- time.end <- rsd <- matrix(nrow=PARAM$n.taxa, ncol=PARAM$n.samples, dimnames=list(rownames(curr.ot), paste0("Site_", 1:PARAM$n.samples)));
	curr.sample_labels <- character(length=PARAM$n.samples); names(curr.sample_labels) <- colnames(sim.ot);
	#Iterate through sites
	for (s in 1:PARAM$n.samples) {
		k <- K[,s];
		#Simulation code modified from Berry & Widder, 2014
		parms <- list(r=r, a=alpha, k=k, kg=PARAM$size.sample);
		#Run simulation
		out <- n.integrate(time=time, init.x=start.ot[,s], model=lvm, parms=parms)
		#Store simulation results
		sim.ot[, s] <- as.numeric(out[nrow(out),2:ncol(out)]);
		#Classify current sample as "oral" or "gut"
		#=> because for completely synthetic data, this is not known a priori
		oral.abd <- sum(sim.ot[otu.pref.sim == "oral", s]);
		gut.abd <- sum(sim.ot[otu.pref.sim == "gut", s]);
		if (oral.abd > gut.abd) {curr.sample_labels[s] <- "oral"} else {curr.sample_labels[s] <- "gut"}
		#Report
		print(paste(Sys.time(), "=> Done with iteration", i, "sample", s, "classified as", curr.sample_labels[s], "with", oral.abd, gut.abd));
	}
	#Filter results to remove all OTUs that do not occur in any sample
	sim.ot[sim.ot < 0.001] <- 0;
	raw.sim.ot <- sim.ot;
	sim.ot <- sim.ot[rowSums(sim.ot) > 0, ];
	keep.otus <- rownames(sim.ot);
	###################################
	#Plot count heatmap of sim.ot
	log.sim.ot <- log10(raw.sim.ot[, c(names(curr.sample_labels)[curr.sample_labels == "oral"], names(curr.sample_labels)[curr.sample_labels == "gut"])]); log.sim.ot[is.infinite(log.sim.ot)] <- 0;
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.sim_ot_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(log.sim.ot, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(0, 4, by=0.01), col=colorRampPalette(brewer.pal(9, "Blues")), trace="none", density.info="none");
	dev.off();
	###################################
	#Calculate SparCC correlations and derived correlation similarities
	sim.sparcc <- sparcc(sim.ot, size.thresh=0, pseudocount=10^-6, nblocks=10, use.cores=20);
	#Filter out OTUs for which no correlation could be computed
	keep.otus.sparcc <- rownames(sim.sparcc)[apply(sim.sparcc, 1, function(x) {!all(is.na(x))})];
	sim.sparcc <- sim.sparcc[keep.otus.sparcc, keep.otus.sparcc];
	sim.ot <- sim.ot[keep.otus.sparcc, ];
	sim.S_sparcc <- 0.5 * (cor(sim.sparcc, method="pearson", use=PARAM$cor.use) + 1);
	#Calculate correlations of interactions in interaction matrix
	sim.S_int <- 0.5 * (cor(alpha[keep.otus.sparcc, keep.otus.sparcc]) + 1);
	sim.S_int[is.na(sim.S_int)] <- 0.5;
	rownames(sim.S_int) <- colnames(sim.S_int) <- rownames(sim.ot);
	###################################
	#Calculate community similarities
	###################################
	t <- 1;
	sim.cs <- list();
	#Jaccard, classical
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, classical";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard", use.cores=10);
	t <- t + 1;
	#Jaccard, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Jaccard, weighted";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="jaccard.abd.frac", use.cores=10);
	t <- t + 1;
	#Bray-Curtis
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Bray-Curtis";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="bray_curtis", use.cores=10);
	t <- t + 1;
	#Morisita-Horn
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "Morisita-Horn";
	sim.cs[[t]]$cs <- community.similarity.par(sim.ot, distance="morisita_horn", use.cores=10);
	t <- t + 1;
	#TINA, unweighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, unweighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.uw.norm", blocksize=10, use.cores=10);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, weighted
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, weighted";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_sparcc, distance="jaccard.corr.w.norm", blocksize=10, use.cores=10);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, unweighted, interaction correlations
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, unweighted, interaction correlations";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_int, distance="jaccard.corr.uw.norm", blocksize=10, use.cores=10);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	t <- t + 1;
	#TINA, weighted, interaction correlations
	sim.cs[[t]] <- list();
	sim.cs[[t]]$name <- "TINA, weighted, interaction correlations";
	sim.cs[[t]]$cs <- community.similarity.corr.par(sim.ot, S=sim.S_int, distance="jaccard.corr.w.norm", blocksize=10, use.cores=10);
	sim.cs[[t]]$cs[sim.cs[[t]]$cs < 0] <- 0;
	###################################
	#Test group separation
	for (t in 1:length(sim.cs)) {
		curr.aov <- adonis(as.dist(sim.cs[[t]]$cs) ~ factor(curr.sample_labels), permutations=2);
		collect.F.sim_5[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.sim_5[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for oral samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[curr.sample_labels == "oral", curr.sample_labels == "oral"])
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), length(which(curr.sample_labels == "oral")), replace=T)), permutations=2);
		collect.F.null_oral.sim_5[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_oral.sim_5[t, i] <- curr.aov$aov.tab$R2[1];
		#Null comparison for gastrointestinal samples
		curr.dist <- as.dist(sim.cs[[t]]$cs[curr.sample_labels == "gut", curr.sample_labels == "gut"]);
		curr.aov <- adonis(curr.dist ~ factor(sample(c("A", "B"), length(which(curr.sample_labels == "gut")), replace=T)), permutations=2);
		collect.F.null_gut.sim_5[t, i] <- curr.aov$aov.tab$F.Model[1];
		collect.R.null_gut.sim_5[t, i] <- curr.aov$aov.tab$R2[1];
	}
	###################################
	#Export heatmap of current interaction matrix
	pdf(file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.alpha_heatmap.iteration_", i, ".pdf"), width=10, height=10, useDingbats=F);
	heatmap.2(alpha, Rowv=F, Colv=F, dendrogram="none", symm=F, na.rm=T, revC=T, breaks=seq(-1,1,by=0.01), col=colorRampPalette(rev(brewer.pal(9, "PRGn"))), trace="none", margins=c(20,20), density.info="none");
	dev.off()
	###################################
	#Export PCoA plot (based on Euclidean distances => i.e., a PCA plot)
	curr.dist <- dist(t(sim.ot^(1/2)));
	#Generate PCoA
	curr.pcoa <- pcoa(as.matrix(curr.dist), rn=colnames(sim.ot));
	#Collect plotting data
	plot.data <- data.frame(
		Group=factor(curr.sample_labels),
		Axis.1=curr.pcoa$vectors[, "Axis.1"],
		Axis.2=curr.pcoa$vectors[, "Axis.2"]
	);
	#Get percent variation explained
	if (curr.pcoa$correction[1] == "none") {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"]/sum(curr.pcoa$values[, "Relative_eig"]), digits=1)} else {plot.var_explained <- round(100*curr.pcoa$values[1:2, "Rel_corr_eig"]/sum(curr.pcoa$values[, "Rel_corr_eig"]), digits=1)}
	#Get convex hulls around points
	plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
	#Get centroids
	plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
	curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
	#Plot ordination
	curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", colour="Group")) +
		#Add group centroids
		geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group), shape=15, size=8) +
		#Add segments of points to centroids
		geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.75) +
		#Add points (per sample)
		geom_point(size=5, alpha=0.6) +
		#Add ellipses (at 95%)
		#stat_ellipse(level=0.95, segments=101, alpha=0.8) +
		#Add convex hull around all points per group
		geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
		#Add x & y axis labels
		xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
		ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep="")) +
		#Set overall plot theme
		theme_bw();
	#Save plot
	ggsave(plot=curr.plot, height=10, width=10, dpi=300, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_5.PCA.iteration_", i, ".pdf"), useDingbats=FALSE);
	###################################
	#Export simulation results
	write.table(collect.F.sim_5, file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.F.tsv"), sep="\t", col.names=NA);
	write.table(collect.R.sim_5, file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.R.tsv"), sep="\t", col.names=NA);
	write.table(collect.F.null_oral.sim_5, file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.F.null_oral.tsv"), sep="\t", col.names=NA);
	write.table(collect.R.null_oral.sim_5, file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.R.null_oral.tsv"), sep="\t", col.names=NA);
	write.table(collect.F.null_gut.sim_5, file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.F.null_gut.tsv"), sep="\t", col.names=NA);
	write.table(collect.R.null_gut.sim_5, file=paste0(PARAM$folder.output, "hmp_samples.simulation_5.R.null_gut.tsv"), sep="\t", col.names=NA);
	###################################
	#Report
	print(paste(Sys.time(), "=> Done with iteration", i));
	###################################
}
###################################
#Prepare data for plotting
df.F.sim_5 <- data.frame(t(collect.F.sim_5)); df.F.sim_5$val_type <- "F"; df.F.sim_5$comparison <- "oral_vs_gut";
df.F.null_oral.sim_5 <- data.frame(t(collect.F.null_oral.sim_5)); df.F.null_oral.sim_5$val_type <- "F"; df.F.null_oral.sim_5$comparison <- "oral_null";
df.F.null_gut.sim_5 <- data.frame(t(collect.F.null_gut.sim_5)); df.F.null_gut.sim_5$val_type <- "F"; df.F.null_gut.sim_5$comparison <- "gut_null";
df.R.sim_5 <- data.frame(t(collect.R.sim_5)); df.R.sim_5$val_type <- "R"; df.R.sim_5$comparison <- "oral_vs_gut";
df.R.null_oral.sim_5 <- data.frame(t(collect.R.null_oral.sim_5)); df.R.null_oral.sim_5$val_type <- "R"; df.R.null_oral.sim_5$comparison <- "oral_null";
df.R.null_gut.sim_5 <- data.frame(t(collect.R.null_gut.sim_5)); df.R.null_gut.sim_5$val_type <- "R"; df.R.null_gut.sim_5$comparison <- "gut_null";
#Melt into long format
plot.df <- rbind(
	melt(df.F.sim_5, id.vars=c("val_type", "comparison")),
	melt(df.F.null_oral.sim_5, id.vars=c("val_type", "comparison")),
	melt(df.F.null_gut.sim_5, id.vars=c("val_type", "comparison")),
	melt(df.R.sim_5, id.vars=c("val_type", "comparison")),
	melt(df.R.null_oral.sim_5, id.vars=c("val_type", "comparison")),
	melt(df.R.null_gut.sim_5, id.vars=c("val_type", "comparison"))
);
###################################
#Plot simulation results (F statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "F" & plot.df$comparison == "oral_vs_gut", ], aes(x=variable, y=value, fill=variable, colour=variable)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitter(), size=5, alpha=0.4) +
	#Scale log10
	scale_y_log10(breaks=c(0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000), limits=PARAM$F.lim) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=10, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_5.F_statistic.pdf"), useDingbats=F);
#Plot simulation results (R2 statistic)
curr.plot <- ggplot(plot.df[plot.df$val_type == "R", ], aes(x=variable, y=value, fill=comparison)) +
	#Add boxplot
	geom_boxplot(alpha=0.7, outlier.colour=NA) +
	#Add jittered points
	geom_point(position=position_jitterdodge(), size=1, alpha=0.4) +
	#Flip coordinates
	coord_flip() +
	#Set overall plot theme
	theme_bw();
#Save
ggsave(curr.plot, width=20, height=10, filename=paste0(PARAM$folder.output, "hmp_samples.simulation_5.R2_statistic.pdf"), useDingbats=F);
################################################################################
################################################################################




q()


















