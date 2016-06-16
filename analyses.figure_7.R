#!/usr/bin/Rscript
################################################################################
#Interaction-adjusted beta diversity (community similarity) estimation.
#
#Prepare analyses shown in Figure 6
#=> analysis of TARA Oceans micro-eukaryotic plankton data
#
#=>	load TARA data for all samples, V9 tag OTU table
#=> filter data down to relevant samples & OTUs
#=>	calculate global cooccurrence network
#=>	for these datasets, estimate community similarity using either traditional or interaction-adjusted formulas
#
#
#2015-10-16
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
PARAM$file.sample_list <- #PUT RAW LIST OF SAMPLES FILE HERE
PARAM$file.sample_metadata <- #PUT RAW SAMPLE METADATA FILE HERE
PARAM$file.otu_table <- #PUT OTU TABLE FILE HERE
PARAM$cor.use <- "na.or.complete";
PARAM$p.adjust.method <- "hochberg";
PARAM$use.cores <- 20;
###################################

###################################
#Set parameters for data processing
###################################
#Minimum OTU size required across all samples
PARAM$thresh.otu_size <- 30;
#Minimum sample size required
PARAM$thresh.sample_size <- 1000;
###################################

###################################
#Import functions
source(PARAM$file.functions)

#Define functions for geodesic distance calculation
deg2rad <- function(deg) return(deg*pi/180)
geodesic <- function(long1, lat1, long2, lat2, R=6371) {acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R}
###################################

###################################
#Read sample raw metadata
#=> as downloaded from http://taraoceans.sb-roscoff.fr/EukDiv/data/Database_W1_Sample_parameters.xls
#=> slightly pre-processed in Excel: spaces replaced by tabs, "mu" (for micron) replaced by "u"
#=> choose PANGAEA_ACCESSION_NUMBER field as rownames (as these correspond to colnames in OTU table)
sample.data.raw <- read.table(PARAM$file.sample_metadata, header=T, sep="\t", row.names=3);
#Annotate filter size fractions as classes
#=>	piconano		=> 0.8-5um
#=>	nano				=> 5-20um
#=>	micro				=> 20-180um
#=>	meso				=> 180-2000um
#=>	"NA"				=> all other size fractions
sample.data.raw$FILTER_SIZE_CLASS <- factor(unlist(lapply(sample.data.raw$FILTER_SIZE_RANGE, function(i) {
	if (i == "0.8-5") {"piconano"} else if (i == "5-20") {"nano"} else if (i == "20-180") {"micro"} else if (i == "180-2000") {"meso"} else {NA}
})));
#Get interaction term of TARA_STATION and SAMPLING_DEPTH
sample.data.raw$STATION_BY_DEPTH <- interaction(sample.data.raw$TARA_STATION_IDENTIFIER, sample.data.raw$SAMPLING_DEPTH);
################################################################################
################################################################################


################################################################################
################################################################################
#Load OTU table
#=> as downloaded from http://taraoceans.sb-roscoff.fr/EukDiv/data/Database_W5_OTU_occurences.tsv.zip
#=> available online also from http://doi.pangaea.de/10.1594/PANGAEA.843022
#=> read raw data
#=> extract OTU table and OTU "metadata"
#=> filter OTU table (by sample and OTU properties)
################################################################################
################################################################################
#Read OTU data
otu.data.raw <- read.table(PARAM$file.otu_table, header=T, sep="\t", row.names=1, quote="");
#Extract OTU table
otu.table.raw <- as.matrix(otu.data.raw[, rownames(sample.data.raw)]);
###################################

###################################
#Filter OTU table
###################################
#Remove samples...
#...which come from stations in North Atlantic (n=1) or the Southern Ocean (n=2)
#...which are singletons
remove.samples <- rownames(sample.data.raw)[c(which(sample.data.raw$STATION_GEOLOC %in% c("North_Atlantic", "Southern_Ocean")), which(sample.data.raw$TARA_STATION_IDENTIFIER %in% names(which(table(sample.data.raw$TARA_STATION_IDENTIFIER) == 1))))];
#remove.samples <- c(rownames(sample.data.raw)[c(which(is.na(sample.data.raw$FILTER_SIZE_CLASS)), which(sample.data.raw$STATION_GEOLOC %in% c("North_Atlantic", "Southern_Ocean")))], "TARA_N000000115", "TARA_N000000238");
#Remove OTUs of global size <= threshold
remove.otus <- rownames(otu.table.raw)[which(rowSums(otu.table.raw[, !colnames(otu.table.raw) %in% remove.samples]) <= PARAM$thresh.otu_size)];
#Get current temporary OTU table
otu.table.tmp <- otu.table.raw[!rownames(otu.table.raw) %in% remove.otus, !colnames(otu.table.raw) %in% remove.samples];
#Prune down sample.data and otu.data
sample.data <- sample.data.raw[!rownames(sample.data.raw) %in% remove.samples, ];
my.sample.data <- do.call(rbind, lapply(unique(sample.data$STATION_BY_DEPTH), function(s) {sample.data[which(sample.data$STATION_BY_DEPTH == s)[1],]}))
rownames(my.sample.data) <- my.sample.data$STATION_BY_DEPTH;
otu.data <- otu.data.raw[!rownames(otu.data.raw) %in% remove.otus, !colnames(otu.data.raw) %in% rownames(sample.data.raw)];
#Sum up counts for different size fractions per sample (stations become "samples")
samples <- unique(as.character(sample.data$STATION_BY_DEPTH));
tmp.sum_ot <- lapply(samples, function (s) {rowSums(otu.table.tmp[, rownames(sample.data)[which(sample.data$STATION_BY_DEPTH == s)]])});
#Coerce into OTU table
#=> sparse matrix, as implemented in the "Matrix" package
ot <- do.call(cbind, tmp.sum_ot);
colnames(ot) <- samples;
rownames(ot) <- otu.data.raw$cid[rownames(otu.data.raw) %in% rownames(ot)];
ot <- Matrix(ot, sparse=T);
#Coerce into phyloseq object
my.ps <- phyloseq(otu_table(as.matrix(ot), taxa_are_rows=T), sample_data(my.sample.data));
###################################
#Save data
save(PARAM, ot, otu.data, sample.data, my.sample.data, file=paste(PARAM$folder.data, "TARA.otu_table.RData", sep=""));
###################################

###################################
#Calculate global SparCC network and derived S matrix
###################################
#SparCC correlation
if ("SparCC" == "SparCC") {
	########################
	#SparCC correlation analysis
	#=> following Friedman & Alm, PLOS CB, 2012
	#=> apparently, the pseudocount has a large influence on the results...
	#=> expects an OTU table in which OTUs are rows, samples are columns
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
	remove.otus <- rowSums(ot) < size.thresh; keep.otus <- ! remove.otus;
	o.t <- ot[! remove.otus, ] + pseudocount;
	otus <- rownames(ot);
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
	save(global.sparcc, file=paste(PARAM$folder.data, "TARA.SparCC_global.RData", sep=""));
	rm(my.combs, dat.split, mat.T, mat.T.desc, var.t, mat.a, omega);
	########################
}
########################
#Calculate derived OTU similarity matrix S
cat(paste("Calculating correlations of SparCC correlations =>", Sys.time()), sep="\n");
tmp.S <- cor.par(global.sparcc, method="pearson", use=PARAM$cor.use, use.cores=40);
cat(paste("Done =>", Sys.time()), sep="\n");
S.sparcc <- 0.5 * (tmp.S + 1);
#Save correlations & tidy up
save(S.sparcc, file=paste(PARAM$folder.data, "TARA.S_SparCC_global.RData", sep=""));
rm(tmp.S);
################################################################################
################################################################################


################################################################################
################################################################################
#Calculate pairwise community (dis)similarities
#=> "classical" indices
#=> SparCC-corrected indices
################################################################################
################################################################################
#Preallocate
get.cs <- list(); t <- 1;
############################
#Classical indices
############################
#Jaccard index, classical, unweighted
get.cs[[t]] <- list(); get.cs[[t]]$name <- "Jaccard, classical"; get.cs[[t]]$call <- "jaccard"; t <- t + 1;
#Jaccard index, classical, weighted
get.cs[[t]] <- list(); get.cs[[t]]$name <- "Jaccard, classical, weighted, no fractions"; get.cs[[t]]$call <- "jaccard.abd"; t <- t + 1;
#Jaccard index, classical, weighted, fraction-based
get.cs[[t]] <- list(); get.cs[[t]]$name <- "Jaccard, classical, weighted"; get.cs[[t]]$call <- "jaccard.abd.frac"; t <- t + 1;
#Jaccard index, classical, weighted, fraction-based, Chao's version
get.cs[[t]] <- list(); get.cs[[t]]$name <- "Jaccard, classical, Chao fraction-based"; get.cs[[t]]$call <- "jaccard.abd.frac_chao"; t <- t + 1;
#Jaccard index, classical, weighted, alternative formula
get.cs[[t]] <- list(); get.cs[[t]]$name <- "Jaccard, classical, alternative formula"; get.cs[[t]]$call <- "jaccard.abd.alt"; t <- t + 1;
#Jaccard index, classical, Chao' abundance correction
get.cs[[t]] <- list(); get.cs[[t]]$name <- "Jaccard, classical, Chao"; get.cs[[t]]$call <- "jaccard.abd.chao"; t <- t + 1;
#Bray-Curtis, classical
get.cs[[t]] <- list(); get.cs[[t]]$name <- "Bray-Curtis"; get.cs[[t]]$call <- "bray_curtis"; t <- t + 1;
#Morisita-Horn, classical
get.cs[[t]] <- list(); get.cs[[t]]$name <- "Morisita-Horn"; get.cs[[t]]$call <- "morisita_horn"; t <- t + 1;
#Iterate through indices and calculate pairwise sample similarities
for (t in seq(1, length(get.cs))) {
	print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
	#Calculate all pairwise distances
	my.cs <- community.similarity.par(ot, distance= get.cs[[t]]$call, use.cores=20);
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
get.cs[[t]] <- list(); get.cs[[t]]$name <- "TINA, unweighted"; get.cs[[t]]$call <- "jaccard.corr.uw.norm"; t <- t + 1;
#TINA, weighted
get.cs[[t]] <- list(); get.cs[[t]]$name <- "TINA, weighted"; get.cs[[t]]$call <- "jaccard.corr.w.norm"; t <- t + 1;
############################
#Iterate through (corrected) indices and calculate pairwise community similarities
for (t in seq(start.t, length(get.cs))) {
	print(paste(Sys.time(), "Starting", get.cs[[t]]$name));
	#Calculate all pairwise similarities
	my.cs <- community.similarity.corr.par(ot, S=S.sparcc, distance=get.cs[[t]]$call, blocksize=1000, use.cores=PARAM$use.cores)
	#Manually correct for rounding errors
	my.cs[my.cs < 0] <- 0;
	#Store current similarities
	get.cs[[t]]$cs <- my.cs;
	#Report
	print(paste(Sys.time(), "Done with", get.cs[[t]]$name));
	rm(my.cs)
}
############################
#Save
save(get.cs, file=paste(PARAM$folder.data, "TARA.get_cs.RData", sep=""));
################################################################################
################################################################################


################################################################################
################################################################################
#Get PCoA plots and PERMANOVA tests
#=> on surface water vs deep chlorophyll maximum
#=> on differences between sampling regions
################################################################################
################################################################################
#Preallocate
all.methods <- unlist(lapply(get.cs, function(i) i$name));
collect.F <- collect.R2 <- collect.P <- matrix(nrow=length(all.methods), ncol=3, dimnames=list(all.methods, c("Region", "Depth", "Region_Depth")));
#Iterate through community similarity metrics
for (t in seq(1, length(get.cs))) {
	curr.cs <- get.cs[[t]]$cs;
	#PCoA on water layer and region
	#curr.pcoa <- pcoa(curr.cs, rn=rownames(my.sample.data));
	curr.pcoa <- ordinate(my.ps, method = "PCoA", distance=as.dist(curr.cs));
	#Prepare data for plotting
	plot.data <- data.frame(my.sample.data, Axis.1=curr.pcoa$vectors[, "Axis.1"], Axis.2=curr.pcoa$vectors[, "Axis.2"], Group=interaction(my.sample.data$STATION_GEOLOC, my.sample.data$SAMPLING_DEPTH));
	plot.var_explained <- round(100*curr.pcoa$values[1:2, "Relative_eig"]/sum(curr.pcoa$values[curr.pcoa$values[, "Relative_eig"] > 0, "Relative_eig"]), digits=1);
	plot.hulls <- ddply(plot.data, "Group", function(df) df[chull(df$Axis.1, df$Axis.2), ]);
	plot.centroids <- ddply(plot.data, "Group", function(df) c(mean(df$Axis.1), mean(df$Axis.2)));
	curr.plot.df <- merge(plot.data, plot.centroids, by="Group");
	#Make plot
	curr.plot <- ggplot(curr.plot.df, aes_string(x="Axis.1", y="Axis.2", color="STATION_GEOLOC", shape="SAMPLING_DEPTH")) +
		geom_point(size=5) +
		geom_polygon(data = plot.hulls, aes(fill=Group), color=NA, alpha = 0.1) +
		geom_segment(data=curr.plot.df, aes(x=Axis.1, y=Axis.2, xend=V1, yend=V2), alpha=0.75) +
		geom_point(data=curr.plot.df, aes(x=V1, y=V2, color=Group), shape=15, size=8) +
		xlab(paste("Axis 1 [", plot.var_explained[1], "%]", sep="")) +
		ylab(paste("Axis 2 [", plot.var_explained[2], "%]", sep=""));
	ggsave(plot=curr.plot,height=20,width=20,dpi=300, filename=paste(PARAM$folder.output, "TARA.", get.cs[[t]]$name, ".PCoA.pdf", sep = ""), useDingbats=FALSE);
	#Calculate PERMANOVA
	curr.adonis <- adonis(curr.cs ~ STATION_GEOLOC*SAMPLING_DEPTH, data=my.sample.data);
	#Store
	get.cs[[t]]$permanova <- curr.adonis;
	#F values
	collect.F[t, "Region"] <- curr.adonis$aov.tab$F.Model[1];
	collect.F[t, "Depth"] <- curr.adonis$aov.tab$F.Model[2];
	collect.F[t, "Region_Depth"] <- curr.adonis$aov.tab$F.Model[3];
	#R2 values
	collect.R2[t, "Region"] <- curr.adonis$aov.tab$R2[1];
	collect.R2[t, "Depth"] <- curr.adonis$aov.tab$R2[2];
	collect.R2[t, "Region_Depth"] <- curr.adonis$aov.tab$R2[3];
	#P value
	collect.P[t, "Region"] <- curr.adonis$aov.tab[1, "Pr(>F)"];
	collect.P[t, "Depth"] <- curr.adonis$aov.tab[2, "Pr(>F)"];
	collect.P[t, "Region_Depth"] <- curr.adonis$aov.tab[3, "Pr(>F)"];
}
write.table(collect.F, file=paste(PARAM$folder.output, "TARA.adonis.F.tsv", sep = ""), col.names=NA, sep="\t");
write.table(collect.R2, file=paste(PARAM$folder.output, "TARA.adonis.R2.tsv", sep = ""), col.names=NA, sep="\t");
write.table(collect., file=paste(PARAM$folder.output, "TARA.adonis.P.tsv", sep = ""), col.names=NA, sep="\t");
################################################################################
################################################################################


################################################################################
################################################################################
#Get biogeographic analyses
#=> calculate pairwise geodesic distances between all samples
#=> match community to geographic distances (correlations)
################################################################################
################################################################################
#Plot all samples on a map of the world
world_map <- map_data("world");
#Empty base plot
p <- ggplot() + coord_fixed();
#Add map to base plot
base_world <- p + geom_polygon(data=world_map, aes(x=long, y=lat, group=group), fill="gray");
#Add data to world map plot
curr.plot <- base_world + geom_point(data=my.sample.data, aes(x=LONGITUDE, y=LATITUDE, color=STATION_GEOLOC), size=5);
ggsave(plot=curr.plot,height=20,width=40,dpi=300, filename=paste(PARAM$folder.output, "TARA.all_samples.world_map.pdf", sep = ""), useDingbats=FALSE);
############################
#Get pairwise geodesic distances between all samples
geodesic.list <- lapply(seq(1, nrow(my.sample.data)), function(r.1) {s.1 <- my.sample.data[r.1,]; lat.1 <- unname(deg2rad(as.numeric(s.1["LATITUDE"]))); lon.1 <- unname(deg2rad(as.numeric(s.1["LONGITUDE"]))); unlist(lapply(seq(1, nrow(my.sample.data)), function(r.2) {s.2 <- my.sample.data[r.2,]; lat.2 <- unname(deg2rad(as.numeric(s.2["LATITUDE"]))); lon.2 <- unname(deg2rad(as.numeric(s.2["LONGITUDE"]))); geodesic(lon.1, lat.1, lon.2, lat.2)}))});
geodesic.dist <- do.call(rbind, geodesic.list);
geo.dist.SUR <- as.dist(geodesic.dist[my.sample.data$SAMPLING_DEPTH == "SUR", my.sample.data$SAMPLING_DEPTH == "SUR"]);
geo.dist.DCM <- as.dist(geodesic.dist[my.sample.data$SAMPLING_DEPTH == "DCM", my.sample.data$SAMPLING_DEPTH == "DCM"]);
#Get pairwise differences in latitude between all samples
geo.dist_lat.list <- lapply(seq(1, nrow(my.sample.data)), function(r.1) {s.1 <- my.sample.data[r.1,]; lat.1 <- unname(deg2rad(as.numeric(s.1["LATITUDE"]))); unlist(lapply(seq(1, nrow(my.sample.data)), function(r.2) {s.2 <- my.sample.data[r.2,]; lat.2 <- unname(deg2rad(as.numeric(s.2["LATITUDE"]))); abs(lat.1 - lat.2)}))});
geo.dist_lat <- do.call(rbind, geo.dist_lat.list);
geo.dist_lat.SUR <- as.dist(geo.dist_lat[my.sample.data$SAMPLING_DEPTH == "SUR", my.sample.data$SAMPLING_DEPTH == "SUR"]);
geo.dist_lat.DCM <- as.dist(geo.dist_lat[my.sample.data$SAMPLING_DEPTH == "DCM", my.sample.data$SAMPLING_DEPTH == "DCM"]);
#Preallocate
collect.cor <- collect.cor.lat <- matrix(nrow=length(all.methods), ncol=2, dimnames=list(all.methods, c("SUR", "DCM")));
collect.cor.to_PCoA <- matrix(nrow=length(all.methods), ncol=6, dimnames=list(all.methods, c("SUR_Absolute_Lat", "DCM_Absolute_Lat", "SUR_Northern_Hemisphere", "DCM_Northern_Hemisphere", "SUR_Southern_Hemisphere", "DCM_Southern_Hemisphere")));
#Iterate through community similarity metrics
for (t in seq(1, length(get.cs))) {
	curr.cs <- get.cs[[t]]$cs;
	curr.cs.SUR <- as.dist(curr.cs[my.sample.data$SAMPLING_DEPTH == "SUR", my.sample.data$SAMPLING_DEPTH == "SUR"]);
	curr.cs.DCM <- as.dist(curr.cs[my.sample.data$SAMPLING_DEPTH == "DCM", my.sample.data$SAMPLING_DEPTH == "DCM"]);
	############################
	#Plot data
	plot.data <- data.frame(Geographic_Distance=c(as.numeric(geo.dist.SUR), as.numeric(geo.dist.DCM)), Community_Distance=c(as.numeric(curr.cs.SUR), as.numeric(curr.cs.DCM)), Depth=c(rep("Surface", length(curr.cs.SUR)), rep("Deep Chlrorophyll Maximum", length(curr.cs.DCM))));
	curr.plot <- ggplot(plot.data, aes(x=Geographic_Distance, y=Community_Distance, color=Depth)) + geom_point(size=3, alpha=0.4) + geom_smooth(method = "lm", se=FALSE, formula = y ~ x);
	ggsave(plot=curr.plot,height=20,width=20,dpi=300, filename=paste(PARAM$folder.output, "TARA.", get.cs[[t]]$name, ".geography_correlations.pdf", sep = ""), useDingbats=FALSE);
	#Get correlation of distances
	collect.cor[t, "SUR"] <- cor(curr.cs.SUR, geo.dist.SUR, method="spearman");
	collect.cor[t, "DCM"] <- cor(curr.cs.DCM, geo.dist.DCM, method="spearman");
	############################
	#Get correlations to latitude only
	plot.data <- data.frame(Geographic_Distance=c(as.numeric(geo.dist_lat.SUR), as.numeric(geo.dist_lat.DCM)), Community_Distance=c(as.numeric(curr.cs.SUR), as.numeric(curr.cs.DCM)), Depth=c(rep("Surface", length(curr.cs.SUR)), rep("Deep Chlrorophyll Maximum", length(curr.cs.DCM))));
	curr.plot <- ggplot(plot.data, aes(x=Geographic_Distance, y=Community_Distance, color=Depth)) + geom_point(size=3, alpha=0.4) + geom_smooth(method = "lm", se=FALSE, formula = y ~ x);
	ggsave(plot=curr.plot,height=20,width=20,dpi=300, filename=paste(PARAM$folder.output, "TARA.", get.cs[[t]]$name, ".geography_correlations_lat.pdf", sep = ""), useDingbats=FALSE);
	#Get correlation of distances
	collect.cor.lat[t, "SUR"] <- cor(curr.cs.SUR, geo.dist_lat.SUR, method="spearman");
	collect.cor.lat[t, "DCM"] <- cor(curr.cs.DCM, geo.dist_lat.DCM, method="spearman");
	############################
	#Plot PCoA first component axis vs latitude
	curr.pcoa <- ordinate(my.ps, method = "PCoA", distance=as.dist(curr.cs));
	#Stitch together in data.frame for plotting
	plot.data <- data.frame(Latitude=my.sample.data$LATITUDE, PCoA_Axis_1=curr.pcoa$vectors[,1], Depth=my.sample.data$SAMPLING_DEPTH, Region=my.sample.data$STATION_GEOLOC, Hemisphere=sapply(my.sample.data$LATITUDE, function(i) {if(i >= 0) {"North"} else {"South"}}));
	#Get current correlations
	cor.abs.SUR <- cor(abs(plot.data$Latitude[plot.data$Depth == "SUR"]), plot.data$PCoA_Axis_1[plot.data$Depth == "SUR"], method="spearman");
	cor.abs.DCM <- cor(abs(plot.data$Latitude[plot.data$Depth == "DCM"]), plot.data$PCoA_Axis_1[plot.data$Depth == "DCM"], method="spearman");
	cor.North.SUR <- cor(plot.data$Latitude[plot.data$Depth == "SUR" & plot.data$Hemisphere == "North"], plot.data$PCoA_Axis_1[plot.data$Depth == "SUR" & plot.data$Hemisphere == "North"], method="spearman");
	cor.North.DCM <- cor(plot.data$Latitude[plot.data$Depth == "DCM" & plot.data$Hemisphere == "North"], plot.data$PCoA_Axis_1[plot.data$Depth == "DCM" & plot.data$Hemisphere == "North"], method="spearman");
	cor.South.SUR <- cor(plot.data$Latitude[plot.data$Depth == "SUR" & plot.data$Hemisphere == "South"], plot.data$PCoA_Axis_1[plot.data$Depth == "SUR" & plot.data$Hemisphere == "South"], method="spearman");
	cor.South.DCM <- cor(plot.data$Latitude[plot.data$Depth == "DCM" & plot.data$Hemisphere == "South"], plot.data$PCoA_Axis_1[plot.data$Depth == "DCM" & plot.data$Hemisphere == "South"], method="spearman");
	collect.cor.to_PCoA[t, ] <- c(cor.abs.SUR, cor.abs.DCM, cor.North.SUR, cor.North.DCM, cor.South.SUR, cor.South.DCM);
	#Get variance explained
	eigvec <- curr.pcoa$values$Relative_eig; percvar <- round(100 * (eigvec[1] / sum(eigvec)), 1); ax.lab <- paste0("Axis 1 [", as.character(percvar), "%]");
	#Plot
	curr.plot <- ggplot(plot.data, aes(x=PCoA_Axis_1, y=Latitude, color=interaction(Depth, Hemisphere))) + geom_point(size=5) + geom_smooth(method = "lm", se=FALSE, formula = y ~ x) + scale_y_continuous(limits=c(-90, 90)) + xlab(ax.lab);
	ggsave(plot=curr.plot,height=20,width=20,dpi=300, filename=paste(PARAM$folder.output, "TARA.", get.cs[[t]]$name, ".PCoA_vs_lat.pdf", sep = ""), useDingbats=FALSE);
}
#Export correlation data
write.table(collect.cor, file=paste(PARAM$folder.output, "TARA.all_samples.geography_correlations.tsv", sep = ""), col.names=NA, sep="\t");
write.table(collect.cor.lat, file=paste(PARAM$folder.output, "TARA.all_samples.latitude_correlations.tsv", sep = ""), col.names=NA, sep="\t");
write.table(collect.cor.to_PCoA, file=paste(PARAM$folder.output, "TARA.all_samples.latitude_correlations_to_PCoA.tsv", sep = ""), col.names=NA, sep="\t");
################################################################################
################################################################################









q()