#!/usr/bin/Rscript
################################################################################
#Interaction-adjusted beta diversity (community similarity) estimation.
#
#Functions to calculate community similarities
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
#Define functions

########################
#SparCC correlation analysis
#=> following Friedman & Alm, PLOS Comp Biol, 2012, http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002687
#=> expects an OTU table in which OTUs are rows, samples are columns
sparcc <- function(otu.table, size.thresh=1, pseudocount=10^-6, nblocks=400, use.cores=detectCores()) {
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
	remove.otus <- rowSums(otu.table) < size.thresh; keep.otus <- ! remove.otus;
	#o.t <- otu.table[1:1000, ] + pseudocount;
	o.t <- otu.table[! remove.otus, ] + pseudocount;
	rm (otu.table)
	otus <- rownames(o.t);
	n.otu <- length(otus);
	########################
	
	########################
	#Preallocate blocks for parallel processing & Aitchinson's T matrix
	#=> based on https://gist.github.com/bobthecat/5024079
	size.split <- floor(n.otu / nblocks);
	if (size.split < 1) {size.split <- 1}
	my.split <- list(); length(my.split) <- nblocks;
	my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split));
	my.split[[nblocks]] <- (size.split*(nblocks-1)):n.otu;
	dat.split <- mclapply(my.split, function(g) {o.t[g,]}, mc.cores=use.cores);
	rm(o.t)
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
	mat.rho <- foreach(i = 1:n.otu, .combine='rbind', .multicombine=T) %dopar% {(omega[i]^2 + omega^2 - mat.T[i,]) / (2 * omega[i] * omega)}
	rownames(mat.rho) <- rownames(ot);
	cat(paste("Done with correlation estimation; returning data matrix =>", Sys.time()), sep="\n");
	########################
	
	#Return
	mat.rho
}
########################

########################
#Parallel computation of correlations
#=> between all columns of the input matrix
#=> following https://gist.github.com/bobthecat/5024079
cor.par <- function(mat, method="pearson", nblocks=400, use="na.or.complete", use.cores=detectCores()) {
	########################
	#Load packages for parallel processing
	require("foreach");
	require("bigmemory");
	library("doMC", quietly=T);
	#Register cluster
	registerDoMC(cores=use.cores);
	########################
	
	########################
	#Preallocate blocks for parallel processing
	#=> based on https://gist.github.com/bobthecat/5024079
	cat(paste("Correlation preallocation steps =>", Sys.time()), sep="\n");
	n.el <- ncol(mat);
	if (n.el < 2*nblocks) {return(cor(mat, method=method, use=use))}
	size.split <- floor(n.el / nblocks);
	if (size.split < 1) {size.split <- 1}
	my.split <- list(); length(my.split) <- nblocks;
	my.split[1:(nblocks-1)] <- split(1:(size.split*(nblocks-1)), rep(1:(nblocks-1), each = size.split));
	my.split[[nblocks]] <- (size.split*(nblocks-1)):n.el;
	dat.split <- mclapply(my.split, function(g) {as.matrix(mat[,g])}, mc.cores=use.cores);
	#Get combinations of splits
	my.combs <- expand.grid(1:length(my.split), 1:length(my.split));
	my.combs <- t(apply(my.combs, 1, sort));
	my.combs <- unique(my.combs);
	#Preallocate correlation matrix as big.matrix ("shared" in memory, so accessible from w/in foreach loop)
	mat.rho <- big.matrix(nrow=n.el, ncol=n.el, dimnames=list(colnames(mat), colnames(mat)), shared=T);
	mat.rho.desc <- describe(mat.rho);
	cat(paste("Done with preallocation steps =>", Sys.time()), sep="\n");
	########################
	
	########################
	#Compute correlation matrix
	#=> iterate through each block combination, calculate matrix
	#		between blocks and store them in the preallocated matrix on both
	#		symmetric sides of the diagonal
	cat(paste("Starting parallel correlation computations =>", Sys.time()), sep="\n");
	results <- foreach(i = 1:nrow(my.combs)) %dopar% {
		#Get current combination
		curr.comb <- my.combs[i, ];
		#Get current data
		g.1 <- my.split[[curr.comb[1]]];
		g.2 <- my.split[[curr.comb[2]]];
		curr.rho <- cor(x=dat.split[[curr.comb[1]]], y=dat.split[[curr.comb[2]]], method=method);
		#Store
		curr.mat.rho <- attach.big.matrix(mat.rho.desc);
		curr.mat.rho[g.1, g.2] <- curr.rho;
		curr.mat.rho[g.2, g.1] <- t(curr.rho);
		#Return
		TRUE
	}
	cat(paste("Done with correlation computation =>", Sys.time()), sep="\n");
	########################
	
	#Return
	as.matrix(mat.rho);
}
########################

########################
#Traditional indices of community similarity
#=> formulated as distances/dissimilarities (not similarities)
#=> calculation in parallel (using mclapply)
community.similarity.par <- function(o.t, distance="jaccard", use.cores=detectCores()) {
	########################
	#Load packages for parallel processing
	require("foreach");
	require("bigmemory");
	library("doMC", quietly=T);
	#Register cluster
	registerDoMC(cores=use.cores);
	########################
	
	########################
	samples <- colnames(o.t);
	#Turn OTU table into list format
	ot.occ.list <- apply((o.t > 0), 2, which);
	ot.count.list <- mclapply(seq(1,length(ot.occ.list)), function(i) {o.t[ot.occ.list[[i]],i]}, mc.cores=use.cores);
	names(ot.occ.list) <- names(ot.count.list) <- samples;
	########################
	
	########################
	#Compute community distance matrix
	########################
	#Classic, unweighted Jaccard
	if (distance == "jaccard") {
		cs.list <- mclapply(ot.occ.list, function(a.list) {unlist(lapply(ot.occ.list, function(b.list) {1 - (length(intersect(a.list, b.list)) / length(union(a.list, b.list)))}))}, mc.cores=use.cores);
	}
	########################
	#Classic, weighted Jaccard
	if (distance == "jaccard.abd") {
		cs.list <- mclapply(ot.count.list, function(a.count) {
			a.N <- sum(a.count); a.list <- names(a.count);
			unlist(lapply(ot.count.list, function(b.count) {ab.shared <- intersect(a.list, names(b.count)); 1 - (sum(a.count[ab.shared]) + sum(b.count[ab.shared])) / (a.N + sum(b.count))}))
		}, mc.cores=use.cores);
	}
	########################
	#Classic, weighted Jaccard, based on fractions (more balancing for uneven sample sizes)
	if (distance == "jaccard.abd.frac") {
		cs.list <- mclapply(ot.count.list, function(a.count) {
			a.N <- sum(a.count); a.list <- names(a.count); a.rel <- a.count / a.N;
			unlist(lapply(ot.count.list, function(b.count) {
				b.N <- sum(b.count); ab.shared <- intersect(a.list, names(b.count)); b.rel <- b.count / b.N;
				1 - (0.5 * (sum(a.rel[ab.shared]) + sum(b.rel[ab.shared])))
			}))
		}, mc.cores=use.cores);
	}
	########################
	#Classic, weighted Jaccard, based on fractions (more balancing for uneven sample sizes), Chao's version
	#=> after Chao et al., 2005, Ecology Letters
	if (distance == "jaccard.abd.frac_chao") {
		cs.list <- mclapply(ot.count.list, function(a.count) {
			a.N <- sum(a.count); a.list <- names(a.count); a.rel <- a.count / a.N;
			unlist(lapply(ot.count.list, function(b.count) {
				b.N <- sum(b.count); ab.shared <- intersect(a.list, names(b.count)); b.rel <- b.count / b.N;
				U <- sum(a.rel[ab.shared]); V <- sum(b.rel[ab.shared]);
				if (U == 0 & V == 0) {return(1)} else {return(1 - ((U*V) / (U + V - (U*V))))}
			}))
		}, mc.cores=use.cores);
	}
	########################
	#Classic, weighted Jaccard, alternative formulation (Bray-Curtis-like)
	if (distance == "jaccard.abd.alt") {
		cs.list <- mclapply(ot.count.list, function(a.count) {
			a.N <- sum(a.count); a.list <- names(a.count);
			unlist(lapply(ot.count.list, function(b.count) {
				b.list <- names(b.count);
				ab.shared <- intersect(a.list, b.list); ab.union <- union(a.list, b.list);
				a.freq <- b.freq <- numeric(length=length(ab.union)); names(a.freq) <- names(b.freq) <- ab.union;
				a.freq[a.list] <- a.count; b.freq[b.list] <- b.count;
				1 - (sum(pmin(a.freq, b.freq)) / sum(pmax(a.freq, b.freq)));
			}))
		}, mc.cores=use.cores);
	}
	########################
	#Classic, weighted Jaccard, Chao's version
	#=> after Chao et al., 2005, Ecology Letters
	if (distance == "jaccard.abd.chao") {
		cs.list <- mclapply(ot.count.list, function(a.count) {
			a.N <- sum(a.count); a.list <- names(a.count);
			unlist(lapply(ot.count.list, function(b.count) {
				b.list <- names(b.count); b.N <- sum(b.count);
				#Get shared taxa & union of taxa
				ab.shared <- intersect(a.list, b.list); ab.union <- sort(union(a.list, b.list));
				a.freq <- b.freq <- numeric(length=length(ab.union)); names(a.freq) <- names(b.freq) <- ab.union;
				a.freq[a.list] <- a.count; b.freq[b.list] <- b.count;
				#If no taxa are shared, return trivial (distance = 1), otherwise compute Chao distance
				if (length(ab.shared) == 0) {d.chao <- 1} else {
					#Which taxa observed in sample a are singletons in sample b and vice versa?
					a_shared.b_singl <- which(a.freq >= 1 & b.freq == 1);
					b_shared.a_singl <- which(a.freq == 1 & b.freq >= 1);
					#How many shared taxa are singletons / doubletons in samples a & b?
					f.a.singl <- length(b_shared.a_singl);
					f.a.doubl <- length(which(a.freq == 2 & b.freq >= 1));
					f.b.singl <- length(a_shared.b_singl);
					f.b.doubl <- length(which(a.freq >= 1 & b.freq == 2));
					if (f.a.doubl == 0) {f.a.doubl = 1}; if (f.b.doubl == 0) {f.b.doubl = 1};
					#Get taxa abundance estimates for samples a & b
					a.adjust <- ((b.N-1)/b.N)*f.b.singl/(2*f.b.doubl)*sum(a.freq[a_shared.b_singl] / a.N);
					a.estimate <- sum(a.freq[ab.shared] / a.N) + a.adjust;
					b.adjust <- ((a.N-1)/a.N)*f.a.singl/(2*f.a.doubl)*sum(b.freq[b_shared.a_singl] / b.N);
					b.estimate <- sum(b.freq[ab.shared] / b.N) + b.adjust;
					if (a.estimate > 1) {a.estimate = 1}; if (b.estimate > 1) {b.estimate == 1};
					d.chao <- 1 - ((a.estimate * b.estimate) / (a.estimate + b.estimate - (a.estimate * b.estimate)));
					if (d.chao < 0) {d.chao = 0}
				}
				#Return
				d.chao
			}))
		}, mc.cores=use.cores);
	}
	########################
	#Classic Bray-Curtis dissimilarity
	if (distance == "bray_curtis") {
		cs.list <- mclapply(ot.count.list, function(a.count) {
			a.N <- sum(a.count); a.list <- names(a.count);
			unlist(lapply(ot.count.list, function(b.count) {
				b.N <- sum(b.count); b.list <- names(b.count);
				ab.shared <- intersect(a.list, b.list); ab.union <- union(a.list, b.list);
				a.freq <- b.freq <- numeric(length=length(ab.union)); names(a.freq) <- names(b.freq) <- ab.union;
				a.freq[a.list] <- a.count; b.freq[b.list] <- b.count;
				1 - (2 * (sum(pmin(a.freq, b.freq)) / (a.N + b.N)));
			}))
		}, mc.cores=use.cores);
	}
	########################
	#Classic Morisita-Horn dissimilarity
	if (distance == "morisita_horn") {
		cs.list <- mclapply(ot.count.list, function(a.count) {
			a.N <- sum(a.count); a.n <- length(a.count); a.list <- names(a.count);
			a.simpson <- sum(a.count ^ 2) / (a.N ^ 2);
			unlist(lapply(ot.count.list, function(b.count) {
				#Get current summary statistics for B
				b.N <- sum(b.count); b.n <- length(b.count); b.list <- names(b.count);
				b.simpson <- sum(b.count ^ 2) / (b.N ^ 2);
				#Get current freq vectors from intersect and union of a & b
				ab.shared <- intersect(a.list, b.list); ab.union <- union(a.list, b.list);
				a.freq <- b.freq <- numeric(length=length(ab.union)); names(a.freq) <- names(b.freq) <- ab.union;
				a.freq[a.list] <- a.count; b.freq[b.list] <- b.count;
				#Get summed abundance product
				ab.prod <- sum(a.freq * b.freq);
				#Return
				1 - ((2 * ab.prod) / ((a.simpson + b.simpson) * a.N * b.N));
			}))
		}, mc.cores=use.cores);
	}
	########################
	
	#Return
	do.call("rbind", cs.list);
}
########################

########################
#Faster UniFrac
#=> see also http://github.com/joey711/phyloseq/issues/524
fasterUniFrac <- function(o.t, tree, distance="UniFrac_unweighted", blocksize=1000, use.cores) {
	#Load packages for parallel processing
	require("foreach");
	require("bigmemory");
	require("Matrix");
	library("doMC", quietly=T);
	#Register cluster
	registerDoMC(cores=use.cores);
	########################
	
	########################
	#Get data
	sp.o.t <- Matrix(o.t, sparse=T); 											#Force sparse <<Matrix>> format
	c.samples <- colnames(sp.o.t);	 											#Assume that OTUs are rows in OTU table
	c.sample.sizes <- colSums(sp.o.t);
	n.sample <- length(c.samples);
	#Match OTU names in OTU table to tree edges
	if( !all(rownames(sp.o.t) == taxa_names(tree)) ){sp.o.t <- sp.o.t[taxa_names(tree), ]}
	#Get N x 2 matrix of pairwise combinations of samples
	spn <- combn(c.samples, 2, simplify=FALSE);
	########################
	
	########################
	#Build pre-requisite matrices
	########################
	# Create a list of descendants, starting from the first internal node (root)
	descList <- prop.part(tree, check.labels = FALSE)
	# Add the terminal edge descendants (tips). By definition, can only have one descendant
	descList <- c(as.list(1:length(tree$tip.label)), descList)
	# Convert `descList` to `edge_array` that matches the order of things in `tree$edge`
	# For each entry in the tree$edge table, sum the descendants for each sample
	# `tree$edge[i, 2]` is the node ID.
	tmp.edge.array <- mclapply(1:nrow(tree$edge), function(i) {colSums(sp.o.t[descList[[tree$edge[i, 2]]], ,drop=F], na.rm=T)}, mc.cores=use.cores);
	edge.array <- do.call("rbind", tmp.edge.array);
	#Get edge occurrence matrix (for unweighted UniFrac)
	edge.occ <- (edge.array > 0) - 0;
	edge.occ.list <- apply((edge.array > 0), 2, which);
	#Get tree in "postorder" order
	z = reorder.phylo(tree, order="postorder");
	#Get "tip ages", or horizontal position of edges in re-ordered phylo object
	tip.ages = ape_node_depth_edge_length(Ntip = length(tree$tip.label), Nnode = tree$Nnode, edge = z$edge, Nedge = nrow(tree$edge)[1], edge.length = z$edge.length);
	#Keep only tips
	tip.ages <- tip.ages[1:length(tree$tip.label)]
	#Rename (when in doubt)
	names(tip.ages) <- z$tip.label;
	#Reorder
	tip.ages <- tip.ages[rownames(sp.o.t)];
	########################
	
	########################
	#Calculate UniFrac distances based on pre-formed arrays
	########################
	#Unweighted UniFrac
	size.chunk <- blocksize;
	if (distance == "UniFrac_unweighted") {
		cs.list <- foreach(i=split(1:length(spn), ceiling(seq_along(1:length(spn))/size.chunk))) %dopar% {
			unlist(lapply(spn[i], function(sp) {
				#Get occurrence vectors for current samples
				a <- edge.occ.list[[sp[1]]]; b <- edge.occ.list[[sp[2]]];
				curr.shared_edges <- intersect(a, b);
				curr.union_edges <- union(a, b);
				#Return: 1 - (shared branch length) / (total branch length)
				#=> formulation as a distance
				1 - sum(tree$edge.length[curr.shared_edges], na.rm=TRUE) / sum(tree$edge.length[curr.union_edges], na.rm=TRUE)
			}));
		}
	}
	########################
	#Weighted UniFrac
	if (distance == "UniFrac_weighted") {
		cs.list <- foreach(i=split(1:length(spn), ceiling(seq_along(1:length(spn))/size.chunk))) %dopar% {
			unlist(lapply(spn[i], function(sp) {
				#Get current samples and total sample sizes
				a <- sp[1]; b <- sp[2];
				a.T <- c.sample.sizes[a]; b.T <- c.sample.sizes[b];
				#Calculate branchweights
				curr.bw <- abs(edge.array[, a]/a.T - edge.array[, b]/b.T);
				#Return: (unshared branch length) / (total branch length)
				#=> formulation as distance
				sum(tree$edge.length * curr.bw, na.rm = TRUE) / sum(tip.ages * (sp.o.t[,a]/a.T + sp.o.t[,b]/b.T), na.rm = TRUE)
			}));
		}
	}
	########################
	
	#Pass output
	#Initialize UniFracMat with NAs
	mat.UF <- matrix(NA_real_, n.sample, n.sample); diag(mat.UF) <- 0;
	rownames(mat.UF) <- colnames(mat.UF) <- c.samples;
	# Matrix-assign lower-triangle of UniFracMat. Then coerce to dist and return.
	mat.idx <- do.call(rbind, spn)[, 2:1];
	#Coerce results into matrices
	mat.UF[mat.idx] <- unlist(cs.list); mat.UF[mat.idx[,2:1]] <- unlist(cs.list);
	mat.UF
}
########################

########################
#Interaction-adjusted indices of community similarity
#=> formulated as distances/dissimilarities
#=> calculated in parallel
community.similarity.corr.par <- function(o.t, S, distance="jaccard.corr.uw", blocksize=1000, use.cores=detectCores()) {
	########################
	#Load packages for parallel processing
	require("foreach");
	require("bigmemory");
	library("doMC", quietly=T);
	#Register cluster
	registerDoMC(cores=use.cores);
	########################
	
	########################
	#Get N x 2 matrix of pairwise combinations of samples
	samples <- colnames(o.t);
	#Turn OTU table into list format
	ot.occ.list <- apply((o.t > 0), 2, which);
	ot.count.list <- mclapply(seq(1,length(ot.occ.list)), function(i) {o.t[ot.occ.list[[i]],i]}, mc.cores=use.cores);
	names(ot.occ.list) <- names(ot.count.list) <- samples;
	########################
	
	########################
	#Unweighted, corrected Jaccard, normalized by sample self-comparisons
	#=> this is for unweighted TINA/PINA
	if (distance == "jaccard.corr.uw.norm") {
		#Pre-calculate sample-wise interaction term sums (normalized by number of OTUs per sample)
		smpl.colsums <- mclapply(ot.occ.list, function(a.list) {colSums(S[a.list, ]) / length(a.list)}, mc.cores=use.cores);
		names(smpl.colsums) <- samples;
		#Pre-calculate sample self-comparisons
		smpl.self <- unlist(lapply(samples, function(a) {a.list <- ot.occ.list[[a]]; sum(smpl.colsums[[a]][a.list]) / length(a.list)}));
		names(smpl.self) <- samples;
		#Calculate similarity index in parallel
		cs.list <- mclapply(samples, function(a) {
			a.sums <- smpl.colsums[[a]]; a.self <- smpl.self[a];
			unlist(lapply(samples, function(b) {
				b.list <- ot.occ.list[[b]]; b.self <- smpl.self[b];
				1 - (sum(a.sums[b.list]) / (length(b.list) * sqrt(a.self*b.self)))
			}));
		}, mc.cores=use.cores);
	}
	########################
	#Weighted, corrected Jaccard, normalized by sample self-comparisons
	#=> this is for weighted TINA/PINA
	if (distance == "jaccard.corr.w.norm") {
		#Pre-calculate relative abundances per sample
		ot.rel_count <- lapply(ot.count.list, function(a.count) {a.count / sum(a.count)});
		#Pre-calculate sample self-comparisons
		smpl.self <- unlist(mclapply(ot.rel_count, function(a.rel) {
			a.list <- names(a.rel); a.S <- S[a.list, a.list];
			sum(a.S * outer(a.rel, a.rel));
		}, mc.cores=use.cores));
		#Iterate through samples in parallel
		cs.list <- mclapply(samples, function(a) {
			#Get current OTU names, relative counts and sample size
			a.rel <- ot.rel_count[[a]]; a.list <- names(a.rel);
			#Get current interaction sub-matrix
			a.S <- a.rel * S[a.list,];
			#Iterate through all potential partner samples
			unlist(lapply(samples, function(b) {
				b.rel <- ot.rel_count[[b]]; b.list <- names(b.rel); curr.S <- a.S[,b.list];
				#Calculate sample-pair weighted interaction term
				1 - (sum(b.rel * t(curr.S)) / sqrt(smpl.self[a] * smpl.self[b]));
			}));
		}, mc.cores=use.cores);
	}
	########################
	
	#Return
	do.call("rbind", cs.list)
}
########################
################################################################################
################################################################################


