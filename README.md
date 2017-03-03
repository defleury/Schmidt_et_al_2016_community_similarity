# Analysis code: A Family of Interaction-Adjusted Indices of Community Similarity

This repository contains the R code for the analyses conducted in the study of Schmidt et al., "A Family of Interaction-Adjusted Indices of Community Similarity". The manuscript is available as a [preprint on biorXiv](http://biorxiv.org/content/early/2016/02/18/040097), through [PubMed Central](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5322292/) and open access via the [ISME Journal](http://www.nature.com/ismej/journal/v11/n3/full/ismej2016139a.html). The full reference is:

Schmidt TSB, Matias Rodrigues FM & von Mering C, A family of interaction-adjusted indices of community similarity, ISME J (2017) 11:791-807, [doi:10.1038/ismej.2016.139](doi:10.1038/ismej.2016.139)

## Study abstract:

"Interactions between taxa are essential drivers of ecological community structure and dynamics, but they are not taken into account by traditional indices of diversity. In this study, we propose a novel family of indices that quantify community similarity in the context of taxa interaction networks. Using publicly available datasets, we assess the performance of two specific indices which are Taxa INteraction-Adjusted (TINA, based on taxa co-occurrence networks), and Phylogenetic INteraction-Adjusted (PINA, based on phylogenetic similarities). TINA and PINA outperformed traditional indices when partitioning human-associated microbial communities according to habitat, even for extremely downsampled datasets, and when organising ocean micro-eukaryotic plankton diversity according to geographical and physicochemical gradients. We argue that interaction-adjusted indices capture novel aspects of diversity outside the scope of traditional approaches, highlighting the biological significance of ecological association networks in the interpretation of community similarity."

## Data availability

The study was conducted based on re-proceesed publicly available data from the [Human Microbiome Project](http://hmpdacc.org) and [TARA Oceans](http://oceans.taraexpeditions.org/en/#). Moreover, this repository contains an additional analysis script on data from the global [Reef Life Survey](http://reeflifesurvey.com). Re-processed data is available via the website [meringlab.org](http://meringlab.org/suppdata/2016-community_similarity/) and in the folder [data](data) in this repository.
In particular, these datasets were used:

* **HMP**: 16S amplicon data, v35, all body sites (raw data available through [hmpdacc.org](http://hmpdacc.org/HMR16S/))
* **TARA Oceans**: 18S rRNA, v9, targeting eukaryotic micro-plankton; [taxa count table](http://doi.pangaea.de/10.1594/PANGAEA.843022) and [raw sample metadata](http://doi.pangaea.de/10.1594/PANGAEA.843017) available online through PANGAEA.
* **Reef Life Survey**: census data on reef-dwelling marine animals, accessed Jan 2016; see the data descriptor in [Nature Scientific Data](http://www.nature.com/articles/sdata20147) for more details.

## Analysis code

The code is organised to re-generate analysis underlying the various figures in the publication. Prior to running any analysis, edit any script to insert the correct paths to data and results folders etc. Moreover, you will need to `source` the script [functions.community_similarity.R](functions.community_similarity.R) (done automatically at the beginning of each script). The community similarity matrices on which most analyses rely are generated in the script [prepare.community_similarity.R](prepare.community_similarity.R).

The code was deposited as-is, so there is no guarantee that scripts will run on any system. In particular, some computations (e.g., SparCC correlation networks) require significant computational resources; parallelization was tested on an Ubuntu system only (it will probably fail on Windows machines). We are currently working to provide more efficient and versatile versions of the TINA/PINA code and will hopefully be able to do so in the near future.

Contact:
Sebastian Schmidt (sebastian.schmidt [at] embl.de)
