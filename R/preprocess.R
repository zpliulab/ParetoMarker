# Main preprocessing adapter, derived from the recovered DE.R and net.R.
# PCA_CMI.R retains the supplied implementation; see docs/IMPLEMENTATION.md.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop("Usage: Rscript R/preprocess.R config/preprocess.example.R")
cfg <- dget(args[1]); if(dir.exists(cfg$output_dir)) stop("Use a fresh output directory")
for(pkg in c("impute","limma","infotheo")) if(!requireNamespace(pkg,quietly=TRUE)) stop(paste("Missing R package",pkg))
source("R/PCA_CMI.R"); set.seed(cfg$seed)
read_csv <- function(p) read.csv(p,check.names=FALSE,stringsAsFactors=FALSE)
d <- read_csv(cfg$expression); genes <- as.character(d[[1]])
if(anyNA(genes)||anyDuplicated(genes)||any(genes=="")) stop("Unique nonempty genes required")
x <- as.matrix(d[,-1,drop=FALSE]); storage.mode(x) <- "double";rownames(x)<-genes
m <- read_csv(cfg$metadata)
if(anyDuplicated(m$sample_id)||anyDuplicated(colnames(x))) stop("Duplicate sample IDs")
i <- match(colnames(x),m$sample_id);if(anyNA(i)) stop("Unmatched metadata")
m <- m[i,,drop=FALSE]
if(any(is.infinite(x))) stop("Infinite expression values")
if(cfg$missing_fraction<0||cfg$missing_fraction>=1) stop("missing_fraction must be in [0,1)")
x <- x[rowMeans(is.na(x))<=cfg$missing_fraction,,drop=FALSE]
if(nrow(x)<=cfg$knn_k) stop("Too few genes for selected KNN k")
if(anyNA(x)) x <- impute::impute.knn(x,k=cfg$knn_k,rng.seed=cfg$seed)$data
x <- round(x,3)
if(any(!is.finite(x))) stop("Nonfinite values after imputation")
minmax <- function(z) {
 z<-as.matrix(z); lo<-apply(z,1,min);span<-apply(z,1,max)-lo
 span[span==0]<-1; (z-lo)/span
}
# Optional limma step is part of candidate/target construction, not downstream validation.
if(cfg$target_mode=="limma") {
 use <- if(is.null(cfg$de_subset_state)) rep(TRUE,nrow(m)) else m$state==cfg$de_subset_state
 mm<-m[use,,drop=FALSE];g<-as.character(mm[[cfg$de_group_column]])
 keep<-g %in% c(cfg$de_reference,cfg$de_comparison);mm<-mm[keep,,drop=FALSE];g<-g[keep]
 xx<-x[,use,drop=FALSE][,keep,drop=FALSE]
 group<-factor(g,levels=c(cfg$de_reference,cfg$de_comparison),labels=c("Reference","Comparison"))
 if(any(table(group)<2)) stop("At least two samples per DE group required")
 if(cfg$de_scale=="minmax") xx<-minmax(xx) else if(cfg$de_scale!="none") stop("Unknown DE scale")
 if(is.null(cfg$de_block_column)) design<-model.matrix(~0+group)
 else {block<-factor(mm[[cfg$de_block_column]]);design<-model.matrix(~0+group+block)}
 if(qr(design)$rank!=ncol(design)||nrow(design)<=ncol(design)) stop("DE design is not estimable")
 contrast<-rep(0,ncol(design));contrast[match("groupComparison",colnames(design))]<-1
 contrast[match("groupReference",colnames(design))]<- -1
 fit<-limma::eBayes(limma::contrasts.fit(limma::lmFit(xx,design),contrast))
 de<-limma::topTable(fit,n=Inf,adjust.method="BH",sort.by="none")
 targets<-rownames(de)[!is.na(de$adj.P.Val)&de$adj.P.Val<cfg$de_fdr&abs(de$logFC)>cfg$de_abs_logfc]
} else if(cfg$target_mode=="supplied") targets<-as.character(read_csv(cfg$targets)$gene)
else stop("target_mode must be supplied or limma")
prior<-read_csv(cfg$network)
if(!all(c("source","target") %in% names(prior))) stop("Network needs source,target")
if(!is.null(cfg$candidate_genes)) x<-x[rownames(x) %in% read_csv(cfg$candidate_genes)$gene,,drop=FALSE]
prior<-unique(prior[,c("source","target")]);prior<-prior[prior$source!=prior$target & prior$source %in% rownames(x) & prior$target %in% rownames(x),,drop=FALSE]
x<-x[rownames(x) %in% union(prior$source,prior$target),,drop=FALSE]
if(nrow(x)<2||nrow(prior)==0) stop("No usable candidate network")
G<-matrix(0,nrow(x),nrow(x));a<-match(prior$source,rownames(x));b<-match(prior$target,rownames(x));G[cbind(a,b)]<-1
refined<-PCA_CMI(as.data.frame(minmax(x)),lamda=cfg$cmi_threshold,order0=cfg$cmi_order,G=G)$G
prior<-prior[(refined*G)[cbind(a,b)]>0,,drop=FALSE]
mi_state <- function(state) {
 z<-x[,m$state==state,drop=FALSE];if(ncol(z)<3) stop("Too few samples per MI state")
 discrete<-infotheo::discretize(as.data.frame(t(minmax(z))),disc="equalwidth",nbins=max(2,floor(sqrt(ncol(z)))))
 infotheo::mutinformation(discrete,method="emp")
}
ma<-mi_state(cfg$state_a);mb<-mi_state(cfg$state_b)
a<-match(prior$source,rownames(x));b<-match(prior$target,rownames(x))
prior$capacity<-abs(ma[cbind(a,b)]-mb[cbind(a,b)])
if(any(!is.finite(prior$capacity))) stop("Invalid MI capacities")
sources<-intersect(as.character(read_csv(cfg$sources)$gene),prior$source)
targets<-intersect(targets,prior$target)
if(length(sources)==0||length(targets)==0) stop("No source or target genes survive; do not relax thresholds after looking at results")
use<-m$state==cfg$optimization_state
labels<-as.character(m$label[use]);if(anyNA(labels)||!all(labels %in% c("0","1"))) stop("Optimization labels must be 0/1")
dir.create(cfg$output_dir,recursive=TRUE)
write.csv(prior,file.path(cfg$output_dir,"weighted_network.csv"),row.names=FALSE)
write.csv(data.frame(gene=sources),file.path(cfg$output_dir,"sources.csv"),row.names=FALSE)
write.csv(data.frame(gene=targets),file.path(cfg$output_dir,"targets.csv"),row.names=FALSE)
write.csv(data.frame(sample_id=m$sample_id[use],t(x[,use,drop=FALSE]),check.names=FALSE),file.path(cfg$output_dir,"expression.csv"),row.names=FALSE)
write.csv(data.frame(sample_id=m$sample_id[use],label=as.integer(labels)),file.path(cfg$output_dir,"labels.csv"),row.names=FALSE)
dput(cfg,file.path(cfg$output_dir,"config_used.R"));capture.output(sessionInfo(),file=file.path(cfg$output_dir,"R_session.txt"))
