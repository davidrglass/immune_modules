##############################################################################
#
# Script: immune_modules.R
# Project: CTCL
# Author: David Glass
# Email: dglass@fredhutch.org
# Date: 8-3-23
#
# This program takes in a table of samples vs features, preprocesses,
# generates immune modules, and runs stats and a random forest classifier
# Identifies clusters of features associated with the response variable
#
##############################################################################

### Instructions:
  ### Starting data is a csv with one or more columns of response variables and
  ### many predictor columns. Example data is provided. Update Inputs below
  ### before starting. You can then run the whole script or if you want to
  ### evaluate outputs step by step, read in all libraries, inputs, and
  ### functions, then run Main line by line. Either way, you will need to
  ### choose the gini cutoff (line 516) based off of the plot (line 515) 
  ###
  ### This has been tested on a matrix of ~30x1500. The
  ### number of rows will probably scale pretty well. The number of columns can
  ### also scale, but you will probably want to adjust hclustTuning to consider
  ### fewer cluster number possibilities (no.clusters argument).

### Pseudocode
  ### start with feature table
  ### specify response factor and non-quantitative variables not used in model
  ### filter high NA features
  ### Impute data for remaining NAs
  ### tune hclust to identify optimal number of clusters
  ### cluster features
  ### Run PCA on each cluster of features to summarize data
  ### tune & run random forest model
  ### visualize and evaluate selected features

##### LIBRARIES ######

require(clusterSim)
require(cluster)
require(ranger)
require(ggbeeswarm)
require(ggplot2)
require(magrittr)
require(data.table)



##### INPUTS #####

### Location of directory for images
images.path <- "images/immune_modules/"
### Location of feature table
example.feature.dat.path <- "tables/example_feature_dat.csv"
### Location to put hclust tuning data generated in this script
hclust.tune.path <- "tables/hclust_tune.csv"
### Non-numeric columns
factors <- c("status", "subject","response", "study")
### The column stratifying patients
response.variable <- "response" 



##### FUNCTIONS #####

filterFeatures <- function(dt=dat,
                           feats=total.features,
                           unique.threshold=15L) {
  # Filters features with many NAs
  # Inputs:
  #   dt - feature dat
  #   feats- numeric features
  #   unique.threshold - minimum number of non-NA observations needed to keep a feature
  # Outputs:
  #   updated dt
  drop.cols <- dt[, lapply(.SD, function(x) sum(!is.na(x))), .SDcols=feats] %>%
    melt() %>%
    .[value<unique.threshold, as.character(variable)]
  
  return(dt[, setdiff(colnames(dt), drop.cols), with=F])
}


imputeFeatures <- function(dt=dat,
                           feats=total.features) {
  # Replaces NAs with mean
  # Inputs:
  #   dt - feature dat
  #   feats- numeric features
  # Outputs:
  #   updated dt
  # Replace NAs with mean
  factors <- setdiff(colnames(dt), feats)
  melted <- melt(dt, id.vars=factors) %>%
    .[, mean:=mean(value, na.rm=T), by=variable]
  melted[is.na(value), value:=mean]
  melted[, mean:=NULL]
  
  updated.dt <- dcast(melted, ...~variable, value.var="value")
  return(updated.dt)
}


hclustTuning <- function(mat=cor.matrix,
                         dm=dist.matrix,
                         no.clusters=seq(5, round(nrow(mat)-5), by=5)) {
  # Quantifies goodness of clustering for number of clusters
  # Inputs:
  #   mat - correlation matrix
  #   dm - distance matrix
  #   no.clusters - sequence of # of clusters to consider
  # Outputs:
  #   scores - data.table with cluster goodness metrics
  scores <- data.table(no.clusters=no.clusters, silhouette=0, db.index=0, correlation=0, diversity=0)
  
  hclust.out <- hclust(d=dm, method="ward.D")
  for (k in no.clusters) {
    print(paste("Number of clusters:", k))
    clusters <- cutree(hclust.out, k=k)
    # silhouette
    sil.sum <- silhouette(clusters, dm) %>%
      summary()
    # DB index
    dbi <- index.DB(dm, clusters, p=2)$DB
    # Intracluster correlation
    means <- c()
    for (i in unique(clusters)) {
      features <- colnames(mat)[clusters==i]
      temp <- mat[features, features]
      if(length(temp)==1) next
      for (feat in features) temp[feat, feat] <- NA # don't bias with correlation with self
      means <- c(means, mean(temp, na.rm=T))
    }
    # Shannon diversity
    cluster.sizes <- table(clusters)
    p <- cluster.sizes/sum(cluster.sizes)
    div <- -sum(p*log(p))
    scores[no.clusters==k, `:=`(silhouette=sil.sum$avg.width,
                                db.index=dbi,
                                correlation=mean(means, na.rm=T),
                                diversity=div)]
    }
  return(scores)
}


visualizeHclustTuning <- function(sc=scores, ip=images.path, hnc=hclust.no.clusters) {
  # Plotting function to visualize goodness of clustering for different number of clusters
  # Inputs:
  #   sc - scores data.table
  #   ip - path to images directory
  #   hnc - integer of selected number of clusters
  # Outputs:
  #   png file
  melted <- melt(sc, id.vars="no.clusters", variable.name="metric", value="score")
  which.max
  ggplot(melted, aes(no.clusters, score)) +
    geom_point() +
    geom_vline(xintercept=hnc) +
    geom_line() +
    theme_bw() +
    facet_wrap(~metric, scales="free")
  ggsave(paste0(ip, "tuning_hclust.png"), width=6, height=6)
}


clusterData <- function(dist=dist.matrix, hnc=hclust.no.clusters, prefix="IM") {
  # Clusters features
  # Inputs:
  #   dist - distance matrix
  #   hnc - integer of selected number of clusters
  # Outputs:
  #   data.table with cluster assigments for each feature
  hclust.out <- hclust(d=dist, method="ward.D")
  clusters <- cutree(hclust.out, k=hnc)
  cluster.dt <- as.data.table(clusters, keep.rownames="feature") %>%
    .[order(clusters)] %>%
    setnames("clusters", "cluster")
  cluster.dt[, cluster:=paste0(prefix, cluster)]
  return(cluster.dt)
}


updateDat <- function(dt=copy(imputed.dat), cd=cluster.dat) {
  # Updates dat with summarized cluster values, which are PC1 of the features in that cluster
  # Inputs:
  #   dt - data.table
  #   clusts - cluster assignments
  #   prefix - character vector to ammend to cluster numbers
  # Outputs:
  #   updated dt
  for (i in unique(cd$cluster)) {
    cluster.features <- cd[cluster==i, feature]
    pca.out <- prcomp(dt[, cluster.features, with=F], center=T, scale.=T)
    dt[, eval(i):=pca.out$x[,1]]
  }
  return(dt)
}


getEffectSize <- function(x, y, hedge=T) {
  # Calculates Hedge's g or Cohen's d effect size
  # Inputs:
  #   x - a numeric vector
  #   y - a numeric vector
  #   hedge - logical, if tru calculate Hedge's g, if false, Cohen's d
  # Ouptut:
  #   the effect size
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]
  x.mean <- mean(x)
  y.mean <- mean(y)
  x.sd <- sd(x)
  y.sd <- sd(y)
  x.n <- length(x) - 1
  y.n <- length(y) - 1
  if (hedge) pooled.sd <- sqrt((x.n*x.sd^2 + y.n*y.sd^2)/(x.n+y.n))
  else pooled.sd <- sqrt((x.sd^2+y.sd^2)/2)
  return((x.mean-y.mean)/pooled.sd)
}


getFeatureStats <- function(dt, cols, comparison=response.variable, paired.test=F) {
  # Gets effect size and wilcoxon p value for each comparison
  # Inputs:
  #   dt - feature data.table
  #   cols - column names to compare
  #   comparison - response variable to compare
  #   paired.test - logical, is it a paired test
  # Outputs:
  #   stats - data.table with stats
  
  if (paired.test) setorder(dt, subject)
  
  labels <- unique(dt[[comparison]])
  stats <- data.table()
  for (col in cols) {
    x <- dt[eval(parse(text=comparison))==labels[2]][[col]]
    y <- dt[eval(parse(text=comparison))==labels[1]][[col]]
    stats <- data.table(feature=col,
                        effect=getEffectSize(x, y),
                        p.value=wilcox.test(x, y, paired=paired.test)$p.value) %>%
      rbind(stats)
  }
  setorder(stats, p.value)
  stats[, q.value:=p.adjust(p.value, method="fdr")]
  return(stats)
}


tuneRF <- function(dt=im.dat, cols=im.features, comparison=response.variable, se=666) {
  # Takes data, runs RF with different settings, and returns hypergrid with oob error for each model
  # Inputs:
  #   dt - data.table
  #   cols - predictive feature
  #   comparison - response column name
  #   se - seed
  # Outputs:
  #   hg - hyper.grid with filled out oob.error column
  hg <- expand.grid(mtry=seq(0.05, 0.5, 0.05),
                    min.node.size=1:5,
                    sample.fraction=seq(0.25, 0.95, 0.05),   
                    replace=c(T, F),
                    oob.error=0) %>%
    as.data.table()
  
  for(i in seq(nrow(hg))) {
    rf <- ranger(y=dt[[comparison]],
                 x=dt[, cols, with=F],
                 mtry=floor(length(cols) * hg$mtry[i]),
                 write.forest=F,
                 min.node.size=hg$min.node.size[i],
                 sample.fraction=hg$sample.fraction[i],
                 probability=T,
                 replace=hg$replace[i],
                 verbose=F,
                 seed=se)
    hg[i, oob.error:=rf$prediction.error]
  }
  hg <- hg[order(oob.error)]
  return(hg)
}


testRF <- function(dt=im.dat, cols=im.features, comparison=response.variable, hg=hyper.grid, seeds=seq(20L)) {
  # Runs tuned RF on different seeds to evaluate accuracy
  # Inputs:
  #   dt - data.table
  #   cols - predictive feature
  #   comparison - response column name
  #   hg - hyper.grid data.table
  #   seeds - sequence of seeds to try
  # Outputs:
  #   results - data.table with seed, oob.error, predictor.rank, predictor.name (for top ten predictors)
  results <- expand.grid(seed=seeds, predictor.rank=1:10) %>%
    as.data.table() %>%
    .[, `:=`(oob.error=0, predictor.name=character())]
  for (i in seeds) {
    rf <- ranger(y=dt[[comparison]],
                 x=dt[, cols, with=F],
                 mtry=floor(length(cols) * hg$mtry[1]),
                 write.forest=F,
                 min.node.size=hg$min.node.size[1],
                 sample.fraction=hg$sample.fraction[1],
                 probability=T,
                 replace=hg$replace[1],
                 importance="impurity",
                 verbose=F,
                 seed=i)
    results[seed==i, oob.error:=rf$prediction.error]
    results[seed==i, predictor.name:=names(sort(rf$variable.importance, decreasing=T))[1:10]]
  }
  return(results)
}


findRepresentativeSeed <- function(accuracy.dt, no.predictors) {
  # Find seed with top X predictors most representative of the cumulative top 10 of all seeds
  # Function is called by runRF()
  # Inputs:
  #   accuracy.dt - data.table with seed, oob.error, predictor.rank, predictor.name 
  #   no.predictors - top X predictors match up
  # Outputs:
  #   se - integer, most representative seed
  top <- accuracy.dt[, sum(10-predictor.rank), by=predictor.name] %>%
    .[order(-V1), predictor.name] %>%
    .[1:no.predictors]
  seeds <- unique(accuracy.dt$seed)
  for (i in no.predictors:10) {
    for (se in seeds) {
      predictors <- accuracy.dt[seed==se & predictor.rank %in% 1:i, predictor.name]
      if(length(intersect(top, predictors))==no.predictors) return(se)
    }
  }
}


runRF <- function(dt=im.dat, cols=im.features, comparison=response.variable, hg=hyper.grid, acc.dt=accuracy.dat, no.pred=5L) {
  # Runs tuned RF on representative seed that has median oob.error
  # Inputs:
  #   dt - data.table
  #   cols - predictive feature
  #   comparison - response column name
  #   hg - hyper.grid data.table
  #   acc.dt - data.table with seed, oob.error, predictor.rank, predictor.name (for top ten predictors)
  # Outputs:
  #   rf - ranger output
  
  ### Select seed with most representative predictors
  final.seed <- findRepresentativeSeed(accuracy.dt=acc.dt, no.predictors=no.pred)
  while(is.null(final.seed)) {
    no.pred <- no.pred - 1
    final.seed <- findRepresentativeSeed(accuracy.dt=acc.dt, no.predictors=no.pred)
  }
  
  rf <- ranger(y=dt[[comparison]],
               x=dt[, cols, with=F],
               mtry=floor(length(cols) * hg$mtry[1]),
               write.forest=T,
               min.node.size=hg$min.node.size[1],
               sample.fraction=hg$sample.fraction[1],
               probability=T,
               replace=hg$replace[1],
               importance="impurity",
               verbose=F,
               seed=final.seed)
  return(rf)
}


plotAccuracy <- function(acc.dt=accuracy.dat, ip=images.path, filename="model_accuracy") {
  # Plots range of OOB errors across seeds on tuned model
  # Inputs:
  #   acc.dt - data.table with seed, oob.error, predictor.rank, predictor.name (for top ten predictors)
  #   ip - path to images directory
  #   filename - name of image
  # Outputs:
  #   png
  acc.dt.lite <- acc.dt[, .(seed, oob.error)] %>%
    unique() %>%
    .[, dummy:="dummy"]
  ggplot(acc.dt.lite, aes(dummy, oob.error)) +
    geom_quasirandom(size=2, width=0.25, color="black") +
    stat_summary(fun=mean, size=0.25, width=0.5, geom="crossbar") +
    theme_bw() +
    labs(x=NULL, y="OOB error") +
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  ggsave(paste0(ip, filename, ".png"), width=1.5, height=2)
}


updateFeatureStats <- function(dt, rf) {
  # Updates feature stats data.table with gini scores from rf model
  # Inputs:
  #   dt -  feature stats data.table
  #   rf - random forest model
  # Outputs:
  #   Updated dt with gini column
  dt <- merge(dt, data.table(feature=names(rf$variable.importance), gini=rf$variable.importance), by="feature") %>%
    .[order(-gini)]
  return(dt)
}


plotGiniDots <- function(rf=rf.out, ip=images.path, co=NULL, filename) {
  # Plots gini vs predictor rank to help decide cutoff for most informative features
  # Inputs:
  #   rf - random forest model
  #   ip - path to images directory
  #   co - gini value cutoff 
  #   filename - images filename
  # Outputs:
  #   eps file
  sorted.predictors <- sort(rf$variable.importance, decreasing=T) %>%
    data.table(module=names(.), gini=., rank=1:length(.))
  ggplot(sorted.predictors, aes(rank, gini)) +
    geom_point() +
    geom_hline(yintercept=co, linetype="dashed", color="red") +
    theme_bw() +
    labs(x="Immune module rank", y="Variable importance") +
    theme(panel.grid.minor=element_blank())
  ggsave(paste0(ip, filename, ".png"), width=3, height=3)
}




plotTopPredictorBars <- function(tp=copy(top.predictors),
                                 fs=feature.stats,
                                 cd=copy(cluster.dat),
                                 ip=images.path,
                                 filename="gini_bars",
                                 w=6,
                                 h=5) {
  # Plots bars of effect size of underlying features in selected immune modules
  # Inputs:
  #   tp - data.table of top predictors
  #   fs - data.table of feature statistics
  #   cd - data.table of cluster names and underlying feature names
  #   ip - path to images directory
  #   filename - filename for image
  #   w - image width
  #   h - image height
  # Outputs:
  #   eps
  setnames(tp, "feature", "cluster")
  merged <- merge(tp, cd, by="cluster") %>%
    merge(fs, by="feature") %>%
    .[, abs.effect:=abs(effect)]
  setorder(merged, gini, abs.effect)
  merged[, feature:=factor(feature, feature)]
  
  ggplot(merged, aes(feature, effect, fill=cluster)) +
    geom_col() +
    theme_bw() +
    theme(panel.grid.major.y=element_blank()) +
    coord_flip() +
    labs(fill=NULL)
  ggsave(paste0(ip, filename, ".png"), width=w, height=h)
} 



##### MAIN #####

if (!dir.exists(images.path)) dir.create(images.path, recursive=T)

dat <- fread(example.feature.dat.path) %>%
  .[, (factors):=lapply(.SD, as.factor), .SDcols=factors]

total.features <- setdiff(colnames(dat), factors)

### Filter and impute NAs
filtered.dat <- filterFeatures()
imputed.dat <- imputeFeatures()

features <- setdiff(colnames(filtered.dat), factors)

### Build distance matrix
cor.matrix <- cor(imputed.dat[, features, with=F], method="spearman")
dist.matrix <- as.dist(1-(cor.matrix)^4)

### Identify optimal number of clusters (go get a coffee maybe)
scores <- hclustTuning()
fwrite(scores, hclust.tune.path)
## Selects highest silhouette score - could be adjusted with other parameters
hclust.no.clusters <- scores[which.max(silhouette), no.clusters] 
visualizeHclustTuning()

### Run optimal clustering
cluster.dat <- clusterData()

### Update dat with cluster PCA values
im.dat <- updateDat()
im.features <- setdiff(colnames(im.dat), c(factors, features))

### Get stats for each feature and each immune module
feature.stats <- getFeatureStats(dt=filtered.dat, cols=features, paired.test=F)
im.feature.stats <- getFeatureStats(dt=im.dat, cols=im.features, paired.test=F)

### Build model
hyper.grid <- tuneRF()
accuracy.dat <- testRF()
rf.out <- runRF()

### Plot accuracy
plotAccuracy()
im.feature.stats <- updateFeatureStats(dt=im.feature.stats, rf=rf.out)

### Select cutoff and identify top immune modules
plotGiniDots(filename="dotplot_gini")
gini.co <- 0.017
plotGiniDots(co=gini.co, filename="dotplot_gini_with_cutoff")
top.predictors <- im.feature.stats[gini>gini.co, .(feature, gini)]
plotTopPredictorBars()
