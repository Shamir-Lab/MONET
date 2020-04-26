SUBTYPES.DATA = list(
  list(name='aml', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='AML'),
  list(name='breast', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='BIC'),
  list(name='colon', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='COAD'),
  list(name='gbm', only.primary=T, is.rna.seq=F, is.mirna.seq=F, display.name='GBM'),
  list(name='kidney', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='KIRC'),
  list(name='liver', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='LIHC'),
  list(name='lung', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='LUSC'),
  list(name='melanoma', only.primary=F, is.rna.seq=T, is.mirna.seq=T, display.name='SKCM'),
  list(name='ovarian', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='OV'),
  list(name='sarcoma', only.primary=T, is.rna.seq=T, is.mirna.seq=T, display.name='SARC'))

ALGORITHMS.DATA = list(
  list(name='nemo', display.name='NEMO', col='black', is_monet=F),
  list(name='snf', display.name='SNF', col='red', is_monet=F),
  list(name='mofa', display.name='MOFA', col='brown', is_monet=F),
  list(name='mdi', display.name='MDI', col='green', is_monet=F),
  list(name='bcc', display.name='BCC', col='purple', is_monet=F),
  list(name='clusternomics', display.name='Clusternomics', col='blue', is_monet=F),
  list(name='twl', display.name='TWL', col='pink', is_monet=F),
  list(name='monet.em.repeats', display.name='Monet', col='grey', is_monet=T),
  list(name='monet.con80.shift2', display.name='Monet', col='grey', is_monet=T),
  list(name='monet.repeats', display.name='Monet_repeats', col='blue', is_monet=T),
  list(name='monet.repeats50', display.name='Monet_repeats100', col='yellow', is_monet=T, is_slow=T))

ALGORITHM.NAMES = sapply(ALGORITHMS.DATA, function(x) x$name)
names(ALGORITHMS.DATA) = ALGORITHM.NAMES

run.monet.repeats50 <- function(...) {
  run.monet.with.params(num.repeats=50, shift=0.2, frac.sampled=0.8, ...)
}

run.monet.repeats <- function(...) {
  run.monet.with.params(num.repeats=1, shift=0.2, frac.sampled=0.8, ...)
}

run.monet.con80.shift2 <- function(...) {
  run.monet.with.params(num.repeats=1, shift=0.2, frac.sampled=0.8, ...)
}

run.monet.em.repeats <- function(...) {
  run.monet.with.params(num.repeats=15, shift=0, ...)
}


#######################
##### Simulations #####
#######################

get.mod.in.omic <- function(omic.dim, nsamples, mod.feats=NULL, feats.per.mod=NULL, sigma=NULL) {
  if (!is.null(sigma)) {
    sigma = diag(omic.dim) * sigma
  }
  if (is.null(mod.feats)) {
    mod.omic = do.call(rbind, lapply(1:nsamples, function(i) {
      mod.omic.mean = rep(0, omic.dim)
      mod.feats = sample(1:omic.dim, feats.per.mod)
      mod.omic.mean[mod.feats] = 1
      mod.omic = rmvnorm(1, mod.omic.mean, sigma=sigma)
      return(mod.omic)
    }))
  } else {
    mod.omic.mean = rep(0, omic.dim)
    mod.omic.mean[mod.feats] = 1
    mod.omic = rmvnorm(nsamples, mod.omic.mean, sigma=sigma)
  }
  return(mod.omic)
}

create.simulation1 <- function() {
  samples.per.mod = 60
  dim1 = 500
  dim2 = 500

  feats.per.mod = 125
  mod1.feats = 1:feats.per.mod
  mod2.feats = 1:feats.per.mod
  mod3.feats1 = (feats.per.mod + 1):(2 * feats.per.mod)
  mod3.feats2 = (feats.per.mod + 1):(2 * feats.per.mod)
  mod4.feats1 = (2 * feats.per.mod + 1):(3 * feats.per.mod)
  mod4.feats2 = (2 * feats.per.mod + 1):(3 * feats.per.mod)
  mod5.feats1 = (3 * feats.per.mod + 1):(4 * feats.per.mod)
  mod5.feats2 = (3 * feats.per.mod + 1):(4 * feats.per.mod)

  mod1.omic1 = get.mod.in.omic(dim1, samples.per.mod, mod1.feats)
  mod1.omic2 = get.mod.in.omic(dim2, samples.per.mod, feats.per.mod=feats.per.mod, sigma=4)

  mod2.omic2 = get.mod.in.omic(dim2, samples.per.mod, mod2.feats)
  mod2.omic1 = get.mod.in.omic(dim1, samples.per.mod, feats.per.mod=feats.per.mod, sigma=4)

  mod3.omic1 = get.mod.in.omic(dim1, samples.per.mod, mod3.feats1)
  mod3.omic2 = get.mod.in.omic(dim2, samples.per.mod, mod3.feats2)

  mod4.omic1 = get.mod.in.omic(dim1, samples.per.mod, mod4.feats1)
  mod4.omic2 = get.mod.in.omic(dim2, samples.per.mod, mod4.feats2)

  mod5.omic1 = get.mod.in.omic(dim1, samples.per.mod, mod5.feats1)
  mod5.omic2 = get.mod.in.omic(dim2, samples.per.mod, mod5.feats2)

  outlier.omic1 = get.mod.in.omic(dim1, 5, feats.per.mod=0, sigma=0.1)
  outlier.omic2 = get.mod.in.omic(dim2, 5, feats.per.mod=0, sigma=0.1)

  omic1 = as.data.frame(t(do.call(rbind, list(mod1.omic1, mod2.omic1, mod3.omic1, mod4.omic1, mod5.omic1, outlier.omic1))))
  omic2 = as.data.frame(t(do.call(rbind, list(mod1.omic2, mod2.omic2, mod3.omic2, mod4.omic2, mod5.omic2, outlier.omic2))))
  sample.names = paste0('sample', 1:ncol(omic1))
  colnames(omic1) = sample.names
  colnames(omic2) = sample.names
  rownames(omic1) = paste0('feat', 1:nrow(omic1))
  rownames(omic2) = paste0('feat', 1:nrow(omic2))
  labels = c(unlist(lapply(1:5, function(i) rep(i, samples.per.mod))), 6:10)
  return(list(list(omic1, omic2), labels))
}

create.simulation2 <- function() {
  num.dim = 500
  feats.per.mod = 100
  feats.per.mod2 = 20
  feats.per.mod3 = 40
  mod1.feats = 1:feats.per.mod
  mod234.feats = (feats.per.mod + 1):(2 * feats.per.mod)
  mod2.feats2 = (feats.per.mod + 1):(2 * feats.per.mod)
  mod3.feats2 = (2 * feats.per.mod + 1):(3 * feats.per.mod)
  mod4.feats2 = (3 * feats.per.mod + 1):(4 * feats.per.mod)

  num.omics = 2
  num.mods = 5
  samples.per.mod = 30
  feats = c(feats.per.mod2, feats.per.mod3)
  all.data = lapply(1:num.omics, function(omic) {
    cur.omic = t(do.call(rbind, lapply(1:num.mods, function(mod) {
      mod.in.omic = get.mod.in.omic(num.dim, samples.per.mod, (1 + (mod - 1) * feats[omic]):(mod * feats[omic]))
      return(mod.in.omic)
    })))
    colnames(cur.omic) = paste0('sample', 1:ncol(cur.omic))
    rownames(cur.omic) = paste0('feat', 1:nrow(cur.omic))
    return(cur.omic)
  })

  mod1.omic1 = get.mod.in.omic(num.dim, samples.per.mod, mod1.feats)
  other.mods.omic1 = get.mod.in.omic(num.dim, samples.per.mod * (num.mods - 1), mod234.feats)
  omic1 = as.data.frame(t(do.call(rbind, list(mod1.omic1, other.mods.omic1))))
  sample.names = paste0('sample', 1:ncol(omic1))
  colnames(omic1) = sample.names
  rownames(omic1) = paste0('feat', 1:nrow(omic1))
  labels = unlist(lapply(1:num.mods, function(i) rep(i, samples.per.mod)))
  return(list(c(list(omic1), all.data), labels))

  sample.names = paste0('sample', 1:ncol(omic1))
  colnames(omic1) = sample.names
  colnames(omic2) = sample.names
  colnames(omic3) = sample.names
  rownames(omic1) = paste0('feat', 1:nrow(omic1))
  rownames(omic2) = paste0('feat', 1:nrow(omic2))
  rownames(omic3) = paste0('feat', 1:nrow(omic3))
  labels = c(rep(1, num.samples.mod1), rep(2, num.samples.mod2), rep(3, num.samples.mod3), rep(4, num.samples.mod4))
  return(list(list(omic1, omic2, omic3), labels))
}

#########################################
###### Code to run other methods ########
#########################################

run.nemo <- function(omics.list, subtype.data, num.clusters=NULL, num.neighbors=NA, num.clusters.per.omic=NULL, dist.without.na=F) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  clustering = nemo.clustering(omics.list, num.clusters, num.neighbors, dist.without.na=dist.without.na)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}
 
run.spectral <- function(omics.list, subtype.data, num.clusters=NULL) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, 
                                 filter.var = T)
  subtype = subtype.data$name
  concat.omics = do.call(rbind, omics.list)
  
  similarity.data = affinityMatrix(dist2(as.matrix(t(concat.omics)),
                                         as.matrix(t(concat.omics))), 20, 0.5)
  if (is.null(num.clusters)) {
    num.clusters = estimateNumberOfClustersGivenGraph(similarity.data, 2:15)[[3]]  
  }
  clustering = spectralClustering(similarity.data, num.clusters)
  names(clustering) = colnames(omics.list[[1]])
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

run.snf <- function(omics.list, subtype.data, num.clusters=NULL, num.clusters.per.omic=NULL) {
  if (length(omics.list) == 1) {
    return(run.spectral(omics.list, subtype.data, num.clusters))
  }
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data)
  subtype = subtype.data$name
  alpha=0.5
  T.val=30
  num.neighbors = round(ncol(omics.list[[1]]) / 10)
  similarity.data = lapply(omics.list, function(x) {affinityMatrix(dist2(as.matrix(t(x)),as.matrix(t(x))), 
                                                                   num.neighbors, alpha)})
  if (length(similarity.data) == 1) {
    W = similarity.data[[1]]
  } else {
    W = SNF(similarity.data, num.neighbors, T.val)  
  }
  
  if (is.null(num.clusters)) {
    num.clusters = estimateNumberOfClustersGivenGraph(W, 2:15)[[3]]  
  }
  clustering = spectralClustering(W, num.clusters)
  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))

}

run.clusternomics <- function(omics.list, subtype.data, num.clusters=NULL, num.clusters.per.omic=NULL) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, normalize=T, filter.var=T)
  omics.list.trans = lapply(omics.list, t)
  if (is.null(num.clusters)) {
    num.clusters = seq(5, 30, by=5)
  }
  if (is.null(num.clusters.per.omic)) {
    num.clusters.per.omic = list(rep(3, length(omics.list)), rep(5, length(omics.list)))
  }

  time.taken.normalization = as.numeric(Sys.time() - start, units='secs')

  all.param.options = expand.grid(1:length(num.clusters), 1:length(num.clusters.per.omic))
  all.rets = mclapply(1:nrow(all.param.options), function(i) {
    set.seed(50 + i)
    cur.num.clusters = num.clusters[all.param.options[i, 1]]
    cur.num.clusters.per.omic = num.clusters.per.omic[[all.param.options[i, 2]]]
    start = Sys.time()
    if (length(cur.num.clusters.per.omic) == 1) {
      num.clusters.per.omic = rep(num.clusters.per.omic, length(omics.list))
    }
    cluster.counts = list(global=cur.num.clusters, context=cur.num.clusters.per.omic)

    # Hack due to a false assertion in clusternomics package
    if (length(omics.list) > 2) {
      cluster.counts = c(cluster.counts, rep('UNUSED', length(omics.list) - 2))
    }
    # parameters are as used in the clusternomics publication.
    results = contextCluster(omics.list.trans, cluster.counts, maxIter=1e4, burnin=5e3, lag=3, dataDistributions='diagNormal', verbose=T)
    cur.clustering = results$samples[[length(results$samples)]]$Global
    dic = results$DIC
    time.taken.per.param = as.numeric(Sys.time() - start, units='secs')
    return(list(clustering=cur.clustering, dic=dic, timing=time.taken.per.param, clusternomics.ret=results))
  }, mc.cores=CONFIG$ncores)

  best.sol.index = which.min(sapply(all.rets, function(x) x$dic))
  wall.timing = max(sapply(all.rets, function(x) x$timing)) + time.taken.normalization
  return(list(clustering=all.rets[[best.sol.index]]$clustering, timing=wall.timing, all.rets=all.rets))
}

run.mdi <- function(omics.list, subtype.data, num.clusters=NULL) {
  start = Sys.time()
  subtype = subtype.data$name
  bin.path = CONFIG$mdi.path
  tmp.mdi.dir = CONFIG$tmp.mdi.dir
  dir.create(tmp.mdi.dir)
  omics.list = log.and.normalize(omics.list, subtype.data, normalize=T, filter.var=T)
  omics.trans = lapply(omics.list, function(omic) t(as.matrix(omic)))
  time.taken = as.numeric(Sys.time() - start, units='secs')
  omic.paths = sapply(1:length(omics.trans), function(i) {
    omic = omics.trans[[i]]
    omic.path = file.path(tmp.mdi.dir, sprintf('%s_%d.csv', subtype, i))
    write.csv(omic, file=omic.path)
    return(omic.path)
  })
  output.path = file.path(tmp.mdi.dir, sprintf('%s_output.csv', subtype))
  timing.path = file.path(tmp.mdi.dir, sprintf('%s_timing.txt', subtype))
  # running mdi through bash to redirect stderr, which contains the timing
  if (!file.exists(output.path)) {
    command = paste('bash -c "', bin.path, 'N', do.call(paste, as.list(omic.paths)), '--benchmark', '--seed 42', '>', output.path, '2>', timing.path, '"')
    command.return = system(command)
    stopifnot(command.return == 0)
  }

  mdi.runtime = strsplit(readLines(timing.path), split=',')[[1]][2]
  parsed.runtime = strsplit(mdi.runtime, split=' ')[[1]]
  stopifnot(parsed.runtime[2] == 'ms')
  mdi.ms = as.numeric(parsed.runtime[1])
  mdi.seconds = mdi.ms / 1000

  start2 = Sys.time()
  mcmc = readMdiMcmcOutput(output.path)
  cpsm = generateConsensusPSM(tail(mcmc))
  if (is.null(num.clusters)) {
    num.clusters = apply(getClustersOccupied(tail(mcmc,101)),2,median)
  } else {
    num.clusters = as.list(rep(num.clusters, length(omics.list)))
  }
  cp = extractPSMClustPartition(cpsm, num.clusters, omics.trans)
  # In order to create a global clustering solution, we look at the cartesian product of omic specific clusters.
  clustering = as.numeric(as.factor(apply(cp, 1, function(x) do.call(paste, c(as.list(x), sep='_')))))
  
  time.taken2 = as.numeric(Sys.time() - start2, units='secs')
  return(list(clustering=clustering, timing=time.taken + time.taken2 + mdi.seconds, mdi.ret=list(mcmc, cp)))
}

run.bcc <- function(omics.list, subtype.data, num.clusters=NULL) {
  start = Sys.time()
  # The original publication used data with a small number of features, so no feature filtering was required.
  # But here we do feature filtering.
  omics.list = log.and.normalize(omics.list, subtype.data, normalize=T, filter.var=T, num.feats=500)

  # If a feature has zero variance within a cluster, bcc crashes. Filtering features with only few different
  # values make this more rare
  #omics.list = filter.features.with.few.different.values(omics.list)
  omics.list = lapply(omics.list, function(omic) omic + matrix(rnorm(ncol(omic) * nrow(omic), mean=0, sd=1e-5), ncol=ncol(omic)))
  time.taken.normalization = as.numeric(Sys.time() - start, units='secs')
  if (is.null(num.clusters)) {
    num.clusters = 2:15
  }
  all.rets = mclapply(num.clusters, function(K) {
    set.seed(50 + K)
    start = Sys.time()
    #IndinAlpha means that each omic has a different tendency to adhere to the global clustering.
    
    bayes.ret = NULL
    for (i in 1:1e2) {
      print(paste('running bcc, K is', K, ' attemp number', i))
      bayes.ret = tryCatch({
        bayesCC(omics.list, K = K, IndivAlpha=T, maxiter=1e4) 
      }, error=function(er) {
        return(NULL)
      })
      if (!is.null(bayes.ret)) {
        break
      }
    }
    if (is.null(bayes.ret)) {
      time.taken.per.param = as.numeric(Sys.time() - start, units='secs')
      return(list(mean.adherence=-Inf, timing=time.taken.per.param))
    }
    cur.clustering = unlist(apply(bayes.ret$Cbest, 1, which.max))
    time.taken.per.param = as.numeric(Sys.time() - start, units='secs')
    mean.adherence = mean((K * bayes.ret$Alpha - 1) / (K - 1))
    return(list(clustering=cur.clustering, timing=time.taken.per.param, bayes.ret=bayes.ret, mean.adherence=mean.adherence))
  }, mc.cores=CONFIG$ncores)

  # should choose solution with highest mean adherence. Note that this approach encourages a
  # small number of clusters with consistency between the different omics.
  best.sol.index = which.max(sapply(1:length(all.rets), function(i) all.rets[[i]]$mean.adherence))
  wall.timing = time.taken.normalization + max(sapply(all.rets, function(x) x$timing))
  return(list(clustering=all.rets[[best.sol.index]]$clustering, timing=wall.timing, all.rets=all.rets))
}

# This is an adaptation of the post_analy_clus method from the TWL package
# so that it will work with a single omic.
my.twl.clust <- function (outpu_new, num_clusts, pdf_path) {
    pdf(pdf_path)
    post_lab <- list()
    for (b in seq_along(outpu_new)) {
	par(cex.lab = 1.6)
        hclus1 <- hclust(as.dist(1 - (outpu_new[[b]])))
	plot(hclus1, cex.lab = 1.6, main = "", labels = FALSE,
	            xlab = "Sample")
	hclus1_rec <- rect.hclust(hclus1, k = num_clusts[b])
        post_lab[[b]] <- matrix(rep("clus_unknown", 2 * (dim(outpu_new[[b]])[1])),
            ncol = 2)
        hclus1_rec <- hclus1_rec[!sapply(hclus1_rec, is.null)]
        for (kk in 1:length(hclus1_rec)) {
            post_lab[[b]][hclus1_rec[[kk]], 2] <- paste0("clus_", kk)
        }
    }
    dev.off()
    return(post_lab[[1]])
}

run.twl <- function(omics.list, subtype.data, num.clusters=NULL, num.clusters.per.omic=NULL) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, normalize=T, filter.var=T, num.feats=500)
  omics.list.trans = lapply(omics.list, t)
  omics.annotation = lapply(1:length(omics.list), function(i) data.table(nam=1:ncol(omics.list[[1]])))
  # we use the default alpha since the size of the datasets we use are in the same
  # order of magnitude as those used in the twl paper. We set beta to maintain 20% label information sharing.
  twl.ret = TWLsample(omics.list.trans, omics.annotation, beta_re=(0.4 / 3) * length(omics.list), num_its=1e4, output_every=50)
  time.taken1 = as.numeric(Sys.time() - start, units='secs')
  posterior.mat = pairwise_clus(twl.ret, BURNIN=2000) 
  time.taken2 = as.numeric(Sys.time() - start, units='secs')

  avg.mat = Reduce(`+`, posterior.mat) / length(posterior.mat)
  if (is.null(num.clusters)) {
    num.cls = SNFtool::estimateNumberOfClustersGivenGraph(avg.mat, 2:15)[[1]]
  } else {
    num.cls = num.clusters
  }
  clustering = as.numeric(as.factor(my.twl.clust(list(avg.mat), num.cls, file.path(CONFIG$twl.tmp.dir, 'twl.pdf'))[,2]))
  names(clustering) = colnames(omics.list[[1]])

  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))

  # Our first implementation, changed due to bad results:
  #if (is.null(num.clusters)) {
  #  num.cls = sapply(posterior.mat, function(post.mat) SNFtool::estimateNumberOfClustersGivenGraph(post.mat, 2:15)[[1]])
  #} else {
  #  num.cls = num.clusters
  #}
  #clustering.per.omic = post_analy_clus(posterior.mat, twl.ret, num.cls, titles=paste0('omic', 1:length(omics.list)), pdf_path=file.path(CONFIG$twl.tmp.dir, 'twl.pdf'))
  #time.taken3 = as.numeric(Sys.time() - start, units='secs')
  #clustering = do.call(paste, c(lapply(clustering.per.omic, function(x) x[[2]]), sep=','))
  #names(clustering) = colnames(omics.list[[1]])
  #time.taken = as.numeric(Sys.time() - start, units='secs')
  #save(twl.ret, posterior.mat, clustering.per.omic, clustering, 
  #     time.taken1, time.taken2, time.taken3, time.taken, 
  #     file=cache.file.path)
}

get.clustering.silhouette <- function(raw.data, clustering) {
  sils = c()
  for (i in 1:length(raw.data)) {
    x = raw.data[[i]]
    distmatrix = dist2(as.matrix(t(x)),as.matrix(t(x)))
    sil = silhouette(clustering, dmatrix = distmatrix)[,3]
    sils = c(sils, mean(sil))
  }
  return(mean(sils))
}

run.mofa <- function(omics.list, subtype.data, num.clusters=NULL, num.clusters.per.omic=NULL) {
  start = Sys.time()
  omics.list = log.and.normalize(omics.list, subtype.data, normalize=T, filter.var=T)
  omics.list = lapply(omics.list, data.matrix)
  names(omics.list) = paste0('omic', 1:length(omics.list))
  mofa.obj = create_mofa(omics.list)
  data.options = get_default_data_options(mofa.obj)
  model.options = get_default_model_options(mofa.obj)
  train.options = get_default_training_options(mofa.obj)
  num.factors = min(15, sapply(omics.list, nrow))
  model.options$num_factors = num.factors
  train.options$seed = 42
  train.options$maxiter = 1e4
  mofa.obj = prepare_mofa(object=mofa.obj, data_options=data.options, 
                    model_options=model.options, training_options=train.options)
  
  mofa.ret = run_mofa(mofa.obj)
  factors = get_factors(mofa.ret, groups='all', factors='all')[[1]]
  if (is.null(num.clusters)) {
    sils = c()
    clustering.per.num.clusters = list()
    for (num.clusters in 2:15) {
      cur.clustering = kmeans(factors, num.clusters, iter.max=1e4, nstart=30)$cluster  
      sil = get.clustering.silhouette(list(t(factors)), cur.clustering)
      sils = c(sils, sil)
      clustering.per.num.clusters[[num.clusters - 1]] = cur.clustering
    }
    clustering = clustering.per.num.clusters[[which.max(sils)]]
  } else {
    clustering = kmeans(factors, num.clusters, iter.max=1e4, nstart=30)$cluster
  }

  time.taken = as.numeric(Sys.time() - start, units='secs')
  return(list(clustering=clustering, timing=time.taken))
}

###################### 
##### MONET CODE ##### 
###################### 

run.monet.with.params <- function(num.repeats=20, ...) {
  if (num.repeats == 1) {
    return(run.monet.with.params.single.run(...))
  } else {
    start = Sys.time()
    all.monet.rets = mclapply(1:num.repeats, function(i) {
      return(run.monet.with.params.single.run(seed=41 + i, ...))
    }, mc.cores=CONFIG$ncores)

    for (i in 1:length(all.monet.rets)) {
      if (class(all.monet.rets[[i]]) != 'list') {
        print(paste('error in run number', i))
      }
    }
    all.monet.rets = all.monet.rets[sapply(all.monet.rets, function(x) class(x) == 'list')]

    all.timings = sapply(all.monet.rets, function(monet.ret) monet.ret$timing)
    all.scores = sapply(all.monet.rets, function(monet.ret) monet.ret$monet.weight)
    best.ret.index = which.max(all.scores)
    monet.ret = all.monet.rets[[best.ret.index]]
    time.taken = as.numeric(Sys.time() - start, units='secs')
    monet.ret$timing = time.taken
    monet.ret$all.timings = all.timings
    monet.ret$all.scores = all.scores
    all.monet.rets = lapply(all.monet.rets, function(x) {x$omic.graphs = NULL;x$all.ems = NULL;return(x)})
    monet.ret$all.monet.rets = all.monet.rets
    return(monet.ret)
  }
}

run.monet.with.params.single.run <- function(omics.list, subtype.data, shift=0.1, frac.sampled=NULL, num.clusters=NULL, 
                                  percentile_shift=NULL, percentile_remove_edge=NULL, nemo.dist.without.na=F, seed=NULL, zero.low.weights=F) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  start = Sys.time()
  time.taken.cons = 0

  all.ems = NULL
  if (is.null(frac.sampled)) {
    em.graph.ret = get.em.graph.per.omic(omics.list)
    omics.dist = em.graph.ret[[1]]
    all.ems = em.graph.ret[[2]]
  } else {
    cons.ret = get.or.create(file.path(CONFIG$results.dir.path, sprintf('nemo_consensus_%f_%s_with_timing', frac.sampled, subtype.data$name)), function() {
      cons.start = Sys.time()
      omics.dist = lapply(omics.list, function(x) get.consensus.clustering(x, nemo.dist.without.na=nemo.dist.without.na, frac.sampled=frac.sampled, num.clusters=num.clusters))
      time.taken.cons = as.numeric(Sys.time() - cons.start, units='secs')
      return(list(omics.dist, time.taken.cons))
    })
    omics.dist = cons.ret[[1]]
    time.taken.cons = cons.ret[[2]]
    start = Sys.time()
  }

  if (!is.null(num.clusters)) {
    omics.dist = lapply(omics.dist, function(x) x - quantile(x, 1 - 1 / num.clusters))
  } else {
    omics.dist = lapply(omics.dist, function(x) x - shift)
  } 

  if (zero.low.weights) {
    for (i in 1:length(omics.dist)) {
      omics.dist[[i]][abs(omics.dist[[i]]) < 1e-5] = 0
    }
  }

  time.taken.normalization = as.numeric(Sys.time() - start, units='secs') + time.taken.cons

  is.any.na = any(sapply(omics.dist, function(x) any(is.na(x))))
  if (is.any.na) {
    pat.names = colnames(omics.dist[[1]])
    clustering = rep(1, length(pat.names))
    names(clustering) = pat.names
    return(list(clustering=clustering, timing=time.taken.normalization, mod.omics=list(), omic.graphs=list()))
  }
  
  tmp.dir = file.path(CONFIG$tmp.monet.dir, uuid::UUIDgenerate())
  dir.create(tmp.dir, showWarnings=F)

  # Empty files in the directory
  init.dir = file.path(tmp.dir, 'init')
  dir.create(init.dir)
  file.remove(file.path(init.dir, list.files(init.dir)))
    
  dist.dir = file.path(tmp.dir, 'dist')
  dir.create(dist.dir)
  file.remove(file.path(dist.dir, list.files(dist.dir)))
  for (i in 1:length(omics.list)) {
    omic = omics.dist[[i]]
    write.csv(omic, file=file.path(dist.dir, paste0('omic.', i, '.csv')))
  }

  file.remove(file.path(tmp.dir, list.files(tmp.dir)))
  for (i in 1:length(omics.list)) {
    omic = omics.list[[i]]
    write.csv(omic, file=file.path(tmp.dir, paste0('omic.', i, '.csv')))
  }
  
  all.sample.names = unique(unlist(lapply(omics.list, colnames)))
  num.samples = length(all.sample.names)
  prev.wd = getwd()
  setwd(CONFIG$monet.python.dir)
  command = sprintf('bash -c "source %s;python %s --dir_path %s --output_dir_path %s --num_seeds 15 --seed_size %d --input_type similarity --min_mod_size %d --max_pats_per_action %d', 
                    CONFIG$python.env.path, CONFIG$monet.main.path, tmp.dir, tmp.dir, 
		    floor(num.samples / 15), max(round(num.samples / 30), 10), ifelse(num.samples < 1000, 10, round(num.samples / 50)))
  
  if (!is.null(percentile_shift)) {
    command = paste0(command, ' --percentile_shift ', as.character(percentile_shift))
  }
  if (!is.null(percentile_remove_edge)) {
    command = paste0(command, ' --percentile_remove_edge ', as.character(percentile_remove_edge))
  }
  if (!is.null(seed)) {
    command = paste0(command, ' --rand_seed ', as.character(seed))
  }
  command = paste0(command, '"')
  
  print('running monet')
  print(command)
  command.return = system(command)
  stopifnot(command.return == 0)
  
  setwd(prev.wd)
  
  monet.secs = as.numeric(readLines(file.path(tmp.dir, 'time')))
  monet.weight = as.numeric(readLines(file.path(tmp.dir, 'weight')))
  monet.iter.times = read.table(file.path(tmp.dir, 'iter_times'))[,1]
  clustering = read.table(file=file.path(tmp.dir, 'pats'), stringsAsFactors=F)
  pat.names = clustering[,1]
  clustering[clustering[,2] == 'None', 2] = '-1'
  clustering = as.numeric(clustering[,2])
  num.outliers = sum(clustering == -1)
  clustering[clustering == -1] = (max(clustering) + 1):(max(clustering) + num.outliers)
  names(clustering) = pat.names
  clustering = clustering[all.sample.names]

  mod.file.content = readLines(file.path(tmp.dir, 'mods'))
  mod.omics = list()
  for (mod.data in strsplit(mod.file.content, ' ')) {
    cur.mod.name = mod.data[1]
    cur.mod.omics = strsplit(mod.data[2], ',')[[1]]
    mod.omics[[cur.mod.name]] = cur.mod.omics
  }
  
  omic.graphs = list()
  for (i in 1:length(omics.list)) {
    cur.omic.graph = as.matrix(read.csv(file.path(tmp.dir, paste0('graph_', i)))[, -1])
    rownames(cur.omic.graph) = all.sample.names
    colnames(cur.omic.graph) = rownames(cur.omic.graph)
    omic.graphs[[i]] = cur.omic.graph
  }

  return(list(clustering=clustering, timing=time.taken.normalization + monet.secs, mod.omics=mod.omics, omic.graphs=omic.graphs, all.ems=all.ems, monet.iter.times=monet.iter.times, monet.weight=monet.weight))
  
}


########################################
###### MONET weighting schemes #########
########################################

get.consensus.clustering <- function(mat, num.repeats = 100, frac.sampled = 0.5, num.clusters=NULL, nemo.dist.without.na=F, seed.offset=100) {
  total.nsamples = ncol(mat)
  nsamples = round(total.nsamples * frac.sampled)
  cons.list = mclapply(1:num.repeats, function(i) {
    set.seed(seed.offset + i)
    print(paste0('in consensus ', i))
    selected.samples = sample(1:total.nsamples, nsamples)
    cur.mat = mat[, selected.samples]
    attr(cur.mat, 'is.seq') = F
    cls = run.nemo(list(cur.mat), 'unused', num.clusters=num.clusters, dist.without.na=nemo.dist.without.na)$clustering

    consensus.mat = outer(cls, cls, function(x, y) x == y)
    prob.in.same.clust = mean(consensus.mat)
    consensus.mat = consensus.mat - prob.in.same.clust
    all.samples.mat = matrix(NA, nrow=total.nsamples, ncol=total.nsamples)
    all.samples.mat[selected.samples, selected.samples] = consensus.mat 
    return(list(samp_mat=all.samples.mat, cls=cls))
  }, mc.cores=CONFIG$ncores)
  consensus = Reduce("+", lapply(cons.list, function(x) {x$samp_mat[is.na(x$samp_mat)] = 0; return(x$samp_mat)})) / 
              Reduce("+", lapply(cons.list, function(x) !is.na(x$samp_mat)))
  colnames(consensus) = colnames(mat)
  rownames(consensus) = colnames(mat)
  return(consensus)
}

get.em.graph.per.omic <- function(raw.data, k=NA, alpha=0.5 ,dist.without.na=F, input.em.rets=NULL) {
  if (is.na(k)) {
    k = as.numeric(lapply(1:length(raw.data), function(i) round(ncol(raw.data[[i]]) / NUM.NEIGHBORS.RATIO)))
  } else if (length(k) == 1) {
    k = rep(k, length(raw.data))
  }
  if (dist.without.na) {
    dist.mats = lapply(1:length(raw.data), function(i) dist2.without.na(as.matrix(t(raw.data[[i]]))))
  } else {
    dist.mats = lapply(1:length(raw.data), function(i) dist2(as.matrix(t(raw.data[[i]])), as.matrix(t(raw.data[[i]]))))
  }
  sim.data = lapply(1:length(dist.mats), function(i) {1 - cor(raw.data[[i]])})

  all.omics.ret = lapply(1:length(sim.data), function(i) {
    cur.sim = sim.data[[i]]
    if (is.null(input.em.rets)) {
      num.clusters = estimateNumberOfClustersGivenGraph(as.matrix(cur.sim), 2:15)[[1]]
    } else {
      num.clusters = input.em.rets[[i]][[2]]
    }
    num.to.sample = nrow(cur.sim)
    chosen = sample(1:nrow(cur.sim), num.to.sample)
    chosen.sims.mat = cur.sim[chosen, chosen]
    chosen.sims = as.vector(chosen.sims.mat[upper.tri(chosen.sims.mat)])
    if (!is.null(input.em.rets)) {
      em.ret = input.em.rets[[i]][[1]]
    } else {
      all.em.rets = lapply(1:20, function(unused) normalmixEM(chosen.sims, k=2, maxit=1e5))
      em.ret = all.em.rets[[which.max(sapply(all.em.rets, function(x) x$loglik))]]
    }
    # calculate probabilities
    prob1 = log(em.ret$lambda[1]) + dnorm(cur.sim, mean=em.ret$mu[1], sd=em.ret$sigma[1], log=T)
    prob2 = log(em.ret$lambda[2]) + dnorm(cur.sim, mean=em.ret$mu[2], sd=em.ret$sigma[2], log=T)
    if (em.ret$mu[1] < em.ret$mu[2]) {
      prob = prob1 - prob2
      shift.by = quantile(prob[upper.tri(prob)], 1 - 1 / num.clusters)
      prob = prob - shift.by
      return(list(prob, em.ret, num.clusters))
    } else {
      prob = prob2 - prob1
      shift.by = quantile(prob[upper.tri(prob)], 1 - 1 / num.clusters)
      prob = prob - shift.by
      return(list(prob, em.ret, num.clusters))
    }
  })
  prob.ratio.data = lapply(all.omics.ret, function(x) x[[1]])
  all.omics.ems = lapply(all.omics.ret, function(x) list(x[[2]], x[[3]]))
  return(list(prob.ratio.data, all.omics.ems))
}

raw.data.to.graph <- function(raw.data, all.ems) {
  dist.mats = lapply(1:length(raw.data), function(i) dist2(as.matrix(t(raw.data[[i]])), as.matrix(t(raw.data[[i]]))))
  sim.data = lapply(1:length(dist.mats), function(i) {1 - cor(raw.data[[i]])})
  prob.ratio.data = lapply(1:length(sim.data), function(i) {
    cur.sim = sim.data[[i]]
    em.ret = all.ems[[i]]
    prob1 = log(em.ret$lambda[1]) + dnorm(cur.sim, mean=em.ret$mu[1], sd=em.ret$sigma[1], log=T)
    prob2 = log(em.ret$lambda[2]) + dnorm(cur.sim, mean=em.ret$mu[2], sd=em.ret$sigma[2], log=T)
    if (em.ret$mu[1] < em.ret$mu[2]) {
      return(list(log(0.9 + 0.1 * exp(prob2 - prob1)) - log(0.05 + 0.95 * exp(prob2 - prob1))))
    } else {
      return(list(log(0.9 + 0.1 * exp(prob1 - prob2)) - log(0.05 + 0.95 * exp(prob1 - prob2))))
    }
  })
  return(prob.ratio.data)
}

###############################
###### MONET analysis #########
###############################

remove.singletons <- function(clustering) {
  module.sizes = table(clustering)
  outlier.clusters = names(module.sizes)[module.sizes == 1]
  clustering = clustering[!(clustering %in% outlier.clusters)]
  pat.names = names(clustering)
  clustering = as.numeric(as.factor(clustering))
  names(clustering) = pat.names
  return(clustering)
}

group.singletons <- function(clustering) {
  module.sizes = table(clustering)
  outlier.clusters = names(module.sizes)[module.sizes == 1]
  clustering[clustering %in% outlier.clusters] = -1
  clustering = clustering + 1
  return(clustering)
}

get.solution.weight <- function(clustering, omic.graphs, mod.omics=NULL) {
  cluster.names = names(table(clustering))
  cls.weights = sapply(cluster.names, function(cls.name) {
    cls.samples = names(clustering)[clustering == cls.name]
    if (is.null(mod.omics)) {
      cur.mod.omics = 1:length(omic.graphs)
    } else {
      cur.mod.omics = as.numeric(mod.omics[[cls.name]])
      if (length(cur.mod.omics) == 0) { # singleton
        return(0)
      }
    }
    mod.weight.per.omic = sapply(cur.mod.omics, function(i) sum(omic.graphs[[i]][cls.samples, cls.samples]))
    cls.mod.weight = sum(mod.weight.per.omic)
  })
  return(cls.weights)
}

save.monet.results.to.file <- function(monet.ret, file.name) {
  con = file(file.name, 'w')
  mod.omics = monet.ret$mod.omics
  omics.strings = sapply(mod.omics, function(omics) paste(omics, collapse=','))
  mod.omics.string = paste(paste(names(mod.omics), omics.strings, sep='->'), collapse='\t')
  writeLines(mod.omics.string, con)
  write.table(monet.ret$clustering, file=con, col.names=F)
  close(con)
}

monet.ret.to.module.membership <- function(monet.ret, omic.graphs=NULL, omic.index=NULL) {
  if (is.null(omic.graphs)) {
    omic.graphs = monet.ret$omic.graphs
  }
  mod.omics = monet.ret$mod.omics
  mod.names = names(mod.omics)
  samples = names(monet.ret$clustering)
  all.module.membership = do.call(rbind, lapply(mod.names, function(mod.name) {
    cur.mod.omics = as.numeric(mod.omics[[mod.name]])
    cur.mod.samples = samples[monet.ret$clustering == mod.name]
    if (is.null(omic.index)) {
      cur.module.membership = Reduce(`+`, lapply(cur.mod.omics, function(i) colSums(omic.graphs[[i]][cur.mod.samples,])))
    } else {
      cur.module.membership = colSums(omic.graphs[[omic.index]][cur.mod.samples,])
    }
    return(cur.module.membership)
  }))
  rownames(all.module.membership) = mod.names
  return(all.module.membership)
}

monet.ret.to.super.graphs <- function(monet.ret) {
  mod.omics = monet.ret$mod.omics
  mod.names = names(mod.omics)
  samples = names(monet.ret$clustering)
  clustering = as.character(monet.ret$clustering)
  num.omics = length(monet.ret$omic.graphs)
  all.super.graphs = lapply(1:num.omics, function(num.omic) {
    omic.super.graph = outer(mod.names, mod.names, FUN = Vectorize(function(name1, name2) {
      mod1.samples = samples[clustering == name1]
      mod2.samples = samples[clustering == name2]
      return(sum(monet.ret$omic.graph[[num.omic]][mod1.samples, mod2.samples]))
    }))
    colnames(omic.super.graph) = mod.names
    rownames(omic.super.graph) = mod.names
    return(omic.super.graph)
  })
  return(all.super.graphs)
}

plot.mod.omics <- function(monet.ret, fig.path, omic.names, mod.omics=NULL, col.set=NULL) {
  if (is.null(mod.omics)) {
    mod.omics = monet.ret$mod.omics
  }
  if (is.null(omic.names)) {
    omic.names = sort(unique(unlist(mod.omics)))
  }
  mod.omics.binary = do.call(rbind, lapply(mod.omics, function(cur.mod.omics) omic.names %in% cur.mod.omics))
  # adding 0 is a trick to make it numeric
  mod.omics.binary.mat = mod.omics.binary * 1:nrow(mod.omics.binary)
  colnames(mod.omics.binary.mat) = omic.names

  if (is.null(col.set)) {
    if (nrow(mod.omics.binary) >= 10) {
      col.set = 'Set3'
    } else {
      col.set = 'Set1'
    }
  }
  cols = c('white', brewer.pal(nrow(mod.omics.binary), col.set))
  breaks = c(0, 0.5 + 0:(nrow(mod.omics.binary)), 1 + nrow(mod.omics.binary.mat))
  pheatmap(mod.omics.binary.mat, breaks=breaks, cluster_rows=F, cluster_cols=F, legend=F, color=cols, filename=fig.path)
}

# General function to perform several analyses on MONET's return value.
# Includes analysis of clinical parameters.
analyze.monet.ret <- function(monet.ret, name, is.survival=F, specific.cluster=NULL, omic.names=NULL, mut.mat=NULL, plot.module.membership=F) {
  fig.dir = file.path(CONFIG$plots.dir.path, sprintf('monet_%s', name))
  dir.create(fig.dir, showWarnings=F)
  clustering = monet.ret$clustering
  orig.clustering = clustering
  clustering = as.numeric(as.factor(clustering))
  names(clustering) = names(orig.clustering)
  clustering = remove.singletons(clustering)
  monet.cls = group.singletons(orig.clustering)
  num.omics = length(monet.ret$omic.graphs)

  module.membership = monet.ret.to.module.membership(monet.ret)
  module.membership.per.omic = lapply(1:num.omics, function(j) monet.ret.to.module.membership(monet.ret, omic.index=j))

  #pheatmap(module.membership, cluster_rows=F, cluster_cols=F)
  pats.ord = names(sort(monet.cls))
  for (j in 0:num.omics) {
    if (j == 0) {
      cur.mod.membership = module.membership
    } else {
      cur.mod.membership = module.membership.per.omic[[j]]
    }

    num.cols = length(table(monet.cls))
    col.pal = brewer.pal(num.cols - 1, 'Set1')
    cols = c('black', col.pal)
    tmp.cls = as.numeric(as.factor(monet.cls))
    names(tmp.cls) = names(monet.cls)
    col.df = as.data.frame(tmp.cls)
    colnames(col.df) = 'clustering'
    max.abs.val = max(abs(cur.mod.membership))
    if (plot.module.membership) {
      pheatmap(cur.mod.membership[, pats.ord], cluster_rows=F, cluster_cols=F, show_colnames=F, breaks=seq(-max.abs.val, max.abs.val, length.out=100),
               annotation_col=col.df, annotation_colors=list(clustering=cols), filename=file.path(fig.dir, sprintf('module_membership_%s.png', j)))
    }
    repr = Rtsne(X=t(module.membership), dims=2, check_duplicates=F)$Y
    png(file.path(fig.dir, sprintf('module_membership_tsne_%s.png', j)))
    plot(repr[,1], repr[,2], pch=19, col=cols[as.numeric(as.factor(monet.cls))], ylab='', xlab='', cex.lab=5, cex.main=5, cex=1.5)
    grid(col='black', lty='dashed')
    dev.off()
  }

  plot.mod.omics(monet.ret, fig.path=file.path(fig.dir, 'mod_omics.png'), omic.names)

  all.super.graphs = monet.ret.to.super.graphs(monet.ret)

  # survival analysis and clinical parameters
  if (is.survival) {
    subtype = name
    surv.data = get.surv.data(clustering, subtype)
    stopifnot(nrow(surv.data) == length(clustering))
    rownames(surv.data) = names(clustering)

    module.membership.bin = module.membership
    max.membership = apply(module.membership, 2, max)
    for (i in 1:ncol(module.membership.bin)) {
      module.membership.bin[,i] = module.membership[,i] == max.membership[i]
    }

    clinical.params = get.clinical.params(subtype)
    patho.stage = sapply(clinical.params[, 'clinical_stage'], patho.stage.to.num)
    names(patho.stage) = rownames(clinical.params)
    patho.stage2 = sapply(clinical.params[, 'clinical_stage'], patho.stage.to.num2)
    names(patho.stage2) = rownames(clinical.params)
    age = as.numeric(clinical.params[, 'age_at_initial_pathologic_diagnosis'])
    is.selected = !(is.na(patho.stage) | is.na(age))
    is.selected2 = !(is.na(patho.stage2) | is.na(age))
    sel.pats = rownames(clinical.params)[is.selected]
    sel.pats2 = rownames(clinical.params)[is.selected2]
    clin.params.filt = data.frame(age=age[is.selected], patho.stage=patho.stage[is.selected])
    clin.params.filt2 = data.frame(age=age[is.selected2], patho.stage=patho.stage2[is.selected2])
    rownames(clin.params.filt) = sel.pats
    rownames(clin.params.filt2) = sel.pats2

    all.cls.names = names(table(clustering))
    for (is.bin in c(F, T)) {
      print('-------- is bin -----------')
      print(is.bin)
      print('---------------------------')
      if (is.bin) {
        data.for.cox = cbind(surv.data, t(module.membership.bin)[rownames(surv.data),])[,c(-1, -4)]
      } else {
        data.for.cox = cbind(surv.data, t(module.membership)[rownames(surv.data),])[,c(-1, -4)]
      }
      for (cur.ind in 0:length(all.cls.names)) {
	cur.data.for.cox = data.for.cox
	if (cur.ind > 0) {
	  specific.cluster = all.cls.names[cur.ind]
          cur.data.for.cox = cur.data.for.cox[,c(1, 2, 2 + cur.ind)]
	  print('LOGRANK WITH PERMUTATION')
	  cls.bin = (clustering == cur.ind) + 1
	  names(cls.bin) = names(clustering)
          survival.file.path = file.path(CONFIG$results.dir.path, 'survival', paste(subtype, 'binary_model', cur.ind, sep='_'))
          logrank.ret = get.or.create(survival.file.path, function() get.cond.perm.surv(cls.bin, subtype))
	  print(logrank.ret)
	}  else {
	  if (is.bin) {
	    cur.data.for.cox = cur.data.for.cox[,1:3]
	    cur.data.for.cox[,3] = as.factor(clustering[rownames(cur.data.for.cox)])
	  }
	}
	print('CURRENT INDEX IS')
	print(cur.ind)
        cox.ret = coxph(Surv(Survival, Death) ~ ., data=cur.data.for.cox)
	print('COX WITHOUT CLINICAL')
        print(summary(cox.ret))
        common.pats = intersect(rownames(clin.params.filt), rownames(cur.data.for.cox))
        common.pats2 = intersect(rownames(clin.params.filt2), rownames(cur.data.for.cox))
        data.for.cox.clin = cbind(cur.data.for.cox[common.pats,], clin.params.filt[common.pats,])
        data.for.cox.clin2 = cbind(cur.data.for.cox[common.pats2,], clin.params.filt2[common.pats2,])
        cox.ret.clin = coxph(Surv(Survival, Death) ~ ., data=data.for.cox.clin)
	print('COX WITH CLINICAL 1')
        print(summary(cox.ret.clin))
        cox.ret.clin2 = coxph(Surv(Survival, Death) ~ ., data=data.for.cox.clin2)
	print('COX WITH CLINICAL 2')
        print(summary(cox.ret.clin2))

	if (cur.ind == 0 & is.bin) {
	  venous = clinical.params[, 'venous_invasion']
	  venous[venous == 'YES'] = 1
	  venous[venous == 'NO'] = 0
	  venous = as.numeric(venous)
	  names(venous) = rownames(clinical.params)
          has.venous = rownames(clinical.params)[!is.na(venous)]
	  common.pats.venous = intersect(has.venous, rownames(data.for.cox.clin))
	  cox.with.venous = cbind(data.for.cox.clin[common.pats.venous,], venous[common.pats.venous])
	  cox.with.venous.rel.cls = cox.with.venous[cox.with.venous[,3] %in% c(2, 4), ]
          cox.ret.venous = coxph(Surv(Survival, Death) ~ ., data=cox.with.venous.rel.cls)
	}
      }
    }


    clin.param.names = c('neoplasm_histologic_grade', 'clinical_stage', 'venous_invasion', 'anatomic_neoplasm_subdivision', 'postoperative_rx_tx')
    for (clin.param.name in clin.param.names) {
      print(paste('chi square on', clin.param.name))
      print(chisq.test(table(clustering, clinical.params[names(clustering), clin.param.name])))
    }
    print('tests for patho stage 1')
    print('chi square on clinical stage')
    print(chisq.test(table(clustering, patho.stage[names(clustering)])))
    print('kruskal wallis on clinical stage')
    print(kruskal.test(as.numeric(patho.stage[names(clustering)]), as.factor(clustering)))
    png(file.path(fig.dir, 'patho_stage_boxplot_1.png'))
    boxplot(split(as.numeric(patho.stage[names(clustering)]), as.factor(clustering)))
    dev.off()

    print('tests for patho stage 1')
    print('chi square on clinical stage')
    print(chisq.test(table(clustering, patho.stage2[names(clustering)])))
    print('kruskal wallis on clinical stage')
    print(kruskal.test(as.numeric(patho.stage2[names(clustering)]), as.factor(clustering)))
    png(file.path(fig.dir, 'patho_stage_boxplot_2.png'))
    boxplot(split(as.numeric(patho.stage2[names(clustering)]), as.factor(clustering)))
    dev.off()

    print('kruskal wallis on age')
    print(kruskal.test(as.numeric(clinical.params[names(clustering), ]$age_at_initial), as.factor(clustering)))
    png(file.path(fig.dir, 'age_at_diagnosis.png'))
    boxplot(split(as.numeric(clinical.params[names(clustering), ]$age_at_initial), as.factor(clustering)))
    dev.off()

    print('chi square on mutations')
    if (!is.null(mut.mat)) {
      rel.pats = intersect(names(clustering), rownames(mut.mat))
      mut.pvals = sapply(colnames(mut.mat), function(mut) { 
        cur.pval = chisq.test(table(mut.mat[rel.pats, mut], clustering[rel.pats]))$p.value
	return(cur.pval)
      })
      adj.mut.pvals = p.adjust(mut.pvals, 'BH')
      print(names(adj.mut.pvals)[adj.mut.pvals < 0.05])
      mut.md.df = data.frame(gene=names(mut.pvals), pval=mut.pvals, qval=adj.mut.pvals, is_known=rep(F, length(mut.pvals)))
      write.csv(mut.md.df, file=file.path(CONFIG$results.dir.path, 'ovarian_mut_df.csv'))
    }

    for (cur.cls in names(table(clustering))) {
      cur.cls.pats = names(clustering)[clustering == cur.cls]
      non.cls.pats = setdiff(names(clustering), cur.cls.pats)
      print(cur.cls)
      print('age')
      print(wilcox.test(clinical.params[cur.cls.pats, ]$age_at_initial, clinical.params[non.cls.pats, ]$age_at_initial))
      print('patho stage 1 wilcox')
      print(wilcox.test(patho.stage[cur.cls.pats], patho.stage[non.cls.pats]))
      print('patho stage 1 chisq')
      print(chisq.test(table(clustering == cur.cls, patho.stage[names(clustering)])))

      print('patho stage 2 wilcox')
      print(wilcox.test(patho.stage2[cur.cls.pats], patho.stage2[non.cls.pats]))
      print('patho stage 2 chisq')
      print(chisq.test(table(clustering == cur.cls, patho.stage2[names(clustering)])))

      print('venous invasion')
      print(chisq.test(table(clustering == cur.cls, clinical.params[names(clustering), 'venous_invasion'])))
    }
    venous.table = table(clustering, clinical.params[names(clustering), 'venous_invasion'])
    write.csv(venous.table, file=file.path(CONFIG$results.dir.path, 'venous_invasion_table.csv'))
  }
}

patho.stage.to.num <- function(patho.stage) {
  if (is.na(patho.stage)) {
    return(NA)
  }
  if (!startsWith(patho.stage, 'Stage ')) {
    return(NA)
  }
  patho.stage = substring(patho.stage, 7, nchar(patho.stage))
  if (patho.stage %in% c('IV', 'IVA', 'IVB')) {
    return(4)
  } else if (patho.stage %in% c('III', 'IIIA', 'IIIB', 'IIIC')) {
    return(3)
  } else if (patho.stage %in% c('II', 'IIA', 'IIB', 'IIC')) {
    return(2)
  } else if (patho.stage %in% c('I', 'IA', 'IB', 'IC')) {
    return(1)
  } else {
    return(NA)
  }
}

patho.stage.to.num2 <- function(patho.stage) {
  if (is.na(patho.stage)) {
    return(NA)
  }
  if (!startsWith(patho.stage, 'Stage ')) {
    return(NA)
  }
  patho.stage = substring(patho.stage, 7, nchar(patho.stage))
  stage.name.to.num = c(I=1, IA=2, IB=3, IC=4, II=5, IIA=6, IIB=7, IIC=8, III=9, IIIA=10, IIIB=11, IIIC=12, IV=13, IVA=14, IVB=15)
  return(stage.name.to.num[patho.stage])
}

############################################################
####### Omic dropping and classification experiments #######
############################################################

run.omic.dropping <- function(is.sim=T) {
  if (is.sim) {
    subtype = 'Digits'
    digits = get.or.create(file.path(CONFIG$ground.truth.datasets.path, 'digits_cached'), get.digits)
    subtype.raw.data = digits[[1]]
    labels = digits[[2]]
    names(labels) = colnames(subtype.raw.data[[1]])
    num.clusters = 10
    current.subtype.data = list(name='Digits', data=digits[[1]], labels=digits[[2]], num.clusters=10)
  } else {
    subtype = 'sarcoma'
    current.subtype.data = SUBTYPES.DATA[[which(sapply(SUBTYPES.DATA, function(x) x$name) == 'sarcoma')]]
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, 
                                           current.subtype.data)
    num.clusters = NULL
  }

  algorithm.name2 = 'monet.repeats'
  clustering.path2 = file.path(CONFIG$results.dir.path,
                            paste(subtype, algorithm.name2, sep='_'))
  algorithm.func.name2 = paste0('run.', algorithm.name2)
  algorithm.func2 = get(algorithm.func.name2)
  algorithm.ret2 = get.or.create(clustering.path2, function() error)
  all.sample.clustering2 = algorithm.ret2$clustering

  adj.to.all.samples = c()
  adj.to.ground.truth = c()

  num.repeats = 10
  adj.to.all.samples.left.out = matrix(NA, nrow=num.repeats, ncol=length(subtype.raw.data))
  adj.to.ground.truth.left.out = matrix(NA, nrow=num.repeats, ncol=length(subtype.raw.data))

  for (i in 1:num.repeats) {
    set.seed(142 + i)
    print(i)
    num.samples = ncol(subtype.raw.data[[1]])

    if (is.sim) {
      samples.per.omic = lapply(1:length(subtype.raw.data), function(j) {
        cur.omic.samples = colnames(subtype.raw.data[[j]])[sample(1:num.samples, 0.8 * num.samples)]
      })
    } else {
      removed.samples = colnames(subtype.raw.data[[1]])[sample(1:num.samples, 0.4 * num.samples)]
      removed.per.omic = split(removed.samples, 1:3)
      samples.per.omic = lapply(1:3, function(j) setdiff(colnames(subtype.raw.data[[j]]), removed.per.omic[[j]]))
    }

    stopifnot(length(table(unlist(samples.per.omic))) == num.samples)
    cur.iteration.data = lapply(1:length(subtype.raw.data), function(j) {
      subtype.raw.data[[j]][,samples.per.omic[[j]]]
    })

    if (!is.sim) {
      cur.iter.subtype.data = current.subtype.data
      # to avoid caching
      cur.iter.subtype.data$name = sprintf('%s_%d', subtype, i)
      cur.iteration.data = set.omics.list.attr(cur.iteration.data, cur.iter.subtype.data)
    } else {
      cur.iter.subtype.data = current.subtype.data
      cur.iter.subtype.data$name = sprintf('%s_%d', subtype, i)
      for (j in 1:length(cur.iteration.data)) {
        attr(cur.iteration.data[[j]], 'is.seq') = F
      }
    }
    
    cur.clustering.path = file.path(CONFIG$results.dir.path, paste(subtype, algorithm.name2, 'rep_vanilla_', i, sep='_'))
    algorithm.ret = get.or.create(cur.clustering.path, function() algorithm.func2(cur.iteration.data, cur.iter.subtype.data, num.clusters=num.clusters))
    new.clustering = algorithm.ret$clustering
    adj.to.all.samples = c(adj.to.all.samples, (adjustedRandIndex(new.clustering, all.sample.clustering2[names(new.clustering)])))
    if (is.sim) {
      adj.to.ground.truth = c(adj.to.ground.truth, (adjustedRandIndex(new.clustering, labels[names(new.clustering)])))
    }

    for (j in 1:length(samples.per.omic)) {
      sel.samples = samples.per.omic[[j]]
      non.sel.samples = setdiff(colnames(subtype.raw.data[[j]]), sel.samples)
      print(j)
      if (is.sim) {
	print('compared to label')
        print(adjustedRandIndex(new.clustering[sel.samples], labels[sel.samples]))
        print(adjustedRandIndex(new.clustering[non.sel.samples], labels[non.sel.samples]))
        adj.to.ground.truth.left.out[i, j] = adjustedRandIndex(new.clustering[non.sel.samples], labels[non.sel.samples])
      }
      print(adjustedRandIndex(new.clustering[sel.samples], all.sample.clustering2[sel.samples]))
      print(adjustedRandIndex(new.clustering[non.sel.samples], all.sample.clustering2[non.sel.samples]))
      adj.to.all.samples.left.out[i, j] = adjustedRandIndex(new.clustering[non.sel.samples], all.sample.clustering2[non.sel.samples])
    }
  }

  png(file.path(CONFIG$plots.dir.path, sprintf('%s_drop_omic_all_samples.png', subtype)))
  boxplot(cbind(adj.to.all.samples.left.out, (adj.to.all.samples)), ylim=c(0, 1), lwd=2)
  dev.off()
  png(file.path(CONFIG$plots.dir.path, sprintf('%s_drop_omic_ground_truth.png', subtype)))
  boxplot(cbind(adj.to.ground.truth.left.out, (adj.to.ground.truth)), ylim=c(0, 1), lwd=2)
  dev.off()
  png(file.path(CONFIG$plots.dir.path, sprintf('%s_drop_omic_both.png', subtype)))
  boxplot(cbind(adj.to.all.samples, (adj.to.ground.truth)), ylim=c(0, 1), lwd=2)
  dev.off()
  return(NULL)

}

run.classification <- function(is.sim=T, subtype='sarcoma') {
  if (is.sim) {
    subtype = 'Digits'
    digits = get.or.create(file.path(CONFIG$ground.truth.datasets.path, 'digits_cached'), get.digits)
    subtype.raw.data = digits[[1]]
    labels = digits[[2]]
    names(labels) = colnames(subtype.raw.data[[1]])
    num.clusters = 10
    labels = digits[[2]]
    names(labels) = colnames(subtype.raw.data[[1]])
    current.subtype.data = list(name='Digits', data=digits[[1]], labels=labels, num.clusters=10)
    for (j in 1:length(subtype.raw.data)) {
      attr(subtype.raw.data[[j]], 'is.seq') = F
    }
  } else {
    current.subtype.data = SUBTYPES.DATA[[which(sapply(SUBTYPES.DATA, function(x) x$name) == subtype)]]
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, 
                                           current.subtype.data)
    num.clusters = NULL
  }

  algorithm.name2 = 'monet.em.repeats'
  all.sample.clustering = get.or.create(file.path(CONFIG$results.dir.path, paste0(subtype, algorithm.name2)), 
                                        function() run.monet.em.repeats(subtype.raw.data, current.subtype.data, num.clusters=num.clusters))$clustering

  num.folds = 10
  all.samples = colnames(subtype.raw.data[[1]])
  set.seed(42)
  samples.per.fold = split(sample(all.samples, length(all.samples)), 1:num.folds)

  adj.to.all.samples.before = c()
  adj.to.ground.truth.before = c()
  adj.to.all.samples.after = c()
  adj.to.ground.truth.after = c()
  num.unclassified = c()

  for (i in 1:num.folds) {
    set.seed(42 + i)
    print(i)
    num.samples = ncol(subtype.raw.data[[1]])
    samples.to.remove = samples.per.fold[[i]]
    samples.to.keep = setdiff(all.samples, samples.to.remove)

    cur.iteration.data = lapply(1:length(subtype.raw.data), function(j) {
      subtype.raw.data[[j]][,samples.to.keep]
    })

    if (!is.sim) {
      cur.iter.subtype.data = current.subtype.data
      # to avoid caching
      cur.iter.subtype.data$name = sprintf('%s_%d', subtype, i)
      cur.iteration.data = set.omics.list.attr(cur.iteration.data, cur.iter.subtype.data)
    } else {
      cur.iter.subtype.data = current.subtype.data
      cur.iter.subtype.data$name = sprintf('%s_%d', subtype, i)
      for (j in 1:length(cur.iteration.data)) {
        attr(cur.iteration.data[[j]], 'is.seq') = F
      }
    }
    cur.iteration.data = log.and.normalize(cur.iteration.data, cur.iter.subtype.data, normalize=F, filter.var=T)
    norm.md = lapply(1:length(cur.iteration.data), function(j) {
      cur.df = cur.iteration.data[[j]]
      cur.means = rowMeans(cur.df, na.rm=T)
      temp = cur.df - cur.means
      cur.sd = apply(temp, 1, sd, na.rm=T)
      return(list(data=temp / cur.sd, means=cur.means, sd=cur.sd))
    })
    cur.iteration.data = lapply(norm.md, function(x) x$data)

    cur.clustering.path = file.path(CONFIG$results.dir.path, paste(subtype, algorithm.name2, 'rep_classification_', i, sep='_'))
    algorithm.ret = get.or.create(cur.clustering.path, function() run.monet.with.params(cur.iteration.data, cur.iter.subtype.data, num.clusters=num.clusters, shift=0, num.repeats=1))
    new.clustering = algorithm.ret$clustering

    # and now do the classification...
    data.for.classification = lapply(1:length(subtype.raw.data), function(j) subtype.raw.data[[j]][rownames(cur.iteration.data[[j]]),])
    if (!is.sim) {
      data.for.classification = set.omics.list.attr(data.for.classification, cur.iter.subtype.data)
    } else {
      for (j in 1:length(data.for.classification)) {
        attr(data.for.classification[[j]], 'is.seq') = F
      }
    }
    data.for.classification = log.and.normalize(data.for.classification, cur.iter.subtype.data, normalize=F, filter.var=T)

    norm.data = lapply(1:length(data.for.classification), function(j) {
      cur.df = data.for.classification[[j]]
      temp = cur.df - norm.md[[j]]$means
      cur.sd = norm.md[[j]]$sd
      return(temp / cur.sd)
    })

    all.graphs = get.em.graph.per.omic(norm.data, input.em.rets=algorithm.ret$all.ems)[[1]]
    if (is.sim) {
      all.graphs = lapply(all.graphs, function(graph) graph - quantile(graph[samples.to.keep, samples.to.keep], 0.9))
    } else {
      all.graphs = lapply(all.graphs, function(graph) graph)
    }
    for (j in 1:length(all.graphs)) {
      diag(all.graphs[[j]]) = 0
    }

    module.membership = monet.ret.to.module.membership(algorithm.ret, omic.graphs=all.graphs)[,samples.to.remove]
    classification.scores = apply(module.membership, 2, max)
    classification = rownames(module.membership)[apply(module.membership, 2, which.max)]
    num.lonely = sum(classification.scores < 0)
    num.unclassified = c(num.unclassified, num.lonely)
    max.cls = max(as.numeric(rownames(module.membership)))
    classification[classification.scores < 0] = (max.cls + 1):(max.cls + num.lonely)
    names(classification) = samples.to.remove
    new.clustering = c(algorithm.ret$clustering, classification)
    adj.to.all.samples.after = c(adj.to.all.samples.after, adjustedRandIndex(new.clustering, all.sample.clustering[names(new.clustering)]))
    adj.to.all.samples.before = c(adj.to.all.samples.before, adjustedRandIndex(algorithm.ret$clustering, all.sample.clustering[names(algorithm.ret$clustering)]))

    if (is.sim) {
      adj.to.ground.truth.after = c(adj.to.ground.truth.after, adjustedRandIndex(new.clustering, labels[names(new.clustering)]))
      adj.to.ground.truth.before = c(adj.to.ground.truth.before, adjustedRandIndex(algorithm.ret$clustering, labels[names(algorithm.ret$clustering)]))
      print(adjustedRandIndex(new.clustering, labels[names(new.clustering)]))
      print(adjustedRandIndex(algorithm.ret$clustering, labels[names(algorithm.ret$clustering)]))
      print(table(new.clustering, labels[names(new.clustering)]))
    }
  }
  fig.name = file.path(CONFIG$plots.dir.path, sprintf('%s_classification.png', subtype))
  print(c('total num unclassified', sum(num.unclassified)))
  png(fig.name)
  if (is.sim) {
    boxplot(cbind(adj.to.all.samples.before, adj.to.all.samples.after,
                  adj.to.ground.truth.before, adj.to.ground.truth.after), ylim=c(0, 1), lwd=2)
  } else {
    boxplot(cbind(adj.to.all.samples.before, adj.to.all.samples.after), ylim=c(0, 1), lwd=2)
  }
  dev.off()
  return(NULL)
}


######################################################
###############Data Processing Functions##############
######################################################

log.and.normalize <- function(omics.data, subtype.data, normalize=T,
                              filter.var=F, num.feats=2000) {
  # filter features with no variance at all
  for (i in 1:length(omics.data)) {
    orig.attr = attr(omics.data[[i]], 'is.seq')
    omics.data[[i]] = omics.data[[i]][apply(omics.data[[i]], 1, function(x) var(x, na.rm=T)) > 0,]
    attr(omics.data[[i]], 'is.seq') = orig.attr
  }
			      
  for (i in 1:length(omics.data)) {
    if (attr(omics.data[[i]], 'is.seq')) {
      omics.data[[i]] = log(1+omics.data[[i]])
    }
  }
  
  if (filter.var) {
    omics.data = lapply(omics.data, function(omic) keep.high.var.features(omic, num.feats))
  }
  
  if (normalize) {
    omics.data = lapply(omics.data, normalize.matrix)    
  }
  
  return(omics.data)
}

filter.features.with.few.different.values <- function(omics.data, threshold=10) {
  filtered.omics = lapply(omics.data, function(omic) {
    features.to.keep = unlist(apply(omic, 1, function(x) length(table(as.numeric(x))) >= threshold))
    return(omic[features.to.keep, ])
  })
  return(filtered.omics)
}

normalize.matrix <- function(data.matrix) {
  temp = data.matrix - rowMeans(data.matrix, na.rm=T)
  should.keep = (apply(temp, 1, sd, na.rm=T) != 0)
  return ((temp / apply(temp, 1, sd, na.rm=T))[should.keep, ])
}

keep.high.var.features <- function(omic, num.features=2000) {
  if (nrow(omic) < num.features) {
    return(omic)
  } else {
    feature.vars = apply(omic, 1, var, na.rm=T)
    threshold = feature.vars[order(feature.vars, decreasing = T)][num.features]
    return(omic[feature.vars >= threshold,])    
  }
}


filter.non.tumor.samples <- function(raw.datum, only.primary=only.primary) {
  # 01 is primary, 06 is metastatic, 03 is blood derived cancer
  if (!only.primary)
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01', '03', '06')])
  else
    return(raw.datum[,substring(colnames(raw.datum), 14, 15) %in% c('01')])
}

get.fixed.names <- function(patient.names, include.type=F) {
  # fix the TCGA names to only include the patient ids
  if (include.type) {
    return(gsub('-', '\\.', toupper(substring(patient.names, 1, 15))))
  } else {
    return(gsub('-', '\\.', toupper(substring(patient.names, 1, 12))))  
  }
}

fix.patient.names <- function(subtype.raw.data, include.type=F) {
  for (i in 1:length(subtype.raw.data)) {
    colnames(subtype.raw.data[[i]]) = get.fixed.names(colnames(subtype.raw.data[[i]]),
                                                      include.type)
  }
  return(subtype.raw.data)
}

set.omics.list.attr <- function(subtype.raw.data, subtype.data) {
  attr(subtype.raw.data[[1]], 'is.seq') = subtype.data$is.rna.seq
  attr(subtype.raw.data[[2]], 'is.seq') = F
  attr(subtype.raw.data[[3]], 'is.seq') = subtype.data$is.mirna.seq
  return(subtype.raw.data)
}

get.raw.data <- function(subtype.name,
                         datasets.path = get.dataset.dir.path(),
                         only.primary=NA,
			 intersect.patients=T) {
  omics.dir = file.path(datasets.path, subtype.name)
  omics.files = list.files(omics.dir)
  omics.files = setdiff(omics.files, c('survival'))  
  raw.data = lapply(file.path(omics.dir, omics.files), read.table)
  
  if (!is.na(only.primary)) {
    raw.data = lapply(raw.data, function(x) filter.non.tumor.samples(x, only.primary = only.primary))
  }
  if (intersect.patients) {
    name.corrected.data = fix.patient.names(raw.data)
    patients.intersection = Reduce(intersect, lapply(name.corrected.data, colnames))
    ret.data = lapply(name.corrected.data, function(datum) datum[,patients.intersection])  
  } else {
    ret.data = raw.data
  }
  
  return(ret.data)
}

##############################################################
#### Code for gene clustering and single cell experiments ####
##############################################################

get.scnmt <- function() {
  scnmt.data = readRDS(CONFIG$scnmt.data.path)
  scnmt.metadata = read.table(CONFIG$scnmt.metadata.path)
  scrna = scnmt.data@ExperimentList[['rna']]
  attr(scrna, 'is.seq') = F
  scmet = scnmt.data@ExperimentList[['met_promoter']]
  attr(scmet, 'is.seq') = F
  all.data = list(scrna, scmet)
  return(list(all.data, scnmt.data))
}

get.gene.microarray.data <- function() {
  data.cache.path = file.path(CONFIG$results.dir.path, 'gene_microarray_data')
  return(get.or.create(data.cache.path, function() {
    microarray.file.path = file.path(CONFIG$microarray.dir.path, 'microarray')
    microarray.data = read.table(microarray.file.path, sep='\t', stringsAsFactors=F)
    rownames(microarray.data) = microarray.data[,1]
    colnames(microarray.data) = microarray.data[1,]
    microarray.data = microarray.data[-c(1, 2, nrow(microarray.data)), -1]
    microarray.data = filter.non.tumor.samples(microarray.data, only.primary=T)
    colnames(microarray.data) = get.fixed.names(colnames(microarray.data))

    microarray.data[microarray.data == 'null'] = NA
    microarray.data = microarray.data[rowSums(is.na(microarray.data)) == 0,]

    exp.data = get.raw.data('breast', only.primary=T, intersect.patients=F)[[1]]
    exp.data = fix.patient.names(list(exp.data))[[1]]
    common.samples = intersect(colnames(microarray.data), colnames(exp.data))
    exp.data = exp.data[,common.samples]
    microarray.data = data.matrix(microarray.data[,common.samples])

    # now intersect the microarrays and exp
    microarray.gene.names = rownames(microarray.data)
    exp.gene.names = sapply(strsplit(rownames(exp.data), '\\|'), function(x) x[1])

    genes.to.remove = names(which(table(exp.gene.names) > 1))
    genes.to.use = setdiff(intersect(exp.gene.names, microarray.gene.names), genes.to.remove)
    exp.data.fil = do.call(rbind, lapply(genes.to.use, function(gene.name) {
      exp.data[exp.gene.names == gene.name, , drop=F]
    }))
    prot.data.fil = do.call(rbind, lapply(genes.to.use, function(gene.name) {
      microarray.data[microarray.gene.names == gene.name, , drop=F]
    }))
    rownames(prot.data.fil) = genes.to.use
    rownames(exp.data.fil) = genes.to.use
    omics.list = list(exp.data.fil, prot.data.fil)
    return(omics.list)
  }))
}

ds.scnmt.data <- function(scnmt.omics) {
  rna = scnmt.omics[[1]]
  methy = scnmt.omics[[2]]
  methy = methy[,colnames(rna)]

  num.nas = colSums(is.na(methy))
  na.thresh = quantile(num.nas, 0.75)
  should.remove = num.nas > na.thresh
  methy.fil = methy[,!should.remove]
  meth.ds = do.call(cbind, lapply(1:ncol(methy.fil), function(i) {
    cur.meth.vector = methy.fil[, i]
    cur.is.na = is.na(cur.meth.vector)
    num.to.na = na.thresh - sum(cur.is.na)
    indices.to.na = sample(which(!cur.is.na), num.to.na)
    cur.meth.vector[indices.to.na] = NA
    return(cur.meth.vector)
  }))
  colnames(meth.ds) = colnames(methy.fil)
  rna.fil = rna[,!should.remove]
  return(list(rna.fil, meth.ds))
}

run.scnmt.analysis <- function() {
  scnmt = get.scnmt()
  scnmt.omics = scnmt[[1]]
  scnmt.data = scnmt[[2]]@colData

  scnmt.omics = ds.scnmt.data(scnmt.omics)

  attr(scnmt.omics[[1]], 'is.seq') = T
  attr(scnmt.omics[[2]], 'is.seq') = F
  num.clusters = NULL

  clustering.path = file.path(CONFIG$results.dir.path, 'monet_scnmt_results')

  scnmt.cell.subtype.data = list(name='scnmt_cell_subtype')
  cell.algorithm.ret = get.or.create(clustering.path, function() run.monet.repeats(scnmt.omics, scnmt.cell.subtype.data, nemo.dist.without.na=T, num.clusters=num.clusters))
  save.monet.results.to.file(cell.algorithm.ret, file.path(CONFIG$results.dir.path, paste('scnmt', 'monet', 'formatted_output', sep='_')))

  nemo.clustering.path = file.path(CONFIG$results.dir.path, 'nemo_scnmt_results')
  nemo.cell.cls = get.or.create(nemo.clustering.path, function() run.nemo(scnmt.omics, 'unused', dist.without.na=T, num.clusters=num.clusters)$clustering)
  
  cell.cls.fil = remove.singletons(cell.algorithm.ret$clustering)
  cell.cls = group.singletons(cell.algorithm.ret$clustering)
  
  cell.tsne.plot.path = file.path(CONFIG$plots.dir.path, 'tsne_scnmt_cells.png')
  cell.tsne.solution.path = file.path(CONFIG$results.dir.path, 'tsne_scnmt_cells_solution')
  plot.clustering.tsne(scnmt.omics, list(monet=cell.cls, nemo=nemo.cell.cls), cell.tsne.plot.path, cell.tsne.solution.path, use.dist=c(F, T), 'Set3')

  methy.tsne.repr = t(get.or.create(paste0(cell.tsne.solution.path, '_2')))
  colnames(methy.tsne.repr) = names(cell.cls)
  attr(methy.tsne.repr, 'is.seq') = F
  methy.nemo.cls = run.nemo(list(methy.tsne.repr), list(name='methy_tsne'), num.clusters=2)$clustering
  
  stage.tbl = table(scnmt.data[names(cell.cls.fil), 'stage'], cell.cls.fil)
  stage.lineage.tbl = table(scnmt.data[names(cell.cls.fil), 'stage_lineage'], cell.cls.fil)
  lineage.10x.tbl = table(scnmt.data[names(cell.cls.fil), 'lineage10x_2'], cell.cls.fil)

  write.csv(stage.tbl, file=file.path(CONFIG$results.dir.path, 'stage_tbl.csv'))
  write.csv(stage.lineage.tbl, file=file.path(CONFIG$results.dir.path, 'stage_lineage_tbl.csv'))
  write.csv(lineage.10x.tbl, file=file.path(CONFIG$results.dir.path, 'lineage_10x_tbl.csv'))

  cls.table = table(cell.algorithm.ret$clustering)
  cls.names = names(cls.table)[cls.table > 1]
  cls.no.outliers = cell.algorithm.ret$clustering[cell.algorithm.ret$clustering %in% cls.names]
  stage.lineage.tbl2 = table(scnmt.data[names(cls.no.outliers), 'stage_lineage'], cls.no.outliers)

  print(stage.tbl)
  print(lineage.10x.tbl)
  print(stage.lineage.tbl)
  print(stage.lineage.tbl2)

  plot.mod.omics(cell.algorithm.ret, file.path(CONFIG$plots.dir.path, 'scnmt_mod_omics.png'), 1:2, col.set='Set3')
}

# run with '_many_genes_no_remove'
run.gene.microarray.analysis <- function(exp.name='') {
  omic.name = paste0('microarray', exp.name)
  omics.list = get.gene.microarray.data()

  attr(omics.list[[1]], 'is.seq') = T
  attr(omics.list[[2]], 'is.seq') = F
  omics.list.orig = omics.list

  omics.list = lapply(log.and.normalize(omics.list, subtype.data, filter.var=T, num.feats=2000), t)
  common.genes = intersect(colnames(omics.list[[1]]), colnames(omics.list[[2]]))
  omics.list[[1]] = omics.list[[1]][,common.genes]
  omics.list[[2]] = omics.list[[2]][,common.genes]
  attr(omics.list[[1]], 'is.seq') = F
  attr(omics.list[[2]], 'is.seq') = F
  
  clustering.path = file.path(CONFIG$results.dir.path, paste0('monet_genes_', omic.name))
  subtype.data = list(name=paste0('gene_subtype_', omic.name))
  algorithm.ret = get.or.create(clustering.path, function() run.monet.repeats(omics.list, subtype.data, num.clusters=num.clusters))
  print(algorithm.ret$mod.omics)

  cls.fil = remove.singletons(algorithm.ret$clustering)
  
  monet.cls = group.singletons(algorithm.ret$clustering)

  cls.table = table(algorithm.ret$clustering)
  cls.names = names(cls.table)[cls.table > 1]
  write.table(colnames(omics.list[[1]]), quote=F, row.names=F, col.names=F, file=file.path(CONFIG$results.dir.path, sprintf('gene_module_bg')))
  write.table(cls.fil, file=file.path(CONFIG$results.dir.path, sprintf('microarray_analysis_gene_modules')))
  for (cur.cls.name in cls.names) {
    cur.genes = names(algorithm.ret$clustering)[algorithm.ret$clustering == cur.cls.name]
    write.table(cur.genes, quote=F, row.names=F, col.names=F, file=file.path(CONFIG$results.dir.path, sprintf('gene_module_%s', cur.cls.name)))
  }
  
  tsne.plot.path = file.path(CONFIG$plots.dir.path, sprintf('tsne_%s_tmp.png', omic.name))
  tsne.solution.path = file.path(CONFIG$results.dir.path, sprintf('tsne_%s_solution_tmp6', omic.name))
  plot.clustering.tsne(omics.list, list(monet=monet.cls), tsne.plot.path, tsne.solution.path)

  omics.norm = log.and.normalize(omics.list, subtype.data, filter.var=F, normalize=T)
  genes.ord = names(sort(algorithm.ret$clustering))

  for (j in 1:length(omics.norm)) {
    cur.omic = omics.norm[[j]]
    cur.omic.cor = cor(cur.omic)
    diag(cur.omic.cor) = NA
    shades = colorRampPalette(c('darkblue', 'blue', 'white', 'red', 'darkred'))(100)

    num.cols = length(table(monet.cls))
    col.pal = brewer.pal(num.cols - 1, 'Set1')
    cols = c('black', col.pal)
    tmp.cls = as.numeric(as.factor(monet.cls))
    names(tmp.cls) = names(monet.cls)
    col.df = as.data.frame(tmp.cls)
    colnames(col.df) = 'clustering'

    png(file.path(CONFIG$plots.dir.path, sprintf('gene_omic_%s.png', j)))
    pheatmap(cur.omic.cor[genes.ord, genes.ord], col=shades, cluster_rows=F, cluster_cols=F, show_colnames=F, show_rownames=F, breaks=seq(-1, 1, length.out=100), annotation_col=col.df, annotation_colors=list(clustering=cols), annotation_names_col=F, annotation_names_row=F, annotation_legend=F)
    dev.off()
  }

  return(algorithm.ret)
}

###################################
##### Code for TCGA benchmark #####
###################################

run.benchmark <- function(subtype.index=NULL) {
  all.timings = matrix(NA, nrow=length(SUBTYPES.DATA), ncol=length(ALGORITHM.NAMES))
  rownames(all.timings) = sapply(SUBTYPES.DATA, function(x) x$name)
  colnames(all.timings) = ALGORITHM.NAMES
  all.logrank.pvalues = all.timings
  all.clinical.params = all.timings
  all.num.clin.params = all.timings
  all.num.clusters = all.timings
  alg.cols = brewer.pal(length(ALGORITHMS.DATA), 'Set3')
  all.clusterings = list()
  monet.num.outliers = c()
  monet.mod.omics = c()
  all.heavy.monet.rets = list()

  if (is.null(subtype.index)) { 
    subtype.indices = 1:length(SUBTYPES.DATA)
  } else {
    subtype.indices = subtype.index
  }

  for (i in subtype.indices) {
    current.subtype.data = SUBTYPES.DATA[[i]]
    subtype = current.subtype.data$name
    #subtype.raw.data = NULL
    subtype.raw.data = get.raw.data(subtype, 
                                    only.primary=current.subtype.data$only.primary)
    
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, 
                                           current.subtype.data)
    
    all.subtype.clusterings = list()
    for (algorithm.name in ALGORITHM.NAMES) {
      set.seed(42)
      print(paste('Running: dataset', subtype, 'algorithm', algorithm.name))
      clustering.path = file.path(CONFIG$results.dir.path,
                                paste(subtype, algorithm.name, sep='_'))

      algorithm.func.name = paste0('run.', algorithm.name)
      algorithm.func = get(algorithm.func.name)
      cur.iteration.data = subtype.raw.data
      algorithm.ret = get.or.create(clustering.path, function() algorithm.func(cur.iteration.data, current.subtype.data))
      pat.names = colnames(cur.iteration.data[[1]])
      #pat.names = names(algorithm.ret$clustering)
      if (algorithm.name == 'monet.repeats50') {
        all.heavy.monet.rets[[subtype]] = algorithm.ret
      }

      clustering = algorithm.ret$clustering
      clustering = as.numeric(as.factor(clustering))
      names(clustering) = pat.names

      timing = algorithm.ret$timing
      all.timings[subtype, algorithm.name] = timing

      survival.dir = file.path(CONFIG$results.dir.path, 'survival')
      dir.create(survival.dir)

      clinical.dir = file.path(CONFIG$results.dir.path, 'clinical')
      dir.create(clinical.dir)

      orig.clustering = clustering

      if (ALGORITHMS.DATA[[algorithm.name]]$is_monet) {
        analyze.monet.ret(algorithm.ret, subtype, omic.names=as.character(1:length(algorithm.ret$omic.graphs)))
	monet.mod.omics = c(monet.mod.omics, sapply(algorithm.ret$mod.omics, function(x) paste(sort(x), collapse='_')))
	monet.num.outliers = c(monet.num.outliers, sum(table(clustering) == 1))
        print(subtype)
	print(algorithm.ret$mod.omics)

        clustering = remove.singletons(clustering)
        survival.file.path = file.path(survival.dir, paste(subtype, algorithm.name, 'no_outliers', sep='_'))
        clinical.file.path = file.path(clinical.dir, paste(subtype, algorithm.name, 'no_outliers', sep='_'))
      } else {
        survival.file.path = file.path(survival.dir, paste(subtype, algorithm.name, sep='_'))
        clinical.file.path = file.path(clinical.dir, paste(subtype, algorithm.name, sep='_'))
      }

      survival.ret = get.or.create(survival.file.path, function() get.cond.perm.surv(clustering, subtype))
      all.logrank.pvalues[subtype, algorithm.name] = survival.ret$pvalue

      if (length(table(clustering)) == 1) {
        all.clinical.params[subtype, algorithm.name] = 0
      } else {
        clinical.ret = get.or.create(clinical.file.path, function() check.clinical.enrichment(clustering, subtype))
        all.clinical.params[subtype, algorithm.name] = sum(clinical.ret * length(clinical.ret) < 0.05)
      }

      all.num.clusters[subtype, algorithm.name] = length(table(clustering))
      all.subtype.clusterings[[algorithm.name]] = orig.clustering

      if (ALGORITHMS.DATA[[algorithm.name]]$is_monet) {
        save.monet.results.to.file(algorithm.ret, file.path(CONFIG$results.dir.path, paste(subtype, algorithm.name, 'formatted_output', sep='_')))
      }

      if (ALGORITHMS.DATA[[algorithm.name]]$is_monet) {

        module.membership = monet.ret.to.module.membership(algorithm.ret)
	cls.names = names(table(clustering))

        monet.cls = group.singletons(orig.clustering)
        monet.cls = monet.cls[pat.names]
        all.subtype.clusterings[[algorithm.name]] = monet.cls
      }
    }



    tsne.solution.path = file.path(CONFIG$results.dir.path, sprintf('tsne_%s', subtype))
    tsne.selected.methods = c('snf', 'mdi', 'monet.em.repeats')
    tsne.plot.path2 = file.path(CONFIG$plots.dir.path, sprintf('tsne_%s_selected_methods.png', subtype))
    plot.clustering.tsne(cur.iteration.data, all.subtype.clusterings[tsne.selected.methods], tsne.plot.path2, tsne.solution.path)

    all.clusterings[[subtype]] = all.subtype.clusterings
  }

  plot.all.pairwise.ari(all.heavy.monet.rets, file.path(CONFIG$plots.dir.path, 'heavy_monet_adj.png'))
  write.csv(all.timings, file=file.path(CONFIG$results.dir.path, 'tcga_timings.csv'))
  write.csv(all.logrank.pvalues, file=file.path(CONFIG$results.dir.path, 'tcga_logrank_pvalues.csv'))
  write.csv(all.num.clusters, file=file.path(CONFIG$results.dir.path, 'tcga_num_clusters.csv'))
  write.csv(all.clinical.params, file=file.path(CONFIG$results.dir.path, 'tcga_num_clinical_params.csv'))

  timing.df = t(data.frame(colMeans(log10(1 + all.timings))))
  colnames(timing.df) = colnames(all.timings)
  plot.bars(timing.df, alg.cols, file.path(CONFIG$plots.dir.path, 'tcga_timings.png'), F)

  survival.df = t(data.frame(colSums(all.logrank.pvalues < 0.05)))
  colnames(survival.df) = colnames(all.timings)
  plot.bars(survival.df, alg.cols, file.path(CONFIG$plots.dir.path, 'tcga_survival.png'), F)

  clinical.df = t(data.frame(colMeans(all.clinical.params > 0)))
  colnames(clinical.df) = colnames(all.clinical.params)
  plot.bars(clinical.df, alg.cols, file.path(CONFIG$plots.dir.path, 'tcga_clinical.png'), F)

  num.cls.df = t(data.frame(colMeans(all.num.clusters)))
  colnames(num.cls.df) = colnames(all.timings)
  plot.bars(num.cls.df, alg.cols, file.path(CONFIG$plots.dir.path, 'tcga_num_cls.png'), F)

  names(monet.mod.omics) = rownames(all.timings)
  plot.bars(table(monet.mod.omics), rep('grey', length(table(monet.mod.omics))), file.path(CONFIG$plots.dir.path, 'tcga_monet_mod_omics.png'), F)

  omic1.mods = grepl('1', monet.mod.omics)
  omic2.mods = grepl('2', monet.mod.omics)
  omic3.mods = grepl('3', monet.mod.omics)

  png(filename=file.path(CONFIG$plots.dir.path, 'tcga_monet_mod_omics_venn.png'))
  draw.triple.venn(sum(omic1.mods), sum(omic2.mods), sum(omic3.mods),
                   sum(omic1.mods & omic2.mods), sum(omic1.mods & omic3.mods), sum(omic2.mods & omic3.mods),
		   sum(omic1.mods & omic2.mods & omic3.mods))
  dev.off()

  num.outliers.df = t(data.frame(monet.num.outliers))
  colnames(num.outliers.df) = rownames(all.timings)
  png(file.path(CONFIG$plots.dir.path, 'tcga_monet_outliers.png'))
  barplot(num.outliers.df, col=rep('grey', length(monet.num.outliers)), las=2)
  dev.off()

}

#####################################
###### Ovarian cancer analysis ######
#####################################

# call with: analyze.cancer.in.depth(), analyze.cancer.in.depth(NULL, 'gbm')
analyze.cancer.in.depth <- function(subtype.raw.data=NULL, subtype='ovarian') {
  return.early = subtype != 'ovarian'
  monet.ret = get.or.create(file.path(CONFIG$results.dir.path, paste(subtype, 'monet.repeats', sep='_')), error)
  current.subtype.data = SUBTYPES.DATA[[which(sapply(SUBTYPES.DATA, function(x) x$name) == subtype)]]

  orig.clustering = monet.ret$clustering
  clustering = remove.singletons(monet.ret$clustering)
  tmp = names(clustering)
  clustering = as.numeric(as.factor(clustering))
  names(clustering) = tmp
  plot.kaplan.meier(clustering, subtype, file.path(CONFIG$plots.dir.path, sprintf('kaplan_meier_%s.png', subtype)), cols=brewer.pal(length(table(clustering)), 'Set1'))

  survival.dir = file.path(CONFIG$results.dir.path, 'survival')
  subtype.logrank.pvalues = sapply(ALGORITHM.NAMES, function(alg.name) {
    survival.file.path = file.path(survival.dir, paste(subtype, alg.name, sep='_'))
    if (ALGORITHMS.DATA[[alg.name]]$is_monet) {
      survival.file.path = paste0(survival.file.path, '_no_outliers')
    }
    pval = get.or.create(survival.file.path, error)$pvalue
    return(pval)
  })

  names(subtype.logrank.pvalues) = sapply(ALGORITHMS.DATA, function(x) x$display.name)
  png(file.path(CONFIG$plots.dir.path, sprintf('%s_logrank_pvals.png', subtype)))
  barplot(subtype.logrank.pvalues, col=brewer.pal(length(ALGORITHM.NAMES), 'Set3'))
  abline(h=0.05, lwd=5)
  dev.off()

  if (is.null(subtype.raw.data)) {
    subtype.raw.data = get.raw.data(subtype, only.primary=current.subtype.data$only.primary)
    subtype.raw.data = set.omics.list.attr(subtype.raw.data, current.subtype.data)
  }

  # cluster using only one omic
  #rna.alg.names = c('snf', 'nemo', 'monet.repeats')
  rna.alg.names = c('snf', 'nemo')
  for (omic.index in 1:length(subtype.raw.data)) {
    omic.name = c('rna', 'methy', 'mirna')[omic.index]
    num.clust.opts = list(list(NULL, 3, 4), list(NULL), list(NULL))[[omic.index]]
    #num.clust.opts = list(list(4), list(4), list(4))[[omic.index]]
    for (num.clust.index in 1:length(num.clust.opts)) {
      cur.num.clust = num.clust.opts[[num.clust.index]]
      print(paste('cur num clust is', cur.num.clust))
      for (rna.alg.name in rna.alg.names) {
        print(rna.alg.name)
        algorithm.func.name = paste0('run.', rna.alg.name)
        algorithm.func = get(algorithm.func.name)
        clustering.path = file.path(CONFIG$results.dir.path, paste(subtype, rna.alg.name, omic.name, cur.num.clust, sep='_'))
        rna.algorithm.ret = get.or.create(clustering.path, function() algorithm.func(subtype.raw.data[omic.index], current.subtype.data, num.clusters=cur.num.clust))

        rna.clustering = rna.algorithm.ret$clustering
        pat.names = names(rna.clustering)
        rna.clustering = as.numeric(as.factor(rna.clustering))
        names(rna.clustering) = pat.names
        survival.dir = file.path(CONFIG$results.dir.path, 'survival')
        survival.file.path = file.path(survival.dir, paste(subtype, rna.alg.name, omic.name, cur.num.clust, sep='_'))
        if (ALGORITHMS.DATA[[rna.alg.name]]$is_monet) {
          rna.clustering = remove.singletons(rna.clustering)
        }
        survival.ret = get.or.create(survival.file.path, function() get.cond.perm.surv(rna.clustering, subtype))
        print(survival.ret$pvalue)
      }
    }
  }

  if (return.early) {
    return(NULL)
  }

  omics.list = log.and.normalize(subtype.raw.data, SUBTYPES.DATA[[9]], filter.var=T, normalize=F)
  cls.names = names(table(clustering))
  pats.ord = c()
  all.cls.vals = list()
  diff.feats.ord = list(list(), list(), list())
  for (i in 1:(length(cls.names) + 2)) {
    if (i <= length(cls.names)) {
      cur.cls.pats = names(clustering)[clustering == cls.names[i]]
      pats.ord = c(pats.ord, cur.cls.pats)
      non.cur.cls.pats = setdiff(names(clustering), cur.cls.pats)
    } else if (i == length(cls.names) + 1) {
      cur.cls.pats = names(clustering)[clustering == cls.names[2]]
      non.cur.cls.pats = names(clustering)[clustering == cls.names[4]]
    } else {
      cur.cls.pats = names(clustering)[clustering == cls.names[4]]
      non.cur.cls.pats = names(clustering)[clustering == cls.names[2]]
    }
    # find differentially expressed features in each omic
    for (j in 1:length(omics.list)) {
      cls.vals = omics.list[[j]][,cur.cls.pats]
      all.cls.vals[[j]] = cls.vals
      non.cls.vals = omics.list[[j]][,non.cur.cls.pats]
      all.pvals = unlist(mclapply(1:nrow(omics.list[[j]]), function(k) {
	tryCatch({
          #return(t.test(cls.vals[k,], non.cls.vals[k,], alternative='greater')$p.value)
          return(t.test(cls.vals[k,], non.cls.vals[k,], alternative='greater')$p.value)
	}, error = function(e) {
	  return(1)
	})
      }, mc.cores=CONFIG$ncores))
      #adj.pvals = p.adjust(all.pvals, 'bonferroni')
      adj.pvals = p.adjust(all.pvals, 'BH')
      names(all.pvals) = rownames(omics.list[[j]])
      names(adj.pvals) = rownames(omics.list[[j]])
      is.high.exp = pmax(rowMeans(cls.vals), rowMeans(non.cls.vals)) > c(2, -Inf, 0.5)[j]
      fc = log(rowMeans(cls.vals) / rowMeans(non.cls.vals))
      is.high.fc = abs(fc) > c(0.5, 0.2, 0.3)[j]

      diff.feats = rownames(omics.list[[j]])[adj.pvals < 0.05 & is.high.fc & is.high.exp]
      print(c(i, j, length(diff.feats)))
      if (j == 1) {
	diff.feats.names = sapply(strsplit(diff.feats, '\\.'), function(x) x[1])
	bg.gene.names = sapply(strsplit(rownames(omics.list[[j]]), '\\.'), function(x) x[1])
        write.table(diff.feats.names, quote=F, row.names=F, col.names=F, file=file.path(CONFIG$results.dir.path, sprintf('diff_feats_%s_%s', i, j)))
        write.table(bg.gene.names, quote=F, row.names=F, col.names=F, file=file.path(CONFIG$results.dir.path, sprintf('bg_feats_%s', j)))
	print(sort(fc[diff.feats], dec=T)[1:50])
	print('------------')
	print(sort(adj.pvals[diff.feats])[1:50])
      }
      if (j == 3 & i == 1) {
	png(file.path(CONFIG$plots.dir.path, 'mirna_feats.png'))
	plot(rowMeans(cls.vals), rowMeans(non.cls.vals), pch=19)
	points(rowMeans(cls.vals)[diff.feats], rowMeans(non.cls.vals)[diff.feats], pch=19, col='red')
	abline(b=1, a=0, col='red')
	grid(col='grey', lwd=0.5)
	dev.off()
        print(diff.feats)

	diff.feats.df = data.frame(name=diff.feats, pval=all.pvals[diff.feats], qval=adj.pvals[diff.feats])
	write.csv(diff.feats.df, file=file.path(CONFIG$results.dir.path, 'diff_mirna.csv'))

	pats.for.mirna.plot = c(cur.cls.pats, non.cur.cls.pats)
	mir1.name = 'hsa.mir.514.1'
	mir2.name = 'hsa.mir.514.3'
	mir1 = unlist(omics.list[[j]][mir1.name, pats.for.mirna.plot])
	mir2 = unlist(omics.list[[j]][mir2.name, pats.for.mirna.plot])
	png(file.path(CONFIG$plots.dir.path, sprintf('%s_%s_comparison.png', mir1.name, mir2.name)))
	plot(mir1[non.cur.cls.pats], mir2[non.cur.cls.pats], pch=19)
	points(mir1[cur.cls.pats], mir2[cur.cls.pats], pch=19, col='red')
	dev.off()

	for (cur.mir in diff.feats) {
	  cur.mir.cls.dens = density(unlist(cls.vals[cur.mir,]))
	  cur.mir.non.cls.dens = density(unlist(non.cls.vals[cur.mir,]))
	  y.max = max(cur.mir.cls.dens[['y']], cur.mir.non.cls.dens[['y']])
	  mir.fig.dir = file.path(CONFIG$plots.dir.path, 'mir_figs')
	  dir.create(mir.fig.dir, showWarnings=F)
	  png(file.path(mir.fig.dir, sprintf('%s_density.png', cur.mir)))
	  plot(cur.mir.cls.dens, col='red', ylim=c(0, y.max), lwd=5)
	  lines(cur.mir.non.cls.dens, ylim=c(0, y.max), lwd=5)
	  dev.off()
	}

      }
      cur.diff.feats.sorted = names(sort(fc[diff.feats], dec=T))
      if (length(cur.diff.feats.sorted) > 50) {
        cur.diff.feats.sorted = cur.diff.feats.sorted[1:50]
      }
      diff.feats.ord[[j]][[i]] = cur.diff.feats.sorted
    }
  }

  clust.for.heatmap = clustering[pats.ord]
  cols = brewer.pal(length(table(clust.for.heatmap)), 'Set1')
  tmp.cls = as.numeric(as.factor(clust.for.heatmap))
  names(tmp.cls) = names(clust.for.heatmap)
  col.df = as.data.frame(tmp.cls)
  colnames(col.df) = 'clustering'
  for (j in 1:length(omics.list)) {
    cur.feats = unlist(diff.feats.ord[[j]])
    cur.feats = cur.feats[!duplicated(cur.feats)]
    png(file.path(CONFIG$plots.dir.path, sprintf('ovarian_image_%s.png', j)))
    norm.omic = normalize.matrix(omics.list[[j]])
    cur.feats = c(cur.feats, setdiff(rownames(omics.list[[j]]), cur.feats))
    norm.omic.lim = pmax(pmin(norm.omic, 3), -3)
    shades = colorRampPalette(c('darkblue', 'blue', 'white', 'red', 'darkred'))(100)
    pheatmap(as.matrix(norm.omic.lim[cur.feats, pats.ord]), breaks=seq(-3, 3, length.out=101), col=shades, 
             cluster_rows=F, cluster_cols=F, show_colnames=F, show_rownames=F, annotation_names_col=F, annotation_names_row=F,
             annotation_col=col.df, annotation_colors=list(clustering=cols))
    dev.off()
  }

  for (j in 1:length(omics.list)) {
    nemo.ret = run.nemo(subtype.raw.data[j], current.subtype.data)
    print(c('old logrank', j))
    print(old.logrank(nemo.ret$clustering, subtype))

    snf.ret = run.spectral(subtype.raw.data[j], current.subtype.data)
    print(c('old logrank', j))
    print(old.logrank(snf.ret$clustering, subtype))
  }

  mut.mat = get.or.create(file.path(CONFIG$results.dir.path, 'ovarian_mut_mat'), mutation.dir.to.mat)
  analyze.monet.ret(monet.ret, subtype, is.survival=T, omic.names=as.character(1:3), mut.mat=mut.mat)
}

mutation.dir.to.mat <- function(mutation.dir=CONFIG$mutation.dir.path) {
  pat.file.names = setdiff(list.files(mutation.dir), 'MANIFEST.txt')
  pat.names = gsub('-', '\\.', toupper(substring(pat.file.names, 1, 12)))
  muts.per.sample = lapply(pat.file.names, function(fname) {
    file.cont = read.table(file.path(mutation.dir, fname), stringsAsFactors=F, sep='\t')
    gene.names = file.cont[-1, 1]
    return(gene.names)
  })
  known.muts = c('TP53', 'BRCA1', 'CSMD3', 'NF1', 'CDK12', 'FAT3', 'GABRA6', 'BRCA2', 'RB1')
  mut.table = table(unlist(muts.per.sample))
  freq.muts = c()
  rel.muts = unique(c(known.muts, freq.muts))
  muts.mat = do.call(rbind, lapply(muts.per.sample, function(x) rel.muts %in% x))
  colnames(muts.mat) = rel.muts
  rownames(muts.mat) = pat.names
  return(muts.mat)
}



###################################
####### Digits dataset code #######
###################################

get.digits.full <- function() {
  digits.dir = CONFIG$digits.path
  all.file.paths = file.path(digits.dir, list.files(digits.dir))
  all.data = lapply(all.file.paths, function(fpath) {
    data.mat = read.table(fpath)
    rownames(data.mat) = paste0('img', 1:nrow(data.mat))
    colnames(data.mat) = paste0('feat', 1:ncol(data.mat))
    ret = as.data.frame(t(data.mat))
    return(ret)
  })
  labs = unlist(lapply(1:10, function(x) rep(x, 200)))
  return(list(all.data, labs))
}

get.digits1000 = function() {get.digits(1000)}
get.digits1500 = function() {get.digits(1500)}

get.digits <- function(num.digits=400) {
  set.seed(42)
  digits.dir = CONFIG$digits.path
  all.file.paths = file.path(digits.dir, list.files(digits.dir))
  all.data = lapply(all.file.paths, function(fpath) {
    data.mat = read.table(fpath)
    rownames(data.mat) = paste0('img', 1:nrow(data.mat))
    colnames(data.mat) = paste0('feat', 1:ncol(data.mat))
    ret = as.data.frame(t(data.mat))
    return(ret)
  })
  labs = unlist(lapply(1:10, function(x) rep(x, 200)))
  indices = sample(1:2000, num.digits)
  all.data = lapply(all.data, function(x) x[,indices])
  labs = labs[indices]
  return(list(all.data, labs))
}

#########################################################
####### Simulations and image dataset experiments #######
#########################################################

run.ground.truth.datasets <- function() {
  digits = get.or.create(file.path(CONFIG$ground.truth.datasets.path, 'digits_cached'), get.digits)
  digits.full = get.or.create(file.path(CONFIG$ground.truth.datasets.path, 'digits_full_cached'), get.digits.full)
  set.seed(42)
  sim1 = get.or.create(file.path(CONFIG$ground.truth.datasets.path, 'sim1_cached'), create.simulation1)
  set.seed(43)
  sim2 = get.or.create(file.path(CONFIG$ground.truth.datasets.path, 'sim2_cached'), create.simulation2)

  ground.truth.md = list(list(name='Sim1', data=sim1[[1]], labels=sim1[[2]], num.clusters=NULL, is.sim=T),
                         list(name='Sim21', data=sim2[[1]][1:2], labels=sim2[[2]], num.clusters=NULL, is.sim=T),
                         list(name='Sim22', data=sim2[[1]], labels=sim2[[2]], num.clusters=NULL, is.sim=T),
                         list(name='Digits', data=digits[[1]], labels=digits[[2]], num.clusters=10, is.sim=F),
                         list(name='Digits_full', data=digits.full[[1]], labels=digits.full[[2]], num.clusters=10, is.sim=F, is.heavy.dataset=T))


  all.adj.rands = matrix(NA, nrow=length(ground.truth.md), ncol=length(ALGORITHM.NAMES) + 1)
  rownames(all.adj.rands) = sapply(ground.truth.md, function(x) x$name)
  monet.no.sing = 'Monet (without lonely)'
  colnames(all.adj.rands) = c(sapply(ALGORITHMS.DATA, function(x) x$display.name), monet.no.sing)
  alg.cols = brewer.pal(length(ALGORITHMS.DATA) + 1, 'Set3')
  all.timings = all.adj.rands[,1:length(ALGORITHMS.DATA)]
  iteration.times = list()
  monet.algorithm.names = ALGORITHM.NAMES[sapply(ALGORITHMS.DATA, function(x) x$is_monet)]

  rownames(all.adj.rands) = sapply(ground.truth.md, function(x) x$name)
  dir.create(CONFIG$plots.dir.path, showWarnings=F)
  for (i in 1:length(ground.truth.md)) {
    current.md = ground.truth.md[[i]]
    name = current.md$name
    labels = current.md$labels
    cur.data = current.md$data
    is.sim = current.md$is.sim
    is.heavy.dataset = current.md$is.heavy.dataset
    if (is.null(is.heavy.dataset)) {
      is.heavy.dataset = F
    }
    cur.num.clusters = current.md$num.clusters
    for (j in 1:length(cur.data)) {
      attr(cur.data[[j]], 'is.seq') = F
    }

    if (is.sim) {
      data.plot.path = file.path(CONFIG$plots.dir.path, sprintf('raw_%s.png', name))
      plot.omics.data(cur.data, paste('Omic', 1:length(cur.data)), labels, data.plot.path)

      data.plot.with.legend.path = file.path(CONFIG$plots.dir.path, sprintf('raw_%s_with_legend.png', name))
      plot.omics.data(cur.data, paste('Omic', 1:length(cur.data)), labels, data.plot.with.legend.path, T)
    }
    
    all.clusterings = list()
    if (is.sim | is.heavy.dataset) {
     cur.algorithm.names = monet.algorithm.names
     slow.algorithm.names = ALGORITHM.NAMES[sapply(ALGORITHMS.DATA, function(x) {
       if (is.null(x$is_slow)) {
         return(F)
       } else {
         return(x$is_slow)
       }
     })]
     if (is.heavy.dataset) {
       cur.algorithm.names = setdiff(cur.algorithm.names, slow.algorithm.names)
     }
    } else {
     cur.algorithm.names = ALGORITHM.NAMES
    }
    for (algorithm.name in cur.algorithm.names) {
      set.seed(42)
      print(paste('Running: dataset', name, 'algorithm', algorithm.name))
      is.monet = ALGORITHMS.DATA[[algorithm.name]]$is_monet
      display.name = ALGORITHMS.DATA[[algorithm.name]]$display.name
      algorithm.func.name = paste0('run.', algorithm.name)
      algorithm.func = get(algorithm.func.name)
      clustering.path = file.path(CONFIG$results.dir.path,
                                paste(name, algorithm.name, sep='_'))

      if (is.heavy.dataset) {
        clustering.ret = get.or.create(clustering.path, function() algorithm.func(cur.data, current.md, num.clusters=cur.num.clusters, zero.low.weights=T))
      } else {
        clustering.ret = get.or.create(clustering.path, function() algorithm.func(cur.data, current.md, num.clusters=cur.num.clusters))
      }
      all.adj.rands[name, display.name] = adjustedRandIndex(clustering.ret$clustering, labels)
      all.timings[name, display.name] = clustering.ret$timing
      all.clusterings[[display.name]] = clustering.ret$clustering
      if (is.monet & !is.sim) {
	names(labels) = names(clustering.ret$clustering)
	cls.no.sing = remove.singletons(clustering.ret$clustering)
        all.adj.rands[name, monet.no.sing] = adjustedRandIndex(cls.no.sing, labels[names(cls.no.sing)])
      }
      if (name == 'Digits' & algorithm.name == 'monet.repeats50') {
        plot.all.pairwise.ari(list(clustering.ret), file.path(CONFIG$plots.dir.path, 'heavy_monet_digits_adj.png'))
        plot.all.pairwise.ari(list(clustering.ret), file.path(CONFIG$plots.dir.path, 'heavy_monet_digits_labels_adj.png'), labels)
      }

      clustering = clustering.ret$clustering
      clustering = as.numeric(as.factor(clustering))
      names(clustering) = colnames(cur.data[[1]])
      print(table(clustering))
      
      if (algorithm.name == 'monet.repeats') {
        iteration.times[[name]] = lapply(clustering.ret$all.monet.rets, function(x) x$monet.iter.times)
      }

      if (is.monet) {
        analyze.monet.ret(clustering.ret, name, omic.names = as.character(1:length(cur.data)))
	# note that some figures were created with singletons and some without, so next few lines might alternate
        monet.cls = group.singletons(clustering.ret$clustering)
        monet.cls = monet.cls[colnames(cur.data[[1]])]

        monet.cls = clustering.ret$clustering[colnames(cur.data[[1]])] + 1
	monet.cls = as.numeric(factor(monet.cls, levels=c('0', monet.cls[!duplicated(monet.cls)])))
	cls.table = table(monet.cls)
	singleton.clusters = names(cls.table)[cls.table == 1]
	monet.cls[monet.cls %in% singleton.clusters] = 0
	names(monet.cls) = names(clustering.ret$clustering[colnames(cur.data[[1]])])

	orig.cls = clustering.ret$clustering
	orig.cls.table = table(orig.cls)
	orig.cls.names = names(orig.cls.table)[orig.cls.table > 1]
	orig.to.new.names = sapply(orig.cls.names, function(orig.name) monet.cls[names(orig.cls)[orig.cls == orig.name]][1])
	names(orig.to.new.names) = orig.cls.names
	new.mod.omics = clustering.ret$mod.omics
	names(new.mod.omics) = orig.to.new.names[names(new.mod.omics)]
	new.mod.omics = new.mod.omics[order(names(new.mod.omics))]

        plot.mod.omics(NULL, file.path(CONFIG$plots.dir.path, sprintf('%s_mod_omics.png', name)), 1:length(cur.data), mod.omics=new.mod.omics)


        all.clusterings[[display.name]] = monet.cls
	print(clustering.ret$mod.omics)

	names(labels) = colnames(cur.data[[1]])
        labels.cls = group.singletons(labels)
        labels.cls = labels.cls[colnames(cur.data[[1]])]
        all.clusterings[['Labels']] = labels.cls
      }
    }

    tsne.plot.path = file.path(CONFIG$plots.dir.path, sprintf('%s.png', name))
    tsne.solution.path = file.path(CONFIG$results.dir.path, sprintf('tsne_%s', name))
    plot.clustering.tsne(cur.data, all.clusterings, tsne.plot.path, tsne.solution.path)
  }

  num.iters = sapply(iteration.times[['Digits']], length)
  mean.time.iters = sapply(iteration.times[['Digits']], mean)
  num.iters.full = sapply(iteration.times[['Digits_full']], length)
  mean.time.iters.full = sapply(iteration.times[['Digits_full']], mean)
  png(file.path(CONFIG$plots.dir.path, 'digits_num_iter.png'))
  boxplot(list(num.iters, num.iters.full))
  dev.off()
  png(file.path(CONFIG$plots.dir.path, 'digits_mean_time_iter.png'))
  boxplot(list(log10(1 + mean.time.iters), log10(1 + mean.time.iters.full)))
  dev.off()
  
  plot.bars(all.adj.rands['Digits',,drop=F], alg.cols=alg.cols, file.path(CONFIG$plots.dir.path, 'ground_truth_adj_rand.png'))
  plot.bars(all.adj.rands['Digits',,drop=F], alg.cols=alg.cols, file.path(CONFIG$plots.dir.path, 'ground_truth_adj_rand_no_legend.png'), F)
  plot.bars(log10(1 + all.timings['Digits', ,drop=F]), alg.cols=alg.cols, file.path(CONFIG$plots.dir.path, 'ground_truth_timings.png'))
  plot.bars(log10(1 + all.timings['Digits', ,drop=F]), alg.cols=alg.cols, file.path(CONFIG$plots.dir.path, 'ground_truth_timings_no_legend.png'), F)
  
}

plot.omics.data <- function(omics.list, omic.names, labels, file.name, is.image.plot=F) {
  png(file.name)
  samp.ord = order(labels)
  max.abs.val = max(sapply(omics.list, function(x) max(abs(x))))
  shades = colorRampPalette(c('darkblue', 'blue', 'white', 'red', 'darkred'))(100)
  par(mfrow=c(length(omics.list), 1), mar=c(2, 2, 2, 2))
  for (i in 1:length(omics.list)) {
    if (is.image.plot) {
      image.plot(as.matrix(omics.list[[i]][, samp.ord]), col=shades, breaks=seq(-max.abs.val, max.abs.val, length.out=length(shades) + 1))
    } else {
      image(as.matrix(omics.list[[i]][, samp.ord]), col=shades, breaks=seq(-max.abs.val, max.abs.val, length.out=length(shades) + 1), axes=F, lab.breaks=F)
      axis(1, at=seq(0, 1, length=6), labels=seq(0, nrow(omics.list[[i]]), length=6))
      axis(2, at=seq(1, 0, length=6), labels=seq(ncol(omics.list[[i]]), 0, length=6))
    }

  }
  dev.off()
}

plot.bars <- function(all.timings, alg.cols, file.name, show.legend=T) {
  png(file.name)
  if (show.legend) {
    barplot(t(all.timings), col=alg.cols, legend=colnames(all.timings), beside=T)
  } else {
    barplot(t(all.timings), col=alg.cols, legend=F, beside=T)
  }
  dev.off()
}


########################################
######### Survival analysis ############
########################################

plot.kaplan.meier <- function(clustering, subtype, file.name, cols=1:length(table(clustering))) {
  surv.data = get.surv.data(clustering, subtype)
  rownames(surv.data) = surv.data[,1]
  surv.data = surv.data[,-1]
  surv.ret = survfit(Surv(Survival, Death) ~ cluster, data=surv.data)
  png(file.name)
  plot(surv.ret, col=cols, lwd=8, legend=T)
  grid(col='black')
  dev.off()
}

old.logrank <- function(clustering, subtype) {
  surv.data = get.surv.data(clustering, subtype)
  rownames(surv.data) = surv.data[,1]
  surv.data = surv.data[,-1]
  surv.ret = survdiff(Surv(Survival, Death) ~ cluster, data=surv.data)
}

get.cond.perm.surv <- function(clustering, subtype) {
  if (length(table(clustering)) == 1) {
    return(list(pvalue=1))
  }
  surv.data = get.surv.data(clustering, subtype)
  rownames(surv.data) = surv.data[,1]
  surv.data = surv.data[,-1]
  get.cond.perm.surv.given.data(surv.data)
}

get.subtype.survival.path <- function(subtype) {
  datasets.path = get.dataset.dir.path()
  survival.file.path = file.path(datasets.path, subtype, 'survival')
  return(survival.file.path)
}


get.surv.data <- function(groups, subtype, survival.file.path) {
  if (missing(survival.file.path)) {
    survival.file.path = get.subtype.survival.path(subtype)
  }
  survival.data = read.table(survival.file.path, header = TRUE)
  patient.names = names(groups)
  patient.names.in.file = as.character(survival.data[, 1])
  patient.names.in.file = toupper(gsub('-', '\\.', substring(patient.names.in.file, 1, 12)))

  stopifnot(all(patient.names %in% patient.names.in.file))
  
  indices = match(patient.names, patient.names.in.file)
  ordered.survival.data = survival.data[indices,]
  ordered.survival.data["cluster"] <- groups
  ordered.survival.data$Survival[is.na(ordered.survival.data$Survival)] = 0
  ordered.survival.data$Death[is.na(ordered.survival.data$Death)] = 0
  return(ordered.survival.data)
  
}

########################################
##### Clinical enrichment code #########
########################################

get.clinical.params.dir <- function() {
  return(CONFIG$clinical.params.dir)
}

get.clinical.params <- function(subtype.name) {
  clinical.data.path = paste(get.clinical.params.dir(), subtype.name, sep = '')
  clinical.params = read.table(clinical.data.path,
                               sep='\t', header=T, row.names = 1, stringsAsFactors = F)
  rownames.with.duplicates = get.fixed.names(rownames(clinical.params))  
  clinical.params = clinical.params[!duplicated(rownames.with.duplicates),]
  rownames(clinical.params) = rownames.with.duplicates[!duplicated(rownames.with.duplicates)]
  return(clinical.params)
}

check.clinical.enrichment <- function(clustering, subtype.name) {
  clinical.params = get.clinical.params(subtype.name)  
  
  clinical.metadata = list(gender='DISCRETE', age_at_initial_pathologic_diagnosis='NUMERIC', 
    leukemia_french_american_british_morphology_code='DISCRETE', PAM50Call_RNAseq='DISCRETE',
    pathologic_M='DISCRETE', pathologic_N='DISCRETE', pathologic_T='DISCRETE', pathologic_stage='DISCRETE')
  
  pvalues = c()
  
  params.being.tested = c()
  
  for (clinical.param in names(clinical.metadata)) {
    
    if (!(clinical.param %in% colnames(clinical.params))) {
      next
    }
    
    clinical.values = clinical.params[names(clustering),clinical.param]
    is.discrete.param = clinical.metadata[clinical.param] == 'DISCRETE'
    is.numeric.param = clinical.metadata[clinical.param] == 'NUMERIC'
    stopifnot(is.discrete.param | is.numeric.param)
    
    # skip parameter if many missing values
    
    if (is.numeric.param) {
      numeric.entries = !is.na(as.numeric(clinical.values))
      if (2 * sum(numeric.entries) < length(clinical.values)) {
        next
      }
    } else {
      not.na.entries = !is.na(clinical.values)
      should.skip = F
      if (2 * sum(not.na.entries) < length(clinical.values)) {
        should.skip = T
      } else if (length(table(clinical.values[not.na.entries])) == 1) {
        should.skip = T
      }
      if (should.skip) {
        next
      }
    }
    
    params.being.tested = c(params.being.tested, clinical.param)
    
    if (is.discrete.param) {
      #clustering.with.clinical = cbind(clustering, clinical.values)
      #tbl = table(as.data.frame(clustering.with.clinical[!is.na(clinical.values),]))
      #test.res = chisq.test(tbl)
      #pvalue = test.res$p.value
      pvalue = get.empirical.clinical(clustering[!is.na(clinical.values)], clinical.values[!is.na(clinical.values)], T)
      
    } else if (is.numeric.param) {
      #test.res = kruskal.test(as.numeric(clinical.values[numeric.entries]),
      #				clustering[numeric.entries])
      #pvalue = test.res$p.value
      pvalue = get.empirical.clinical(clustering[numeric.entries], as.numeric(clinical.values[numeric.entries]), F)
    }
    
    pvalues = c(pvalues, pvalue)
    
  }
  names(pvalues) = params.being.tested
  return(pvalues)
}

get.empirical.clinical <- function(clustering, clinical.values, is.chisq) {
  set.seed(42)
  if (is.chisq) {
      clustering.with.clinical = cbind(clustering, clinical.values)
      tbl = table(as.data.frame(clustering.with.clinical))
      test.res = chisq.test(tbl)
  } else {
    test.res = kruskal.test(as.numeric(clinical.values), clustering)
  }
  orig.pvalue = test.res$p.value
  num.iter = 1000
  total.num.iters = 0
  total.num.extreme = 0
  should.continue = T
  
  while (should.continue) {
    print('another iteration in empirical clinical')
    perm.pvalues = as.numeric(mclapply(1:num.iter, function(i) {
      cur.clustering = sample(clustering)
      names(cur.clustering) = names(clustering)
    
      if (is.chisq) {
        clustering.with.clinical = cbind(cur.clustering, clinical.values)
        tbl = table(as.data.frame(clustering.with.clinical))
        test.res = chisq.test(tbl)
      } else {
        test.res = kruskal.test(as.numeric(clinical.values), cur.clustering)
      }
      cur.pvalue = test.res$p.value
      return(cur.pvalue)
    }, mc.cores=CONFIG$ncores)) 
    total.num.iters = total.num.iters + num.iter
    total.num.extreme = total.num.extreme + sum(perm.pvalues <= orig.pvalue)
    
    binom.ret = binom.test(total.num.extreme, total.num.iters)
    cur.pvalue = binom.ret$estimate
    cur.conf.int = binom.ret$conf.int
    
    sig.threshold = 0.05
    is.threshold.in.conf = cur.conf.int[1] < sig.threshold & cur.conf.int[2] > sig.threshold
    if (!is.threshold.in.conf | total.num.iters > 1e5) {
      should.continue = F
    }
  }
  return(cur.pvalue)
}

#################################
###### Utils and config #########
#################################

plot.clustering.tsne <- function(omics.list, all.clusterings, plot.file.name, tsne.file.name, use.dist=NULL, col.set=NULL) {
  if (is.null(use.dist)) {
    use.dist = rep(F, length(omics.list))
  }
  omics.list = log.and.normalize(omics.list, 'UNUSED', normalize=T, filter.var=T)
  if (!is.null(plot.file.name)) {
    png(plot.file.name, height=480 * length(all.clusterings), width=480 * length(omics.list))
    par(mfrow=c(length(all.clusterings), length(omics.list)), mar=c(2, 7, 4, 2))
  } else {
    par(mfrow=c(length(all.clusterings), length(omics.list)), mar=c(2, 2, 2, 2))
  }
  for (i in 1:length(all.clusterings)) {
    clustering = all.clusterings[[i]]
    num.cols = length(table(clustering))

    if (is.null(col.set)) {
      if (num.cols >= 10) {
        col.set = 'Set3'
      } else {
        col.set = 'Set1'
      }
    }
    if (min(clustering) == 0) {
      cols = c('black', brewer.pal(num.cols - 1, col.set))
    } else {
      cols = c(brewer.pal(num.cols, col.set))
    }
    for (j in 1:length(omics.list)) {
      omic = omics.list[[j]]
      if (use.dist[j]) {
	omic.dist = 1 - cor(omic, use='pairwise.complete.obs')
        repr = get.or.create(paste0(tsne.file.name, '_', j), function() Rtsne(X=omic.dist, dims=2, check_duplicates=F, is_distance=T)$Y)
      } else {
        repr = get.or.create(paste0(tsne.file.name, '_', j), function() Rtsne(X=t(omic), dims=2, check_duplicates=F)$Y)
      }
      ylabel = ifelse(j == 1, names(all.clusterings)[[i]], '')
      title = ''
      plot(repr[,1], repr[,2], pch=21, bg=cols[as.numeric(as.factor(clustering))], ylab=ylabel, xlab='', main=title, cex.lab=5, cex.main=5, cex=3)
    }
  }
  if (!is.null(plot.file.name)) {
    dev.off()
  }
}

plot.all.pairwise.ari <- function(monet.ret.list, fig.path, labels=NULL) {
  all.adj.per.dataset = lapply(1:length(monet.ret.list), function(i) {
    monet.ret = monet.ret.list[[i]]
    all.clusterings = lapply(monet.ret$all.monet.rets, function(x) x$clustering)
    all.adj = c()
    if (is.null(labels)) {
      for (i in 2:length(all.clusterings)) {
        for (j in 1:(i - 1)) {
          all.adj = c(all.adj, adjustedRandIndex(all.clusterings[[i]], all.clusterings[[j]]))
        }
      }
    } else {
      for (i in 1:length(all.clusterings)) {
        all.adj = c(all.adj, adjustedRandIndex(all.clusterings[[i]][names(labels)], labels))
      }
    }
    return(all.adj)
  })
  names(all.adj.per.dataset) = names(monet.ret.list)
  png(fig.path)
  boxplot(all.adj.per.dataset, ylim=c(0, 1))
  dev.off()
}

get.dataset.dir.path <- function() {
  return(CONFIG$dataset.dir.path)
}

get.or.create <- function(data.path, func=NULL) {
  if (file.exists(data.path)) {
    var.name = load(data.path)
    obj = get(var.name)
    return(obj)
  } else {
    obj = func()
    save(obj, file=data.path)
    return(obj)
  }
}

load.all <- function() {
  # libraries for R.3.5.3
  library('MOFA2')
  library('RColorBrewer')
  library('fields')
  library('mclust')
  library('SNFtool')
  library('cluster')
  library('survival')
  library('parallel')

  library('clusternomics')
  library('uuid')
  library('fields')
  library('VennDiagram')
  library('pheatmap')
  library('mclust')
  library('igraph')
  library('umap')
  library('RColorBrewer')
  library('Rtsne')
  library('bayesCC')
  library('SNFtool')
  library('parallel')
  library('survival')
  library('mcclust')
  library('mixtools')
  library('fossil')
  library('parallel')
  library('ggplot2')
  library('rmatio')
  library('twl')
  
  source(CONFIG$nemo.path)
  source(CONFIG$logrank.code.path)
  source(CONFIG$mdi.script.path)
}

CONFIG = list(# number of cores to use
              ncores=20, 
	      # path to a script activating MONET's python virtual environment (in the experiments here MONET was run in a virtualenv)
              python.env.path='/path/to/bin/activate', 
	      # path to MONET's main python file (main_bash.py)
              monet.main.path='/path/to/code/multiomic_dim_reduction/main_bash.py',
	      # path to a directory containing all clinical data for clinical parameter enrichment calculations.
	      # See supplementary information to see how to obtain all used clinical data.
	      clinical.params.dir='/path/to/benchmark/clinical/',
	      # path to a directory where temporary files will be written for MONET
              tmp.monet.dir='/path/to/tmp_monet_dir/',
	      # path to NEMO's R file (alternatively, NEMO's package can be installed from github)
              nemo.path='/path/to/nemo/NEMO.R',
	      # path to mdi executable file (the file we run was called mdipp)
              mdi.path='/path/to/mdipp_bin/bin/mdipp',
	      # path to directory for mdi temporaray files
              tmp.mdi.dir='/path/to/tmp_mdi_dir/',
	      # path to R script analysing mdi's output
              mdi.script.path='/path/to/mdipp-1.0.1.tar/mdipp-1.0.1/scripts/analysis.R',
	      # path for twl's temp files
              twl.tmp.dir ='/path/to/tmp',
	      # path to code for permutation based survival analysis (obtained from https://github.com/Shamir-Lab/Logrank-Inaccuracies)
              logrank.code.path='/path/to/perm_logrank/cond_logrank.R',
	      # path to all TCGA datasets.
	      # See supplementary information to see how to obtain all these datasets.
              dataset.dir.path='/path/to/datasets/',
	      # path to direcotry where all results are saved
              results.dir.path='/path/to/results/',
	      # path to direcotry where all plots are saved
              plots.dir.path='/path/to/plots/',
	      # path to direcotry containing all TCGA ovarian cancer mutations. 
              mutation.dir.path='/path/to/ovarian_mutations/',
	      # path to digits dataset.
	      # See supplementary information to see how to obtain the dataset.
              digits.path='/path/to/mv_datasets/handwritten/',
	      # path to single-cell dataset (this is rds file that can be loaded with readRDS).
	      # See supplementary information to see how to obtain the dataset (file name is scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds).
              scnmt.data.path='/path/to/mv_datasets/scnmt/scnmtseq_gastrulation_mae_826-cells_orderedFeatures.rds',
	      # path to single-cell metadata file.
	      # See supplementary information to see how to obtain the metadata file name is sample-metadata_matching-cells.txt).
              scnmt.metadata.path='/path/to/mv_datasets/scnmt/sample-metadata_matching-cells.txt',
	      # path to a directory where the microarray data is saved.
	      # See supplementary information to see how to obtain the data.
              microarray.dir.path='/path/to/mv_datasets/protein/',
	      # path to all datasets where the ground truth is known (simulations and digits dataset).
              ground.truth.datasets.path='/path/to/mv_datasets/',
	      # path to the directory containing MONET's python code
	      monet.python.dir='/path/to/code/multiomic_dim_reduction/')


