
# R functions for processing bobaseq data.
# Derived from Btcompl.R, but more modular

# The key functions are:
# findCodesFiles -- find the .codes file for a sample
# getCounts -- read in the codes files and make a table of counts per sample
# time0Counts -- aggregate counts across Time0 sample(s)
# logRatios -- compute log ratios
# highBarcodes -- identify statistically significant hits
# confirmedByOverlap -- determine which hits overlap
# highBarcodeReplicates -- identify barcodes that are significant in both replicates
# topHitPerRegion -- cluster the hits from hiRep into regions
# regionsToCandidates -- all the candidate proteins for each region
# candidateSimilarity -- candidate proteins that are supported by similarity
# topCandAddInserts, topCandAddSim -- compute metrics for each candidate
# topProtPerRegion -- select the best candidate protein for each region
# show -- show the log ratio for each insert across a region
#   (can also average across replicates)

# See fit2.R for an example of how to use these to load bobaseq data,
# compute fitness values for each barcode, identify statistically
# significant barcodes, cluster the hits into regions, and select
# the best candidate protein behind each region's benefit.

# Bugs: If two proteins in a region are similar to each other, this is
#  incorrectly considered as support for that region.

# The samples table must include set, sampleId, and index, and,
# each set must be a subdirectory (or symbolic link) in the base directory
findCodesFiles = function(samples, baseDir="/usr2/people/mprice/workbook/htcompl") {
  mapply(findCodesFile, as.character(samples$set), samples$sampleId, samples$index, baseDir);
}

# The samples table must include name and file; maps must include barcode
getCounts = function(samples, maps) {
  counts = data.frame(barcode = maps$barcode, stringsAsFactors=FALSE);
  for (i in 1:nrow(samples)) {
    cat("Loading", samples$name[i], "\n");
    counts[[ samples$name[i] ]] = readCountsFile(samples$file[i], counts$barcode);
  }
  return(counts);
}

# Compute some statistics for each sample, like
# total usable counts (in the correct libraries)
# fraction of reads in the libraries
# across libraries that should be present, the minimum fraction of reads
#
# samples must include name and lib; maps must include barcode and lib
# strain numbers include only strains with >1 count
totalCounts = function(counts, samples, maps) {
  out = data.frame(tot=rep(NA, nrow(samples)), totLib=NA, nStrains=NA, minNStrainsLib=NA, minStrainsLib="",
    stringsAsFactors=F);
  for(i in 1:nrow(samples)) {
    cat("sample", samples$name[i], "\n");
    countsThis = counts[[ as.character(samples$name[i]) ]];
    lib = samples$lib[i];
    libs = libToLibs(lib);
    u = counts$barcode %in% maps$barcode[maps$lib %in% libs];
    stopifnot(any(u));
    out$tot[i] = sum(countsThis);
    out$totLib[i] = sum(countsThis[u])
    out$nStrains[i] = sum(countsThis > 1);
    nLibStrains = sapply(libs, function(lib)
       sum(counts$barcode[countsThis > 1] %in% maps$barcode[maps$lib==lib]));
    names(nLibStrains) = libs;
    out$minNStrainsLib[i] = min(nLibStrains);
    out$minStrainsLib[i] = names(nLibStrains)[which.min(nLibStrains)];
  }
  out$fLib = out$totLib/out$tot;
  out$minStrainFrac = out$minNStrainsLib / out$nStrains;
  return(out);
}

# samples must include name (matching column names in counts), t0set, and used
time0Counts = function(samples, counts, maps) {
  countTime0 = data.frame(barcode=counts$barcode, stringsAsFactors=F);
  t0s = subset(samples, Group=="Time0" & used);
  for (t0set in unique(t0s$t0set)) {
    t0names = t0s$name[t0s$t0set == t0set];
    libs = t0s$lib[t0s$t0set == t0set];
    if (!all(libs == libs[1])) stop("Mismatching libs for t0set = ", t0set);
    countTime0[[t0set]] = rowSums(counts[,t0names,drop=F]);
  }
  return(countTime0);
}

# samples must include name, Group, lib, used
logRatios = function(samples, counts, countTime0, maps) {
  lr = data.frame(barcode = counts$barcode, stringsAsFactors=F);
  for (i in which(samples$Group != "Time0" & samples$used)) {
    n = as.character(samples$name[i]);
    t0set = as.character(samples$t0set[i]);
    if (!t0set %in% names(countTime0)) {
      cat("Skipping sample with no matching Time0: ",
           n, " ", samples$Group[i], " ", samples$desc[i], " ", t0set, "\n");
      next;
    }
    t0lib = samples$lib[samples$Group=="Time0" & samples$t0set == t0set][1];
    t0libs = libToLibs(t0lib);
    libs = libToLibs(samples$lib[i]);

    t0 = countTime0[[t0set]];
    if(t0lib != samples$lib[i]) {
      cat("Warning, lib for", t0set, " does not match lib for experiment ", samples$name[i], "\n");
      stopifnot(all(libs %in% t0libs));
    }
    # only consider barcodes that should have been in the experimental sample
    # this allows for cases where the Time0s were pooled across libraries
    # but the experimental samples were not
    t0 = ifelse(maps$lib %in% libs, t0, 0);
    tot0 = sum(t0);
    countThis = counts[[n]];
    countThis = ifelse(maps$lib %in% libs, countThis, 0);
    totThis = sum(countThis);
    stopifnot(tot0 > 0 && totThis > 0);
    # Balance the pseudocounts so that a count of 0 vs. 0 gives a log ratio of 0
    psi = sqrt(totThis/tot0);

    u = counts$barcode %in% maps$barcode[maps$lib %in% libs];
    lr[[n]] = ifelse(u, log2( (psi + countThis) / (1/psi + t0) ) - log2(totThis/tot0), NA);
  }
  return(lr);
}

highBarcodes = function(samples, counts, countTime0, lr, maps,
                        minLogRatio=5, minZ=4) {
  hiIndex = which(lr[,-1] >= minLogRatio, arr.ind=T);
  hiRow = hiIndex[,1];
  hiCol = hiIndex[,2];
  # Need to use counts[,-1] so that all types are uniform, otherwise
  # counts[hiCounts] has problems.
  lrColToCountsCol = match(names(lr)[-1], names(counts)[-1]);
  hiCountsIndex = as.matrix(data.frame(row=hiRow, col=lrColToCountsCol[hiCol]));
  hi = data.frame(barcode=lr$barcode[hiRow], name=names(lr)[-1][hiCol],
                  lr=lr[,-1][hiIndex],
                  n = counts[,-1][hiCountsIndex]);
  hi$t0set = samples$t0set[match(hi$name, samples$name)];
  time0CountsIndex = as.matrix(data.frame(row=hiRow, col=match(hi$t0set, names(countTime0)[-1])));
  hi$n0 = countTime0[,-1,drop=F][time0CountsIndex];

  # Add total counts and experiment metadata
  # (Do not include experiment lib, as that may be broader than the fragment lib added below)
  hi = merge(hi, samples[,c("name","background","Group","desc","totLib")]);
  tot0Totals = colSums(countTime0[,-1,drop=F]);
  hi = merge(hi, data.frame(t0set=names(tot0Totals), tot0=tot0Totals));

  # Add fragment metadata
  hi = merge(hi, maps[,c("barcode", "lib", "contig", "start", "end")]);
  hi = hi[order(hi$name, -hi$lr),];

  # Compute z
  psi = sqrt(hi$totLib/hi$tot0);
  hi$stdev = sqrt(1/(psi + hi$n) + 1/(1/psi + hi$n0)) / log(2);
  hi$z = hi$lr/hi$stdev;
  hi = subset(hi, z >= minZ);
  return(hi);
}

# For each entry in hi that is not unique (by t0set/background/Group/desc/lib/barcode)
# returns the average (including across non-high replicates, if there are multiple replicates).
# The input must include confirmedByOverlap.
highBarcodesReplicates = function(samples, lr, hi) {
  d = subset(hi, !is.unique(paste(t0set, background, Group, desc, lib, barcode)));
  if (nrow(d)==0) stop("No high barcodes replicate");

  n = aggregate(d[,"lr",drop=F], d[,words("t0set background Group desc lib barcode")], length);
  names(n)[names(n)=="lr"] = "nHi";
  mu = aggregate(d[,"lr",drop=F], d[,words("t0set background Group desc lib barcode")], mean);
  names(mu)[names(mu)=="lr"] = "avgHi";
  nConf = aggregate(d[,"confirmedByOverlap",drop=F], d[,words("t0set background Group desc lib barcode")], sum);
  out = merge(merge(n, mu), nConf);

  # compute the actual averages
  s = subset(samples, name %in% names(lr));
  s = split(s$name, paste(s$t0set, s$background, s$Group, s$desc));
  # some entries in this t0set could be NA if they did not include all libraries
  # so count non-NA entries to get #replicates
  lrN = lapply(s, function(snames) rowSums(!is.na(lr[,snames,drop=F])));
  iBarcode = match(out$barcode, lr$barcode);
  out$nRep = mapply(function(iBarcode, repsName) lrN[[repsName]][iBarcode],
    iBarcode, paste(out$t0set, out$background, out$Group, out$desc) );
  lrAvg = lapply(s, function(snames) rowMeans(lr[,snames,drop=F], na.rm=T));
  out$avg = mapply(function(iBarcode, repsName) lrAvg[[repsName]][iBarcode],
    iBarcode, paste(out$t0set, out$background, out$Group, out$desc) );
  out = merge(out, unique(hi[,c("barcode","lib","contig","start","end")]));
  out = out[order(out$t0set,out$background,out$Group,out$desc,-out$avg,out$barcode),];
  return(out);
}

# Report if a hit overlaps with another hit
# hits must include barcode, lib, contig, start, end
# matchBy might be the experiment or the Group/desc combination
confirmedByOverlap = function(hits, matchBy) {
  d = cbind(hits[,c("barcode","lib","contig","start","end")], matchBy=matchBy);
  cmp = merge(d, d, by=c("matchBy","lib","contig"), suffixes=c("","2"));
  cmp = subset(cmp, barcode != barcode2);
  cmp$overlap = overlaps(cmp$start, cmp$end, cmp$start2, cmp$end2);
  paste(d$barcode, d$matchBy) %in% paste(cmp$barcode, cmp$matchBy)[cmp$overlap];
}

# hi must include t0set, background, lib, contig, start, end, barcode
# cond is an arbitrary string to group the observations by
#  (but, they are always grouped by t0set, background, lib, and contig as wel)
# For regions with only 1 barcode, returns that 1
topHitPerRegion = function(hi, cond, score) {
  d = cbind(hi[,c("t0set","background","lib","barcode","contig","start","end")],
            row=1:nrow(hi), cond=cond, score=score);
  d = merge(d, d, by=c("t0set", "background", "lib", "cond", "contig"), suffixes=c("","2"));
  d = subset(d, barcode != barcode2 & overlaps(start, end, start2, end2));
  d = d[order(d$t0set, d$background, d$cond, -d$score),];
  # Will keep the best barcode in a region only.
  d = subset(d, score >= score2);
  # If x > y > z, then both x > y and x > z and y > z are in there;
  # removing duplicated barcode2 will remove y > z
  d = subset(d, !duplicated(paste(background, cond, barcode2)));
  d = subset(d, !duplicated(paste(background, cond, barcode)));

  return(hi[1:nrow(hi) %in% d$row | !confirmedByOverlap(hi, paste(hi$t0set,hi$background,hi$lib,cond)),]);
}

# Given a table with barcode and t0set, background, Group, desc, makes a table
# that lists the proteins in those barcodes and scores their average fitness
# in those conditions. Also adds the iRegion field.
regionsToCandidates = function(regions, samples, lr, p) {
  stopifnot(all(words("barcode t0set background Group desc") %in% names(regions)));
  regions$iRegion = 1:nrow(regions);
  cand = merge(regions, p[,c("barcode","locus_tag")]);
  # Compute average fitness across all barcodes containing this locus_tag
  cand$row = 1:nrow(cand);
  d = merge(cand[,words("t0set background Group desc row locus_tag")], p[,c("barcode","locus_tag")]);
  d$iBarcode = match(d$barcode, lr$barcode);
  d = split(d, d$row);
  avgs = sapply(d, function(x) {
    n = intersect(names(lr),
                  samples$name[samples$t0set==x$t0set[1] & samples$background==x$background[1]
                            & samples$Group == x$Group[1] & samples$desc == x$desc[1]]);
    if(length(n)==0) { print(x); stop("No data for these\n"); }
    # Some values could be NA due to t0set/lib mismatch
    mean(unlist(lr[x$iBarcode, n]), na.rm=T);
  });
  avgs = data.frame(row = names(d), avgProt = avgs);
  cand = merge(cand, avgs);
  cand$row = NULL;
  cand = cand[order(cand$t0set, cand$background, cand$Group, cand$desc, cand$locus_tag),];
  return(cand);
}

# Given a table of regions with candidates, with fields barcode, locus_tag, t0set, background, Group, desc,
# make a table of all pairs of similar proteins, with suffixes 1:2
#   Note -- t0set is not used to match on. Not sure that is ideal.
candidateSimilarity = function(topCand, psim) {
  stopifnot(all(words("barcode locus_tag t0set background Group desc") %in% names(topCand)));
  stopifnot(all(words("locus_tag1 locus_tag2 identity") %in% names(psim)));
  psimUse = subset(psim, locus_tag1 %in% topCand$locus_tag & locus_tag2 %in% topCand$locus_tag);
  sim = merge(psim[,words("locus_tag1 locus_tag2 identity")], topCand, by.x="locus_tag1", by.y="locus_tag");
  sim = merge(sim, topCand,
    by.x=words("locus_tag2 background t0set Group desc"),
    by.y=words("locus_tag background t0set Group desc"), suffixes=1:2);
  return(sim);
}

# For each candidate, compute the number of inserts that are in hi for this set of replicates
# (background/t0set/Group/desc).
# Returns a new data frame with the additional field nHiInserts.
#    (This can be run with either per-sample hi or per-replicate-group hiRep)
topCandAddInserts = function(topCand, hi) {
  d = merge(without(topCand, words("start end")), prot[,words("locus_tag start end")]);
  d = merge(d, hi, by=words("background t0set Group desc lib contig"));
  # Protein is contained within the insert
  d = subset(d, start.x >= start.y & end.x <= end.y);
  d$nHiInserts = 1;
  d = aggregate(d[,"nHiInserts",drop=F], d[,words("background t0set Group desc locus_tag")], sum);
  return(merge(topCand, d));
}

# For each candidate, compute the number of similar proteins that are candidates.
# Returns a new data frame with the additional field nSim
topCandAddSim = function(topCand, topCandSim) {
  d = aggregate(data.frame(nSim = rep(1, nrow(topCandSim))),
    topCandSim[,words("background t0set Group desc locus_tag1")], sum);
  names(d)[names(d)=="locus_tag1"] = "locus_tag";
  d = merge(topCand, d, all.x=T); # 884 rows, ok
  d$nSim[is.na(d$nSim)] = 0;
  return(merge(topCand, d));
}

# Select best protein per "region" (defined by barcode, t0set, background, Group, desc),
# using a combination of nHiInserts, nSim, and avgProt
# Renames desc to expDesc and adds the protId field and protDesc fields
topProtPerRegion = function(topCand, prot) {
  stopifnot(all(words("barcode locus_tag t0set background Group desc nHiInserts nSim avgProt") %in% names(topCand)));
  all = topCand;
  names(all)[names(all)=="desc"] = "expDesc";
  all = all[order(all$t0set, all$background, all$Group, all$expDesc,
                  all$barcode, -all$nHiInserts, -all$nSim, -all$avgProt),];
  out = subset(all, !duplicated(paste(t0set,background,Group,expDesc,barcode)));
  second = subset(all, duplicated(paste(t0set,background,Group,expDesc,barcode)));
  second = subset(second, !duplicated(paste(t0set,background,Group,expDesc,barcode)));
  out = merge(out,
    second[,words("t0set background Group expDesc barcode locus_tag nHiInserts nSim")],
    by=words("t0set background Group expDesc barcode"),
    all.x=T, suffixes=c("","2"));
  out = merge(out, prot[,c("locus_tag","protId","desc")]);
  names(out)[names(out)=="desc"] = "protDesc";
  return(out);
}


# Given a vector of locus tags, return a table of locus => cluster number
clusterLoci = function(loci, psim) {
  loci = unique(loci);
  if (length(loci) == 0) return(NULL);
  lsim = subset(psim, locus_tag1 %in% loci & locus_tag2 %in% loci);
  clust = 1:length(loci);
  lsim$i = match(lsim$locus_tag1, loci);
  lsim$j = match(lsim$locus_tag2, loci);
  iValues = sort(unique(lsim$i));
  lsim = split(lsim, lsim$i);
  for(i in iValues) {
    j = lsim[[as.character(i)]]$j;
    j = j[j < i];
    if (length(j) > 0) {
      # Put i in j[1]'s cluster
      clust[i] = clust[j[1]];
      # Merge clusters across the hits if necessary
      clusts = unique(clust[j]);
      if (length(clusts) > 1) clust[clust %in% clusts] = clust[j[1]];
    }
  }
  return(data.frame(locus_tag=loci, cluster=clust));
}

# Compare two experiments or replicates
# Use add=TRUE to add points to an existing plot
cmpPlot = function(name1, name2, lr, hi,
               xlab=NULL, add=FALSE, legend=!add, ylab=NULL, cex=1,
               minFit=3, ...) {
  stopifnot(all(name1 %in% names(lr)));
  stopifnot(all(name2 %in% names(lr)));
  x = rowMeans(lr[, name1, drop=F]);
  y = rowMeans(lr[, name2, drop=F]);
  i = which(x >= minFit | y >= minFit);
  if (length(i) == 0) {
    cat("All log ratios are close to zero, skipping plot for ", name1, " vs. ", name2, "\n");
    return;
  }

  if(is.null(hi$confirmedByOverlap)) hi$confirmedByOverlap = FALSE;

  if(is.null(xlab)) xlab = paste("Fitness in", paste(name1, collapse=", "));
  if(is.null(ylab)) ylab = paste("Fitness in ", paste(name2, collapse=", "));
  xymax = max(c(x,y), na.rm=T)+1;
  # in hi (significant) and confirmed by overlap is level 3; in hi is level 2
  conf = ifelse(lr$barcode[i] %in% hi$barcode[hi$name %in% c(name1,name2) & hi$confirmedByOverlap], 3,
           ifelse(lr$barcode[i] %in% hi$barcode[hi$name %in% c(name1,name2)], 2, 1));
  if (!add) plot(c(0,xymax), c(0,xymax),
                 xlab=xlab, ylab=ylab,
                 xaxs="i", yaxs="i", col=NA, bty="l", ...);
  points(pmax(0, x[i]), pmax(y[i],0),
       col=c("darkgrey","mediumblue","darkgreen")[conf],
       pch=c(1,4,5)[conf],
       xpd=T, cex=cex);
  if(!add) {
    abline(h=5, lty=2);
    abline(v=5, lty=2);
  }
  if(legend) {
    text(xymax, 0 + 2.5*strheight("A"), "significant", adj=c(1,0), col="mediumblue");
    text(xymax, 0 + 1.25*strheight("A"), "and confirmed", adj=c(1,0), col="darkgreen");
  }
}

# shows log-ratios for a region of a genome; use cond and optionally background to choose the
# experiment(s) (which are averaged) and use locus, contig/at, or barcode
# to choose the region.
# (locus may be a list of nearby locus tags, not just 1)
# Use rainbow to color code every barcode in the region (to allow easier comparison
#     across plots; every barcode will always get the same color, but they are not unique)
show = function(cond,
     lr=get("lr",env=globalenv()),
     samples=get("samples",env=globalenv()),
     prot=get("prot",env=globalenv()),
     maps=get("maps",env=globalenv()),
     locus=NULL, contig=NULL, barcode=NULL, at=NULL, width=8*1000,
     ymin = -2, ymax=18, background="Bt", t0set=NULL, rainbow=FALSE,
     col.default="black", col.fg="darkgreen", col.bg="darkred",
     arrow.begin=0.03, arrow.end=0.06, arrow.nonsig.lty=NULL,
     extraLabels=NULL, ylab=NULL,
     main="", debug=FALSE) {
  if (!is.null(locus) && length(locus) == 1 && grepl("^[ACGT]+", locus) && !locus %in% prot$locus_tag) {
    barcode = locus;
    locus = NULL;
  }
  if(!is.null(locus)) {
    genes = subset(prot, locus_tag %in% locus);
    if (nrow(genes) < length(locus)) stop("Cannot find prot named ", locus);
    stopifnot(all(genes$contig == genes$contig[1]));
    contig = genes$contig[1];
    at = mean( (genes$start+genes$end)/2 );
    cat("Focusing on ", genes$locus_tag, genes$desc, "\n");
  } else if (!is.null(barcode)) {
    b = maps[maps$barcode==barcode,];
    if(nrow(b) != 1) stop("No mapped barcode ", barcode);
    contig = b$contig;
    at = (b$start+b$end)/2;
    cat("Focusing on barcode ", b$barcode, "\n");
  }
  
  if (is.null(contig) || is.null(at)) stop("Must specify locus, barcode, or conting and at");
  left = at - width/2;
  if(left < 1) left = 1;
  right = at + width/2;
  i = which(maps$contig == contig & maps$end >= left & maps$start <= right);

  if (length(cond) > 1) {
    d = samples[samples$name %in% cond,];
  } else {
    d = samples[samples$name == cond | samples$desc == cond,];
    if (nrow(d) == 0) d = samples[grepl(cond, samples$desc, ignore.case=T),];
  }
  d = d[d$background == background,];
  if (nrow(d) == 0) stop("No experiments match ", cond, " for background ", background);
  if (!is.null(t0set)) {
    d = d[d$t0set %in% t0set,];
    if (nrow(d)==0) stop("No experiments in t0set ", t0set, " for this condition");
  }
  # Which sample are relevant to this contig?
  if (!is.null(contig)) {
    lib = unique(prot$lib[prot$contig == contig]);
    d$u = sapply(d$lib, function(x) lib %in% libToLibs(x));
    if (!any(d$u)) stop("No experiments match ", cond , " and library ", lib);
    d = d[d$u,];
  } else { # give up
    d$u = TRUE;
  }
  n = intersect(names(lr), d$name[d$u]);
  cat("Showing (average of) ", n, "\n");
  writeDelim(subset(samples, name %in% n)[, c("background", "t0set", "lib", "Group", "name", "sampleId", "desc")]);
  y = rowMeans(lr[, n, drop=F]);
  y = pmin(ymax, pmax(ymin, y));

  if (length(i) == 0) stop("No inserts overlap ", contig, " : ", left, " to ", right);
  col = rep(col.default, length(i));
  lwd = rep(1, length(i));
  if (!is.null(locus)) {
    within = mapply(function(s,e) all(s <= genes$start & e >= genes$end), maps$start[i], maps$end[i]);
#    if (!is.null(col.lo) && background == "Bt") {
#      n = samples$name[samples$used & samples$background=="Bt" & samples$Group != "control"];
#      lo = rowSums(counts[i,n]) < 10;
#      cat("lo\n"); print(table(lo));
#      cat("Both ", sum(lo & within), "\n");
#      col = ifelse(within, ifelse(lo, col.lo, col.fg), col.bg);
#      cat("col","\n"); print(table(col));
#    } else {
    col = ifelse(within, col.fg, col.bg);
  }

  if (!is.null(barcode)) {
    col = rep("darkgrey", length(i));
    col[maps$barcode[i] == barcode] = col.fg;
    lwd[maps$barcode[i] == barcode] = 2;
  }
  lty = 1;
  if (rainbow) {
    rainbowColors = c(1,2,3,4,6,8)
    col = 1 + (i %% length(rainbowColors));
    lty = 1 + (i %% 5);
  } else if (!is.null(arrow.nonsig.lty)) {
    nonsig = !lr$barcode[i] %in% hi$barcode[hi$name %in% n];
    lty = ifelse(nonsig, arrow.nonsig.lty, 1);
  }
    
  if (any(andNoNA(y[i] > 5))) {
    hibc = as.character(lr$barcode[i[ andNoNA(y[i] >= 5) ]]);
    cat("High barcodes: ", hibc, "\n");
    cat("library: ", unique(maps$lib[maps$barcode %in% hibc]), "\n");
  }
  condShow = if(background=="Bt") paste(cond, collapse=", ") else paste(background, cond, collapse=", ");
  if (is.null(ylab)) ylab = paste("Fitness in", condShow);
  main = sprintf("%s (%d exps.)", main, length(n));
  plot(c(left, right), c(ymin, ymax),
       pch=NA,
       xlab=paste("Position on", contig),
       ylab=ylab,
       main=main,
       bty="l", xaxs="i");
  beg = ifelse(maps$strand[i]=="+", maps$start[i], maps$end[i]);
  end = ifelse(maps$strand[i]=="+", maps$end[i], maps$start[i]);
  # angled arrowhead at end
  arrows(beg, y[i], end, y[i], code=2, length=arrow.end, col=col, lty=lty, lwd=lwd);
  # flat mark at beginning
  arrows(beg, y[i], end, y[i], code=1, length=arrow.begin, col=col, lty=lty, lwd=lwd, angle=90);

  i = which(prot$contig==contig & prot$end >= left & prot$start <= right);
  if (length(i) == 0) {
    cat("No protein-coding genes overlap ", contig, " : ", left, " to ", right);
  } else {
    top = ymax - 1.5 * strheight("A");
    d = prot[i,];
    d$col = "black";
    d$always = 0;
    if(!is.null(locus)) {
      d$col[d$locus_tag %in% locus] = "darkgreen";
      d$always[d$locus_tag %in% locus] = 1;
    }
    # remove text that extends past the edges of the plot
    d$label = sub(".*_", "_", d$locus_tag);
    d$width = strwidth(d$label);
    d$mid = (d$start + d$end)/2;
    d$label[ d$mid - d$width/2 < left | d$mid + d$width/2 > right ] = "";
    d$width = strwidth(d$label);
    d$below = 1:nrow(d) %% 2 == 0;
    d$y = top + ifelse(d$below, -1, 1) * strheight("A")/4;
    # remove overlapping labels
    d = d[order(d$below,d$mid),];
    for (i in 1:nrow(d)) {
      if (i > 1) {
        prev = d[i-1,];
        here = d[i,];
        if (prev$label != "" && here$label != "" 
            && prev$below == here$below
            && abs(prev$mid - here$mid) < (prev$width + here$width + strwidth("A"))/2) {
          if (here$always && prev$always) {
            d$label[i] = "";
          } else if (here$always) {
            d$label[i-1] = "";
          } else {
            d$label[i] = "";
          }
        }
      }
    }
    arrows(ifelse(d$strand=="+", d$start, d$end), d$y,
           ifelse(d$strand=="+", d$end, d$start), d$y,
           code = 2, length=0.05, col=d$col);
    text(d$mid,
         ifelse(d$below, top - 1.5*strheight("A"), top + 1.5*strheight("A")),
          d$label,
         adj=c(0.5, 0.5), col=d$col);
    text(left + strwidth("A")/2, top - 3*strheight("A"), paste("prefix: ", sub("_.*", "", d$locus_tag[1])),
         adj=c(0,1));
    if(!is.null(extraLabels)) {
      d2 = merge(d[,c("mid","locus_tag", "col", "below")], extraLabels, by="locus_tag");
      if (nrow(d2) > 0)
        text(d2$mid, top + ifelse(d2$below, 1.5, 2.7)*strheight("A"),
             d2$label, adj=c(0.5,0.5), col=d2$col, xpd=T, font=2);
    }
    cat("Protein-coding genes:\n");
    writeDelim(d[order(d$locus_tag),c("locus_tag", "protId", "start","end","strand","desc")]);
  }
}

# Given a library specifier like "Btheta" or "6" or "Pjohn,Pmerdae",
# return the list of relevant libraries
libToLibs = function(lib) {
  if (lib=="6") return(c("Bcaccae","Bfrag","Bsa","Buni","Pjohn","Pmerdae"));
  if (lib=="BtBuniPj") return(c("Btheta","Buni","Pjohn"));
  if (grepl(",", lib)) lib = words(lib, ",");
  lib[lib=="Bc"] = "Bcaccae";
  lib[lib=="Bf"] = "Bfrag";
  lib[lib=="Pj"] = "Pjohn";
  lib[lib=="Pm"] = "Pmerdae";
  lib[lib=="Bt"] = "Btheta";
  return(lib);
}

findCodesFile = function(set, sampleId, index, baseDir) {
  dir = paste0(baseDir, "/", set);
  if (sampleId=="") sampleId = index;
  patterns = c(paste0(sampleId, "[_.].*codes$"),
               paste0(index,"[_.].*codes$"),
               paste0(".*_",index,"_.*codes$"));
  for (pattern in patterns) {
    file = list.files(path=dir, pattern=pattern);
    if (length(file) == 1) return(paste0(dir,"/",file));
  }
  stop(sprintf("Error finding .codes file for set %s sampleId %s index %s\nusing patterns %s",
                                       set, sampleId, index, paste(patterns, collapse=" ")));
}

# prot must include locus_tag and protId (the identifiers in the BLAST hits)
# blastFile shuold be the all-vs-all blastp results in tab-delimited format
# lenFile should include the fields id and len (for length in amino acids)
proteinSimilarity = function(prot, blastFile, lenFile, minIdentity = 40, minCoverage = 0.75) {
  sim = read.delim(blastFile, as.is=T, header=F,
     col.names=c("query","subject","identity","alen","mm","gap","qbeg","qend","sbeg","send","eval","bits"));
  # Leave self-hits in for onw so that we can get cases
  # where locus tags differ but sequences are identical
  sim = unique(subset(sim, identity >= minIdentity));
  faaLen = read.delim(lenFile, as.is=T, header=F, col.names=c("id","len"));
  sim = merge(merge(sim, faaLen, by.x="query", by.y="id"),
              faaLen, by.x="subject", by.y="id", suffixes=c(".query",".subject"));
  sim = unique(sim);
  sim$cov.q = (sim$qend-sim$qbeg+1)/sim$len.query;
  sim$cov.s = (sim$send-sim$sbeg+1)/sim$len.subject;
  sim = subset(sim, cov.q >= minCoverage & cov.s >= minCoverage);
  sim$cov = pmin(sim$cov.q, sim$cov.s);
  psim = data.frame(protId1=sim$subject, protId2=sim$query, identity=sim$identity, cov=sim$cov);
  psim = merge(merge(psim, prot, by.x="protId1", by.y="protId"),
              prot, by.x="protId2", by.y="protId", suffixes=1:2);
  subset(psim, locus_tag1 != locus_tag2);
}

# Expects the file to have a barcode column and a second (arbitrarily named)
# column with the counts. Returns a vector of counts (for the mapped barcodes only)
readCountsFile = function(file, mappedBarcodes) {
    d = read.delim(file, as.is=T);
    row = match(mappedBarcodes, d$barcode);
    ifelse(is.na(row), 0, d[row,2]);
}

overlaps = function(beg1, end1, beg2, end2) beg1 <= end2 & end1 >= beg2;

writeDelim = function(table,file=stdout(),report=FALSE,col.names=T,...) {
	write.table(table,file,sep="\t",quote=FALSE,row.names=FALSE,col.names=col.names,...);
	if(report) cat("Wrote",nrow(table),"rows","to",file,"\n")
}

words = function(s, by=" ", ...) strsplit(s[1], by, ...)[[1]];

andNoNA = function(x) ifelse(is.na(x), FALSE, x);

na0 = function(x) ifelse(is.na(x),0,x);

is.unique = function(x) {
  counts = as.data.frame.table(table(x), dnn="x");
  return(x %in% counts$x[counts$Freq==1]);
}

without = function(list,columns=list$without) {
  for (i in columns) list[[i]] = NULL;
  return(list);
}
