args <- commandArgs(trailingOnly=T)
n_args <- length(args)

# Possible values for string matching sex
male_keywords <- c("Male","MALE","male","m","M","1")
female_keywords <- c("Female","female","FEMALE","F","f","2")

f_gammai_matrix <- as.character(args[1])
f_chrom_abs <- as.character(args[2])
f_qcmeta <- as.character(args[3])
f_windows <- as.character(args[4])
chrY_dup_factor <- as.numeric(args[5])
ab_covar <- as.character(args[6])
plot_flag <- as.logical(args[7])
output_dir <- as.character(args[8])

fit.chrY <- FALSE
if (chrY_dup_factor < 0) fit.chrY <- TRUE

qcmeta <- as.data.frame(read.table(f_qcmeta,header=T,row.names=NULL,sep="\t"))
last_covar_col <- length(colnames(qcmeta))-1
sex_colname <- tail(colnames(qcmeta), n=1)
if (! grepl("gender|sex", sex_colname, ignore.case=T) ) {
    sex_colname <- ""
    last_covar_col <- last_covar_col + 1
}
qc_covar_cols <- colnames(qcmeta)[2:last_covar_col]

# Check for updated/independent sex assignment file
if (n_args > 8) {
    f_sex <- as.character(args[9]);
    ks_sex <- as.data.frame(read.table(f_sex,header=T,row.names=NULL,sep="\t"))
    new_sex_colname <- colnames(ks_sex)[2]
    pre_len <- nrow(qcmeta)
    qcmeta <- merge(qcmeta, ks_sex, by="SampleID")
    post_len <- nrow(qcmeta)
    if (pre_len != post_len) {
        print(paste0("The sample lists in ",f_qcmeta," and ",f_sex,
            " do not match. Proceeding with only the intersection."));
    }
    if (paste0(new_sex_colname,".x") %in% colnames(qcmeta)) {
        qcmeta[[new_sex_colname]] <- qcmeta[[paste0(new_sex_colname,".y")]]
        qcmeta <- qcmeta[,!(names(qcmeta) %in% c(paste0(new_sex_colname,".x"),
                    paste0(new_sex_colname,".y")))]
    }
    sex_colname <- new_sex_colname
}

if (sex_colname == "") {
    print(paste0("Sex column not recognized in input data. Expecting the last column ",
        "to represent sex, with column name containing \"gender\" or \"sex\". ",
        "Alternatively, a separate sex mapping file can be provided. Quitting"));
    quit()
}

# Create matrix to hold D-statistics for pairwise chromosome GC-content comparisons
windows <- as.data.frame(read.table(f_windows,header=F,row.names=NULL,sep="\t"))
lastcol <- length(colnames(windows))
colnames(windows)[1:3] <- c("CHR","START","END")
colnames(windows)[(lastcol-2):lastcol] <- c("GCFRAC","MAPPABILITY","FLAG")

ds <- matrix(0, nrow=24, ncol=24)
dimnames(ds) <- list(as.vector(unique(windows$CHR)), as.vector(unique(windows$CHR)))
for (i in 1:nrow(ds)) {
  for (j in 1:nrow(ds)) {
    col.i <- colnames(ds)[i]; col.j <- colnames(ds)[j];
    ks.out <- ks.test(windows$GCFRAC[windows$CHR == col.i], windows$GCFRAC[windows$CHR == col.j])
    ds[i,j] <- ks.out$statistic
}}

ds.df <- as.data.frame(ds)
chr.adjusts <- data.frame(chrom=character(), adjustment=character())
for (i in 1:nrow(ds.df)) {
  row.i <- row.names(ds.df)[i];
  adjust <- paste0("chr",row.i,"\n","~ (",
                paste(colnames(sort(ds.df[row.i,])[c(2,3)]),sep="",collapse=", "),
                ")");
  chr.adjusts <- rbind(chr.adjusts, data.frame(chrom=paste0("chr",row.i), adjustment=adjust))
}

# Read coverage (gamma_i) matrix and merge on QC meta-data (including sex)
coverage <- as.data.frame(read.table(f_gammai_matrix,sep="\t",header=T,row.names=NULL))
pre_len <- nrow(coverage)
coverage <- merge(coverage, qcmeta, by="SampleID")
post_len <- nrow(coverage)
if (pre_len != post_len) {
  print(paste0("The sample lists in the QC meta-data (",f_qcmeta,") and gamma_i matrix (",
                f_gammai_matrix,") do not match. Proceeding with only the intersection."));
}

# Fit gamma_i regression models
for (i in 1:nrow(ds.df)) {
  row.i <- row.names(ds.df)[i]
  if (row.i == "Y" & !fit.chrY) next
  
  chr.i <- paste("chr",row.i,sep="")
  qc.covars <- paste0("+",paste(qc_covar_cols,sep="",collapse="+"))
  if (row.i == "X") {
    qc.covars <- paste0(qc.covars,"+",sex_colname);
  }
  if (row.i == "Y") {
    fit.chr <- lm(as.formula(paste0(chr.i," ~ chr",
                    paste(colnames(sort(ds.df[row.i,])[c(2,3)]),sep="",collapse="+chr"),
                    qc.covars)),data=coverage[coverage[[sex_colname]] %in% male_keywords,])
  } else {
    fit.chr <- lm(as.formula(paste0(chr.i," ~ chr",
                    paste(colnames(sort(ds.df[row.i,])[c(2,3)]),sep="",collapse="+chr"),
                    qc.covars)),data=coverage)

    coverage[[paste(chr.i,".hat",sep="")]] <- hatvalues(fit.chr)
    coverage[[paste(chr.i,".hat.np",sep="")]] <- coverage[[paste(chr.i,".hat",sep="")]]*
                nrow(coverage)/(length(qc_covar_cols)+3) # QC metrics, sex, gamma_i, +1
  }

  # For large sample sizes, refit models after removing high-leverage samples. Whether or not
  # this is necessary depends entirely on the sample set. It effectively reduces standard
  # error estimates for well-behaved samples.
  # Hard-coding for now to be >1k samples. Note hat/hat.np values are not re-calculated here.
  if (nrow(coverage) > 1000 & row.i != "Y") { 
    coverage.tmp <- subset(coverage, coverage[[paste(chr.i,".hat.np",sep="")]] < 5)
    fit.chr <- lm(as.formula(paste0(chr.i," ~ chr",
                    paste(colnames(sort(ds.df[row.i,])[c(2,3)]),sep="",collapse="+chr"),
                    qc.covars)),data=coverage.tmp)
  }
  # Apply new model to full data set
  predict.chr <- predict.lm(fit.chr,
                    newdata=coverage,
                    interval=c("prediction"),
                    level=0.95,
                    se.fit=T)
  coverage[[paste(chr.i,".fit",sep="")]] <- predict.chr$fit[,1]
  coverage[[paste(chr.i,".lwr",sep="")]] <- predict.chr$fit[,2]
  coverage[[paste(chr.i,".upr",sep="")]] <- predict.chr$fit[,3]
  coverage[[paste(chr.i,".se.fit",sep="")]] <- predict.chr$se.fit
  coverage[[paste(chr.i,".res",sep="")]] <- coverage[[chr.i]] - 
        coverage[[paste0(chr.i,".fit")]]
  coverage[[paste(chr.i,".p",sep="")]] <- 2*pnorm(-abs(coverage[[paste0(chr.i,".res")]]/
        (coverage[[paste0(chr.i,".se.fit")]]*sqrt(nrow(coverage)))))
  coverage[[paste(chr.i,".q",sep="")]] <- p.adjust(coverage[[paste0(chr.i,".p")]],
        method="fdr")
  
  # If fitting chrY, we only train on male samples but predict on all. Predictions for
  # female samples will be meaningless and skew the output, NA all the relevant values
  if (row.i == "Y") {
    coverage[coverage[[sex_colname]] %in% female_keywords,
             c("chrY.fit","chrY.lwr","chrY.upr",
                "chrY.se.fit","chrY.res","chrY.p","chrY.q")] <- rep(NA,7)
    coverage$chrY.p <- 2*pnorm(-abs(coverage$chrY.res / (coverage$chrY.se.fit * sqrt(
            nrow(coverage[coverage[[sex_colname]] %in% male_keywords,])
        ))))
    coverage$chrY.q <- p.adjust(coverage$chrY.p, method="fdr")
  }
}

write.table(coverage,
  paste0(output_dir,"/gammai_regressed.txt"),
  sep="\t",
  quote=F,
  row.names=F,
  col.names=T)

# Calculate percentiles of y-hat (Update plot function to reflect these cutoffs)
hat.np.all <- c(
  coverage$chr1.hat.np, coverage$chr2.hat.np, coverage$chr3.hat.np,
  coverage$chr4.hat.np, coverage$chr5.hat.np, coverage$chr6.hat.np,
  coverage$chr7.hat.np, coverage$chr8.hat.np, coverage$chr9.hat.np,
  coverage$chr10.hat.np, coverage$chr11.hat.np, coverage$chr12.hat.np,
  coverage$chr13.hat.np, coverage$chr14.hat.np, coverage$chr15.hat.np,
  coverage$chr16.hat.np, coverage$chr17.hat.np, coverage$chr18.hat.np,
  coverage$chr19.hat.np, coverage$chr20.hat.np, coverage$chr21.hat.np,
  coverage$chr22.hat.np, coverage$chrX.hat.np)

# Leverage cutoffs for plot.aneuploidy function
cutoffs <- data.frame(
  Percentile=c("99","99.5"),
  hat_val=quantile(hat.np.all,
  prob=c(.99,.995)))

# Build data frame of putatively anomalous samples via coverage
cov.anomalies <- data.frame(
  SampleID = character(),
  Chromosome = character(),
  Fraction = numeric(),
  p = numeric(),
  q = numeric(),
  hatnp=numeric())

for (i in 1:nrow(ds.df)) {
  chri <- paste0("chr",row.names(ds.df)[i]);
  if (chri == "chrY" & !fit.chrY) next

  sample_set <- as.vector(coverage$SampleID[coverage[[paste0(chri,".p")]] <= 0.25 &
                                            !is.na(coverage[[paste0(chri,".p")]])]);
  if (length(sample_set) > 0) {
    fractions <- coverage[coverage$SampleID %in% sample_set, paste0(chri,".res")] /
      (0.5*coverage[coverage$SampleID %in% sample_set, paste0(chri,".fit")]);
    ps <- coverage[coverage$SampleID %in% sample_set, paste0(chri,".p")]
    qs <- coverage[coverage$SampleID %in% sample_set, paste0(chri,".q")]
    hatnps <- rep(NA, length(sample_set))
    if (chri != "chrY")
        hatnps <- coverage[coverage$SampleID %in% sample_set, paste0(chri,".hat.np")]
    
    cov.anomalies <- rbind(
      cov.anomalies,
      data.frame(SampleID=sample_set, 
                  Chromosome=rep(chri, length(sample_set)),
                  Fraction=fractions,
                  p=ps,
                  q=qs,
                  hatnp=hatnps));
  } else {
    print("No coverage anomalies identified!");
  }
}
# Add chrY dups using median method
if (!fit.chrY) {
  chrY_median <- median(coverage[coverage[[sex_colname]] %in% male_keywords,"chrY"])
  
  ydup_samples <- as.vector(coverage$SampleID[coverage$chrY > 
      chrY_median * chrY_dup_factor])
  fractions <- ((coverage[coverage$SampleID %in% ydup_samples, "chrY"] / chrY_median) - 1.0)*2
  ps <- rep(0.0,length(ydup_samples))
  qs <- ps
  hatnps <- ps
  cov.anomalies <- rbind(
      cov.anomalies,
      data.frame(SampleID=ydup_samples,
                  Chromosome=rep("chrY", length(ydup_samples)),
                  Fraction=fractions,
                  p=ps,
                  q=qs,
                  hatnp=hatnps));
}
# Finalize
cov.anomalies <- merge(cov.anomalies, coverage[,c("SampleID",sex_colname)], by="SampleID")
cov.anomalies$Fraction[cov.anomalies$Chromosome %in% c("chrX","chrY") & 
                       cov.anomalies[[sex_colname]] %in% male_keywords] <-
    cov.anomalies$Fraction[cov.anomalies$Chromosome %in% c("chrX","chrY") &
                           cov.anomalies[[sex_colname]] %in% male_keywords] / 2;

cov.anomalies$Tier <- 4
cov.anomalies$Tier <- cov.anomalies$Tier - ifelse(cov.anomalies$Fraction > 0.1 | cov.anomalies$p < 0.05,1,0)
cov.anomalies$Tier <- cov.anomalies$Tier - ifelse(cov.anomalies$q < 0.05,1,0)
cov.anomalies$Tier <- cov.anomalies$Tier - ifelse(abs(cov.anomalies$Fraction) > 0.15,1,0)

write.table( format(subset(cov.anomalies[with(cov.anomalies, order(SampleID)),],Tier <= 3), digits=6), 
    paste0(output_dir,"/coverage_anomalies.txt"),
    sep="\t",quote=F,row.names=F,col.names=T)

# Analyzed chromosome-wide median het allele balance stats
chrom_het_abs <- as.data.frame(read.table(
    f_chrom_abs,
    sep="\t",
    header=F,
    row.names=NULL,
    col.names=c("Chr","AB","NHetSNP","lnP","SampleID")))

# Drop the preliminary p-value (lnP column), we are going to compute it more accurately
chrom_het_abs <- chrom_het_abs[,(names(chrom_het_abs) != "lnP")]

# Merge QC meta-data for regression covariate and sex assignment
chrom_het_abs <- merge(chrom_het_abs, qcmeta[,c("SampleID",ab_covar,sex_colname)], by="SampleID")
colnames(chrom_het_abs)[names(chrom_het_abs) == ab_covar] <- "QC_covar"
colnames(chrom_het_abs)[names(chrom_het_abs) == sex_colname] <- "Sex"

ab.anomalies <- data.frame(
    SampleID=character(),
    Chr=character(),
    AB=numeric(),
    NHetSNP=numeric(),
    QC_covar=numeric(),
    Sex=character(),
    AB.fit=numeric(),
    AB.res=numeric(),
    AB.p=numeric(),
    AB.q=numeric())

for (chri in row.names(ds.df)) {
  chr.subset <- chrom_het_abs[chrom_het_abs$Chr == chri,]
  if (chri == "X") {
    chr.subset <- chr.subset[chr.subset$Sex %in% female_keywords |
                            (chr.subset$NHetSNP > 75 & chr.subset$AB > 0.15),]
  }
  fit.chr <- lm("AB ~ QC_covar", data=chr.subset)
  chr.subset$fit <- fit.chr$fitted.values
  chr.subset$res <- residuals(fit.chr)
  chr.subset$AB.p <- pnorm(chr.subset$res/summary(fit.chr)$sigma)
  chr.subset$AB.q <- p.adjust(chr.subset$AB.p, method="fdr")

  ab.anomalies <- rbind(
    ab.anomalies,
    chr.subset[chr.subset$AB.p < 0.05 | 
        (chr.subset$Chr == "X" & chr.subset$Sex %in% male_keywords),
        c('SampleID','Chr','AB','QC_covar','Sex','fit','res','AB.p','AB.q')])
}
ab.bonf.corrector <- length(subset(chrom_het_abs,
                                    Chr != "X" |
                                    Sex %in% female_keywords |
                                    (NHetSNP > 75 & AB > 0.15) )[,1])
ab.anomalies$Tier <- ifelse(ab.anomalies$AB.p < (0.05/ab.bonf.corrector) |
                            (ab.anomalies$Chr == "X" &
                             ab.anomalies$Sex %in% male_keywords),
                            1, ifelse(ab.anomalies$AB.q < 0.05, 2, 3))

write.table(format(ab.anomalies[with(ab.anomalies, order(SampleID)),], digits=6),
    paste0(output_dir,"/chrom_ab_anomalies.txt"),
    sep="\t",quote=F,row.names=F,col.names=T)

if (!plot_flag) {
    quit()
}

# PLOTTING
suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(gridBase)
    library(gridExtra)
    library(stringr)
    library(reshape2)
})

# Helper function to plot individual sample coverage distributions
# Requires that the following variables are pre-defined:
# - covearge (prefit gammai matrix)
# - chr.adjusts (covariate chromosomes for each chromosome)
# - cutoffs (hat.np 99 and 99.5th percentiles for leverage flagging)
# - output_dir (plots will be saved to output_dir/plots/)
karyoscan.plot.save <- function(sampleid) {
  tmp <- melt(subset(coverage, SampleID == sampleid));
  tmp1 <- data.frame(str_split_fixed(tmp$variable,"[.]",2));
  colnames(tmp1) <- c("chrom","variable.split");
  tmp <- cbind(tmp, tmp1);
  tmpd <- dcast(data = tmp,
                formula = chrom ~ variable.split,
                fun.aggregate = sum,
                value.var = "value");
  sex <- subset(coverage, SampleID == sampleid)[['KARYOTYPIC_SEX']];
  if (sex == "Male") {
    tmpd$res[tmpd$chrom == "chrX"] <- tmpd$res[tmpd$chrom == "chrX"]/2 
  };
  tmpdm <- merge(tmpd, chr.adjusts, by="chrom");
  colnames(tmpdm)[2] <- "value";
  tmpdm$chrom <- factor(tmpdm$chrom, levels=c("chr1","chr2","chr3","chr4","chr5","chr6",
                    "chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15",
                    "chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX"));
  tmpdm$sig <- ifelse(tmpdm$q <= 0.05, "**", ifelse(tmpdm$p <= .05, "*",""));

  cutoff_995 <- cutoffs$hat_val[cutoffs$Percentile == "99.5"]
  cutoff_99 <- cutoffs$hat_val[cutoffs$Percentile == "99"]
  tmpdm$lev <- ifelse(tmpdm$hat.np >= cutoff_995, "!!",
                 ifelse(tmpdm$hat.np >= cutoff_99, "!", ""))

  title <- paste(sampleid, " (", sex, ")", sep="");
  p1 <- ggplot(tmpdm[tmpdm$fit > 0.001,], aes(x=chrom,y=((res)/(0.5*fit))))
  p1 <- p1 + geom_point(aes(shape=lev),size=3) 
  p1 <- p1 + geom_errorbar(aes(ymin=(lwr-fit)/(0.5*fit), ymax=(upr-fit)/(0.5*fit)),width=0.25) 
  p1 <- p1 + geom_hline(yintercept=1,lty="dashed") 
  p1 <- p1 + geom_hline(yintercept=-1, lty="dashed") 
  p1 <- p1 + xlab("") + ylab("Normalized Coverage Ratio") 
  p1 <- p1 + theme(axis.text = element_text(size=12,colour="black"),
    text = element_text(size=15,colour="black"), axis.text.x = element_text(angle=90, vjust=.5)) 
  p1 <- p1 + ggtitle(title) 
  p1 <- p1 + geom_text(aes("chr1", 1.05, label="Estimated Ratio of Trisomy", hjust=0), colour="black") 
  p1 <- p1 + geom_text(aes("chr1", -1.05, label="Estimated Ratio of Monosomy", hjust=0), colour="black") 
  p1 <- p1 + geom_text(aes(label=sig, size=sig, vjust=((ceiling(res)+1)%%2)*1.7)) 
  p1 <- p1 + scale_shape_manual(name="Leverage",labels=c("Non-Outlier","99th Percentile",
    "99.5th Percentile"),values=c(16,1,13),limits=c("","!","!!")) 
  p1 <- p1 + scale_x_discrete(labels=tmpdm[order(tmpdm$chrom),'adjustment']) 
  p1 <- p1 + scale_size_manual("Significance",labels=c("p < .05","q < .05"),values=c(8,8,8),
    breaks=c("*","**"),limits=c("","*","**"));

  p1.grob <- ggplotGrob(p1);
  p1.grob$grobs[[8]]$grobs[[2]]$grobs[[4]]$label <- "*";
  p1.grob$grobs[[8]]$grobs[[2]]$grobs[[6]]$label <- "**";
  p1.grob.2 <- arrangeGrob(p1.grob)
  ggsave(paste0(output_dir, "/plots/", sampleid,".read_coverage.png"), p1.grob.2, units="in", width=16, height=10);
}

# Plot each sample
for (sample in coverage$SampleID) {
    karyoscan.plot.save(sample);
}
