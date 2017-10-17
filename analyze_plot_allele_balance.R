infile <- as.character(commandArgs(trailingOnly=T)[1]);
samplename <- as.character(unlist(strsplit(as.character(tail(unlist(strsplit(infile,"/")),n=1)),"[.]"))[[1]])
rohfile <- as.character(commandArgs(trailingOnly=T)[2]);
centrofile <- as.character(commandArgs(trailingOnly=T)[3]);
homAB.threshold <- as.numeric(commandArgs(trailingOnly=T)[4]);
plot_flag <- as.logical(commandArgs(trailingOnly=T)[5]);

df <- as.data.frame(read.table(
    infile,
    col.names=c("Chr", "Pos", "AB"),
    colClasses=c("factor", "numeric", "numeric"),
    sep="\t",
    header=F,
    row.names=NULL
));

df$Chr <- factor(df$Chr, levels=c(as.character(seq(1,22)),'X','Y'))

# Compute running median with window size of 19 SNPs (only SNPs with AB > homAB.threshold)
# and merge as new column onto df
warn_lvl <- getOption("warn")
options(warn=-1)
df.tmp <- df[0,]
for (i in c(as.character(seq(1,22)),"X","Y")) {
    df.tmp.chr <- df[df$AB > homAB.threshold & df$Chr == i,];
    df.tmp.chr$runmed <- runmed(df.tmp.chr$AB, k=19, endrule="constant");
    df.tmp <- rbind(df.tmp, df.tmp.chr);
}
options(warn=warn_lvl)
df <- merge(df, df.tmp, by=c('Chr','Pos','AB'), all.x=T);

# Compute chromosome-wide median (standard) allele balances (ABs)
chr.AB.medians.meds <- aggregate(AB ~ Chr, data=df[df$AB > homAB.threshold & !is.na(df$runmed),], median)
chr.AB.medians.lens <- setNames(aggregate(AB ~ Chr, data=df[df$AB > homAB.threshold & !is.na(df$runmed),], length),c("Chr","NHetSNPs"))
chr.AB.medians <- merge(chr.AB.medians.meds, chr.AB.medians.lens, by="Chr")

# Read in centromere location info (defined as the nearest exon boundary to the centromere)
centromeres <- as.data.frame(read.table(
    centrofile,
    sep="\t",
    header=F,
    row.names=NULL,
    col.names=c("Chr","Start","End","ChrEnd","ExpectedMedAB"),
    colClasses=c("factor","numeric","numeric","numeric","numeric")
));

# Search for regions of major deviations in running median vs. chromosome-wide median
df.devs <- data.frame(Chr=character(), Start=numeric(), End=numeric(), Nhet=numeric(), Area=numeric(), AB=numeric());
cur.chr <- df.tmp$Chr[1];
cur.chr.median <- chr.AB.medians$AB[chr.AB.medians$Chr == cur.chr];
cur.cent.start <- centromeres$Start[centromeres$Chr == cur.chr];
cur.cent.end <- centromeres$End[centromeres$Chr == cur.chr];
cur.chr.end <- centromeres$ChrEnd[centromeres$Chr == cur.chr];
cur.start <- -1;
cur.end <- -1;
cur.hets <- 0;
cur.area <- -1;
cur.min.ab <- 0.5;  # Minimum of the smoothed (running median) allele balance
last.runmed <- 0;   # This is the pre-computed running median allele balance.
for (i in 1:nrow(df.tmp)) {
    snp <- df.tmp[i,];
    if (snp$Chr != cur.chr | snp$runmed >= cur.chr.median | i == nrow(df.tmp)) { # Triggers for end of active anomalous AB event
        if (cur.start >= 0) { # Check occurrence of event
            # Normalize area for callable length of chromosome (x) & chromosome median (y)
            cur.area <- cur.area / ((cur.chr.end - (cur.cent.end - cur.cent.start)) * cur.chr.median);
            df.devs <- rbind(df.devs, data.frame(Chr=as.character(cur.chr), Start=as.numeric(cur.start), End=as.numeric(cur.end), Nhet=as.numeric(cur.hets), Area=as.numeric(cur.area), AB=as.numeric(cur.min.ab)));
            cur.start <- -1;
            cur.end <- -1;
            cur.hets <- 0;
            cur.area <- -1;
            cur.avg.ab <- 0.5;
            last.runmed <- 0;
        }
        if (snp$Chr != cur.chr) {
            cur.chr <- snp$Chr;
            cur.chr.median <- chr.AB.medians$AB[chr.AB.medians$Chr == cur.chr];
            cur.chr.end <- centromeres$ChrEnd[centromeres$Chr == cur.chr];
            cur.cent.start <- centromeres$Start[centromeres$Chr == cur.chr];
            cur.cent.end <- centromeres$End[centromeres$Chr == cur.chr];
        }
    }
    if (snp$runmed < cur.chr.median) { # Flag for new or continuing event
        if (cur.start < 0) { # New event
            cur.start <- snp$Pos;
            cur.end <- snp$Pos;
            cur.hets <- 1;
            cur.area <- cur.chr.median - snp$runmed;
            cur.min.ab <- snp$runmed;
            last.runmed <- snp$runmed;
        } else { # Continuing event
            new.dist <- snp$Pos - cur.end;
            if (cur.end <= cur.cent.start & snp$Pos >= cur.cent.end ) {  # Bridging the centromere
                new.dist <- new.dist - (cur.cent.end - cur.cent.start);
            }
            new.area <- (((cur.chr.median - snp$runmed) + (cur.chr.median - last.runmed)) / 2) * new.dist; # ((y1+y2)/2) * delta(x);
            cur.area <- cur.area + new.area;
            if (snp$runmed < cur.min.ab) {
                cur.min.ab <- snp$runmed;
            }
            cur.end <- snp$Pos;
            cur.hets <- cur.hets + 1;
            last.runmed <- snp$runmed;
        }
    }
}

roh <- as.data.frame(read.table(
    rohfile,
    col.names=c("Chr","Start","End","Length","NSNPs","MaxConf","Sample"),
    colClasses=c("factor", "numeric", "numeric", "numeric", "numeric", "numeric", "character"),
    sep="\t",
    header=F,
    row.names=NULL
));

if (plot_flag) {
    suppressPackageStartupMessages({
        library(ggplot2)
    })
    ggplot_version <- as.numeric(substr(packageVersion("ggplot2"),0,1))
    
    abplot <- ggplot(df, aes(x=Pos, y=AB))
    abplot <- abplot + geom_point(size=0.75)
    abplot <- abplot + facet_wrap(~Chr,scales="free_x")
    if (ggplot_version > 1) {
        abplot <- abplot + geom_hline(data=aggregate(AB ~ Chr, df[df$AB > homAB.threshold,], median), aes(yintercept=AB), color='red',lwd=.75)
        abplot <- abplot + geom_line(data=df[!is.na(df$runmed),], aes(y=runmed, group=Chr), color='blue', linetype='longdash',lwd=.75)
    } else {
        abplot <- abplot + geom_line(data=df[df$AB > homAB.threshold,], stat='hline', yintercept=median, color='red',lwd=.75) 
        abplot <- abplot + geom_line(data=df[!is.na(df$runmed),], stat='hline', aes(yintercept=runmed, group=Chr), color='blue', linetype='longdash',lwd=.75)
    }
    abplot <- abplot + geom_rect(data=centromeres, aes(group=Chr, xmax=Start+1, xmin=0, ymin=ExpectedMedAB, ymax=0.5), fill="gold", alpha=.5, inherit.aes=FALSE) 
    abplot <- abplot + geom_rect(data=centromeres, aes(group=Chr, xmin=End-1, xmax=ChrEnd, ymin=ExpectedMedAB, ymax=0.5), fill="gold", alpha=.5, inherit.aes=FALSE)
    if (nrow(roh) > 0) {
        abplot <- abplot + geom_rect(data=roh, aes(xmin=Start, xmax=End, ymin=-0.025, ymax=-0.01, group=Chr), fill="limegreen", inherit.aes=FALSE)
    }
    abplot <- abplot + ylab("SNP Allele Balance [0, 0.5]")
    abplot <- abplot + ggtitle(samplename)
    
    plot.filename <- paste(substr(infile, 1, nchar(infile)-4), ".png", sep="")
    ggsave(plot.filename, abplot, units="in", width=12*1.618, height=12, type="cairo")
}

chrmedians.filename <- paste(substr(infile, 1, nchar(infile)-4), ".chrmedians.txt", sep="")
chr.AB.medians$Sample <- rep(samplename, nrow(chr.AB.medians))
write.table(chr.AB.medians, chrmedians.filename, quote=F, sep="\t",col.names=F,row.names=F)

eventstats.filename <- paste(substr(infile, 1, nchar(infile)-4), ".eventstats.txt", sep="")
df.tmp <- df.devs[df.devs$Area >= 0.01,]
df.tmp$Sample <- rep(samplename, nrow(df.tmp))
write.table(df.tmp, eventstats.filename, quote=F, sep="\t",col.names=F,row.names=F)
