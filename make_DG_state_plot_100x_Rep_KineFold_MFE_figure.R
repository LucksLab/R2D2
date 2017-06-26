require(data.table)
library(Hmisc)

# Make a DG state plot from the DG state .dump file
#
# Usage: R < make_DG_state_plot_Rep_KineFold_MFE_figure.R --no-save <outfile> <Cotranscriptional DG dump file with 100 repititions> <Directory containing DG dump files from KineFold> <title of plot> <smallest length to plot> <largest length to plot>
#
# Author: Angela M Yu, 2017
# Version: 0.0.1
#
# Copyright (C) 2017  Julius B. Lucks, Angela M Yu.
# All rights reserved.
# Distributed under the terms of the GNU General Public License, see 'LICENSE'.

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
outfile <- args[1]
co.dumpfiles <- strsplit(args[2], ",")[[1]]
kf.dumpfile.dir <- args[3]
EQ.MFE.dumpfile.dir <- args[4]
MFE.dumpfile.dir <- args[5]
title <- args[6]
start <- as.numeric(args[7])
end <- as.numeric(args[8])

# Specify output to a pdf
pdf(outfile, width=7, height=3)
#par(mar=c(3, 3, 1.5, 0.25)+0.1)

back.col <- c("#7F7F7F2A", "#7F7F7F3A") #gray50 w/ transparancy -> gray 234/255
co.col <- c("#4E2A842A", "#4E2A84AA", "#4E2A843A")  #ED1C25

co.data <- read.table(co.dumpfiles[1], header=TRUE, sep="\t")
co.data <- co.data[order(co.data$nt),]
for(cdf in co.dumpfiles[-1]){
 co.data.next <- read.table(cdf, header=TRUE, sep="\t")
 co.data <- merge(co.data, co.data.next, by=c("nt", 'DG'), all=TRUE)
 co.data[is.na(co.data)] <- 0
}

# not necessary to call unique(), but in case of not unique dump file
unique.data.co <- unique(co.data)
unique.data.co$nt <- unique.data.co$nt - 14

# reduce to only plot between start and end inclusive
if(length(args) > 6){
 unique.data.co <- unique.data.co[unique.data.co$nt <= end & unique.data.co$nt >= start,]
}
unique.data.co <- unique.data.co[order(unique.data.co$nt),]

# set data tables and unique.points
unique.points.co <- as.data.table(unique.data.co[,c(1,2)])

# plot all points
plot(unique.points.co, cex=0.25, main=title, ylab="", xlab="", cex.lab=1, cex.axis=1, cex.main=1, col="#FFFFFF00", xaxs="i", bty='l', xlim=c(min(unique.points.co$nt - 0.2), max(unique.points.co$nt) + 0.2))
minor.tick(nx=2, ny=2, tick.ratio=0.5)
title(xlab="RNA Length (nt)", ylab=expression(paste(Delta, " G (kcal/mol)", sep="")), line=2)

# shade range of DG's
# Cotranscriptional data
nts <- as.vector(unique(unique.points.co$nt))
max.dg <- unique.points.co[,.SD[which.max(DG)],by=nt]
min.dg.co <- unique.points.co[,.SD[which.min(DG)],by=nt]
polygon(c(nts, rev(nts)), c(min.dg.co$DG, rev(max.dg$DG)), col=back.col[1], border=back.col[2], lwd=0.1)

# plot R2D2 cotranscriptional pathways
for(it in 3:ncol(unique.data.co)){
 unique.data.co.it <- unique.data.co[which(unique.data.co[,it] == 1),]
 nts <- unique(unique.data.co.it[,1])

 for(i in 2:length(nts)){
  if(nts[i] - nts[i-1] == 1){
   for(prev in which(unique.data.co.it[,1] == nts[i-1])){
    for(curr in which(unique.data.co.it[,1] == nts[i])){
     lines(c(nts[i-1], nts[i]), c(unique.data.co.it[prev,2], unique.data.co.it[curr,2]), col=co.col[2], lwd=0.1)
    }
   }
  }
 }
} 

# plot KineFold folding pathway
kf.dumpfiles <- list.files(kf.dumpfile.dir, pattern="\\.dump$", full.names=TRUE)
for(kf.dumpfile in kf.dumpfiles){
 kf.data <- read.table(kf.dumpfile, header=TRUE, sep="\t")

 lines(kf.data$nt, kf.data$DG, col="red", lwd=0.1)
}

# plot EQ SHAPE directed MFE line
EQ.MFE.dumpfiles <- list.files(EQ.MFE.dumpfile.dir, pattern="\\.dump$", full.names=TRUE)
for(EQ.MFE.dumpfile in EQ.MFE.dumpfiles){
 EQ.MFE.data <- read.table(EQ.MFE.dumpfile, header=TRUE, sep="\t")
 #lines(EQ.MFE.data$nt, EQ.MFE.data$DG, col="navy", lwd=0.6)
}

# plot SHAPE directed MFE line(s)
MFE.dumpfiles <- list.files(MFE.dumpfile.dir, pattern="\\.dump$", full.names=TRUE)
for(MFE.dumpfile in MFE.dumpfiles){
 MFE.data <- read.table(MFE.dumpfile, header=TRUE, sep="\t")
 lines(MFE.data$nt, MFE.data$DG, col="darkslategray", lwd=0.6)
}


dev.off()

