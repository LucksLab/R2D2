require(data.table)

# Make a DG state plot from the DG state .dump file
#
# Usage: R < make_DG_state_plot_100x_EQ_Rep_figure.R --no-save <outfile> <Equilibrium refolded DG dump file with 100 repititions> <Cotranscriptional DG dump file with 100 repititions> <DG dump file from KineFold> <title of plot> <smallest length to plot> <largest length to plot>
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
eq.dumpfile <- args[2]
co.dumpfiles <- strsplit(args[3], ",")[[1]]
kf.dumpfile <- args[4]
title <- args[5]
start <- as.numeric(args[6])
end <- as.numeric(args[7])

# Specify output to a pdf
pdf(outfile, width=7, height=3)

eq.col <- c("#838B8B5A", "#1DADADAA", "#838B8BAA")
co.col <- c("#8B3E2F5A", "#FF3030AA", "#8B3E2FAA")

eq.data <- read.table(eq.dumpfile, header=TRUE, sep="\t")
co.data <- read.table(co.dumpfiles[1], header=TRUE, sep="\t")
for(cdf in co.dumpfiles){
 co.data.next <- read.table(cdf, header=TRUE, sep="\t")
 co.data <- merge(co.data, co.data.next, by=c("nt", 'DG'), all=TRUE)
 co.data[is.na(co.data)] <- 0
}
kf.data <- read.table(kf.dumpfile, header=TRUE, sep="\t")

# not necessary to call unique(), but in case of not unique dump file
unique.data.eq <- unique(eq.data)
unique.data.co <- unique(co.data)
unique.data.co$nt <- unique.data.co$nt - 14

# reduce to only plot between start and end inclusive
if(length(args) > 5){
 unique.data.co <- unique.data.co[unique.data.co$nt <= end & unique.data.co$nt >= start,]
 unique.data.eq <- unique.data.eq[unique.data.eq$nt <= end & unique.data.eq$nt >= start,]
}

# set data tables and unique.points
unique.points.co <- as.data.table(unique.data.co[,c(1,2)])
unique.points.eq <- as.data.table(unique.data.eq[,c(1,2)])
unique.points <- rbind(unique.points.co, unique.points.eq)

# plot all points
plot(unique.points, cex=0.25, main=title, ylab="", xlab="", cex.lab=1, cex.axis=1, cex.main=1, col="#FFFFFF00", xaxs="i", bty='l')
title(xlab="RNA Length (nt)", ylab=expression(paste(Delta, " G (kcal/mol)", sep="")), line=2)

# shade range of DG's
# Cotranscriptional data
nts <- as.vector(unique(unique.points.co$nt))
max.dg <- unique.points.co[,.SD[which.max(DG)],by=nt]
min.dg.co <- unique.points.co[,.SD[which.min(DG)],by=nt]
polygon(c(nts, rev(nts)), c(min.dg.co$DG, rev(max.dg$DG)), col=co.col[1], border=co.col[3], lwd=0.1)
# Equilibrium refolded data
nts <- as.vector(unique(unique.points.eq$nt))
max.dg <- unique.points.eq[,.SD[which.max(DG)],by=nt]
min.dg.eq <- unique.points.eq[,.SD[which.min(DG)],by=nt]
polygon(c(nts, rev(nts)), c(min.dg.eq$DG, rev(max.dg$DG)), col=eq.col[1], border=eq.col[3], lwd=0.1)

# plot MFE line
min.dg <- rbind(min.dg.co, min.dg.eq)[,.SD[which.min(DG)],by=nt]
lines(min.dg, col="black", lwd=0.6)

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

# plot R2D2 equilibrium refolded pathways
for(it in 3:ncol(unique.data.eq)){
 unique.data.eq.it <- unique.data.eq[which(unique.data.eq[,it] == 1),]
 nts <- unique(unique.data.eq[,1])

 for(i in 2:length(nts)){
  if(nts[i] - nts[i-1] == 1){
   for(prev in which(unique.data.eq.it[,1] == nts[i-1])){
    for(curr in which(unique.data.eq.it[,1] == nts[i])){
     lines(c(nts[i-1], nts[i]), c(unique.data.eq.it[prev,2], unique.data.eq.it[curr,2]), col=eq.col[2], lwd=0.1)
    }
   }
  }
 }
}

# plot KineFold folding pathway
lines(kf.data$nt, kf.data$DG, col="olivedrab", lwd=0.1)
lines(kf.data$nt, kf.data$KineFold_DG, col="olivedrab1", lwd=0.1)


dev.off()

