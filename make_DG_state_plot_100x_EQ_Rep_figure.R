require(data.table)
library(Hmisc)
#library("shape")

# Make a DG state plot from the DG state .dump file
#
# Usage: R < make_DG_state_plot_100x_EQ_Rep_figure.R --no-save <outfile> <Equilibrium refolded DG dump file with 100 repititions> <Cotranscriptional DG dump file with 100 repititions> <title of plot> <smallest length to plot> <largest length to plot> <Vertical lines x positions>
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
title <- args[4]
start <- as.numeric(args[5])
end <- as.numeric(args[6])
lines.x <- as.numeric(strsplit(args[7], ",")[[1]])

# Specify output to a pdf
pdf(outfile, width=7.007874, height=3)

#rgb(t(col2rgb("turquoise3")), maxColorValue=255)
back.col <- c("#7F7F7F2A", "#7F7F7F3A") #gray50 w/ transparancy -> gray 234/255
co.col <- c("#4E2A842A", "#4E2A84AA", "#4E2A843A")  #ED1C25
eq.col <- c("#00CED12A", "#00CED1AA", "#00CED13A") #darkturquiose

eq.data <- read.table(eq.dumpfile, header=TRUE, sep="\t")
co.data <- read.table(co.dumpfiles[1], header=TRUE, sep="\t")
co.data <- co.data[order(co.data$nt),]
for(cdf in co.dumpfiles[-1]){
 co.data.next <- read.table(cdf, header=TRUE, sep="\t")
 co.data <- merge(co.data, co.data.next, by=c("nt", 'DG'), all=TRUE)
 co.data[is.na(co.data)] <- 0
}

# not necessary to call unique(), but in case of not unique dump file
unique.data.eq <- unique(eq.data)
unique.data.co <- unique(co.data)
unique.data.co$nt <- unique.data.co$nt - 14
# reduce to only plot between start and end inclusive
if(length(args) > 4){
 unique.data.co <- unique.data.co[unique.data.co$nt <= end & unique.data.co$nt >= start,]
 unique.data.eq <- unique.data.eq[unique.data.eq$nt <= end & unique.data.eq$nt >= start,]
}

# set data tables and unique.points
unique.points.co <- as.data.table(unique.data.co[,c(1,2)])
unique.points.eq <- as.data.table(unique.data.eq[,c(1,2)])
unique.points <- rbind(unique.points.co, unique.points.eq)
head(unique.points)

# plot all points
plot(unique.points, cex=0.25, main=title, ylab="", xlab="", cex.lab=1, cex.axis=1, cex.main=1, col="#FFFFFF00", xaxs="i", bty='l', xlim=c(min(unique.points$nt - 0.2), max(unique.points$nt) + 0.2))
minor.tick(nx=2, ny=2, tick.ratio=0.5)
title(xlab="RNA Length (nt)", ylab=expression(paste(Delta, " G (kcal/mol)", sep="")), line=2)

# shade range of DG's
# Equilibrium refolded data
nts <- as.vector(unique(unique.points.eq$nt))
max.dg.eq <- unique.points.eq[,.SD[which.max(DG)],by=nt]
min.dg.eq <- unique.points.eq[,.SD[which.min(DG)],by=nt]
# Cotranscriptional data
nts <- as.vector(unique(unique.points.co$nt))
max.dg.co <- unique.points.co[,.SD[which.max(DG)],by=nt]
min.dg.co <- unique.points.co[,.SD[which.min(DG)],by=nt]
# Merge both sets and shade
min.dg <- rbind(min.dg.co, min.dg.eq)[,.SD[which.min(DG)],by=nt]
max.dg <- rbind(max.dg.co, max.dg.eq)[,.SD[which.max(DG)],by=nt]
polygon(c(nts, rev(nts)), c(min.dg$DG, rev(max.dg$DG)), col=back.col[1], border=back.col[2], lwd=0.1)

# plot MFE line
lines(min.dg, col="black", lwd=0.6)

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

# place vertical lines at specified lengths
for(l in lines.x){
 abline(v=l, lty=2)
}
dev.off()
