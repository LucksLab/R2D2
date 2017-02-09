require(data.table)
library("gplots")

# Make a DG state plot from the DG state .dump file
#
# Usage: R < make_DG_state_plot_100x_figure.R --no-save <outfile> <comma-separated list of DG dump files> <title of plot> <smallest length to plot> <largest length to plot>
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
outfile = args[1]
dumpfile_list = args[2]
title <- args[3]
start <- as.numeric(args[4])
end <- as.numeric(args[5])

# Specify output to a pdf
pdf(outfile, width=7, height=3)

# Create different colors for each iteration of R2D2
colorp <- colorpanel(100, "#00ffff", "#ffbf00", "#7f0000")

dumpfiles = strsplit(dumpfile_list, ",")[[1]]
for(df in 1:length(dumpfiles)){
 print(dumpfiles[df])
 in.data <- read.table(dumpfiles[df], header=TRUE, sep="\t")

 # not necessary to call unique(), but in case of not unique dump file
 unique.data <- unique(in.data)
 if(length(args) > 3){
  unique.data <- unique.data[unique.data$nt <= end & unique.data$nt >= start,]
 }

 unique.data.states.dt <- as.data.table(unique.data)
 nts <- as.vector(unique(unique.data[,1]))

 if(df == 1){
  plot(unique.data[,c(1,2)], cex=0.25, main=title, ylab="", xlab="", cex.lab=1, cex.axis=1, cex.main=1, col="#FFFFFF00", xaxs="i", bty='l')
  title(xlab="RNA Length (nt)", ylab=expression(paste(Delta, " G (kcal/mol)", sep="")), line=2)
  max.dg <- unique.data.states.dt[,.SD[which.max(DG)],by=nt]
  min.dg <- unique.data.states.dt[,.SD[which.min(DG)],by=nt]
  #text(max.dg[,nt], max.dg[,DG] + 1, labels=max.dg[,nt], col="cadetblue")
  print(c(nts, rev(nts)))
  print(c(min.dg$DG, rev(max.dg$DG)))
  polygon(c(nts, rev(nts)), c(min.dg$DG, rev(max.dg$DG)), col="#838B8B96", border=NA)  # azure4
  lines(min.dg, col="black", lwd=1.5)
 }
 rm(unique.data.states.dt)

 for(it in 3:102){
  unique.data.pass <- unique.data[which(unique.data[,it] == 1),]
  nts <- unique(unique.data.pass[,1])

  #points(unique.data.pass[,1], unique.data.pass[,2], col=colorp[it-2], cex=0.5)
  for(i in 2:length(nts)){
   if(nts[i] - nts[i-1] == 1){
    for(prev in which(unique.data.pass[,1] == nts[i-1])){
     for(curr in which(unique.data.pass[,1] == nts[i])){
      lines(c(nts[i-1], nts[i]), c(unique.data.pass[prev,2], unique.data.pass[curr,2]), col=colorp[it-2], lwd=0.1)
     }
    }
   }
  }
 } 
}

dev.off()

