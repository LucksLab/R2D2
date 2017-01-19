require(data.table)
library("gplots")

# Make a DG state plot from the DG state .dump file
#
# Usage: R < make_DG_state_plot_100x.R --no-save <outfile> <comma-separated list of DG dump files> <smallest length to plot> <largest length to plot>
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
start <- as.numeric(args[3])
end <- as.numeric(args[4])

# Specify output to a pdf
pdf(outfile, height=24, width=48)

# Create different colors for each iteration of R2D2
colorp <- colorpanel(100, "#00ffff", "#ffbf00", "#7f0000")

dumpfiles = strsplit(dumpfile_list, ",")[[1]]
for(df in 1:length(dumpfiles)){
 print(dumpfiles[df])
 in.data <- read.table(dumpfiles[df], header=TRUE, sep="\t")

 # not necessary to call unique(), but in case of not unique dump file
 unique.data <- unique(in.data)
 if(length(args) > 2){
  unique.data <- unique.data[unique.data$nt <= end & unique.data$nt >= start,]
 }

 unique.data.states.dt <- as.data.table(unique.data)

 if(df == 1){
  plot(unique.data[,c(1,2)], cex=0.5, main="States Plot", ylab="Delta G", xlab="nt length", cex.lab=1.5, cex.axis=1.5, cex.main=2)
  max.dg <- unique.data.states.dt[,.SD[which.max(DG)],by=nt]
  min.dg <- unique.data.states.dt[,.SD[which.min(DG)],by=nt]
  text(max.dg[,nt], max.dg[,DG] + 1, labels=max.dg[,nt], col="cadetblue")
 }
 rm(unique.data.states.dt)

 colorpanel
 for(it in 3:102){
  unique.data.pass <- unique.data[which(unique.data[,it] == 1),]
  nts <- unique(unique.data.pass[,1])

  points(unique.data.pass[,1], unique.data.pass[,2], col=colorp[it-2])
  for(i in 2:length(nts)){
   if(nts[i] - nts[i-1] == 1){
    for(prev in which(unique.data.pass[,1] == nts[i-1])){
     for(curr in which(unique.data.pass[,1] == nts[i])){
      lines(c(nts[i-1], nts[i]), c(unique.data.pass[prev,2], unique.data.pass[curr,2]), col=colorp[it-2])
     }
    }
   }
  }
 } 
}

dev.off()

