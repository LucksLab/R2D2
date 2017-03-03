require(data.table)

# Plot SI Sample Size Determination
#
# Usage: R < make_SI_sample_size_determination.R --no-save <outfile> <Rep1> <Rep2> <Rep3> <EQ>
#
# Author: Angela M Yu, 2016-2017
# Version: 0.0.1
#
# Copyright (C) 2016, 2017  Julius B. Lucks, Angela M Yu.
# All rights reserved.
# Distributed under the terms of the GNU General Public License, see 'LICENSE'.

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
print(args)
outfile <- args[1]
rep1 <- t(read.table(args[2], header=FALSE)[,-1])
rep2 <- t(read.table(args[3], header=FALSE)[,-1])
rep3 <- t(read.table(args[4], header=FALSE)[,-1])
eq <- t(read.table(args[5], header=FALSE)[,-1])

# Sample size is actually 3x amount listed because we sample with 3 different methods
rep1[,1] <- rep1[,1]*3
rep2[,1] <- rep2[,1]*3
rep3[,1] <- rep3[,1]*3
eq[,1] <- eq[,1]*3

all.points <- rbind(rep1, rep2, rep3, eq)

# Specify output to a pdf
pdf(outfile, width=7)

colors <- c("red", "pink", "orange", "blue")
par(mfrow=c(1,2))
plot(all.points[,1], all.points[,2], col=rep(colors, each=8), xlab="Sample Size", ylab="Number Unique Structures")
lines(rep1[,c(1,2)], col=colors[1])
lines(rep2[,c(1,2)], col=colors[2])
lines(rep3[,c(1,2)], col=colors[3])
lines(eq[,c(1,2)], col=colors[4])
legend("topleft", legend=c("SRP Replicate 1", "SRP Replicate 2", "SRP Replicate 3", "SRP Equilibrium Refolded"), lty=c(1,1), col=colors, cex=0.6)

plot(all.points[,1], all.points[,3], xlab="Sample Size", ylab="Minimum Distance in Sample", col=rep(colors, each=8))
lines(rep1[,c(1,3)], col=colors[1])
lines(rep2[,c(1,3)], col=colors[2])
lines(rep3[,c(1,3)], col=colors[3])
lines(eq[,c(1,3)], col=colors[4])
legend("topright", legend=c("SRP Replicate 1", "SRP Replicate 2", "SRP Replicate 3", "SRP Equilibrium Refolded"), lty=c(1,1), col=colors, cex=0.6)


dev.off()

