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
method.colors <- c("green", "purple", "orange", "black")
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


par(mfrow=c(2,1))
all_methods <- c(rep1[,4], rep1[,5], rep1[,6], rep1[,2])
plot(rep(rep1[,1], 4), all_methods, col=rep(method.colors, each=8), xlab="Sample Size", ylab="Number Unique Structures", main="Replicate 1")
lines(rep1[,c(1,4)], col=method.colors[1])
lines(rep1[,c(1,5)], col=method.colors[2])
lines(rep1[,c(1,6)], col=method.colors[3])
lines(rep1[,c(1,2)], col=method.colors[4])
legend("topleft", legend=c("No SHAPE", "SHAPE", "Hard Constrained", "All"), lty=c(1,1), col=method.colors, cex=0.6)

all_methods <- c(rep2[,4], rep2[,5], rep2[,6], rep2[,2])
plot(rep(rep2[,1], 4), all_methods, col=rep(method.colors, each=8), xlab="Sample Size", ylab="Number Unique Structures", main="Replicate 2")
lines(rep2[,c(1,4)], col=method.colors[1])
lines(rep2[,c(1,5)], col=method.colors[2])
lines(rep2[,c(1,6)], col=method.colors[3])
lines(rep2[,c(1,2)], col=method.colors[4])
legend("topleft", legend=c("No SHAPE", "SHAPE", "Hard Constrained", "All"), lty=c(1,1), col=method.colors, cex=0.6)

all_methods <- c(rep3[,4], rep3[,5], rep3[,6], rep3[,2])
plot(rep(rep3[,1], 4), all_methods, col=rep(method.colors, each=8), xlab="Sample Size", ylab="Number Unique Structures", main="Replicate 3")
lines(rep3[,c(1,4)], col=method.colors[1])
lines(rep3[,c(1,5)], col=method.colors[2])
lines(rep3[,c(1,6)], col=method.colors[3])
lines(rep3[,c(1,2)], col=method.colors[4])
legend("topleft", legend=c("No SHAPE", "SHAPE", "Hard Constrained", "All"), lty=c(1,1), col=method.colors, cex=0.6)

all_methods <- c(eq[,4], eq[,5], eq[,6], eq[,2])
plot(rep(eq[,1], 4), all_methods, col=rep(method.colors, each=8), xlab="Sample Size", ylab="Number Unique Structures", main="Equilibrium Refolded")
lines(eq[,c(1,4)], col=method.colors[1])
lines(eq[,c(1,5)], col=method.colors[2])
lines(eq[,c(1,6)], col=method.colors[3])
lines(eq[,c(1,2)], col=method.colors[4])
legend("topleft", legend=c("No SHAPE", "SHAPE", "Hard Constrained", "All"), lty=c(1,1), col=method.colors, cex=0.6)



dev.off()

