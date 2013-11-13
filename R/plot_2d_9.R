# plot info about fits from 9 sources
"model2d_9_plots" <- function(fitsum, L, source=9) {

	fsr <- fitsum
	corrs <- vector("list", 9)

	for (i in 1:9) { # get correlations
		corrs[[i]] <- fsr[grep(paste0("corrSigma_f\\[.*.,",i,",",source,"\\]"), rownames(fsr)),1]
	}

	# plot correlation vs time
	pdf(paste0("pdf/2d_9/corr_L",L,"_source",source,".pdf"))
		par(bty="l")
		plot(uf,corrs[[1]],
			xlab="Grid Cells per Cycle",ylab="Correlation",
			main=paste0("Correlation vs Grid Cells per Cycle, L=",L,"\n","DIC=",round(DIC,1),", pD=",round(pD,1)),
			type="l",lwd=2,col="blue",ylim=c(-0.5,1), xaxt='n',xlim=rev(range(uf)))   # GCM

			axis(1, at=f.samet[seq.samet], labels=p.samet[seq.samet])

			# RCMs
			for (i in 2:7) lines(uf,corrs[[i]],lwd=0.5,col="red")

			# data
			for (i in 8:9) lines(uf,corrs[[i]],lwd=1,col="black")

			abline(h=0, lty=3)

			if (L > 1) {
				points(knots, rep(-0.5, length(knots)), col="red", pch=4, cex=1.5)
			}

			legend("topleft",c("GCM","RCM","Data"),ncol=1,inset=0.05,col=c("blue","red","black"),lty=c(1,1,1))

	graphics.off()

	if (source == 8 || source == 9) {
		# plot conditional correlation of RCMs given GCM
		cond.corrs <- matrix(0, nrow=6, ncol=length(uf))

		for (csource in 2:7) {
			cond.corrs[csource-1,] <- unlist(sapply(1:length(uf), function(my_f) {
				cmat <- matrix(fsr[grep( paste0("corrSigma_f\\[",my_f,","), rownames(fsr)),1],nrow=9,ncol=9)

				# compute conditional corr mat
				cmat.cond <- cmat[c(csource,source),c(csource,source)] - tcrossprod(cmat[c(csource,source),1])

				cmat.cond[1,2]
			}))
		}

		print(cond.corrs[,1:10])

		pdf(paste0("pdf/2d_9/cond_corr_L",L,"_source",source,".pdf"))
			par(bty="l")
			plot(uf,cond.corrs[1,],
				xlab="Grid Cells per Cycle",ylab="Conditional Correlation",
				main=paste0("Conditional Correlation given GCM vs Grid Cells per Cycle\nL=",L,", DIC=",round(DIC,1),", pD=",round(pD,1)),
				type="l",lwd=0.5,col="red",ylim=c(-0.5,1), xaxt='n',xlim=rev(range(uf)))   # first RCM

				axis(1, at=f.samet[seq.samet], labels=p.samet[seq.samet])

				for (i in 2:6) {
					lines(uf,cond.corrs[i,],lwd=0.5,col="red")       # ith RCM
				}

				abline(h=0, lty=3)

				if (L > 1) {
					points(knots, rep(-0.5, length(knots)), col="red", pch=4, cex=1.5)
				}

				#legend("topleft",c("GCM given RCM","RCM given GCM"),ncol=1,inset=0.05,col=c("blue","red"),lty=c(1,1))

		graphics.off()
	}

}

#d.samet <- D[row(D)==col(D)]
#N.samet <- length(d.samet)
#seq.samet <- ( c(1, 20, 30, 40, 50, 60, 70, 80, 100) )

d.samet <- D[row(D)==col(D)]
f.samet <- F[row(F)==col(F)]
p.samet <- round(P[row(P)==col(P)],1)
N.samet <- length(f.samet)
seq.samet <- c(2, 15, 30, 40, 50, 60, 70, 80, 90, min(n1,n2))

if (TRUE) {
	# linear
	for (L in c(5,10,15,20)) {
		load(paste0("fitsums/fitsum_linL",L,".RData"))
		for (s in 8:9) model2d_9_plots(fitsum, L, s)
	}
}

if (FALSE) {
	# b-spline
	for (L in c(5,10,15,20)) {
		load(paste0("fitsums/fitsum_bsL",L,".RData"))
		for (s in 8:9) model2d_9_plots(fitsum, L, s)
	}
}
