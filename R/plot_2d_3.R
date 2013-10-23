# plot info about fits from 3 sources
"model2d_3_plots" <- function(fitsum, L) {

	fsr <- fitsum
	corr_1_3 <- fsr[grep("corrSigma_d\\[.*.,1,3\\]", rownames(fsr)),]
	corr_2_3 <- fsr[grep("corrSigma_d\\[.*.,2,3\\]", rownames(fsr)),]

	# plot correlation vs time
	pdf(paste0("pdf/2d_3/corr_L",L,".pdf"))
		par(bty="l")
		plot(ud,corr_1_3[,1],
			xlab="Grid Cells Apart",ylab="Correlation",
			main=paste0("Correlation vs Grid Cell Distance, L=",L,"\n","DIC: ",round(DIC,1),", pD=",round(pD,1)),
			type="l",lwd=2,col="blue",ylim=c(-0.5,1), xaxt='n',xlim=rev(range(ud)))   # GCM

			axis(1, at=d.samet[seq.samet], labels=rev(seq.samet))

			lines(ud,corr_2_3[,1],lwd=0.5,col="red")       # RCM

			# 95% intervals

			# GCM
			lines(ud,corr_1_3[,4],lwd=0.5,col="blue",lty=2)
			lines(ud,corr_1_3[,8],lwd=0.5,col="blue",lty=2)

			# RCM
			lines(ud,corr_2_3[,4],lwd=0.5,col="red",lty=2)
			lines(ud,corr_2_3[,8],lwd=0.5,col="red",lty=2)

			abline(h=0, lty=3)

			if (L > 1) {
				points(knots, rep(-0.5, length(knots)), col="red", pch=4, cex=1.5)
			}

			legend("topleft",c("GCM","RCM"),ncol=1,inset=0.05,col=c("blue","red"),lty=c(1,1))

	graphics.off()

	# plot conditional correlation
	cond.corrs <- matrix(0, nrow=2, ncol=length(ud))
	for (source in c(1,2)) {

		if (source == 1) {
			o_source <- 2
		} else {
			o_source <- 1
		}

		cond.corrs[source,] <- unlist(sapply(1:length(ud), function(my_d) {
			cmat <- matrix(fsr[grep( paste0("corrSigma_d\\[",my_d,","), rownames(fsr)),1],nrow=3,ncol=3)

			# compute conditional corr mat
			cmat.cond <- cmat[c(source,3),c(source,3)] - tcrossprod(cmat[c(source,3),o_source])

			cmat.cond[1,2]
		}))
	}
	print(cond.corrs[,1:10])

	pdf(paste0("pdf/2d_3/cond_corr_L",L,".pdf"))
		par(bty="l")
		plot(ud,cond.corrs[1,],
			xlab="Grid Cells Apart",ylab="Conditional Correlation",
			main=paste0("Conditional Corr vs Grid Cell Distance, L=",L,"\n","DIC: ",round(DIC,1),", pD=",round(pD,1)),
			type="l",lwd=2,col="blue",ylim=c(-0.5,1), xaxt='n',xlim=rev(range(ud)))   # GCM

			axis(1, at=d.samet[seq.samet], labels=rev(seq.samet))

			lines(ud,cond.corrs[2,],lwd=0.5,col="red")       # RCM

			abline(h=0, lty=3)

			if (L > 1) {
				points(knots, rep(-0.5, length(knots)), col="red", pch=4, cex=1.5)
			}

			legend("topleft",c("GCM given RCM","RCM given GCM"),ncol=1,inset=0.05,col=c("blue","red"),lty=c(1,1))

	graphics.off()

}

d.samet <- D[row(D)==col(D)]
N.samet <- length(d.samet)
seq.samet <- ( c(1, 20, 30, 40, 50, 60, 70, 80, 100) )

if (TRUE) {
	# linear
	for (L in c(5,10,15,25)) {
		load(paste0("fitsums/fitsum_linL",L,".RData"))
		model2d_3_plots(fitsum, L)
	}
}

if (FALSE) {
	# b-spline
	for (L in c(5,10,15,25,30)) {
		load(paste0("fitsums/fitsum_bsL",L,".RData"))
		model2d_3_plots(fitsum, L)
	}
}

if (FALSE) {
	for (L in c(4,7,11,14,17,20,24,27,30,35)) { #,20,24,27,30)) {
		load(paste0("fitsums/fitsumL",L,".RData"))
		model2d_3_plots(fitsum, L)
	}
}

done
load("fitsums/fitsumL1.RData"); model2d_3_plots(fitsum, 1)
load("fitsums/fitsumL2.RData"); model2d_3_plots(fitsum, 2)
load("fitsums/fitsumL3.RData"); model2d_3_plots(fitsum, 3)
load("fitsums/fitsumL4.RData"); model2d_3_plots(fitsum, 4)
load("fitsums/fitsumL6.RData"); model2d_3_plots(fitsum, 6)
load("fitsums/fitsumL8.RData"); model2d_3_plots(fitsum, 8)
load("fitsums/fitsumL10.RData"); model2d_3_plots(fitsum, 10)
#load("fitsums/fitsumL11.RData"); model2d_3_plots(fitsum, 11)
#load("fitsums/fitsumL12.RData"); model2d_3_plots(fitsum, 12)
#load("fitsums/fitsumL13.RData"); model2d_3_plots(fitsum, 13)
#load("fitsums/fitsumL14.RData"); model2d_3_plots(fitsum, 14)
#load("fitsums/fitsumL15.RData"); model2d_3_plots(fitsum, 15)
load("fitsums/fitsumL16.RData"); model2d_3_plots(fitsum, 16)
load("fitsums/fitsumL17.RData"); model2d_3_plots(fitsum, 17)
load("fitsums/fitsumL18.RData"); model2d_3_plots(fitsum, 18)
load("fitsums/fitsumL19.RData"); model2d_3_plots(fitsum, 19)
load("fitsums/fitsumL20.RData"); model2d_3_plots(fitsum, 20)
