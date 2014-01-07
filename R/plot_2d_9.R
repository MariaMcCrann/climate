# plot info about fits from 9 sources
"model2d_9_plots" <- function(fitsum, L, source=9) {

	fsr <- fitsum
	sds <- vector("list", 9)
	corrs <- vector("list", 9)

	max_sd <- 0
	min_sd <- 0

	for (i in 1:9) { # get correlations: median, 2.5%, 97.5%
		sds[[i]] <- fsr[grep(paste0("sdSigma_f\\[.*.,",i,"\\]"), rownames(fsr)),c(6,4,8)]
		max_sd <- ifelse(max_sd > max(sds[[i]][,1]), max_sd, max(sds[[i]][,1]))
		min_sd <- ifelse(min_sd > min(sds[[i]][,1]), min_sd, min(sds[[i]][,1]))
		corrs[[i]] <- fsr[grep(paste0("corrSigma_f\\[.*.,",i,",",source,"\\]"), rownames(fsr)),c(6,4,8)]
	}

	if (WHICH_CDAT == "ST") { m <- "Summer Temp"
	} else if (WHICH_CDAT == "WT") { m <- "Winter Temp"
	} else if (WHICH_CDAT == "SP") { m <- "Summer Precip"
	} else if (WHICH_CDAT == "WP") { m <- "Winter Precip"
	}

print(c(min_sd,max_sd))

	# plot SD over time
	pdf(paste0("pdf/2d_9/sd_L",L,"_",WHICH_CDAT,"_source",source,".pdf"), height=7/3, width=7)
		par(mar=c(2,2,3,1))
		par(bty="l")
		plot(uf,sds[[1]][,1],
			xlab="Grid Cells per Cycle",ylab="SD",
			main=paste0(m, " SD vs Grid Cells per Cycle, L=",L,"\n","DIC=",round(DIC,1),", pD=",round(pD,1)),
			type="l",lwd=2,col="blue",ylim=c(0,max_sd), xaxt='n',xlim=rev(range(uf)))   # BC

			axis(1, at=f.samet[seq.samet], labels=p.samet[seq.samet])
			for (i in 2:7) lines(uf,sds[[i]][,1],lwd=0.5,col="red")
			lines(uf,sds[[8]][,1],lwd=1,col="gray")
			lines(uf,sds[[9]][,1],lwd=1,col="black")
			abline(h=0, lty=3);
			if (L > 1) points(knots, rep(0, length(knots)), col="red", pch=4, cex=1.5)
			if (WHICH_CDAT == "ST")
				legend("top",c("UDEL","CRU","BC","RCM"),ncol=4,inset=0.01,col=c("gray","black","blue","red"),lty=c(1,1,1,1),lwd=1,cex=1,bty="n")
	graphics.off()

	# plot correlation vs time
	pdf(paste0("pdf/2d_9/corr_L",L,"_",WHICH_CDAT,"_source",source,".pdf"), height=7/3, width=7)
		par(mar=c(2,2,3,1))
		par(bty="l")
		plot(uf,corrs[[1]][,1],
			xlab="Grid Cells per Cycle",ylab="Correlation",
			main=paste0(m, " Correlation vs Grid Cells per Cycle, L=",L,"\n","DIC=",round(DIC,1),", pD=",round(pD,1)),
			type="l",lwd=2,col="blue",ylim=c(-0.5,1.25), xaxt='n',xlim=rev(range(uf)))   # BC

			#lines(uf,corrs[[1]][,2],lwd=1,col="blue",lty=2)
			#lines(uf,corrs[[1]][,3],lwd=1,col="blue",lty=2)

			axis(1, at=f.samet[seq.samet], labels=p.samet[seq.samet])

			# RCMs
			for (i in 2:7) lines(uf,corrs[[i]][,1],lwd=0.5,col="red")
			#for (i in 2:7) lines(uf,corrs[[i]][,2],lwd=0.15,col="red",lty=2)
			#for (i in 2:7) lines(uf,corrs[[i]][,3],lwd=0.15,col="red",lty=2)

			# data
			if (source == 8) lines(uf,corrs[[9]][,1],lwd=1,col="black")
			if (source == 9) lines(uf,corrs[[8]][,1],lwd=1,col="black")
			#for (i in 8:9) lines(uf,corrs[[i]][,2],lwd=0.25,col="black",lty=2)
			#for (i in 8:9) lines(uf,corrs[[i]][,3],lwd=0.25,col="black",lty=2)

			abline(h=1.0, lty=3); abline(h=0.5, lty=3); abline(h=0.0, lty=3)

			if (L > 1) points(knots, rep(-0.4, length(knots)), col="red", pch=4, cex=1.5)

			#legend("topleft",c("BC","RCM","Data"),ncol=1,inset=0.05,col=c("blue","red","black"),lty=c(1,1,1))
			if (WHICH_CDAT == "ST") legend("top",c("UDEL","BC","RCM"),ncol=3,inset=-0.10,col=c("black","blue","red"),lty=c(1,1,1),lwd=1,cex=1,bty="n")

	graphics.off()

	if (source == 8 || source == 9) {
		# plot conditional correlation of RCMs given BC
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

		pdf(paste0("pdf/2d_9/cond_corr_L",L,"_",WHICH_CDAT,"_source",source,".pdf"), height=7/3, width=7)
			par(mar=c(2,2,3,1))
			par(bty="l")
			plot(uf,cond.corrs[1,],
				xlab="Grid Cells per Cycle",ylab="Conditional Correlation",
				main=paste0(m, " Conditional Correlation given BC vs Grid Cells per Cycle\nL=",L,", DIC=",round(DIC,1),", pD=",round(pD,1)),
				type="l",lwd=0.5,col="red",ylim=c(-0.5,1.25), xaxt='n',xlim=rev(range(uf)))   # first RCM

				axis(1, at=f.samet[seq.samet], labels=p.samet[seq.samet])

				for (i in 2:6) {
					lines(uf,cond.corrs[i,],lwd=0.5,col="red")       # ith RCM
				}

				abline(h=1.0, lty=3); abline(h=0.5, lty=3); abline(h=0.0, lty=3)

				if (L > 1) {
					points(knots, rep(-0.4, length(knots)), col="red", pch=4, cex=1.5)
				}

				#legend("topleft",c("BC given RCM","RCM given BC"),ncol=1,inset=0.05,col=c("blue","red"),lty=c(1,1))
				#if (WHICH_CDAT == "ST") legend("top",c("UDEL","BC","RCM"),ncol=3,inset=-0.10,col=c("black","blue","red"),lty=c(1,1,1),lwd=1,cex=1,bty="n")

		graphics.off()
	}

}

if (!("WHICH_CDAT" %in% ls())) {
	stop("Missing data type")
}

if (FALSE) {
	# linear
	for (L in c(5,10,15,20)) {
		load(paste0("fitsums/fitsum_linL",L,".RData"))
		for (s in 8:9) model2d_9_plots(fitsum, L, s)
	}
}

if (FALSE) {
	# b-spline
	for (L in c(5,10)) { #,15,20)) {
		load(paste0("fitsums/fitsum_bsL",L,"_",WHICH_CDAT,".RData"))
		for (s in 8:9) model2d_9_plots(fitsum, L, s)
	}
}

# ST = 10, WT = 10, SP = 15, WP = 15
WHICH_CDAT <- "ST"; load("fitsums/fitsum_bsL10_ST.RData"); for (s in 8:9) model2d_9_plots(fitsum, L, s)
WHICH_CDAT <- "WT"; load("fitsums/fitsum_bsL10_WT.RData"); for (s in 8:9) model2d_9_plots(fitsum, L, s)
WHICH_CDAT <- "SP"; load("fitsums/fitsum_bsL10_SP.RData"); for (s in 8:9) model2d_9_plots(fitsum, L, s)
WHICH_CDAT <- "WP"; load("fitsums/fitsum_bsL10_WP.RData"); for (s in 8:9) model2d_9_plots(fitsum, L, s)
