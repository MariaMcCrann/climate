# plot info about fits
"model2d_plots" <- function(n, fitsum, L, data.index=9) {
	# plot correlation vs frequency
	pdf(paste0("pdf/2d/corr_L",L,"_D",data.index,".pdf"))
		par(bty="l")
		plot(1:n,fitsum$summary[grep(paste0("corrSigma_d\\[.*.,1,",data.index,"]"), rownames(fitsum$summary)),1],
		     xlab="Frequency (T)",ylab="Correlation",main=paste0("Correlation vs Frequency, L=",L),
		     type="l",lwd=2,col="blue",ylim=c(-1,1))   # GCM

#		for (i in 2:7)
		i <- 2
			lines(1:n,fitsum$summary[grep(paste0("corrSigma_d\\[.*.,",i,",",data.index,"]"), rownames(fitsum$summary)),1],lwd=0.5,col="red")   # RCMs
#		for (i in 8:9)
#			if (i != data.index)
		i <- 3
				lines(1:n,fitsum$summary[grep(paste0("corrSigma_d\\[.*.,",i,",",data.index,"]"), rownames(fitsum$summary)),1],lwd=1,col="black") # data
	graphics.off()
}
