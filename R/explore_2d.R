library(Cairo)

if (!("cdata" %in% ls())) {
	source("R/load_data2.R")
}

i.data <- 9
i.bc <- 1
i.rcm <- 5

# plot original data
if (FALSE) {
#r <- range(cdata[,,c(i.data,i.bc,i.rcm)])
r <- range(cdata[,,1:9])
#pdf("pdf/2d/data.pdf",height=2)
Cairo("pdf/2d/data.png",type="png",pointsize=20,width=1024,height=900)
  par(bty="l")
  par(mfrow=c(3,3))
	par(mar=c(2,2,3,5))

  image.plot(tlon,tlat,cdata[,,8],xlab="",ylab="",axes=F,main="UDEL Data",zlim=r); map("world",add=T)
  image.plot(tlon,tlat,cdata[,,9],xlab="",ylab="",axes=F,main="CRU Data",zlim=r); map("world",add=T)
  image.plot(tlon,tlat,cdata[,,1],xlab="",ylab="",axes=F,main="NCEP BC",zlim=r); map("world",add=T)

  image.plot(tlon,tlat,cdata[,,2],xlab="",ylab="",axes=F,main="CRCM RCM",zlim=r);   map("world",add=T)
  image.plot(tlon,tlat,cdata[,,3],xlab="",ylab="",axes=F,main="ECP2 RCM",zlim=r);   map("world",add=T)
  image.plot(tlon,tlat,cdata[,,4],xlab="",ylab="",axes=F,main="HRM3 RCM",zlim=r);   map("world",add=T)

  image.plot(tlon,tlat,cdata[,,5],xlab="",ylab="",axes=F,main="MM5I RCM",zlim=r);   map("world",add=T)
  image.plot(tlon,tlat,cdata[,,6],xlab="",ylab="",axes=F,main="RCM3 RCM",zlim=r);   map("world",add=T)
  image.plot(tlon,tlat,cdata[,,7],xlab="",ylab="",axes=F,main="WRFG RCM",zlim=r);   map("world",add=T)
graphics.off()}

# plot Zs
if (FALSE) {pdf("pdf/2d/data_trans.pdf")
  par(bty="l")
  par(mfrow=c(1,3))
  image.plot(tlon,tlat,Z[,,i.data],xlab="",ylab="",axes=F,main="Data"); map("world",add=T)
  image.plot(tlon,tlat,Z[,,i.bc],xlab="",ylab="",axes=F,main="BC");   map("world",add=T)
  image.plot(tlon,tlat,Z[,,i.rcm],xlab="",ylab="",axes=F,main="RCM");   map("world",add=T)
graphics.off()}

# smooth correlations
if (FALSE) {
if (TRUE) {
	#knots <- seq(min(d), max(d), len=80)
	inc <- 0.10
	knots <- seq(min(d), max(d), len=1000)
	cor.dat <- sapply(1:length(knots), function(k) {
		cor(z[d >= (knots[k]-inc)&d <= (knots[k]+inc),c(i.data,i.bc,i.rcm)])[1,2:3]
	})
	pdf("pdf/2d/data_z_corr.pdf", height=7/2)
		plot(knots, cor.dat[1,], col="blue", type="l", ylab="Correlation", xlab="Grid Cells per Cycle", xlim=rev(range(knots)))#, xaxt="n")
		#axis(1, at=d.samet[seq.samet], labels=period.samet)
		points(knots, cor.dat[2,], col="red", type="l")
		abline(h=0, lty=3)
		legend("topleft",c("BC","RCM"),ncol=1,inset=0.05,col=c("blue","red"),lty=c(1,1))
	graphics.off()
}

if (FALSE) {
	inc <- 0.05
	knots <- seq(min(F), max(F), len=1000)
	cor.dat <- sapply(1:length(knots), function(k) {
		cor(z[F >= (knots[k]-inc)&F <= (knots[k]+inc),c(i.data,i.bc,i.rcm)])[1,2:3]
	})
	pdf("pdf/2d/data_z_corr.pdf", height=7/2)
		plot(knots, cor.dat[1,], col="blue", type="l", ylab="Correlation", xlab="Grid Cells per Cycle", xlim=rev(range(knots)), xaxt="n")
		axis(1, at=d.samet[seq.samet], labels=period.samet)
		points(knots, cor.dat[2,], col="red", type="l")
		abline(h=0, lty=3)
		legend("topleft",c("BC","RCM"),ncol=1,inset=0.05,col=c("blue","red"),lty=c(1,1))
	graphics.off()
}

}

# smooth correlations for all models
if (TRUE) {
	#knots <- seq(min(d), max(d), len=80)
	inc <- 0.005
	knots <- seq(min(f), max(f), len=(max(f)-min(f))/inc)
	cor.dat8 <- sapply(1:length(knots), function(k) { cor(z[f >= (knots[k]-inc)&f <= (knots[k]+inc),])[8,1:9] })
	cor.dat9 <- sapply(1:length(knots), function(k) { cor(z[f >= (knots[k]-inc)&f <= (knots[k]+inc),])[9,1:9] })

	# data source 8
	pdf("pdf/2d/data_z_corr_d8.pdf", height=7/2)
		# BC
		plot(knots, cor.dat8[1,], col="blue", type="l", main="Empirical Marginal Correlations for UDEL Data",
		     ylab="Correlation", xlab="Grid Cells per Cycle", xaxt="n", xlim=rev(range(knots)), ylim=c(-.3,1.25))
			axis(1, at=f.samet[seq.samet], labels=p.samet[seq.samet])

		# RCMs
		for (i in 2:7)
			points(knots, cor.dat8[i,], col="red", type="l")

		# other data
		points(knots, cor.dat8[9,], col="black", type="l")

		abline(h=0, lty=3)
		legend("top",c("BC","RCM","CRU"),ncol=3,inset=0.02,col=c("blue","red","black"),lty=c(1,1,1),lwd=2)
	graphics.off()

	# data source 9
	pdf("pdf/2d/data_z_corr_d9.pdf", height=7/2)
		# BC
		plot(knots, cor.dat9[1,], col="blue", type="l", main="Empirical Marginal Correlations for CRU Data",
		     ylab="Correlation", xlab="Grid Cells per Cycle", xaxt="n", xlim=rev(range(knots)), ylim=c(-.3,1.25))
			axis(1, at=f.samet[seq.samet], labels=p.samet[seq.samet])

		# RCMs
		for (i in 2:7)
			points(knots, cor.dat9[i,], col="red", type="l")

		# other data
		points(knots, cor.dat9[8,], col="black", type="l")

		abline(h=0, lty=3)
		legend("top",c("BC","RCM","UDEL"),ncol=3,inset=0.02,col=c("blue","red","black"),lty=c(1,1,1),lwd=2)
	graphics.off()

}
