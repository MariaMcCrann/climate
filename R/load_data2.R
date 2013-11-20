source("R/DCTfx.R")

if (!("SunTemp" %in% ls())) {
cat("Loading data\n")

load("data/temp.save")
load("data/prec.save")

#Used to remove rows and columns that are all ocean.
cut<-function(Y,good1,good2){
  Y<-Y[good1,];Y<-Y[,good2]
Y}

good1<-1:128 > 11 & 1:128 < 128-14
good2<-1:112 >  0 & 1:112 < 112-42

lon<-cut(summer.temp$lon,good1,good2)
lat<-cut(summer.temp$lat,good1,good2)

n1<-nrow(lon)
n2<-ncol(lat)
nm<-9
models<-c("NCEP","CRCM","ECP2","HRM3","MM5I","RCM3","WRFG","UDEL","CRU")
type<-c("BC","RCM","RCM","RCM","RCM","RCM","RCM","Data","Data")

##################################################
#Put the data in array form
##################################################

SumTemp<-WinTemp<-SumPrec<-WinPrec<-array(0,c(n1,n2,9))

SumTemp[,,1]<-cut(summer.temp$NCEP,good1,good2)
SumTemp[,,2]<-cut(summer.temp$CRCM,good1,good2)
SumTemp[,,3]<-cut(summer.temp$ECP2,good1,good2)
SumTemp[,,4]<-cut(summer.temp$HRM3,good1,good2)
SumTemp[,,5]<-cut(summer.temp$MM5I,good1,good2)
SumTemp[,,6]<-cut(summer.temp$RCM3,good1,good2)
SumTemp[,,7]<-cut(summer.temp$WRFG,good1,good2)
SumTemp[,,8]<-cut(summer.temp$UDEL,good1,good2)
SumTemp[,,9]<-cut(summer.temp$CRU,good1,good2)

WinTemp[,,1]<-cut(winter.temp$NCEP,good1,good2)
WinTemp[,,2]<-cut(winter.temp$CRCM,good1,good2)
WinTemp[,,3]<-cut(winter.temp$ECP2,good1,good2)
WinTemp[,,4]<-cut(winter.temp$HRM3,good1,good2)
WinTemp[,,5]<-cut(winter.temp$MM5I,good1,good2)
WinTemp[,,6]<-cut(winter.temp$RCM3,good1,good2)
WinTemp[,,7]<-cut(winter.temp$WRFG,good1,good2)
WinTemp[,,8]<-cut(winter.temp$UDEL,good1,good2)
WinTemp[,,9]<-cut(winter.temp$CRU,good1,good2)

SumPrec[,,1]<-cut(summer.prec$NCEP,good1,good2)
SumPrec[,,2]<-cut(summer.prec$CRCM,good1,good2)
SumPrec[,,3]<-cut(summer.prec$ECP2,good1,good2)
SumPrec[,,4]<-cut(summer.prec$HRM3,good1,good2)
SumPrec[,,5]<-cut(summer.prec$MM5I,good1,good2)
SumPrec[,,6]<-cut(summer.prec$RCM3,good1,good2)
SumPrec[,,7]<-cut(summer.prec$WRFG,good1,good2)
SumPrec[,,8]<-cut(summer.prec$UDEL,good1,good2)
SumPrec[,,9]<-cut(summer.prec$CRU,good1,good2)

WinPrec[,,1]<-cut(winter.prec$NCEP,good1,good2)
WinPrec[,,2]<-cut(winter.prec$CRCM,good1,good2)
WinPrec[,,3]<-cut(winter.prec$ECP2,good1,good2)
WinPrec[,,4]<-cut(winter.prec$HRM3,good1,good2)
WinPrec[,,5]<-cut(winter.prec$MM5I,good1,good2)
WinPrec[,,6]<-cut(winter.prec$RCM3,good1,good2)
WinPrec[,,7]<-cut(winter.prec$WRFG,good1,good2)
WinPrec[,,8]<-cut(winter.prec$UDEL,good1,good2)
WinPrec[,,9]<-cut(winter.prec$CRU,good1,good2)

rm(summer.temp,winter.temp,summer.prec,winter.prec,cut)
} # end load data

if (TRUE) {
#Perform the discete cosine transformation

DCT1<-DCT1D(n1)
M1<-DCT1$M
period1<-DCT1$period
cycle1<-DCT1$cycles
d1<-DCT1$d
P1<-matrix(period1,n1,n2,byrow=FALSE)
C1<-matrix(cycle1,n1,n2,byrow=FALSE)
D1<-matrix(d1,n1,n2,byrow=FALSE)
rm(DCT1)

DCT2<-DCT1D(n2)
M2<-DCT2$M
period2<-DCT2$period
period2[1]<-2*period2[2]
cycle2<-DCT2$cycles
d2<-DCT2$d
P2<-matrix(period2,n1,n2,byrow=TRUE)
C2<-matrix(cycle2,n1,n2,byrow=TRUE)
D2<-matrix(d2,n1,n2,byrow=TRUE)
D<-D1+D2
rm(DCT2)

SumTempTrans<-WinTempTrans<-SumPrecTrans<-WinPrecTrans<-array(0,c(n1,n2,9))

for (j in 1:9) {
	SumTempTrans[,,j]<-M1%*%SumTemp[,,j]%*%t(M2)
	WinTempTrans[,,j]<-M1%*%WinTemp[,,j]%*%t(M2)
	SumPrecTrans[,,j]<-M1%*%SumPrec[,,j]%*%t(M2)
	WinPrecTrans[,,j]<-M1%*%WinPrec[,,j]%*%t(M2)
}


# A function of plots:
splot<-function(lon,lat,y,map=TRUE,main=""){
	require(fields)
	require(maps)
	image.plot(lon,lat,y,xlab="",ylab="",axes=F,main=main)
	if(map){map("world",add=T)} 
}

if (FALSE) {
	cat("Filtering\n")
	# Example of filtering and plotting
	Z  <-WinTempTrans[,,9]
	Zl <-ifelse(D<0.25,Z,0)
	Zm <-ifelse(D>0.25 & D<1,Z,0)
	Zh <-ifelse(D>1,Z,0)
	Yl <-t(M1)%*%Zl%*%M2
	Ym <-t(M1)%*%Zm%*%M2
	Yh <-t(M1)%*%Zh%*%M2
	pdf("filtering.pdf")
		par(mfrow=c(2,2))
		splot(lon,lat,WinTemp[,,9],main="Data")
		splot(lon,lat,Yl,main="low res")
		splot(lon,lat,Ym,main="medium res")
		splot(lon,lat,Yh,main="high res")
	graphics.off()
done
}


#rm(good1,good2,j,DCT1D)
#save.image("RCM.RData")
} # end trans test

# which data is of interest?
if (!("WHICH_CDAT" %in% ls()) || WHICH_CDAT == "ST") {
	cdata <- SumTemp
} else if (WHICH_CDAT == "WT") {
	cdata <- WinTemp
} else if (WHICH_CDAT == "SP") {
	cdata <- SumPrec
} else if (WHICH_CDAT == "WP") {
	cdata <- WinPrec
}
tlon  <- lon
tlat  <- lat
n1    <- nrow(cdata)
n2    <- ncol(cdata)

# get Zs
icar1 <- DCT1D(n1)
icar2 <- DCT1D(n2)

#M <- kronecker(icar1$M, icar2$M)
#d <- diag(kronecker.spam(as.spam(diag(icar1$d)), as.spam(diag(n2))) + kronecker.spam(as.spam(diag(n1)), as.spam(diag(icar2$d))))
D1 <- matrix(icar1$d,n1,n2,byrow=FALSE)
D2 <- matrix(icar2$d,n1,n2,byrow=TRUE)
D  <- D1 + D2
T1 <- matrix(1:n1,n1,n2,byrow=FALSE)-1
T2 <- matrix(1:n2,n1,n2,byrow=TRUE)-1
T  <- T1+T2
#F  <- T1/(2*n1) + T2/(2*n2)
F  <- sqrt( (T1/(2*n1))^2 + (T2/(2*n2))^2 )
P  <- 1/F

Z <- array(NA, dim=c(n1,n2,9))
sapply(1:9, function(i) {
  Z[,,i] <<- icar1$M %*% cdata[,,i] %*% t(icar2$M)
})

d <- as.vector(D)
z <- matrix(NA, nrow=n1*n2, ncol=9)
sapply(1:9, function(i) {
  z[,i] <<- as.vector(Z[,,i])
})
zstar <- sqrt(d) * z
t <- as.vector(T) + 1 # start at 1
f <- as.vector(F)
p <- as.vector(P)
