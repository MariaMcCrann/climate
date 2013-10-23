library(compiler)

DCT1D<-function(k){  # intrinsic
  #Spectral rep of an AR1 inverse covariance
  #t(M)%*%M = M%*%t(M) = diag(k)
  #Precision = t(M)%*%diag(d)%*%M

  M<-diag(k)
  for(j in 1:k){
     M[1,j]<-1/sqrt(k)
     for(i in 2:k){
       M[i,j]<-sqrt(2/k)*cos(pi*(i-1)*(j-1/2)/k)
  }}
  i<-1:k
  d<-2*(1-cos(pi*(i-1)/k))
  cycles<-(i-1)/2
  period<-2*k/(i-1)
list(M=M,d=d,cycles=cycles,period=period)}

"dct" <- function(n, y) {
	a <- c(1/sqrt(n), rep(sqrt(2/n), n-1))
	s <- t <- 0:(n-1)

	z <- a * sapply(t, function(time) {
		sum( cos(time * pi * (s+0.5)/n)*y )
	})
}

"dct2d" <- cmpfun(function(n1, n2, y) {
	a <- c(1/sqrt(n1), rep(sqrt(2/n1), n1-1))
	b <- c(1/sqrt(n2), rep(sqrt(2/n2), n2-1))
	s1 <- t1 <- 0:(n1-1)
	s2 <- t2 <- 0:(n2-1)

	z <- matrix(NA, nrow=n1, ncol=n2)
	sapply(0:(n1-1), function(t1) {
		sapply(0:(n2-1), function(t2) {
			z[t1+1,t2+1] <- a[t1+1]*b[t2+1]*sum(sapply(0:(n1-1), function(s1) {
				sapply(0:(n2-1), function(s2) {
					cos(t1 * pi * (s1+0.5)/n1) * cos(t2 * pi * (s2+0.5)/n2) * y[s1+1,s2+1]
				})
			}))
		})
	})

	z
})

if (FALSE) {
DCT1D<-function(k){
  #Spectral rep of an AR1 inverse covariance
  #t(M)%*%M = M%*%t(M) = diag(k)
  #Precision = t(M)%*%diag(d)%*%M

  M<-diag(k)
  for(j in 1:k){
     M[1,j]<-1/sqrt(k)
     for(i in 2:k){
       M[i,j]<-sqrt(2/k)*cos(pi*(i-1)*(j-1/2)/k)
  }}
  i<-1:k-1
  d<-2*(1-cos(pi*i/k))
  cycles<-i/2
  period<-2*k/i
list(M=M,d=d,cycles=cycles,period=period)}
}


#Some example to understand these operations:
if(FALSE){

  #Plots of the basis functions for a simple case
  k<-25
  M<-DCT1D(k)$M
  d<-DCT1D(k)$d
  p<-DCT1D(k)$period
  c<-DCT1D(k)$cycles
  j<-4

  library(fields)
  plot(M[j,],type="l",main=paste("# cycles =",c[j]));abline(v=p[j]+.5)
  image.plot(M%*%t(M))
  image.plot(t(M)%*%M)
  image.plot(t(M)%*%diag(d)%*%M)


  #Filtering example
  library(fields)
  n1<-20
  n2<-50
  Y<-matrix(0.1*(1:n1-n1/10)^2,n1,n2,byrow=F)+
     matrix(5*sin(1*1:n2),n1,n2,byrow=T)+
     rnorm(n1*n2)
  image.plot(Y)

  DCT1<-DCT1D(n1)
  M1<-DCT1$M
  period1<-DCT1$period
  period1[1]<-2*period1[2]
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
  rm(DCT2)


  Z<-M1%*%Y%*%t(M2)
  Zl<-ifelse(C1+C2<5,Z,0)
  Zh<-ifelse(C1+C2>=5,Z,0)
  Yl<-t(M1)%*%Zl%*%M2
  Yh<-t(M1)%*%Zh%*%M2

  par(mfrow=c(2,2))
  image.plot(Y,main="Raw data")
  image.plot(cycle1,cycle2,Z,main="Transformed data")
  image.plot(Yl,main="low res")
  image.plot(Yh,main="high res")
}




