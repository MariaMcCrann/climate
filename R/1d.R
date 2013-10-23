# shows that Q = M D M' in one dimension

DCT1D_i<-function(k){  # intrinsic
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
  d<-2*(1-cos(pi*(i-0)/k))
  cycles<-i/2
  period<-2*k/i
list(M=M,d=d,cycles=cycles,period=period)}

DCT1D_r<-function(k, rho){  # with rho
  #Spectral rep of an AR1 inverse covariance
  #t(M)%*%M = M%*%t(M) = diag(k)
  #Precision = t(M)%*%diag(d)%*%M

  M<-diag(k)
  for(j in 1:k){
     M[1,j]<-1/sqrt(k)
     for(i in 2:k){
       M[i,j]<- sqrt(rho)*sqrt(2/k)*cos(pi*(i-1)*(j-1/2)/k)
  }}
  i<-1:k-1
  d<-2*(1-cos(pi*i/k))
  #d<-2-2*cos(pi*i/k)
  cycles<-i/2
  period<-2*k/i
list(M=M,d=d,cycles=cycles,period=period)}

icar <- DCT1D_i(4)
print(icar$M); print(icar$d); print(round(t(icar$M) %*% diag(icar$d) %*% icar$M,4))

done

if (FALSE) {
Q <- round(t(icar$M) %*% diag(icar$d) %*% icar$M,4)

print(round( (icar$M) %*% Q %*% t(icar$M), 4))

Q[Q==-1] <- -0.5
print(round(Q,4))

print(round( (icar$M) %*% Q %*% t(icar$M), 4))
done

C <- sqrt(diag(0.35/diag(Q)))

print(round( (icar$M) %*% C %*% Q %*% C %*% t(icar$M), 4))
done

}

car <- DCT1D_r(10, 0.5)
print(car$M); print(car$d); print(round(t(car$M) %*% diag(car$d) %*% car$M,4))
