####
alphamatrix<-function(n){
  A<-matrix(NA,ncol = (n-1),nrow = (n-1))
  p<-(n-1)*n/2

  A[row(A)-col(A) >= 0] = c(1:(n-1),-(n:p))
  A[row(A)-col(A) < 0] = n:p

  return(A)
}

#### Determine zero or non-yero output
zerooneoutput<-function(j,rho,A){

  v<-j*rho

  eval<-function(x){
    return(sum(sign(x)*v[abs(x)]))
  }

  if(all(apply(X = A,MARGIN = 1,FUN = eval) == 0)){
    return(1)
  }else{
    return(0)
  }
}


#### Check non-zero's



#### b's


btcol<-function(col){

  t<-col[1]
  a<-col[2]
  jt<-col[3]

  if(jt == 0){
    return(besselI(abs(a), nu = 0))
  }else{
    if(a < 0){
      return(2*(-1)^jt*besselI(abs(a),nu = jt))
    }else{
      return(2*besselI(abs(a),nu = jt))
    }
  }
}


btprod<-function(t,a,Jstar){
  return(prod(apply(X = cbind(1:ncol(Jstar),a,Jstar[t,]),MARGIN = 1,FUN = btcol)))
}




# DENSITIES

drayl3D<-function(dK,Ccomp,lim=3){

  n=3
  p=3

  #### Setup of order Matrix

  permat<-expand.grid(rep(list(c(-1,1)), p))

  #### Exponential Product/ sum loops

  # set up j's, j_1=1,...,lim, others always smaller than the precursor

  A = alphamatrix(n = n)
  J=as.matrix(expand.grid(rep(list(seq(0,lim)),p)))
  J=unique(t(apply(J,1,sort,decreasing=T)))
  Jstar<-(J-cbind(J[,-1],0))[,n:1]

  retcalculator<-function(l,Jstar,permat){

    eval<-function(i){
      return(zerooneoutput(j = Jstar[l,],rho = permat[i,],A = A))
    }

    return(sum(sapply(X = 1:nrow(permat),FUN = eval)))
  }

  retjslst<-sapply(X = 1:nrow(Jstar),FUN = retcalculator,Jstar = Jstar,permat = permat)
  nonzeroes<-which(retjslst != 0)

  func<-function(r){

    #### a's

    a<-(-1/dK)*c(r[1]*r[2]*Ccomp[1,2],r[1]*r[3]*Ccomp[1,3],r[2]*r[3]*Ccomp[2,3])
    inside_sum<-sum(sapply(X = nonzeroes,FUN = btprod,Jstar = Jstar, a =a)*retjslst[nonzeroes])

    return(prod(r)/sqrt(dK)*exp(-1/(2*dK)*sum(diag(Ccomp)*r*r))/(2^p)*inside_sum)
  }

  return(func)
}



drayl4D<-function(dK,Ccomp,lim=3){

  n=4
  p=6

  #### Setup of order Matrix

  permat<-expand.grid(rep(list(c(-1,1)), p))

  #### Exponential Product/ sum loops

  # set up j's, j_1=1,...,lim, others always smaller than the precursor

  A = alphamatrix(n = n)
  Jstar=as.matrix(expand.grid(rep(list(seq(0,floor(lim/3))),p)))
  fr<-which(rowSums(Jstar[,1:3]) %% 2 == 0)
  ar2<-which(rowSums(Jstar)!= 2)
  ar5<-which(rowSums(Jstar)!= 5)
  Jstar=Jstar[Reduce(intersect, list(fr,ar2,ar5)),]


  retcalculator<-function(l,Jstar,permat){

    eval<-function(i){
      return(zerooneoutput(j = Jstar[l,],rho = permat[i,],A = A))
    }

    return(sum(sapply(X = 1:nrow(permat),FUN = eval)))
  }

    retjslst<-sapply(X = 1:nrow(Jstar),FUN = retcalculator,Jstar = Jstar,permat = permat)


  nonzeroes<-which(retjslst != 0)
  sqdK<-1/sqrt(dK)/2^p

  func<-function(r){

    #### a's

    a<-(-1/dK)*c(r[1]*r[2]*Ccomp[1,2],r[1]*r[3]*Ccomp[1,3],r[1]*r[4]*Ccomp[1,4],r[2]*r[3]*Ccomp[2,3],r[2]*r[4]*Ccomp[2,4],r[3]*r[4]*Ccomp[3,4])
    inside_sum<-sum(sapply(X = nonzeroes,FUN = btprod,Jstar = Jstar, a =a)*retjslst[nonzeroes])

    return(prod(r)*sqdK*exp(-1/(2*dK)*sum(diag(Ccomp)*r*r))*inside_sum)}

  return(func)
}







drayl_int3D<-function(r,omega,sigma,cor,method = "mixed"){

  # scalars
  g1<-r[1]/sigma[1]^2*exp(-1/(2*det(omega))*(r[1]^2*(1-cor[3]^2)/sigma[1]^2))
  g2<-r[2]/sigma[2]^2*exp(-1/(2*det(omega))*(r[2]^2*(1-cor[2]^2)/sigma[2]^2))
  g3<-r[3]/sigma[3]^2*exp(-1/(2*det(omega))*(r[3]^2*(1-cor[1]^2)/sigma[3]^2))

  scalar<-g1*g2*g3

  L1<-r[2]*r[3]/(sigma[2]*sigma[3]*det(omega))*(cor[3]-cor[1]*cor[2])
  L2<-r[1]*r[3]/(sigma[1]*sigma[3]*det(omega))*(cor[2]-cor[1]*cor[3])
  L3<-r[1]*r[2]/(sigma[1]*sigma[2]*det(omega))*(cor[1]-cor[2]*cor[3])

  integrand<-function(x){
    return(exp(L2*cos(x))*besselI(x = sqrt(L1^2+L3^2+2*L1*L3*cos(x)),nu = 0)/(pi*det(omega)))
  }

  if(method == "Kronrod" | method == "Simpson" | method == "Clenshaw"){
    return(scalar*pracma::integral(fun = integrand,xmin = 0,xmax = 3.14159265,method = method))
  }
  if(method == "Romberg" | method == "TOMS614"){
    return(scalar*rmutil::int(f = integrand,a = 0,b = 3.14159265,type = method))
  }
  if(method == "mixed"){
    return(scalar*stats::integrate(f = integrand,lower = 0,upper = 3.14159265)$value)
  }


}






drayl_int4D<-function(r,omega,sigma,cor,method = "Quadrature"){

  phi1<-matrix(c(1,cor[4],cor[5],cor[4],1,cor[6],cor[5],cor[6],1),ncol = 3)
  phi2<-matrix(c(1,cor[2],cor[3],cor[2],1,cor[6],cor[3],cor[6],1),ncol = 3)
  phi3<-matrix(c(1,cor[1],cor[3],cor[1],1,cor[5],cor[3],cor[5],1),ncol = 3)
  phi4<-matrix(c(1,cor[1],cor[2],cor[1],1,cor[4],cor[2],cor[4],1),ncol = 3)

  detphi<-c(det(phi1),det(phi2),det(phi3),det(phi4))

  scalar<-prod(r)/(4*pi^2*prod(sigma^2)*det(omega))*exp(-1/(2*det(omega))*sum(r^2/sigma^2*detphi))

  M12<-t(matrix(c(cor[1],cor[4],cor[5],cor[2],1,cor[6],cor[3],cor[6],1),ncol = 3))
  M13<-t(matrix(c(cor[2],cor[4],cor[6],cor[1],1,cor[5],cor[3],cor[5],1),ncol = 3))
  M14<-t(matrix(c(cor[3],cor[5],cor[6],cor[1],1,cor[4],cor[2],cor[4],1),ncol = 3))
  M23<-t(matrix(c(cor[4],cor[2],cor[6],cor[1],1,cor[3],cor[5],cor[3],1),ncol = 3))
  M24<-t(matrix(c(cor[5],cor[3],cor[6],cor[1],1,cor[2],cor[4],cor[2],1),ncol = 3))
  M34<-t(matrix(c(cor[6],cor[3],cor[5],cor[2],1,cor[1],cor[4],cor[1],1),ncol = 3))

  # Final L- Coefficient setup
  L12<-r[1]*r[2]/(sigma[1]*sigma[2]*det(omega))*det(M12)
  L13<-r[1]*r[3]/(sigma[1]*sigma[3]*det(omega))*det(M13)
  L14<-r[1]*r[4]/(sigma[1]*sigma[4]*det(omega))*det(M14)
  L23<-r[2]*r[3]/(sigma[2]*sigma[3]*det(omega))*det(M23)
  L34<-r[3]*r[4]/(sigma[3]*sigma[4]*det(omega))*det(M34)
  L24<-r[2]*r[4]/(sigma[2]*sigma[4]*det(omega))*det(M24)



  if(method == "Romberg"){
    # Integrand function
    integrandR<-function(x,y){
      return(exp(L13*cos(x) + L14*cos(y)+ L34*cos(y-x))*besselI(sqrt(L12^2 + L23^2 + L24^2 + 2*L12*L23*cos(x) + 2*L12*L24*cos(y) + 2*L23*L24*cos(x-y)),nu = 0))
    }
    return(scalar*rmutil::int2(f = integrandR,a = c(0,0),b = c(1,1)*2*3.1415926))
  }
  if(method == "Cubature"){
    # Integrand function
    integrandC<-function(t){
      return(exp(L13*cos(t[1]) + L14*cos(t[2])+ L34*cos(t[2]-t[1]))*besselI(sqrt(L12^2 + L23^2 + L24^2 + 2*L12*L23*cos(t[1]) + 2*L12*L24*cos(t[2]) + 2*L23*L24*cos(t[1]-t[2])),nu = 0))
    }
    return(scalar*cubature::pcubature(f = integrandC,lowerLimit = c(0,0),upperLimit = c(1,1)*2*3.1415926)$integral)}
  if(method == "Quadrature"){
    # Integrand function
    integrandQ<-function(x,y){
      return(exp(L13*cos(x) + L14*cos(y)+ L34*cos(y-x))*besselI(sqrt(L12^2 + L23^2 + L24^2 + 2*L12*L23*cos(x) + 2*L12*L24*cos(y) + 2*L23*L24*cos(x-y)),nu = 0))
    }
    return(scalar*pracma::quad2d(f = integrandQ,xa = 0,xb = 2*3.1415926,ya = 0,yb = 2*3.1415926))}

}



