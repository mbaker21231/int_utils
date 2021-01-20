mata:
mata clear 
mata set matastrict on
real matrix int_add(real matrix X,real matrix Y,|real scalar digits)
{
	real matrix Z
	real scalar i

	if (args()==2) digits=1e-8
	Z=X+Y
	for (i=1;i<=cols(Z);i=i+2) {
		Z[,i]=r_down(Z[,i],digits)
		Z[,i+1]=r_up(Z[,i+1],digits)
							   }
	return(Z)
}
real matrix int_sub(real matrix A,real matrix B,|real scalar digits)
{
	real matrix X,Y,Z
	real scalar i

	if (args()==2) digits=1e-8
	Z=J(rows(A),cols(A),.)
	for (i=1;i<=cols(A);i=i+2) {
		X=A[.,i::i+1]
		Y=B[.,i::i+1]
		Z[.,i::i+1]=r_down(X[.,1]-Y[.,2],digits),r_up(X[.,2]-Y[.,1],digits)
								}
	return(Z)
}
real matrix r_up(real matrix X,real scalar digits)
{
	real matrix Xh
	real scalar tola

	tola=.1/digits
	Xh=X:*tola
	Xh=ceil(Xh)
	Xh=Xh/tola
	return(Xh)
}
real matrix r_down(real matrix X,real scalar digits)
{
	real matrix Xh
	real scalar tola
	
	tola=.1/digits
	Xh=X:*tola
	Xh=floor(Xh)
	Xh=Xh/tola
	return(Xh)
}
real matrix int_mult(real matrix X,real matrix Y,|real scalar digits)
{
	if (args()==2) digits=1e-8
	return(r_down(rowmin((X[.,1]:*Y[.,1],X[.,1]:*Y[.,2],X[.,2]:*Y[.,1],X[.,2]:*Y[.,2])),digits),
		r_up(rowmax((X[.,1]:*Y[.,1],X[.,1]:*Y[.,2],X[.,2]:*Y[.,1],X[.,2]:*Y[.,2])),digits))
}
real matrix int_pow(real matrix X,real scalar pow,|real scalar digits)
{
	real matrix Marker1,Marker2,X1
	real scalar foo
	
	if (args()==2) digits=1e-8

	if (pow>0) {
		X1=X:^pow
		Marker1=(X[,2]):<0
		Marker2=(X[,2]:>=0):*(X[,1]:<0)
		X1=abs(X):^pow
		X1[,1]=X1[,1]:*(1:-Marker2)
		X1=X1:*(1:-Marker1):+Marker1:*X1[,2::1]
	}
	else if (pow==0) X1=J(rows(X),2,0)
	else {
		Marker2=X:<=0
		X1=X[,2::1]
		X1=X1:^pow
		foo=runiform(1,1)
		X1=X1:*(1:-Marker2):+foo:*Marker2
		_editvalue(X1,foo,.)
	}
	X1=r_down(X1[,1],digits),r_up(X1[,2],digits)
	return(X1)
}
real matrix int_matmult(real matrix X,real matrix Y,|real scalar digits)
{
	real matrix Z,P
	real scalar i,j

	if (args()==2) digits=1e-8
	Z=J(rows(X),cols(Y),.)
	for (i=1;i<=rows(X);i++) {
		P=rowshape(X[i,.],cols(X)/2)
	for (j=1;j<=cols(Y)/2;j=j+1) {
		Z[i,(2*j-1)::2*j]=colsum(int_mult(P,Y[.,2*j-1::2*j],digits))
								}
							}
	return(Z)
}

real matrix int_hun(real matrix X, real matrix Y,|real scalar digits)
{
	if (args()==2) digits=1e-8
	return(r_down(rowmin((X[.,1],Y[.,1])),digits),r_up(rowmax((X[.,2],Y[.,2])),digits))
}
real matrix int_rowhun(real matrix X,real matrix Y,|real scalar digits)
{
	real matrix XY
	real scalar i
	
	if (args()==2) digits=1e-8
	XY=int_hun(X,Y,digits)
	if (cols(X)==2) return(XY)
	for (i=3;i<=cols(X);i=i+2) XY=XY,int_hun(X[,i::i+1],Y[,i::i+1],digits)
	return(XY)
}
real matrix int_int(real matrix X,real matrix Y)
{
	real scalar ind

	ind=rowmax((X[,1],Y[,1])):<=rowmin((X[,2],Y[,2]))
	_editvalue(ind,0,.)
	return(ind:*(rowmax((X[,1],Y[,1])),rowmin((X[,2],Y[,2]))))
}
real matrix int_rowmult(real matrix A,|real scalar digits) 
{
	real matrix R
	real scalar i

	if (args()==1) digits=1e-8
	R=J(rows(A),2,1)
	for (i=1;i<=cols(A)-1;i=i+2) {
		R=int_mult(R,A[.,i::i+1],digits)
								}
	return(R)
}
real matrix int_rowadd(real matrix A,|real scalar digits)
{
	real matrix R
	real scalar i

	if (args()==1) digits=1e-8	
	R=J(rows(A),2,0)
	for (i=1;i<=cols(A)-1;i=i+2) {
		R=int_add(R,A[.,i::i+1],digits)
					}
	return(R)
}
real matrix int_mid(real matrix P,|real scalar digits)
{
	real scalar i
	real matrix Pn,ave
	
	if (args()==1) digits=1e-8

	Pn=J(rows(P),0,.)
	for (i=1;i<=cols(P);i=i+2) {
		ave=(P[,i]:+P[,i+1])/2
		Pn=Pn,r_down(ave,digits),r_up(ave,digits)
								}
return(Pn)
}
real matrix mid(real matrix P)
{
	real scalar i
	real matrix Pn

	Pn=J(rows(P),0,.)
	for (i=1;i<cols(P);i=i+2) Pn=Pn,(P[,i]:+P[,i+1]):/2
	return(Pn)
}
real matrix int_mince(real matrix P,real scalar j,real scalar n,|real scalar digits)
{
	real matrix pl,pu,IL,ctr,pln,pun,NP
	
	if (args()==3) digits=1e-8
	
	pl=P[.,2*j-1]
	pu=P[.,2*j]
	IL=pu-pl
	ctr=J(rows(pl),1,(0::(n-1))/n)
	
	pl=pl#J(n,1,1)
	IL=IL#J(n,1,1)
	
	pln=pl:+ctr:*IL
	pun=pln+IL/n

	NP=P#J(n,1,1)
	NP[.,2*j-1]=r_down(pln,digits)
	NP[.,2*j]=r_up(pun,digits)
	return(NP)
}
real matrix int_mesh(real scalar n,real scalar t)
{
	real matrix I,K,Knew
	real scalar i

	I=J(0,2,.)
	for (i=1;i<=t;i++) {
		I=I\((i-1)/t,i/t)
			   }
	K=I
	for (i=2;i<=n;i++) {
		Knew=J(t,1,1)#K,I#J(rows(K),1,1)
		K=Knew
			  }
	return(K)
}
real matrix int_widths(real matrix P) 
{
	real matrix W
	real scalar i

	W=J(rows(P),0,.)
	for (i=1;i<=cols(P);i=i+2) W=W,P[.,i+1]-P[.,i]
	return(W)
}
real matrix int_rowint(real matrix R1,real matrix R2)
{
	real matrix Int
	real scalar i

	Int=J(rows(R1),0,.)
	for (i=1;i<cols(R1);i=i+2) {
		Int=Int,int_int(R1[,i::i+1],R2[,i::i+1])
				  }
	return(Int)	
}
real matrix int_rowhull(real matrix R1,|real matrix R2)
{
	real matrix Hull
	real scalar i

	if (args()==2) {
		Hull=J(rows(R1),0,.)
		for (i=1;i<cols(R1);i=i+2) Hull=Hull,rowmin(R1[,i],R2[i]),rowmax(R1[,i+1],R2[i+1])
		       }
	if (args()==1) {
		Hull=J(1,0,.)
		for (i=1;i<cols(R1);i=i+2) Hull=Hull,colmin(R1[,i]),colmax(R1[,i+1])
		       }
	return(Hull)
}
real matrix int_collect(real matrix P)
{
	real matrix Po,Pn,Test,Rec,Rec2
	real scalar i

	if (rows(P)==1) return(P)
	if (rows(P)==0) return(J(0,cols(P),.))

	Po=P
	Pn=P[1,.]

		for (i=2;i<=rows(Po);i++) {
			Test=int_rowint(Pn,Po[i,.]#J(rows(Pn),1,1))
			Rec=select(Pn,rownonmissing(Test):>0)
			if (rows(Rec)==0) Pn=Pn \ Po[i,.]
			else {
				Rec2=select(Pn,rownonmissing(Test):==0)
				Pn=Rec2 \ int_rowhull(Rec \ Po[i,.])
			     }
				          } 
	   	return(Pn)
}
real matrix int_rowmatmult(real matrix A,real matrix B,|real scalar digits)
{
	real matrix C,Bt,NE
	real scalar i,j,k,dim

	if (args()==2) digits=1e-8
	
	C=J(rows(A),0,.)
	dim=sqrt(cols(A)/2)
	Bt=int_rowtranspose(B)
	for (i=1;i<cols(A);i=i+2*dim) {
	for (j=1;j<cols(A);j=j+2*dim) {
		NE=J(rows(A),2,0)
		for (k=1;k<2*dim;k=k+2) {
			NE=int_add(NE,int_mult(A[,i+k-1::i+k],Bt[,j+k-1::j+k],digits),digits)
					}
		C=C,NE
				     }
				     }
	return(C)
}
real matrix int_rowmatvecmult(real matrix A,real matrix B,|real scalar digits)
{
	real matrix C,NE
	real scalar dim,k,i
	
	if (args()==2) digits=1e-8
	
	C=J(rows(A),0,.)
	dim=sqrt(cols(A)/2)
	for (i=1;i<cols(A);i=i+2*dim) {
		NE=J(rows(A),2,0)
		for (k=1;k<2*dim;k=k+2) {
			NE=int_add(NE,int_mult(A[,i+k-1::i+k],B[,k::k+1],digits),digits)
				      }
		C=C,NE
		                       }

	return(C)
}
real matrix int_rowtranspose(real scalar A)
{
	real matrix PV
	real scalar dim,i,j

	PV=J(1,0,.)
	dim=sqrt(cols(A)/2)
	for (i=1;i<=dim;i++) {	
	for (j=1;j<=dim;j++) {
	PV=PV,(2*(i+(j-1)*dim)-1),2*(i+(j-1)*dim)
	}
	}
	return(A[,PV])
}
real matrix radius(real matrix P)
{
	real matrix R
	real scalar i

	R=J(rows(P),0,.)
	for (i=1;i<cols(P);i=i+2) R=R,(P[,i+1]:-P[,i])/2	
	return(R)           
}
real matrix int_transpose(real matrix pI) 
{
	real scalar i,k
	real matrix pHolder
	
	k=1
	pHolder=J(cols(pI)/2,rows(pI)*2,.)
	for (i=1;i<=cols(pI);i=i+2) {
		pHolder[k,]=rowshape(pI[,i::i+1],1)
		k++
	}

	return(pHolder)
}
mata mlib create lint_utils, dir(PERSONAL) replace
mata mlib add lint_utils *()
mata mlib index
end
