/* literature

Willim J. Thompson, <a href="https://doi.org/10.1002/9783527617821">Angular Momentum</a>, Wiley 1994, Eq 7.51

A 6j coefficient is invariant under interchange of any two columns
and under interchange of the upper and lower arguments in each of any two columns

G Junker <a href="https://doi.org/10.1088/0305-4470/26/7/021">Explicit evaluation of coupling coefficients for the most degendrate representations of SO(n)</a>, J Phys A 26 (1993) 1649
*/

/*  return the square of Delta(a,b,c)
*/
GreekD(a,b,c) = 
{
	if( a+b-c < 0 || denominator(a+b-c) != 1,
		return(0)  ;
	) ;
	if( a+c-b < 0 || denominator(a+c-b) != 1,
		return(0)  ;
	) ;
	if( b+c-a < 0 || denominator(b+c-a) != 1,
		return(0)  ;
	) ;
	if( a+b+c+1 < 0 || denominator(a+b+c+1) != 1,
		return(0)  ;
	) ;
	return((a+b-c)!*(a+c-b)!*(b+c-a)!/(a+b+c+1)!) ;
}

/* return the Clebsch-Gordon coefficient, including the sign
*/
cgord(j1,j2,m1,m2,J,M) = 
{
	local(k,kmin,kmax,resul,f,sg) ;
	if(m1+m2 != M, return(0) );
	if(J>abs(j1+j2), return(0) );
	if(J<abs(j1-j2), return(0) );
	kmax=j2+J+m1 ;
	kmax=min(kmax,J-j1+j2) ;
	kmax=min(kmax,J+M) ;
	kmin=m1-j1 ;
	kmin=max(kmin,j2+M-j1) ;
	kmin=max(kmin,0) ;
	resul=0 ;
	for(k=kmin,kmax,
		resul +=(-1)^(k+j2+m2)*(j2+J+m1-k)!*(j1-m1+k)!/k!/(J-j1+j2-k)!
			/(J+M-k)!/(k+j1-j2-M)!
	) ;
	if(resul==0,
		return(0) ;
	) ;
	/*
	* if(2*J+1<0, return(0) );
	* if(J<abs(M), return(0) );
	* if(j1<abs(m1), return(0) );
	* if(j2<abs(m2), return(0) );
	*/
	resul *= sqrt(2*J+1) ;
	f = (J+j1-j2)!*(J-j1+j2)!*(j1+j2-J)!*(J+M)!*(J-M)!/(j1+j2+J+1)!/(j1-m1)!
		/(j1+m1)!/(j2-m2)!/(j2+m2)! ;
	resul *= sqrt(f) ;
	return(resul) ;
}

coupl3j(j1,j2,j3,m1,m2,m3) =
{
	/* sum rule of the magnetic quantum numbers must be satisfied
	*/
	if( m1+m2+m3 != 0,
		return(0)
	) ;
	/* the three ji-mi must be integer
	*/
	if( denominator(j1-m1) != 1,
		return(0) 
	) ;
	if( denominator(j2-m2) != 1,
		return(0) 
	) ;
	if( denominator(j3-m3) != 1,
		return(0) 
	) ;
	/* the absolute values of the three mi must be <= ji
	*/
	*/
	if( abs(m1) > j1,
		return(0) 
	) ;
	if( abs(m2) > j2,
		return(0) 
	) ;
	if( abs(m3) > j3,
		return(0) 
	) ;
	return ( (-1)^(j1-j2-m3)/sqrt(2*j3+1)*cgord(j1,j2,m1,m2,j3,-m3)) ;
}

/*
* J. Rasch and A. C. H. Yu, <a href="https://doi.org/10.1137/S1064827503422932">Efficient storage scheme for precalculated Wigner 3j, 6j and Gaunt Coefficeints</a>, SIAM J Sci Comput 25 (4) (2003) 1416-1428
* 
* A. R. Edmonds, Drehimpulse in der Quatenmechanik, BI Hochschultaschenbucher 53/53a (1964) eq. (6.37)
*
* A. R. Edmonds, <a href="https://cds.cern.ch/record/212213/">Angular momentum in quantum mechanics</a>, (1957) eq. (4.6.1)
*/
coupl6jRasch(j1,j2,j3,j4,j5,j6)=
{
	local(D,kmin,kmax) ;
	D=GreekD(j1,j2,j3)*GreekD(j3,j4,j5)*GreekD(j1,j5,j6)*GreekD(j2,j4,j6) ;
	if ( D == 0,
		return(0)
	) ;

	/*
	* if( denominator(a+b+c+d) != 1
	* 	|| denominator(e+f-a-d) != 1
	* 	|| denominator(e+f-b-c) != 1
	* 	|| denominator(a+b-e) != 1
	* 	|| denominator(c+d-e) != 1
	* 	|| denominator(a+c-f) != 1
	* 	|| denominator(b+d-f) != 1,
	* 	return(0)
	*)  ;
	*/

	kmax=j1+j2+j4+j5 ;
	kmin=j1+j2+j3 ;

	kmin=max(kmin,j1+j5+j6) ;
	kmin=max(kmin,j2+j4+j6) ;
	kmin=max(kmin,j3+j4+j5) ;

	kmax=min(kmax,j1+j3+j4+j5) ;
	kmax=min(kmax,j2+j3+j5+j6) ;

	resul=0 ;
	for(k=kmin,kmax,
		resul +=(-1)^k*(1+k)!/(k-j1-j2-j3)!/(k-j1-j5-j6)!/(k-j2-j4-j6)!
		/(k-j3-j4-j5)!/(j1+j2+j4+j5-k)!/(j1+j3+j4+j6-k)!/(j2+j3+j5+j6-k)! ;
	) ;
	if(resul==0,
		return(0) ;
	) ;
	sg=sign(resul) ;
	resul = resul*resul ;
	resul *= D ;
	return(sg*resul) ;
}

coupl6jAlt(a,b,e,d,c,f) =
{
	local(resul,epsilon,Gamma,Phi) ;
	resul=0 ;
	/* We use the sum over (2c+1), whence the sum over Gamma is already done, and the
	* other coefficients have to adapt. We obviously chose a Gamma that is always valid, +- 1/2 for half-integer c,
	* 0 for integer c.
	*/
	if ( denominator(c) == 1,
		Gamma=0,
		Gamma=1/2
	) ;
	for(alpha=-a,a,
		for(beta=-b,b,
			epsilon=alpha+beta ;
			delta=Gamma-epsilon ;
			Phi=beta+delta ;
			if( Gamma==alpha+Phi,
				resul += (-1)^(f-e-alpha-delta)*(2*c+1)
					*coupl3j(a,b,e,alpha,beta,-epsilon)
					*coupl3j(e,d,c,epsilon,delta,-Gamma)
					*coupl3j(b,d,f,beta,delta,-Phi)
					*coupl3j(a,f,c,alpha,Phi,-Gamma) ;
			) ;
		) ;
	) ;
	return( resul*(-1)^(a+b+c+d)) ;
}

/** Richard N. Zare, Angular Momentum, Wiley (1988), Eq. (4.14)
*/
coupl6jZare(j1,j2,j3,j4,j5,j6) =
{
	local(resul,m3,m5,m6) ;
	/* below we multiply with 2j5+1 which implies that the corresponding m5 is fixed (because it's sum is already done)
	* A valid value of m5 depends on whether j5 is integer or half-integer
	*/
	resul=0 ;
	if ( denominator(j5) == 1,
		m5=0,
		m5=1/2
	) ;
	for(m1=-j1,j1,
		for(m2=-j2,j2,
			m3=m1+m2 ;
			m4= -m5-m3 ;
			m6=m2+m4 ;
			if( m5+m1+m6==0,
				resul += (-1)^(-m1-m4)*(2*j5+1)
					*coupl3j(j1,j2,j3,m1,m2,-m3)
					*coupl3j(j4,j5,j3,m4,m5,m3)
					*coupl3j(j2,j4,j6,m2,m4,-m6)
					*coupl3j(j5,j1,j6,m5,m1,m6) ;
			) ;
		) ;
	) ;
	return( resul*(-1)^(j1+j2-j3+j4+j5+j6)) ;
}

/** Richard N. Zare, Angular Momentum, Wiley (1988), eq. 4.13)
*/
coupl6jZareBase(j1,j2,j3,j4,j5,j6) =
{
	local(resul,m3,m6,m4) ;
	resul=0 ;
	for(m1=-j1,j1,
		for(m2=-j2,j2,
			m3=-m1-m2 ;
			for(m5=-j5,j5,
				m6=m5-m1 ;
				m4=m6-m2 ;
				if( m5+m3-m4==0,
					resul += (-1)^(j4-m4+j5-m5+j6-m6)
						*coupl3j(j1,j2,j3,m1,m2,m3)
						*coupl3j(j1,j5,j6,m1,-m5,m6)
						*coupl3j(j4,j2,j6,m4,m2,-m6)
						*coupl3j(j4,j5,j3,-m4,m5,m3) ;
				) ;
			) ;
		) ;
	) ;
	return( resul) ;
}

/**
* Return the product of the 4 Delta-coefficients.
*/
tstDelta(a,b,e,d,c,f) =
{
	return(
	GreekD(a,b,e)*GreekD(a,c,f)*GreekD(b,d,f)*GreekD(c,d,e) 
	) ;
}

/*
* return the square of the 6j coefficient, with a sign
*/
coupl6j(a,b,e,d,c,f) = 
{
	local(D,kmin,kmax) ;
	D=GreekD(a,b,e)*GreekD(a,c,f)*GreekD(b,d,f)*GreekD(c,d,e) ;
	if ( D == 0,
		return(0)
	) ;

	if( denominator(a+b+c+d) != 1
		|| denominator(e+f-a-d) != 1
		|| denominator(e+f-b-c) != 1
		|| denominator(a+b-e) != 1
		|| denominator(c+d-e) != 1
		|| denominator(a+c-f) != 1
		|| denominator(b+d-f) != 1,
		return(0)
	) ;

	kmax=a+b+c+d+1 ;
	kmin=0 ;
	kmin=max(kmin,-e-f+a+d) ;
	kmin=max(kmin,-e-f+b+c) ;
	kmax=min(kmax,a+b-e) ;
	kmax=min(kmax,c+d-e) ;
	kmax=min(kmax,a+c-f) ;
	kmax=min(kmax,b+d-f) ;

	resul=0 ;
	for(k=kmin,kmax,
		resul +=(-1)^k*(a+b+c+d+1-k)!/k!/(e+f-a-d+k)!/(e+f-b-c+k)!/(a+b-e-k)!/(c+d-e-k)!/(a+c-f-k)!/(b+d-f-k)! ;
	) ;
	if(resul==0,
		return(0) ;
	) ;
	resul *= (-1)^(a+b+c+d) ;
	sg=sign(resul) ;
	resul = resul*resul ;
	resul *= D ;
	return(sg*resul) ;
}

/**
* Print a human-readable version of the square root.
* This essentially pulls out a common rational square from the square root.
* @param sg The sign of the entire expression
* @param truecoe the value underneath the square root
*/
prsqroo(sg,truecoe)=
{
	local(num,numf,numff,den,denf,denff,pre) ;
	num=numerator(truecoe) ;
	den=denominator(truecoe) ;
	denf=factor(den) ;
	pre=1 ;
	if(num != 1,
		numf=factor(num) ;
		numff=(matsize(numf))[1] ;
		for(i=1,numff,
			while( numf[i,2] >= 2,
				pre *= numf[i,1] ;
				num /= numf[i,1]^2 ;
				numf[i,2] -= 2 ;
			) ;
		) ;
	\\ print(num," ",numf," ",numff) ;
	) ;
	if(den != 1,
		denf=factor(den) ;
		denff=(matsize(denf))[1] ;
		for(i=1,denff,
			while( denf[i,2] >= 2,
				pre /= denf[i,1] ;
				den /= denf[i,1]^2 ;
				denf[i,2] -= 2 ;
			) ;
		) ;
	\\ print(den," ",denf," ",denff) ;
	) ;
	print1(sg,"[",truecoe,"]^(1/2)");
	if( pre!=1,
		if( num != den,
			print1(" = ",sg,"(",pre,")[",num/den,"]^(1/2)") ,
			print1(" = ",sg,pre)
		) ;
	) ;
	print() ;
}

/* print a human-readable version of the 6j symbol
* @param a first numerator half-integer
* @param b second numerator half-integer
* @param e third numerator half-integer
* @param d first denominator half-integer
* @param c second denominator half-integer
* @param f third denominator half-integer
* @param coe The 6j coefficient
*/
pr6j(a,b,e,d,c,f,coe)=
{
	local(sg,truecoe) ;
	sg=sign(coe) ;
	truecoe=abs(coe) ;
	if(sg>0,
		print1("{",a," ",b," ",e," ; ",d," ",c," ",f"} = ") ;
		if( truecoe != 1,
			prsqroo("",truecoe),
			print("1") ;
		) ;
		,
		if(coe==0,
			print("{",a," ",b," ",e," ; ",d," ",c," ",f"} = 0") ;
			,
			print1("{",a," ",b," ",e," ; ",d," ",c," ",f"} = ") ;
			if( truecoe != 1,
				prsqroo("-",truecoe),
				print("-1") ;
			) ;
		) ;
	) ;
}

/* test a sum rule on the 6j coefficients.
* See e.g. (4.17) in the book by Zare
*/
tstSum(a,b,c,d,f,g)=
{
	local(tst,tstAl,tstZa,testZaBc1,tstRa,c2,sg,efa) ;
	tst=0 ;
	tstAl=0 ;
	tstZa=0 ;
	tstZaB=0 ;
	tstRa=0 ;
	for(e=-a-b,a+b,
		efa=e ;
		c1 =coupl6j(a,b,e,d,c,f) ;
		c2 = coupl6j(a,b,e,d,c,g) ;
		sg=sign(c1*c2) ;
		c1=abs(c1) ;
		c2=abs(c2) ;
		tst += (2*efa+1)*(2*f+1)*sg*sqrt(c1)*sqrt(c2) ;
		tstAl += (2*efa+1)*(2*f+1)*coupl6jAlt(a,b,e,d,c,f)*coupl6jAlt(a,b,e,d,c,g) ;
		tstZa += (2*efa+1)*(2*f+1)*coupl6jZare(a,b,e,d,c,f)*coupl6jZare(a,b,e,d,c,g) ;
		tstZaB += (2*efa+1)*(2*f+1)*coupl6jZareBase(a,b,e,d,c,f)*coupl6jZareBase(a,b,e,d,c,g) ;

		c1 =coupl6jRasch(a,b,e,d,c,f) ;
		c2 = coupl6jRasch(a,b,e,d,c,g) ;
		sg=sign(c1*c2) ;
		c1=abs(c1) ;
		c2=abs(c2) ;
		tstRa += (2*efa+1)*(2*f+1)*sg*sqrt(c1)*sqrt(c2) ;
		\\print("e ",e," adds ",(2*e+1)*(2*f+1)*sg*sqrt(c1)*sqrt(c2) ) ;
	) ;
	if( f != g,
		print("tst ",tst," tstAl ",tstAl," tstZa ",tstZa," tstZaB ",tstZaB," tstRa ",tstRa," should be 0 ",f," " ,g) ,
		print("tst ",tst," tstAl ",tstAl," tstZa ",tstZa," tstZaB ",tstZaB," tstRa ",tstRa," should be 1 ",f," " ,g)
	) ;
}

/**
* Tabulate 6j coefficients
* @param amax A integer or half integer positive vlue whish is the
*  upper limit of the table for the 3 "numerator" elements of the 6j symbols.
*/
main(amax) =
{
	/* tests of the sum rule.. debugging
	*	tstSum(1,1,1,1,1,1) ;
	*	tstSum(1,1,1,1/2,1/2,1/2) ;
	*	tstSum(1,1,1/2,1/2,1/2,1/2) ;
	*/
	/* print a header line to demonstrate the format
	*/
	print("{a b e ; d c f} = {j1 j2 j3 ; j4 j5 j6} = ...") ;

	/* 6-fold loop over the parameters a to f in 1/2 steps
	*/
	forstep(a=1/2,amax,1/2,
		forstep(b=1/2,a,1/2,
			forstep(e=1/2,b,1/2,
				forstep(d=1/2,a,1/2,
					forstep(c=1/2,b,1/2,
						forstep(f=1/2,e,1/2,
							/* print only non-zero coefficients
							*/
							if( tstDelta(a,b,e,d,c,f) != 0,
								coef2=coupl6j(a,b,e,d,c,f) ;
								pr6j(a,b,e,d,c,f,coef2) ;
								/* test
								* coef2=coupl6jAlt(a,b,e,d,c,f) ;
								* print("or ",coef2) ;
								*/
							) ;
						) ;
					) ;
				) ;
			) ;
		) ;
	) ;

}

main(11/2) ;
