/* Literature on the subject

Willim J. Thompson, <a href="https://doi.org/10.1002/9783527617821">Angular Momentum</a>, Wiley 1994, Eq 7.51

Milton Abramowitz, Irena A. Stegun, Handbook of mathematical functions, chapter 27.9.1

K. Hagiwara et al, <a href="https://doi.org/10.1103/PhysRevD.66.010001">Review of Particle Properities</a>, Phys Rev D 66 (2002) 010001

K Schulten and R G Gordon, <a href="https://doi.org/10.1016/0010-4655(76)90058-8">Recursive evaluation fo 3j and 6j coefficients</a>, Comp Phys Comm 11 (2) (1976) 269-278

<a href="https://en.wikipedia.org/wiki/Table_of_Clebsch-Gordan_coefficients">Wikipedia</a>

P. McNamee, Frank Chilton, <a href="">Tables of Clebsch-Gordan Coefficients of SU3</a>, Rev Mod Phys 36 (1964) 1005
Note that these here are for the SO3 group, not SU3.
*/

/* compute the square of the CG coefficient, 
* including a the sign, so the CG coefficient is the
* sign times the square root of the absolute value of this
* @param j1 first spin orbital momentum number
* @param j2 second spin orbital momentum number
* @param m1 first spin third (magnetic) component
* @param m2 first spin third (magnetic) component
* @param J Coupled spin orbital momentum
* @param M Coupled spin third (magnetic) component
*/
cgord2(j1,j2,m1,m2,J,M) = 
{
	local(k,kmin,kmax,resul,f,sg) ;
	\\print("j1 ",j1," j2 ",j2, " m1 ",m1," m2 ",m2," J ",J," M ",M) ;
	/* there are three trivial cases where the coefficient is
	* zero: magnetic quantum numbers don't match, or the total
	* J is not in the range |j1-j2|<=J<=j1+j2
	*/
	if(m1+m2 != M, return(0) );
	if(J>abs(j1+j2), return(0) );
	if(J<abs(j1-j2), return(0) );
	kmax=j2+J+m1 ;
	kmax=min(kmax,J-j1+j2) ;
	kmax=min(kmax,J+M) ;
	\\print(" kmax ",kmax) ;
	kmin=m1-j1 ;
	kmin=max(kmin,j2+M-j1) ;
	kmin=max(kmin,0) ;
	\\print("kmin ",kmin," 0") ;
	resul=0 ;
	for(k=kmin,kmax,
		/* debugging
		print("j2+J+m1-k ",j2+J+m1-k) ;
		print("j1-m1+k ",j2-m1+k) ;
		print("k ",k) ;
		print("J-j1+j2-k ",J-j1+j2-k) ;
		print("J+M-k ",J+M-k) ;
		print("k+j1-j2-M ",k+j1-j2-M) ;
		*/
		resul +=(-1)^(k+j2+m2)*(j2+J+m1-k)!*(j1-m1+k)!/k!/(J-j1+j2-k)!
			/(J+M-k)!/(k+j1-j2-M)!
	) ;
	\\print("kmin ",kmin," kmax ",kmax," resul ",resul) ; \\ debugging
	if(resul==0,
		return(0) ;
	) ;
	sg=sign(resul) ;
	resul = resul*resul ;
	resul *= 2*J+1 ;
	/* debugging
		print("J+j1-j2 ",J+j1-j2) ;
		print("J-j1+j2 ",J-j1+j2) ;
		print("j1+j2-J ",j1+j2-J) ;
		print("J+M ",J+M) ;
		print("J-M ",J-M) ;
		print("j1+j2+J+1 ",j1+j2+J+1) ;
		print("j1-m1 ",j1-m1) ;
		print("j1+m1 ",j1+m1) ;
		print("j2-m2 ",j2-m2) ;
		print("j2+m2 ",j2+m2) ;
	*/
	f = (J+j1-j2)!*(J-j1+j2)!*(j1+j2-J)!*(J+M)!*(J-M)!/(j1+j2+J+1)!/(j1-m1)!
		/(j1+m1)!/(j2-m2)!/(j2+m2)! ;
	resul *= f ;
	\\print("si ",sg) ;
	return(sg*resul) ;
}

/* Print a human-readable representation of a CG coefficient
* @param j1 first spin orbital momentum number
* @param j2 second spin orbital momentum number
* @param m1 first spin third (magnetic) component
* @param m2 first spin third (magnetic) component
* @param J Coupled spin orbital momentum
* @param M Coupled spin third (magnetic) component
* @param cg The signed square of the CG coefficient
*/
prcg(j1,j2,m1,m2,J,M,cg)=
{
	local(sg,truecg) ;
	sg=sign(cg) ;
	truecg=abs(cg) ;
	/* two branches: sign was positive or zero-negative
	*/
	if(sg>0,
		\\print("<",j1," ",j2,",",m1," ",m2,"|",J," ",M,">= [",cg,"]^(1/2)= ",cg^(1/2)) ;
		/* start of the line are the specs of the individual and total spin
		*/
		print1("<",j1," ",j2," ; ",m1," ",m2," | ",J," ",M,"> = ") ;
		/* if the CG coefficient has not absolute value 1,
		* print a representation as an exact square root
		*/
		if( abs(cg) != 1,
			print("[",cg,"]^(1/2) = ",cg^(1/2)) ,
			print("1 = ",cg^(1/2)) ;
		) ;
		,
		if(cg==0,
			/* case where the total spin is zero
			*/
			print("<",j1," ",j2," ; ",m1," ",m2," | ",J," ",M,"> = ",0) ;
			,
			/* case where the sign was negative
			*/
			\\print("<",j1," ",j2,",",m1," ",m2,"|",J," ",M,">= -[",truecg,"]^(1/2)= -",truecg^(1/2)) ;
			print1("<",j1," ",j2," ; ",m1," ",m2," | ",J," ",M,"> = ") ;
			if( truecg != 1,
				print("-[",truecg,"]^(1/2) = -",truecg^(1/2)) ,
				print("-1 = -",truecg^(1/2)) ;
			) ;
		) ;
	) ;
}

/* print a human-readable table of all results up to some maximum
*/
w3(j1max)=
{
	/* a sort of header that indicates the format of the
	* lines that follow
	*/
	print("<",j1," ",j2," ; ",m1," ",m2," | ",J," ",M,"> =...") ;

	/* outer loop is for first main quantum numbers 1/2 to 9/2.
	* for larger tables set the 9/2 to something larger
	*/
	forstep(j1=1/2,j1max,1/2,
		/* because the coefficients are symmetri w.r.t. swapping j1 <->j2,
		* it only makes sense to let j2 run up to j1.
		*/
		forstep(j2=1/2,j1,1/2,
			/* loops over the magnetic quanutm number give nonzero
			* results only if in the range of plus-minus main quantum number
			*/
			for(m1=-j1,j1,
				for(m2=-j2,j2,
					M=m1+m2 ;
					Jmin=abs(j1-j2) ;
					Jmin=max(Jmin,abs(M)) ;
					for(J=Jmin,abs(j1+j2),
						cg2=cgord2(j1,j2,m1,m2,J,M) ;
						prcg(j1,j2,m1,m2,J,M,cg2) ;
					) ;
				) ;
			) ;
		) ;
	) ;
}


AtlLis()=
{
	forstep(j1=1/2,9/2,1/2,
		forstep(j2=1/2,j1,1/2,
			for(m1=-j1,j1,
				for(m2=-j2,j2,
					M=m1+m2 ;
					Jmin=abs(j1-j2) ;
					Jmin=max(Jmin,abs(M)) ;
					for(J=Jmin,abs(j1+j2),
						cg2=cgord2(j1,j2,m1,m2,J,M) ;
						print1(sqrt(abs(cg2))" Glebsch_Gordan<"j1","j2";"m1","m2"|"J","M">=sqrt("cg2")\n") ;
					) ;
				) ;
			) ;
		) ;
	) ;
}

\\AtlLis for constants
/* print a table with a maximum j1 of 9/2
*/
w3(9/2) ;
