\\Testing with modular form functions

\\default(parisize,4000000000)
\\read("/home/trevor/Documents/REU/Ek.gp");

\\Checks if M1*M2 = M3*M4
compare4MF(M1,M2,M3,M4) =
{
   my(E1,E2);
   E1 = mfmul(M1,M2);
   E2 = mfmul(M3,M4);
   mfisequal(E1,E2);
}

\\Checks if M1*M2 = M3
compare3MF(M1,M2,M3) =
{
   my(E1);
   E1 = mfmul(M1,M2);
   mfisequal(E1,M3);
}

\\Checks combinations of two-products of Eisenstien series for weights minW to maxW
checkProd2(minW,maxW) =
{
   maxW = floor(maxW/2);
   print("j\tk\tl\tm\ti\t isEqual");
   my(k,Mj,Mk,Ml,Mm,Mi,totalS,result);
   for(i = minW, maxW,
      for(j = 1, floor(i/2),
	 k = i-j;
	 Mj = mfEk(2*j);
	 Mk = mfEk(2*k);
	 for(l = 1, j-1,
	    m = i-l;
	    Ml = mfEk(2*l);
	    Mm = mfEk(2*m);
	    result = compare4MF(Mj,Mk,Ml,Mm);
	    if(result == 1, totalS++;print(2*j "\t" 2*k "\t" 2*l "\t" 2*m "\t" 2*i "\t" result););
	 );
	 \\Case where m = 0
	 Mi = mfEk(2*i);
	 result = compare3MF(Mj,Mk,Mi);
	 if(result == 1, totalS++;print(2*j "\t" 2*k "\t\t\t" 2*i "\t" result););
      );      
   );
   print("Total: "totalS);
}

\\ Checks all combinations of the three-product (or less) of Eisenstien series for weights minW to max 
checkProd3(minW,maxW) =
{
   \\Range of weights
   minW = floor(minW/2);
   maxW = floor(maxW/2);

   \\Checking if M_1*M_j*M_k = M_l*M_n*M_m for all combinations of i,j,k,l,m,n
   print("j\tk\ti\tl\tn\tm\t Weight\t isEqual");
   my(m,i,Mi,Mj,Mk,Ml,Mm,Mn,totalS,result);
   for(weight = minW, maxW,
      for(j = 1, floor(weight/2),
	 for(k = 1, weight-j,
	    i = weight-j-k;
	    Mi = mfEk(2*i);
	    Mj = mfEk(2*j);
	    Mk = mfEk(2*k);
	    Mk = mfmul(Mk,Mj);
	    for(l = 1, j-1,
	       for(n = 1, weight-l,  
		  m = weight-n;
		  Ml = mfEk(2*l);
		  Mm = mfEk(2*m);
		  Mn = mfEk(2*n);
		  Mn = mfmul(Mn,Ml);
		  if((i == 0) && (m != 0), result = compare3MF(Mk,Mn,Mm), if((i != 0) && (m ==0), result = compare3MF(Mk,Mn,Mi), if((m == 0) && (i == 0), result = mfisequal(Mn,Mk), result = compare4MF(Mk,Mi,Mn,Mm););););
		  if(result == 1, totalS++;print(2*j "\t" 2*k "\t" 2*i "\t" 2*l "\t" 2*n "\t" 2*m "\t" 2*weight "\t" result););
	       );
	    );
	 );
      );
   );
   print("Total: "totalS);
}

\\Computes the rational part of zeta(n)
Zeta(n) =
{
   my(result);
   n = n/2;
   result = (-1)^(n+1)*bernfrac(2*n)*2^(2*n)/(2*(2*n)!);
   result;
}

\\Writing a weight k Eisenstien series as a linear combination of E_4E_6 (See Wikipedia)
ek(k) =
{
   if(k%2 != 0,
      print("ERROR: not an even weight\n");
      return(1);
     );
   
   if(mapisdefined(Ek,k),
      return(mapget(Ek,k)),
      k = (k-4)/2;
      return(dk(k)/((2*k+3)*k!*2*Zeta(2*k+4)));
     );
}

dk(k) =
{
   my(c,n,s);
   \\Terminating cases
   if(mapisdefined(Ek,2*k+4),
      return((2*k+3)*k!*2*Zeta(2*k+4)*mapget(Ek,2*k+4));,
      \\Recurrence
      n = k-2;
      c = (3*n+6)/(2*n+9);
      s = sum(i=0,n, binomial(n,i)*dk(i)*dk(n-i));
      return(c*s);
     );
}

gk(k) =
{
   my(n,c,s);
   if(mapisdefined(Ek,k), return(2*Zeta(k)*mapget(Ek,k)));
   if(k%6 != 2, return(2*Zeta(k)*ek(k)));
   n = (k-2)/6;
   c = ((4*n+1)!)/((6*n+1)*((2*n)!^2));
   s = sum(i=1, n,
	     binomial(2*n,2*i-1)/binomial(6*n,2*n+2*i-1)*gk(2*n+2*i)*gk(4*n-2*i+2)
       );
   return(c*s);
}

faster_ek(k) =
{
   return(gk(k)/(2*Zeta(k)));
}

\\Computes the dimension of the Eisenstien subspace
dimk(k) =
{
   my(result);
   if(k%12 == 2, result = floor(k/12),
      result = floor(k/12)+1);
   result;
}

\\Returns the coordinate vector (row) of E_k wrt E_4E_6
ekcoef(k) =
{
   my(B1,coef,d,k12,a,b);
   if(mapisdefined(Ek,k), coef = mapget(Ek,k),
      coef = faster_ek(k);
      );
   d = dimk(k);
   \\Need to remove the extra E4 and E6 terms
   \\ 12c+4a+6b = k a<= 2 and b <=1;
   k12 = k%12;
   if(k12 == 0, a=0;b=0;,
      if(k12 == 2, a=2;b=1;,
	 if(k12 == 4, a=1;b=0;,
	    if(k12 == 6, a=0;b=1;,
	       if(k12 == 8, a=2;b=0;,
		  if(k12 == 10, a=1;b=1;
		    );
		 );
	      );
	   );
	);
     );
   coef = coef/(E4^a*E6^b);
   coef = subst(coef,E6,1);
   B1 = vector(d,i,polcoef(coef,3*(d-i)));
   B1;
}

addtoEk(k) =
{
   my(F,image);
   if(mapisdefined(Ek,k) == 0,
      image = faster_ek(k);
      mapput(Ek,k,image);
      F = fileopen("Ek.gp","w");
      filewrite(F,"Ek = "Ek";");
      fileclose(F);
      );
}

\\Make the Change of basis matrix A
ekA(k) =
{
   my(A,d,N);
   N = 1728;
   d = dimk(k);
   A = matrix(d,d,i,j,
	      if(i>j, 0,
		 (-1)^(i+1)*N^(i-1)*binomial(j-1,i-1)
		);
	     );
   
}

\\Compute the coordinate vector (column) of E_k wrt E_4*\delta
ZeroPoly(k) =
{
   my(B2);
   B2 = ekA(k)*ekcoef(k)~;
   B2;
}

\\mod p irreducibility test
irre(c) =
{
   my(result,f,temp);
   forprime(p = 2, 1000000,
	    result = true;
	    temp = c%p;
	    f = Pol(temp);
	    for(q = 0, p-1,
	       if(subst(f,x,q)%p == 0, result = false; break;);
	    );
	    if(result == true,
	       if(p != 999983,
		  print("p: "p);
		  break,
		  result == false;
		  print("ERROR: Reached max prime");
		  );
	      );
	   );
}

\\Checking the zero polynomial for irreducibility
zeroisirre(k) =
{
   my(c,f,N,temp);
   N = numerator(bernfrac(k)/(2*k));
   c = ZeroPoly(k);
   c = c*N;
   if(irre(c) == false, print("\Phi_k is Reducible"));
}

Amu(k,mu) =
{
   sin(2*Pi*(ceil(k/4)+mu)/k);
}   

\\Remark: if mu == 1, Amu-1 terms have to be ignored (t1,t4 ish) (
imaginaryzeroEk(k,mu) =
{
   \\should have 1 <= mu <= N
   my(Theta,k0,N,s,gamk,del,t1,t2,t3,ct4,t4);
   if(k%12 != 2,
      N = (k-(k%12))/12;
      s = k%12;,
      N = (k-14)/12;
      s = 14;
      );
   k0 = ceil(k/4);
   if(s == 6 || s == 10 || s == 14, gamk = 1/2, gamk = 0);
   if(mu == 1, del = 1, del =0);
   if(mu != 1,
      t2 = -1*(mu-1+gamk)*Amu(k,mu-1);
      t4 = intnum(x = -1/2,1/2, log(abs(mfeval(mfinit(mfEk(k)),mfEk(k),x+I*Amu(k,mu))))-log(abs(mfeval(mfinit(mfEk(k)),mfEk(k),x+I*Amu(k,mu-1)))));,
      t2 = 0;
      t4 = intnum(x = -1/2,1/2, log(abs(mfeval(mfinit(mfEk(k)),mfEk(k),x+I*Amu(k,mu))));
		 );
      );
   t1 = (mu+gamk)*Amu(k,mu);
   t3 = -1*gamk*del;
   ct4 = 1/(4*Pi);
   

   Theta = t1+t2+t3+(ct4*t4);
   
   Theta = asin(Theta);
   Theta = Pi-Theta;
   Theta = Theta/(2*Pi);
   Theta;
}


\\Computes all zeros of Ek or if none-zero the mu-th component
zeroEk(k,{mu = 0}) =
{
   my(N,c,z);
   c = 1;
      
   if(mu == 0,
      if((k)%12 != 2,
	 N = (k-((k)%12))/12;,
	 N = (k-14)/12;,
	);,
      N = mu;
      c = mu;
     );
   if(mu == 0,
      z = vector(N,i,exp(2*Pi*I*imaginaryzeroEk(k,i)));,
      z = vector(1,i,exp(2*Pi*I*imaginaryzeroEk(k,mu)));
     );
   z;
}

temp() =
{
   for(i = 1,10,
      print(ellj(zeroEk(12*i+2,1)[1]));
      \\print(ellj(zeroEk(12*i+4,1)[1]));
      print(ellj(zeroEk(12*i+6,1)[1]));
      \\print(ellj(zeroEk(12*i+8,1)[1]));
      print(ellj(zeroEk(12*i+10,1)[1]));
   );
   
}

\\mu == 1
integralPart(k) =
{
   my(t4);
   t4 = 1/(4*Pi)*intnum(x = -1/2,1/2, log(abs(mfeval(mfinit(mfEk(k)),mfEk(k),x+I*Amu(k,1))));
	       );
   t4;
}

nuk(k,{mu = 0}) =
{
   my(N,z);
   
   if(mu == 0 ,
      if((k)%12 != 2,
	 N = (k-((k)%12))/12;,
	 N = (k-14)/12;
	 );
	 z = vector(N,i,sin(2*Pi*imaginaryzeroEk(k,i)));,
	 z = vector(1,i,sin(2*Pi*imaginaryzeroEk(k,mu)));
     );
   z;
}

temp2(k,{numc = 10}) =
{
   my(f,g);
   f = mfcoefs(mfEk(k),numc);
   g = mfcoefs(mfEk(k+4),numc);
   for(i = 1,numc, print(f[i] > g[i]));

}
