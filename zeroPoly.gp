read("Ek.gp");

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