\\default(parisize,4000000000) \\Uncomment to allocate necessary RAM, only needed to be called once per session

\\Returns A_\mu as defined in Kohnen
Amu(k,mu) =
{
   sin(2*Pi*(ceil(k/4)+mu)/k);
}   

\\Returns \theta_mu as defined in Kohnen
\\Run with mu=1 to recreate data for Tables 1 and 2
theta_mu(k,mu) =
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
   
   Theta = (Pi-asin(t1+t2+t3+(ct4*t4)))/(2*Pi); 
   Theta;
}


\\Computes all zeros of Ek or the mu-th zero
\\Run with mu = 1 to recreate data for Tables 1 and 2 (returns vector instead of arr element)
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
      z = vector(N,i,exp(2*Pi*I*theta_mu(k,i)));,
      z = vector(1,i,exp(2*Pi*I*theta_mu(k,mu)));
     );
   z;
}

alphakn(k,n) =
{
   my(res);
   if(k%4 == 0,
      res = (1/2 + (2*n-1)/k)*Pi;,
      if(k%4 == 2,
	 res = (1/2 + (2*n)/k)*Pi;
	);
      );
   res;	
}


\\Can be used to recreate data for Tables 3 and 4
range(k) =
{
   my(t,c,z);
   c = 1;
   t = Pi/(c*k^2);
   z = vector(2);
   z[1] = alphakn(k,1)-t;
   z[2] = alphakn(k,1)+t;
   z;

}