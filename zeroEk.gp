\\Testing with modular form functions

\\default(parisize,4000000000)

\\Returns A_\mu as defined in Kohnen
Amu(k,mu) =
{
   sin(2*Pi*(ceil(k/4)+mu)/k);
}   

\\Returns \theta_mu as defined in Kohnen
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
