default(parisize,4000000000)
COS = vector(83,i,read(Strexpand("data/Cos/"i"_zeros.gp"))); \\ Map of Maps defining zero locations for cos and Ek
EKZEROS = vector(83,i,read(Strexpand("data/Ek/"i"_zeros.gp")));

/*
* Functions for populating COS and EKZEROS
*/

addtoCOS(k) =
{
   my(N,F,image);
   N = numZeros(k);
   if(!mapisdefined(COS,k),
      image = vector(N,i,alphakn(k,i));
      mapput(COS,k,image);
      F = fileopen("cos.gp","w");
      filewrite(F,"COS = "COS";");
      fileclose(F);
      );
}

addtoEKZEROS(k) =
{
   my(F,image);
   if(!mapisdefined(EKZEROS,k),
      image = zeroEk(k,1,1);
      mapput(EKZEROS,k,image);
      F = fileopen("ekzeros.gp","w");
      filewrite(F,"EKZEROS = "EKZEROS";");
      fileclose(F);
      );
}

deletefromEKZEROS(k) =
{
   my(F);
   if(mapisdefined(EKZEROS,k),
      mapdelete(EKZEROS,k);
      F = fileopen("ekzeros.gp","w");
      filewrite(F,"EKZEROS = "EKZEROS";");
      fileclose(F);
      );
}

/*
* Computes all zeros of Ek or the mu-th zero
* Run with mu = 1 to recreate data for Tables 1 and 2 "Monomials of EIsenstien Series"
* set flag == 1 to return vector of theta values instead of zeros
*/
zeroEk(k,{mu = 0},{flag = 0}) =
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
      if(flag == 0,
	 z = vector(N,i,exp(I*faster_theta(k,i)));,
	 z = vector(N,i,faster_theta(k,i));
	);,
      if(flag == 0,
	 z = vector(1,i,exp(I*faster_theta(k,mu)));,
	 z = vector(1,i,faster_theta(k,mu));
	);
     );
   z;
}

\\Can be used to recreate data for Tables 3 and 4 in "Monomials of Eisenstien Series"
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

alphadiff(k,{mu = 1}) =
{
   my(res);
   res = mapget(EKZEROS,k)[mu] - mapget(COS,k)[mu];
   res;
}

numZeros(k) =
{
   my(N);
   if((k)%12 != 2,
      N = (k-((k)%12))/12;,
      N = (k-14)/12;,
     );
   N;
}

/*
* The Bisection Method for finding Ek zeros (fast)
*/

Fk(k,ang) =
{
   my(res);
   res = 2*cos(k*ang/2) + (1/(2*cos(ang/2)))^k + (1/(2*sin(ang/2)))^k;
   res;
}

alphakn(k,{j=1}) =
{
   my(res);
   if(k%4 == 0,
      res = (1/2 + ((2*j-1)/k))*Pi,
      res = (1/2 + ((2*j)/k))*Pi
     );
   res;
}

faster_theta(k,{mu = 1},{loops = 1000}) = 
{
   my(res,a1,a2,temp,loop);
   \\Initialize the bounds
   a1 = alphakn(k,mu);
   a2 = a1 + (6*Pi/(k*(k+12)));
   a1 = a1 - (6*Pi/(k*(k+12)));
   loop = 1;
   while(loop <= loops,
	 temp = (a2+a1)/2;
	 if(Fk(k,(a1))*Fk(k,temp) < 0,
	    a2 = temp,
	    a1 = temp;
	   );
	 loop++;
	);  
   res = (a2+a1)/2;
   res;
}

addEkData(filename,k,{mu = 1}) =
{
   my(F,image,mapname);
   iferr(mapname = read("data/"filename".gp"),E,mapname = Map()); \\Case of non-existent file
   if(type(mapname) == "t_INT",mapname = Map()); \\Case of empty file
   if(!mapisdefined(mapname,k),
      mapput(mapname,k,zeroEk(k,mu,1));
      F = fileopen("data/"filename".gp","w");
      filewrite(F,mapname);
      fileclose(F);
      );
}

addCosData(filename,k,{mu = 1}) =
{
   my(F,image,mapname);
   iferr(mapname = read("data/"filename".gp"),E,mapname = Map()); \\Case of non-existent file
   if(type(mapname) == "t_INT",mapname = Map()); \\Case of empty file
   if(!mapisdefined(mapname,k),
      mapput(mapname,k,alphakn(k,mu));
      F = fileopen("data/"filename".gp","w");
      filewrite(F,mapname);
      fileclose(F);
      );
}

/*
* Functions for verifying the interlacing of zeros 
 */
checkInterlacingCasesCOS(k) =
{
   interlacingCOS(k,2);
   interlacingCOS(k,4);
   interlacingCOS(k,6);
   interlacingCOS(k,8);
   interlacingCOS(k,10);
   interlacingCOS(k,12);
   interlacingCOS(k,14);
   interlacingCOS(k,16);
   interlacingCOS(k,18);
   interlacingCOS(k,22);
}

checkInterlacingCasesEK(k) =
{
   interlacingEK(k,2);
   interlacingEK(k,4);
   interlacingEK(k,6);
   interlacingEK(k,8);
   interlacingEK(k,10);
   interlacingEK(k,12);
   interlacingEK(k,14);
   interlacingEK(k,16);
   interlacingEK(k,18);
   interlacingEK(k,22);
}

interlacingCOS(k,a) = \\ Creates a table to manually check for interlacing WARNING: table directory not included in repo, must manually add
{
   if(numZeros(k+a)-numZeros(k) > 1, return());
   my(N,z1,z2);
   N = numZeros(k+a);
   write(Strexpand("data/tables/interlacing/cos/k+"a".txt"),"\nFor "k" and "k+a);
   write(Strexpand("data/tables/interlacing/cos/k+"a".txt"),Strexpand("i-th zero \ta_(i,k+"a") \ta_(i,k) \tdifference"));
   for(i = 1,N,
      if(mapisdefined(COS[i],k), z1 = mapget(COS[i],k), z1 = 0);
      if(mapisdefined(COS[i],k+a), z2 = mapget(COS[i],k+a), z2 = 0);
      write(Strexpand("data/tables/interlacing/cos/k+"a".txt"),i"\t\t"z2"\t"z1"\t"z1-z2); 
   );
}

interlacingEK(k,a) = \\ Creates a table to manually check for interlacing
{
   if(numZeros(k+a)-numZeros(k) > 1, return());
   my(N,z1,z2);
   N = numZeros(k+a);
   write(Strexpand("data/tables/interlacing/Ek/k+"a".txt"),"\nFor "k" and "k+a);
   write(Strexpand("data/tables/interlacing/Ek/k+"a".txt"),Strexpand("i-th zero \ta_(i,k+"a") \ta_(i,k) \tdifference"));
   for(i = 1,N,
      if(mapisdefined(EKZEROS[i],k), z1 = mapget(EKZEROS[i],k), z1 = 0);
      if(mapisdefined(EKZEROS[i],k+a), z2 = mapget(EKZEROS[i],k+a), z2 = 0);
      write(Strexpand("data/tables/interlacing/Ek/k+"a".txt"),i"\t\t"z2"\t"z1"\t"z1-z2); 
   );
}

isInterlacingCOS(k,a) = \\ automatically checks for interlacing
{
   if(numZeros(k) == 0, return(1));
   if(numZeros(k+a)-numZeros(k) > 1, return(-1));
   my(N,z1,z2,SIGN);
   N = numZeros(k+a);
   for(i = 1,N,
      if((i == N) && (numZeros(k+a) > numZeros(k)),
	 if(mapisdefined(COS[i-1],k), z1 = mapget(COS[i-1],k), z1 = 0);
	 if(mapisdefined(COS[i],k+a), z2 = mapget(COS[i],k+a), z2 = 0);
	 if(z1-z2 > 0, return(-1), return(1));
	);	
      if(mapisdefined(COS[i],k), z1 = mapget(COS[i],k), z1 = 0);
      if(mapisdefined(COS[i],k+a), z2 = mapget(COS[i],k+a), z2 = 0);
      if(i == 1, SIGN = sign(z1-z2));
      if(SIGN == -1 && numZeros(k+a) > numZeros(k), return(-1));
      if(SIGN != sign(z1-z2), return(-1));
   );
   return(1);
}

isInterlacingCOS(k,a) = \\ automatically checks for interlacing
{
   if(numZeros(k) == 0, return(1));
   if(numZeros(k+a)-numZeros(k) > 1, return(-1));
   my(N,z1,z2,SIGN);
   N = numZeros(k+a);
   for(i = 1,N,
      if((i == N) && (numZeros(k+a) > numZeros(k)),
	 if(mapisdefined(COS[i-1],k), z1 = mapget(COS[i-1],k), z1 = 0);
	 if(mapisdefined(COS[i],k+a), z2 = mapget(COS[i],k+a), z2 = 0);
	 if(z1-z2 > 0, return(-1), return(1));
	);	
      if(mapisdefined(COS[i],k), z1 = mapget(COS[i],k), z1 = 0);
      if(mapisdefined(COS[i],k+a), z2 = mapget(COS[i],k+a), z2 = 0);
      if(i == 1, SIGN = sign(z1-z2));
      if(SIGN == -1 && numZeros(k+a) > numZeros(k), return(-1));
      if(SIGN != sign(z1-z2), return(-1));
   );
   return(1);
}

compareZeros(k,zero) = 
{
   my(v,res1,res2);
   res1 = mfeval(mfinit(mfEk(k)),mfEk(k),exp(I*mapget(EKZEROS[1],k)[1]));
   res2 =  mfeval(mfinit(E),E,exp(I*zero));
   v = vector(2);
   v[1] = res1;
   v[2] = res2;
   v;
}

/*
* OBSOLETE FUNCTIONS: Old functions used for implementing Kohnen's method of finding zeros of Ek. Replaced by Bisection method which was extracted from Nozaki's paper.
* this method is more accurate for small k. Use "compareZeros" for comparing methods.
*/

\\Returns A_\mu as defined in Kohnen
Amu(k,mu) =
{
   sin(2*Pi*(ceil(k/4)+mu)/k);
}   

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
