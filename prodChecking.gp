\\Computational verification of known relationships between Eisenstein sereis and initial evidence of results in paper.

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