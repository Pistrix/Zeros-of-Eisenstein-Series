read("lfundata.gp"); \\Dirichlet coefficients of L-functions and epsilon values taken from the LMFDB

\\Sanity check, equivalent to EVAL(12,l)
Delta_case(k,l) =
{
   my(L,res);
   L = lfunmf(mfinit(mfDelta()),mfDelta());
   res = (l*(k-l)!*lfun(L,k-l))/(bernfrac(l)*bernfrac(k-l)*(2*Pi*I)^(k-l));
   res;
}

\\Expression in final theorem, defined for k = 12,16,18,20,22,26
EVAL(k1,l) =
{
   my(L,res,an,eps,k);
   an = mapget(Dirichlet,k1);
   k = floor(k1);
   eps = mapget(EPS,k);
   L = lfuncreate([an,1,[0,1],k,1,eps]);
   res = (l*(k-l)!*lfun(L,k-l))/(bernfrac(l)*bernfrac(k-l)*(2*Pi*I)^(k-l));
   res;
}

\\Check all possible l1,l2 for a given k
checkEval(k) =
{
   cases = 0;
   for(l1 = 4,floor(k/2),
      for(l2 = l1+2,k-l1-2,
	 print("& " "("l1","l2") & (" EVAL(k,l1) ", "EVAL(k,l2)") \\\\" );
	 l2++;
      );
      l1++;
   );
}

\\Recreates data for Table 6
gammal(k1) =
{
   my(k);
   k = floor(k1);
   for(l = 4,k-6,
      print("& $"l"$ & $"EVAL(k1,l)"$ \\\\"
	   );
      l++;
   );
}