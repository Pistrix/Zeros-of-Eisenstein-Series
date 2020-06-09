# Zeros-of-Eisenstein-Series
Computations of Eisenstein Series and Related Objects.

All scripts are Pari/GP scripts, see https://pari.math.u-bordeaux.fr/pub/pari/manuals/2.11.1/users.pdf for more details on scripting with Pari/GP.

The repository contains six files:
  - this README
  - prodChecking.gp
  - zeroPoly.gp
  - zeroEk.gp
  - eigenform.gp
  - lfundata.gp
  
prodChecking.gp uses some basic functions of the modular forms package to verify the known relationships between 1-dimensional spaces of Eisenstein series. Included also are some functions for testing if the three product of Eisenstein series ever equates to another 3 or less product. This data was motivating evidence towards the theorem proved in the main paper.

zeroPoly.gp focuses on the zero polynomial of Eisentein series. Included are functions for computing the zero polynomial via a change of basis on the representation of Ek by E4^3 and E6^2. These finding were not dicussed in the main paper but discussion on them can be found at the following: https://jstor.org/stable/23807544. Furthermore, one can find function for computing the coefficients of Ek when expressed in terms of E4 and E6. The algorithm automatically stores results in Ek.gp for reducing redundancy. Improvements on computing the representation Ek were made by implementing the recurrence relation given in:

LACUNARY RECURRENCES FOR EISENSTEIN SERIES by MICHAEL H. MERTENS and LARRY ROLEN.

zeroEk.gp contains the functions necessary for reproducing the results given at the end of the paper. When recomputing these zeroes, it is reccomended to increase the precision of Pari/GP to at least 100, for accurate results. 

eigenform.gp houses functions needed for computing the L-values in Table 6. EVAL does the actual computation while the other functions print the results in latex format. checkEval can be used to see all possible cases of l_1 and l_2, however, it does not exclude the cases arising from trivial identites. 

lfundata.gp stores data pulled from the LMFDB needed to run the computations of eigenform.gp, i.e Dirichlet coefficients of normalized eigenforms of a given weight.

Further details on the functions and running them can be found in the comments of the code.

WARNING: Functions to not check weight. User might get unexpected results from odd or small weights.
