/*
 * =============================================================
 * findimps3D.c 
  
 * takes two input matrices X1,X2 ond two vectors t1,t2
 * 
 * equivalent to: 
 *     Dist=distance(X1.'*X2);
 *     imp1=find(Dist<repmat(t1,N1,1))';
 *     imp2=find(Dist<repmat(t2,N1,1))';
 *     [a,b]=ind2sub([N1,N2],[imp1;imp2]);
 *
  * =============================================================
 */


/* $Revision: 1.3 $ */

#include "mex.h"

/* If you are using a compiler that equates NaN to zero, you must
 * compile this example using the flag -DNAN_EQUALS_ZERO. For 
 * example:
 *
 *     mex -DNAN_EQUALS_ZERO findnz.c  
 *
 * This will correctly define the IsNonZero macro for your
   compiler. */

#if NAN_EQUALS_ZERO
#define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d) != 0.0)
#endif


double square(double x) { return(x*x);}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  /* Declare variables. */ 
  int N1,N2, o1,o2,i1,i2;
  int c,r,size1,size2;
  double *po,*X1,*X2,*C,*sx1,*sx2,*Thresh1,*Thresh2;
  int m,p,n, oi;
  char *chn="N"; 
  char *chn2="T";
  double minustwo=-2.0,one=1.0, zero=0.0;



  /* Check for proper number of input and output arguments. */    
  if (nrhs != 4)
    mexErrMsgTxt("Four input arguments required.");
   

  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }

  /* Check data type of input argument. */
  if (!(mxIsDouble(prhs[0]))) {
   mexErrMsgTxt("Input array must be of type double.");
  }

    
  /* Get the number of elements in the input argument. */
/*  N1 = mxGetNumberOfElements(prhs[0]);/*


  /* Get the data. */
  X1  = mxGetPr(prhs[0]);
  X2  = mxGetPr(prhs[1]);
  Thresh1=mxGetPr(prhs[2]);
  Thresh2=mxGetPr(prhs[3]);

  p = mxGetM(prhs[0]);
  m = mxGetN(prhs[0]);
  n = mxGetN(prhs[1]);
  
 if(p!=mxGetM(prhs[1])) mexErrMsgTxt("Inner dimensions must agree!\n");
 if(m!=mxGetN(prhs[2])) mexErrMsgTxt("Threshold1 must be of same length as first input has columns!\n");
 if(n!=mxGetN(prhs[3])) mexErrMsgTxt("Threshold2 must be of same length as second input has columns!\n");
  
  C=(double *)mxCalloc(m*n*2,sizeof(double));

  dgemm_ (chn2,chn,&m, &n, &p, &minustwo, X1, &p, X2, &p, &zero, C, &m);

  sx1=mxCalloc(m,sizeof(double));
  i1=0;i2=0;
  for(c=0;c<m;c++){
   for(r=0;r<p;r++) {
     sx1[i2]=sx1[i2]+square(X1[i1]);   
     i1=i1+1;
    }
   i2=i2+1;
  }

  sx2=mxCalloc(n,sizeof(double));
  i1=0;i2=0;
  for(c=0;c<n;c++){
   for(r=0;r<p;r++) {
     sx2[i2]=sx2[i2]+square(X2[i1]);   
     i1=i1+1;
    }
   i2=i2+1;
  }

  plhs[0]= mxCreateDoubleMatrix(2,m*n,mxREAL);
  po=mxGetPr(plhs[0]);
  
  i1=0;oi=0;
 for(c=0;c<n;c++) {
   for(r=0;r<m;r++){
    if(C[i1]+sx1[r]+sx2[c]<Thresh2[c] || C[i1]+sx1[r]+sx2[c]<Thresh1[r]) {po[oi]=r+1;po[oi+1]=c+1;oi=oi+2;};
    i1=i1+1;
   }
  }
 /* 
  mxSetN(plhs[0],ceil(oi/2));  
 */

}



