addpath(pwd);
addpath([pwd '/mexfunctions']);
addpath([pwd '/helperfunctions']);
cd mexfunctions
fprintf('Compiling all mex functions ...\n');
mex SOD.c  -lmwlapack 
mex count.c  -lmwlapack 
mex mulv.c  -lmwlapack 
mex SODE2.c	-lmwlapack -lmwblas
mex findimps3.c	 -lmwblas -lmwlapack
mex sd.c  -lmwlapack 
mex SODd.c		  -lmwlapack 
mex findimps3D.c	  -lmwlapack -lmwblas
mex sumiflessh2.c  -lmwlapack 
mex addh.c		  -lmwlapack 
mex findlessh.c	  -lmwlapack 
mex sumiflessv2.c  -lmwlapack 
mex addv.c		  -lmwlapack 
mex mink.c  -lmwlapack 
mex cdist.c		  -lmwlapack 
mex mulh.c  -lmwlapack 
mex findimps3Dac.c  -lmwlapack
mex SODW.c -lmwlapack
addpath(pwd);
cd ..

cd mtrees
mex buildmtreec.cpp
mex findknnmtree.cpp  
mex findNimtree.cpp  
addpath(pwd);
cd ..


%%mex distance.c  %% de-comment this line for maximum speed on 64 processors
fprintf('done\n\n');





