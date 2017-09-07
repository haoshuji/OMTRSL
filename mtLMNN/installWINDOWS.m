addpath(pwd);
addpath([pwd '/mexfunctions']);
addpath([pwd '/helperfunctions']);
cd mexfunctions
fprintf('Compiling all mex functions ...\n');
mex SOD.c
mex count.c
mex mulv.c
mex SODE2.c
mex findimps3.c
mex sd.c
mex SODd.c
mex findimps3D.c
mex sumiflessh2.c
mex addh.c
mex findlessh.c
mex sumiflessv2.c
mex addv.c
mex mink.c
mex cdist.c
mex mulh.c
mex findimps3Dac.c
mex SODW.c
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





