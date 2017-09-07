function [L,dd]=decompose(M);
% function [L,dd]=decompose(M);
%
% decomposees the (positive semidefinite) matrix M into L and dd such that M=L'*L
% where the ROWS of L are the eigenvectors of M
% sorted with decreasing eigenvalues (which are stored in dd)
% 
% copyright by Kilian Q. Weinberger, 2006
%

[L,dd]=eig(M);
L=real(transpose(L*sqrt(max(dd,0))));
dd=diag(dd);
[temp,ii]=sort(-dd);
L=L(ii,:);
dd=dd(ii);
