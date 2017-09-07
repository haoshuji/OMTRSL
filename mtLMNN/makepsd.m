function [M,L]=makepsd(Q);
    % decompose Q
    [L,dd]=eig(Q);
    dd=real(diag(dd));
    L=real(L);
    % reassemble Q (ignore negative eigenvalues)
    j=find(dd<1e-10);
    dd(j)=0;
    [temp,ii]=sort(-dd);
    L=L(:,ii);
    dd=dd(ii);
    
    L=(L*diag(sqrt(dd)))';
    M=L'*L;

