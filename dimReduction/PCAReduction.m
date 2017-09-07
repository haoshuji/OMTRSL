function X_PCA = PCAReduction(X, variance_to_achieve,file_name)
    
%%
% X: M*N, M: number of instance, N: number of dimensions
    % variance_to_achieve, eg. 90%,95%
    %PCA01.m
    %PCA01 Do PCA on a set of input data vectors using Matlab's princomp().
    %Works well on signals with near zero mean (test files EMG_set1.txt, 
    %EMG_set1L1.txt, EMG_set1L2.txt; see Excel file).
    %Fails to properly reconstruct non-zero mean signals (e.g. EMG_set2.txt). 
    %Matlab's PCA function subtracts mean from each trace first.  
    %Inspection  of reconstruction plots suggests that saving 
    %and later adding back mean of each trace to the reconstructed signal will 
    %not be sufficient to get a good reconstruction. 
    %W. Rose 10/22/2006

%     clear all; close all;

    %get the data
%     file_name = input('Enter name of text file containing data (.txt will be added): ','s');
%     X=load([file_name,'.txt']);  
    [M,N]=size(X);
    fprintf('%d Instances; %d Dimensions\n',M,N);

    %% analyze the data
%     [coefs,scores,variances] = princomp(X); 
    [coefs,scores,variances] = pca(X); 
        %coefs=[NxN]. Columns of coefs=PCs.
        %scores=[MxN].  Each row of scores=PC weights for one data record.
        %Col j of scores = weights of jth PC for all data records.
        %Varianes=[Nx1]=variance acctd for by each PC.
    pc_exp = 100*variances/sum(variances); %percent of variance explained by each PC
    %find k according to 90%, 95%, and 1% criteria
    pc_exp_cum(1)=pc_exp(1);
    for i=2:N
        pc_exp_cum(i)=pc_exp_cum(i-1)+pc_exp(i);
    end;
    i=1;
    while pc_exp_cum(i)<variance_to_achieve
        i=i+1;
    end
    kVariance=i;
    fprintf('%d dimension to keep\n',kVariance);
    X_PCA = X*coefs(:,1:kVariance);
    
%     i=1;
%     while pc_exp(i)>1
%         i=i+1;
%     end
%     k01=i-1;
%     kmax=max([kVariance, k01]);
%    
%     %flip the "upside down" PCs and their scores
% %     for j=1:N
% %         if max(coefs(:,j))<abs(min(coefs(:,j)))
% %             coefs(:,j)=-coefs(:,j);
% %             scores(:,j)=-scores(:,j);
% %         end;
% %     end;
%    
% %% reconstruct the data using PCs 1 thru k90
%     xrecon=zeros(3,N);
%     for j=1:3
%         for k=1:kVariance
%             xrecon(j,:)=xrecon(j,:)+scores(j,k)*coefs(:,k)';
%         end;
%     end;
% 
%     %% plot and print some results
%     fprintf('k(%f %%)=%d, k(1%%)=%d\n',variance_to_achieve, kVariance, k01);
%     fprintf('Percent of variance explained by PCs 1-%d:\n',kmax+2);
%     for i=1:kmax+2; fprintf('%3.1f  ',pc_exp(i)); end;
%     fprintf('\n');
%     
%     t=(0:N-1);
%    
%     h = figure;
%     plot(t,coefs(:,1),'-ro',t,coefs(:,2),'-gx',t,coefs(:,3),'-b+',t,coefs(:,4),'-c^',t,coefs(:,5),'-mv',t,coefs(:,6),'-y.');
%     title(['PRINCIPAL COMPONENTS: input=',file_name,'.txt']);
%     xlabel('TIME'); ylabel('AMPLITUDE'); grid on; 
%     legend('PC1','PC2','PC3','PC4','PC5','PC6');
%     print(h,'-depsc',strcat(file_name,'_Coefficients.eps'));    
%     saveas(h,strcat(file_name,'-Coefficients.fig'));
%     close(h);
%     
%     h = figure;
%     plot(t,X(1,:),'ro',t,xrecon(1,:),'-r',t,X(2,:),'gx',t,xrecon(2,:),'-g',t,X(3,:),'bs',t,xrecon(3,:),'-b');
%     title(['PCA Reconstruction using PC1 thru PC',num2str(kVariance),': input=',file_name,'.mat']);
%     xlabel('TIME'); ylabel('AMPLITUDE'); grid on; 
%     legend('X1','X1r','X2','X2r','X3','X3r');
%     print(h,'-depsc',strcat(file_name,'_Reconstruction.eps'));    
%     saveas(h,strcat(file_name,'_Reconstruction.fig'));
%     close(h);
%     
%     h = figure;
%     pareto(pc_exp);
%     grid on; axis([0.5 kmax+1.5 0 100]);
%     title(['PCA: input=',file_name,'.mat']);
%     xlabel('Principal Component');
%     ylabel('Variance Explained (%)');
%     print(h,'-depsc',strcat(file_name,'_Variance.eps'));    
%     saveas(h,strcat(file_name,'_Variance.fig'));
%     close(h);
% 
%     %save results to disk
%     outfilename=[file_name,'_PCA'];
%     save(outfilename,'kVariance','k01','coefs','scores','variances');
%     fprintf(['Variables saved in ',outfilename,'.mat.\n']);

    
