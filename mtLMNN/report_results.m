function report_results(errs,errs_te)
close all

% All validation Errors 
teErr = errs(2:2:end,:);
% Train Errors
trErr = errs(1:2:end,:);
% Test Errors
teErr2 = errs_te(2:2:end,:);


M0te = teErr(3:2:end,:);
M0te2 = teErr2(3:2:end,:);
M0tr = trErr(3:2:end,:);

MLte = teErr(4:2:end,:);
MLte2 = teErr2(4:2:end,:);
MLtr = trErr(4:2:end,:);

% fprintf('Best of %d  of parameters \n',size(M0te,1) ); 

% Pick parameter that gave minimum mean error
[minAvg argAvg] = min(mean(MLte,2));
[minAvg0 argAvg0] = min(mean(M0te,2));
m1 = mean(MLte2,2);
minAvg2 = m1(argAvg);
m2 = mean(M0te2,2);
minAvg02 = m2(argAvg0);

% Validation set minimums
minML = MLte(argAvg,:);
minM0 = M0te(argAvg0,:);

% Test set minimums
minML2 = MLte2(argAvg,:);
minM02 = M0te2(argAvg0,:);


fprintf('\\begin{tabular} { |l | c | c | c |c| c| c | c|} \n')
fprintf('\\hline\n')
fprintf('Task & lmnn-val & Final-val & Early-val && lmnn-te & Final-te & Early-te \\\\ \n')
fprintf('\\hline\n')
for i = 1:length(minML)
fprintf('%d & %2.3f & %2.3f & %2.3f && %2.3f & %2.3f & %2.3f \\\\ \n',i,teErr(2,i),minM0(i), minML(i),teErr2(2,i),minM02(i), minML2(i))
end
fprintf('\\hline\n')
fprintf('\\hline\n')
fprintf('Best Avg & %2.3f & %2.3f & %2.3f && %2.3f & %2.3f & %2.3f \\\\ \n',mean(teErr(2,:)),minAvg0, minAvg,mean(teErr2(2,:)),minAvg02, minAvg2)
fprintf('\\hline\n')
fprintf('\\end{tabular}\n')

