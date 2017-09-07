% This script load 'task' variables from files task_set%d files 
% Example task variable can be found in task_set1
% max_rep = # of different task_sets, i.e, number of distinct tr/va/te cuts


setpaths;
% Repeat 'max_rep' times with different validation sets
% (In this implementation, all these different cuts are loaded from files)
max_rep = 1; 

% single task max iterations
maxiter_st = 500;
% multi task max iterations
maxiter_mt = 500;
Kg=3;
% reuglarization constants
reg0 = [1];
reg = [500];

for r = 1:max_rep
 clear task;
 
 load(sprintf('task_set%d',r), 'task');
 clear mtD errs* mtwt D2


 for i=1:length(task)
    %% run LMNN
    val = length(task(i).yv)./(length(task(i).yv) + length(task(i).y));
    ntr = length(task(i).y);
    [L,D(i)]=lmnn2([task(i).x task(i).xv],[task(i).y task(i).yv],'quiet',0,'validation',val,'obj',0,'ntr',ntr,'maxiter',maxiter_st);
    task(i).L_for_nn = D(i).bestL;
    % task(i).L_for_nn = eye(size(D(i).bestL));
 end;

 for i=1:length(reg0)
  for j = 1:length(reg)
    %% run mtLMNN
    mtD{i,j}=MTlmnn_sp(task,Kg,'quiet',0,'regweight',reg(j),'regweight0',reg0(i),'weight1',0.5,'maxiter',maxiter_mt);
  end;
end;

 run_knns_REG; % get errs
 
 % collect the errors in cells
 Error_val{i} = errs_2;
 Error_test{i} = errs_2_te;
end;

% Now average the errors and print in latex format
AvgErr_val = Error_val{1};
AvgErr_test = Error_test{1};
for i=2:max_rep
  AvgErr_val = AvgErr_val + Error_val{i};
  AvgErr_test = AvgErr_test + Error_test{i};
end
AvgErr_val = AvgErr_val ./ max_rep;
AvgErr_test = AvgErr_test ./ max_rep;
% Read off errors. If more than one parameter setting is tested, use 
% validation set to read off errors of the best setting
report_results(AvgErr_val,AvgErr_test);

