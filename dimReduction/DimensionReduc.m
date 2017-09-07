data_name = 'news20';
data_dir = '../datasets/';
data_full_name = [data_dir,data_name,'.mat'];
load(data_full_name);
num_ins_task = [];
all_X = [];
variance_to_keep = 95 ;
for i=1:length(task)
    num_ins_task(end+1) = size(task(i).x,2);
    all_X = [all_X,task(i).x];
end;
clear task;
% all_X = all_X';
all_X_After_PCA = PCAReduction(all_X', variance_to_keep,data_name);
[n,d] = size(all_X_After_PCA);
fprintf('Instances:%d,Dimensions:%d\n',n,d);

for i=1:length(task)
    task2(i).x = all_X_After_PCA(1:num_ins_task(i),:)';
    task2(i).y = task(i).y;
end;
clear task;
task = task2;
save(sprintf('%s%s%d_%d.mat',data_dir,data_name,variance_to_keep,d),'task');