clear
input_dir = './datasets/';
output_dir = './datasets/';
test_ratio = 0.2;
val_ratio = 0.2;
num_folds = 5;
data_names ={'face' };%'news20' 'isolet'}; % {'isolet' 'tic'}; %isolet
% X: each column corresponds to an instance
for id=1:length(data_names)
    fprintf('%s\n',data_names{id});
    clear task;
    clear data;
    val_ifo_ind = [2,3,4,5,1];
    fprintf('%s%s.mat',input_dir,data_names{id});
    load(sprintf('%s%s.mat',input_dir,data_names{id}));
    if strcmp(data_names{id},'isolet') 	
        num_tasks = length(task);
        for it = 1:num_tasks
            task_all(it).X = [task(it).x task(it).xv task(it).xTe];
            task_all(it).Y = [task(it).y task(it).yv task(it).yTe];            
        end;
    elseif strcmp(data_names{id},'landmine')
        num_tasks = length(label);
        fprintf('There are %d tasks totally\n',num_tasks);
        for it=1:num_tasks
            task_all(it).X = cell2mat(feature(it))';
            task_all(it).Y = cell2mat(label(it))';
            if size(task_all(it).X,2) ~= length(task_all(it).Y)
                error('size of feature X and label Y does not equal\n');
            end;
        end;
    elseif strcmp(data_names{id},'underwatermines')
        all_feasible_tasks=[1,2];
        num_tasks = length(all_feasible_tasks);
        fprintf('There are %d tasks totally\n',num_tasks);        
        for it=1:num_tasks
            task_all(it).X = cell2mat(xEachTask_Part(all_feasible_tasks(it)));
            task_all(it).Y = cell2mat(yEachTask_Part(all_feasible_tasks(it)))';
            if size(task_all(it).X,2) ~= length(task_all(it).Y)
                error('size of feature X and label Y does not equal\n');
            end;
        end;
    elseif strcmp(data_names{id},'tic')
        num_tasks = 4;
        task_all(1).X = task1_x'; task_all(1).Y = task1_y;
        task_all(2).X = task2_x'; task_all(2).Y = task2_y;
        task_all(3).X = task3_x'; task_all(3).Y = task3_y;
        task_all(4).X = task4_x'; task_all(4).Y = task4_y;
%         task_all(5).X = task5_x'; task_all(5).Y = task5_y;
%         task_all(6).X = task6_x'; task_all(6).Y = task6_y;        
    elseif strcmp(data_names{id},'email')
        load(sprintf('%s%s.mat',input_dir,data_names{id}));
        num_tasks = length(task_all);
    elseif strcmp(data_names{id},'imagenet')
        load(sprintf('%s%s.mat',input_dir,data_names{id}),'task');
        num_tasks = length(task);
        for idx_task = 1:num_tasks
            task_all(idx_task).X = task(idx_task).x;
            task_all(idx_task).Y = task(idx_task).y;
        end;
        clear task;
    else
        load(sprintf('%s%s.mat',input_dir,data_names{id}),'task');
        num_tasks = length(task);
        for idx_task = 1:num_tasks
            task_all(idx_task).X = task(idx_task).X;
            task_all(idx_task).Y = task(idx_task).Y;
        end;
        clear task;
    end;    
    
    fileID = fopen(sprintf('%s_info.txt',data_names{id}),'w');
    fprintf(fileID,'task &\t \\# classes &\t \\#Dim &\t \\#tr &\t \\#va &\t \\#te &\t \\# tri \\\\ \n');
    for it = 1:num_tasks
        fprintf(fileID,'%d &\t',it);
        unique_y = unique(task_all(it).Y);
        for i_un=1:length(unique_y)
            if sum(unique_y(i_un)==task_all(it).Y(:)) < 10
                Ind = find(task_all(it).Y == unique_y(i_un));
                task_all(it).Y(Ind) = [];
                task_all(it).X(:,Ind) = [];
            elseif sum(unique_y(i_un)==task_all(it).Y(:)) > 100
                Ind = find(task_all(it).Y == unique_y(i_un));
                task_all(it).Y(Ind(100+1:end)) = [];
                task_all(it).X(:,Ind(100+1:end)) = [];
            end;            
        end;
        num_all_ins = length(task_all(it).Y);
        d = size(task_all(it).X,1);
        fprintf(fileID,'%d &\t %d &\t %d &\t %d &\t %d &\t %d \\\\ \n',...
            length(unique(task_all(it).Y)), d, ...
            round(num_all_ins*0.6), round(num_all_ins*0.2),...
            round(num_all_ins*0.2), 10*(round(num_all_ins*0.6)));
        
        unique_y = unique(task_all(it).Y);
        for i_un = 1:length(unique_y)
            fprintf('%d\t',sum(unique_y(i_un)==task_all(it).Y(:)));
        end;
        
        if length(unique_y) >= 2
            task_all2(it).X = task_all(it).X;
            task_all2(it).Y = task_all(it).Y;
        end;           
    end;
    
    fprintf('\n');
    for it=1:num_tasks
        idx_perm = randperm(length(task_all2(it).Y));
        task_all2(it).X = task_all(it).X(:,idx_perm);
        task_all2(it).Y = task_all(it).Y(idx_perm);
    end;
    
    
    num_tasks = length(task_all2);
    last_len_unique_y = 100;
    for it = 1:num_tasks
        unique_y = unique(task_all2(it).Y);
        for i_un = 1:length(unique_y)
            Ind = find(task_all2(it).Y == unique_y(i_un));
            task_all2(it).Y(Ind) = last_len_unique_y + i_un;
        end;
        last_len_unique_y = last_len_unique_y + length(unique_y);
    end;
    
    norm_methods = {'orignal' 'length1' 'gaussian' 'scale01' };
%     norm_methods = {'length1'};
    for idx_norm = 1:length(norm_methods)
        fprintf('%s normalization\n',norm_methods{idx_norm});
        
        for it = 1:num_tasks
            task_all_pro(it).Y = task_all2(it).Y;

            if strcmp(norm_methods{idx_norm},'length1')% for each instance
                task_all_pro(it).X = normc(task_all2(it).X);
            elseif strcmp(norm_methods{idx_norm},'gaussian')% for each feature
                for i=1:size(task_all2(it).X,1)
                    sd=std(task_all2(it).X(i,:));
                    if sd~=0.0,
                        mu=mean(task_all2(it).X(i,:));
                        task_all_pro(it).X(i,:)=(task_all2(it).X(i,:)-mu)/sd;
                    end        
                end
            elseif strcmp(norm_methods{idx_norm},'scale01')
                for idx_fea = 1:size(task_all2(it).X,1)
                    fea_min = min(task_all2(it).X(idx_fea,:));
                    fea_max = max(task_all2(it).X(idx_fea,:));
                    task_all_pro(it).X(i,:)=(task_all2(it).X(idx_fea,:)-fea_min)/(fea_max - fea_min);
                end;
            elseif strcmp(norm_methods{idx_norm},'scaleminus1')
                for idx_fea = 1:size(task_all2(it).X,1)
                    fea_min = min(task_all2(it).X(idx_fea,:));
                    fea_max = max(task_all2(it).X(idx_fea,:));
                    task_all_pro(it).X(i,:)=(task_all2(it).X(idx_fea,:)-fea_min)/(fea_max - fea_min);
                    task_all_pro(it).X(i,:)=-1 + task_all2(it).X(idx_fea,:)*2;
                end;
            elseif strcmp(norm_methods{idx_norm},'orignal')
                task_all_pro(it).X = task_all2(it).X;
            else
                error('Unknown normalization method');    
            end;
            %% permutation
%             num_ins = length(task_all_pro(it).Y);
%             perm_idx = randperm(num_ins);
%             task_all_pro(it).Y = task_all_pro(it).Y(perm_idx);
%             task_all_pro(it).X = task_all_pro(it).X(:,perm_idx);
            
            %% crowwvalind
            [test,sort_idx] = sort(task_all_pro(it).Y);
            task_all_pro(it).Y = task_all_pro(it).Y(sort_idx);
            task_all_pro(it).X = task_all_pro(it).X(:,sort_idx);
            
            uni_y = unique(task_all_pro(it).Y);
            for i_un = 1:length(uni_y)
                task_all_pro(it).class(i_un).num_ins = sum(task_all_pro(it).Y == uni_y(i_un));
                task_all_pro(it).class(i_un).crs_idx= crossvalind('Kfold', task_all_pro(it).class(i_un).num_ins, num_folds);      
            end;
                 
        end;
        
        for ifo=1:num_folds
            for it=1:num_tasks
                uni_y = unique(task_all_pro(it).Y);
                crs_idx_start = 0;
                
                x = [];y=[];xv=[];yv=[];xt=[];yt=[];
                for i_un = 1:length(uni_y)
                    test_ind = find(task_all_pro(it).class(i_un).crs_idx == ifo) ;
                    val_ind = find(task_all_pro(it).class(i_un).crs_idx == val_ifo_ind(ifo)) ;                    
                    
                    fullidx = [1:1:task_all_pro(it).class(i_un).num_ins];                    
                    train_ind = setdiff(fullidx,[test_ind',val_ind']);
                    
                    test_ind = test_ind + crs_idx_start;
                    val_ind = val_ind + crs_idx_start;
                    train_ind = (train_ind + crs_idx_start)';
                    
                    x = [x, task_all_pro(it).X(:,train_ind)];
                    y = [y, task_all_pro(it).Y(train_ind)];
                    xv = [xv, task_all_pro(it).X(:,val_ind)];
                    yv = [yv, task_all_pro(it).Y(val_ind)];
                    xt = [xt, task_all_pro(it).X(:,test_ind)];
                    yt = [yt, task_all_pro(it).Y(test_ind)];
                    
                    crs_idx_start = crs_idx_start + task_all_pro(it).class(i_un).num_ins;
                    
                end;
                
                data(ifo).task(it).x = x;
                data(ifo).task(it).y = y;
                data(ifo).task(it).xv = xv;
                data(ifo).task(it).yv = yv;
                data(ifo).task(it).xt = xt;
                data(ifo).task(it).yt = yt;
            end;		
        end;

        save(sprintf('%s%s_%s_all_folds.mat',output_dir,data_names{id},norm_methods{idx_norm}),'data');%,'-v7.3'
    end;
end;