function triplets = GenerateTriplets(X, Y,it,options)
% generate triplets based on Y and it

triplets = []; % x, x_1, x_2, y, it


num_instances = length(Y);
unique_Y = unique(Y);

% num_fail_triplets = 0;

for iu = 1:length(unique_Y)
    idx_unique_Y(iu).idx = find(Y(:) == unique_Y(iu));
%     
%     n_tmp = sum(Y(:)==unique_Y(iu));    
% %     fprintf('%d instance for class %d\n',n_tmp,unique_Y(iu));
%     if n_tmp < 3
%         continue;
%     else
%         num_fail_triplets = num_fail_triplets + nchoosek(n_tmp,3);        
% %         fprintf('%d infeasible instance for class %d\n',nchoosek(n_tmp,3),unique_Y(iu));
%     end;
end;
% num_ok_triplets = nchoosek(num_instances,3) - num_fail_triplets;
% 
% % num_to_generate = num_ok_triplets * ratio;

% num_to_generate = num_instances * ratio;
num_to_generate = num_instances * options.tri_time;
if it == -1
    num_to_generate = options.global_trip_num;
end;

hard = options.tri.hard;
hard_triplets_margin = 1;

num_to_generate = round(num_to_generate);
% fprintf('Num of triplets to generate is: %d\n', num_to_generate);
num_easy_trip = 0;

% ind_Y = 1:1:num_instances;
% all_tri = combnk(ind_Y, 3); %never use this, time-consuming
% num_all_tri = size(all_tri,1);
num_generated = 0;

while num_generated < num_to_generate
    same_class = randint(length(unique_Y));
    idx_i1 = randint(length(idx_unique_Y(same_class).idx));
    idx_i2 = randint(length(idx_unique_Y(same_class).idx));
    while idx_i2 == idx_i1
        idx_i2 = randint(length(idx_unique_Y(same_class).idx));
    end;
    i1 = idx_unique_Y(same_class).idx(idx_i1);
    i2 = idx_unique_Y(same_class).idx(idx_i2);
    
    diff_class = randint(length(unique_Y));
    while same_class == diff_class
        diff_class = randint(length(unique_Y));
    end;
    idx_ins_diff_class = randint(length(idx_unique_Y(diff_class).idx));
    i3 = idx_unique_Y(diff_class).idx(idx_ins_diff_class);
    
    if (Y(i1) == Y(i2))  && (Y(i2) == Y(i3))
        fprintf('Error in generate i1,i2,i3\n');
        continue;
    else
        if hard
            if X(:,i1)'*X(:,i2) > (X(:,i1)'*X(:,i3) - hard_triplets_margin)
                num_easy_trip = num_easy_trip + 1;
                continue;
            end;
        end;
        if rand <= 0.5            
            triplets(end+1,:) = [i1,i2,i3,1,it];            
        else           
            triplets(end+1,:) = [i1,i3,i2,-1,it];            
        end;	    
    end;  
    num_generated = num_generated + 1;
end;