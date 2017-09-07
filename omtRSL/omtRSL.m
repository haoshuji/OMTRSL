function [M_omtRSL, que] = omtRSL(task, triplets_all_tasks, options)
%online multiple tasks Relative Similarity Learning


if strcmp(options.rel,'cov')
    C = options.omtRSL_CP_C;
elseif strcmp(options.rel,'fixed')
    C = options.omtRSL_FP_C;
end;
K = options.K;
b = options.omtRSL_FP_b;

if strcmp(options.query,'active') || strcmp(options.query, 'random')
    delta = options.delta;
end;

delta = options.delta;

num_queries = 0;
NT = options.NT;

update = true;
A = [];

d = size(task(1).x,1);
% M_omtRSL = [];
% X = [];
all_rand = [];


for idx_task =1:NT
	M_omtRSL(idx_task).M = eye(d);		
end;

if strcmp(options.rel, 'cov')
    A_tmp = [];
    for i=1:NT
        A_tmp = [A_tmp,vec(M_omtRSL(i).M)];
    end;
    A = cov(A_tmp);
elseif strcmp(options.rel, 'fixed')	
    A = eye(NT);
    A(~~eye(NT)) = (b+K)/((1+b)*K);
	A(~eye(NT)) = (b)/((1+b)*K);
else
    error('Unknown rel strategy');
end;

if NT < 1
	error('number of tasks should be at least more than 1\n');
end;

stream = RandStream('mt19937ar','Seed',200);
all_pt = [];
que_pt = [];
idx_que_pt = 1;
num_update = 0;
num_tri = size(triplets_all_tasks,1);
for ind_tri = 1:num_tri
    tri_cell = num2cell(triplets_all_tasks(ind_tri,:));
	[xi, xi1, xi2, y, it_tri] = tri_cell{:};
    
	X_t = task(it_tri).x(:,xi) * ( task(it_tri).x(:,xi1)-task(it_tri).x(:,xi2) )';
    
%     X_t = X_t / norm(X_t,'fro');
    
%     X_t_tmp = X_t / norm(X_t,'fro');
%     p_t_tmp = trace(M_omtRSL(it_tri).M*X_t_tmp');
	
    p_t = trace(M_omtRSL(it_tri).M*X_t');

    all_pt(end+1) = p_t;
    query = false;
	if strcmp(options.query,'passive')
		query = true;
	elseif strcmp(options.query, 'active') % to-be-done
        rand_float = rand(stream);
        all_rand(end+1) = rand_float;
        if rand_float < delta/(delta + abs(p_t))
            query = true;
        end;
    elseif strcmp(options.query, 'random') 
        if rand(stream) < delta
            query = true;
        end;
    elseif strcmp(options.query, 'active_WM')
    	p_t = 0;
    	for idx_task=1:NT
    		p_t = p_t + trace(M_omtRSL(idx_task).M*X_t');
    	end;
        p_t = p_t / NT;
    	rand_float = rand;
%         all_rand(end+1) = rand_float;
        if rand_float < delta/(delta + abs(p_t))
            query = true;
        end;
    else
        error('Unknown query strategy\n');
	end;
	
	if query        
        que_pt(idx_que_pt,1) = ind_tri;
        que_pt(idx_que_pt,2) = p_t;
        
        idx_que_pt = idx_que_pt+1;
		num_queries = num_queries + 1;
		ell = max(0, 1-y*p_t); 
		if ell > 0.0
            num_update = num_update + 1;
			tau = min(C, ell/norm(X_t,'fro')^2);			
			for idx_task = 1:NT
				M_omtRSL(idx_task).M = M_omtRSL(idx_task).M + A(idx_task,it_tri) * tau * y * X_t;
			end;
		end;
	end;
       
    if strcmp(options.rel, 'fixed')
        continue;
    elseif strcmp(options.rel, 'cov')          
        if ~rem(ind_tri,30)
            A_tmp = [];
            for i=1:NT
                A_tmp = [A_tmp,vec(M_omtRSL(i).M)];
            end;
            A = cov(A_tmp);
        end;
% 		d=dbstack;
%         error('%s:%s:%d To-be-done\n', d.file,d.name,d.line);
	else 
        error('unkonw rel update rule');
    end;
end;

que = num_queries / size(triplets_all_tasks,1);
% fprintf('num_update/num_tri = %f\n',num_update/num_tri);
% h = figure('visible','off');
% subplot(1,2,1);
% scatter([1:1:size(triplets_all_tasks)],all_pt)
% subplot(1,2,2);
% scatter(que_pt(:,1),que_pt(:,2));
% print(h,'-depsc',strcat('./tmp/',options.data_name,'_pt',options.query,num2str(delta),'.eps'));
% close(h);
% fprintf('%s_%s_delta=%f\n',options.rel,options.query,delta);