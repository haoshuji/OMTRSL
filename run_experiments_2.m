clear;
% clc;
setpaths;
cur_path = pwd;
data_dir = strcat(cur_path,'/datasets');
base_results_dir = strcat(cur_path,'/results');

%% all options
options.use_cov = false;
options.topk = [1,2,3,4,5,6,7,8,9,10];
NTK = length(options.topk);

options.validation_measure = (10);
options.b_candidate = [1,2,4,8,16,32,64,128];
options.C_candidate = [0.001,0.01,0.1,1,10,100,1000];
options.reg0_candidate = [1,100,300,400,500]; 
options.reg_candidate = [1,100,300,400,500];

options.triplets_using_existed = false; 
options.tri.hard= false;  
options.tri_time = 5;

options.mtLMNN_using_existed = false; 
options.mtLMNN_reg_search = false;
options.maxitaer_st = 500; 
options.maxitaer_mt = 500;
options.Kg = 1; 
options.mtLMNN.best_reg0 = 1; 
options.mtLMNN.best_reg = 500;

options.ostRSL_global_using_existed = false;   
options.ostRSL_global_C_search = true;    
options.ostRSL_global_C = 1; 
options.ostRSL_global_instances_using_existed = false;

options.ostRSL.exi = false;    
options.ostRSL_C_search = true;     
options.ostRSL_C = 0.001;

options.omtRSL_FP_using_existed = false;  
options.omtRSL_FP_bC_search = 1;  
options.omtRSL_FP_b = 2; 
options.omtRSL_FP_C =0.001;

options.run_omtRSL_CP = false;
options.omtRSL_CP_using_existed = false;  
options.omtRSL_CP_C_search = false;   
options.omtRSL_CP_C = 1;

options.omtRSL_AR_using_existed = false;
all_iq = (1:1:20);
num_que = length(all_iq);
options.delta = -1;
options.query = 'passive'; 
options.rel = 'fixed';
options.delta_start = -15; 
option.delta_multiplier = 1.2;

data_names = { 'face'};%'isolet' 'letter'
norm_methods = {'orignal' 'length1' 'gaussian' 'scale01' };
% norm_methods= {'length1'}; 


%% all datasets
for idx_data=1:length(data_names)
    data_name = data_names{idx_data};
    for id = 1:length(norm_methods)
        full_name = strcat(data_name,'_',norm_methods{id});
        fprintf('Experiments on %s dataset\n',full_name);   
%         current_date = datestr(now,'mmmm-dd-yyyy-HH');
        current_date =  '30-Jan-2016';

        results_dir = strcat(base_results_dir,'/',full_name,'/',current_date);
        
        if ~exist(results_dir, 'dir')  
            mkdir(results_dir);    
        end;

        load(sprintf('%s/%s_all_folds.mat',data_dir,full_name),'data');    

        NF = length(data);        NT = length(data(1).task);
        options.NT = NT;        options.data_name = full_name;        options.K = NT;
        
        %% allocating results space        
        %%time, que: NF size; others: NT * NTK * NF size
        time_mtLMNN = zeros(NF,1);  
        AP_tr_mtLMNN = zeros(NT,NTK,NF); mAP_tr_mtLMNN = zeros(NT,NTK,NF);
        AP_te_mtLMNN = zeros(NT,NTK,NF); mAP_te_mtLMNN = zeros(NT,NTK,NF);

        time_ostRSL = zeros(NF,1);  
        AP_tr_ostRSL = zeros(NT,NTK,NF); mAP_tr_ostRSL = zeros(NT,NTK,NF);
        AP_te_ostRSL = zeros(NT,NTK,NF); mAP_te_ostRSL = zeros(NT,NTK,NF);

        time_ostRSL_global = zeros(NF,1); 
        AP_tr_ostRSL_global = zeros(NTK,NF); mAP_tr_ostRSL_global = zeros(NTK,NF);  
        AP_te_ostRSL_global = zeros(NTK,NF); mAP_te_ostRSL_global = zeros(NTK,NF);

        time_omtRSL_FP = zeros(NF,1);  
        AP_tr_omtRSL_FP = zeros(NT,NTK,NF); mAP_tr_omtRSL_FP = zeros(NT,NTK,NF); 
        AP_te_omtRSL_FP = zeros(NT,NTK,NF); mAP_te_omtRSL_FP  = zeros(NT,NTK,NF); 

        time_omtRSL_CP = zeros(NF,1);  
        AP_tr_omtRSL_CP = zeros(NT,NTK,NF); mAP_tr_omtRSL_CP = zeros(NT,NTK,NF); 
        AP_te_omtRSL_CP = zeros(NT,NTK,NF); mAP_te_omtRSL_CP   = zeros(NT,NTK,NF);  

        %%time, que: num_que *   NF, other: num_que * num_task * NTK * num_fol 
        time_omtRSL_FA = zeros(num_que,NF); que_omtRSL_FA  = zeros(num_que,NF);    
        AP_tr_omtRSL_FA = zeros(num_que,NT,NTK,NF);   mAP_tr_omtRSL_FA = zeros(num_que,NT,NTK,NF);    
        AP_te_omtRSL_FA = zeros(num_que,NT,NTK,NF);   mAP_te_omtRSL_FA = zeros(num_que,NT,NTK,NF);

        time_omtRSL_FR = zeros(num_que,NF); que_omtRSL_FR = zeros(num_que,NF);    
        AP_tr_omtRSL_FR = zeros(num_que,NT,NTK,NF);   mAP_tr_omtRSL_FR = zeros(num_que,NT,NTK,NF);
        AP_te_omtRSL_FR = zeros(num_que,NT,NTK,NF);   mAP_te_omtRSL_FR = zeros(num_que,NT,NTK,NF);
                
        %% generate global dataset for ostRSL_global
        if options.ostRSL_global_instances_using_existed
            load(sprintf('%s/%s_global.mat',results_dir,full_name));
        else
            global_ins = struct;
            for ifo = 1:NF 
                global_ins(ifo).x =[]; global_ins(ifo).y =[];
                global_ins(ifo).xv =[]; global_ins(ifo).yv =[];
                global_ins(ifo).xt =[]; global_ins(ifo).yt =[];
            end;
            for ifo=1:NF
                for ita =1:NT
                    global_ins(ifo).x = [global_ins(ifo).x, data(ifo).task(ita).x];
                    global_ins(ifo).y = [global_ins(ifo).y, data(ifo).task(ita).y];

                    global_ins(ifo).xv = [global_ins(ifo).xv, data(ifo).task(ita).xv];
                    global_ins(ifo).yv = [global_ins(ifo).yv, data(ifo).task(ita).yv];

                    global_ins(ifo).xt = [global_ins(ifo).xt, data(ifo).task(ita).xt];
                    global_ins(ifo).yt = [global_ins(ifo).yt, data(ifo).task(ita).yt];
                end;
            end;
            save(sprintf('%s/%s_global.mat',results_dir,full_name),'global_ins');
        end;
        
        %% generate triplets
        if options.triplets_using_existed
            load(sprintf('%s/%s_triplets',results_dir,full_name));
        else         
            fprintf('Generating triplets\n');     
            for ifo = 1:NF
%                 fprintf('%d-fold\n',ifo);                                
                tri_tmp = [];
                for ita=1:NT                    
                    triTF(ifo).task(ita).trip = GenerateTriplets(data(ifo).task(ita).x,data(ifo).task(ita).y,ita,options);                     
                    perm_ind = randperm(size(triTF(ifo).task(ita).trip,1));
                    triTF(ifo).task(ita).trip  = triTF(ifo).task(ita).trip(perm_ind,:);
                    tri_tmp =[tri_tmp;triTF(ifo).task(ita).trip];                        
%                     fprintf('%d-task:%d\t',ita,size(triTF(ifo).task(ita).trip,1));
                end;
                triF(ifo).trip = tri_tmp; 	num_trip_all_tasks = size(triF(ifo).trip,1);
%                 fprintf('\nAll task:%d triplets\t',num_trip_all_tasks);                
                options.global_trip_num = num_trip_all_tasks;
                triGlobal(ifo).trip = GenerateTriplets(global_ins(ifo).x, global_ins(ifo).y, -1, options);
%                 fprintf('Global trip:%d triplets\n',size(triGlobal(ifo).trip,1)); 
            end;
            save(sprintf('%s/%s_triplets.mat',results_dir,full_name),'triTF','triF','triGlobal');
        end;
        
        %% mtLMNN    
        if options.mtLMNN_reg_search
            fprintf('\nFind best reg for mtlnn\n');
            mtlmnn_best_AP_va = 0;
            all_AP_va_mtLMNN = zeros(length(options.reg0_candidate),length(options.reg_candidate));
            for i=1:length(options.reg0_candidate)
                for j = 1:length(options.reg_candidate)
                    mtlmnn_AP_va_tmp = 0;
                    for ifo=1:NF               
                        mtD_tmp=MTlmnn_sp(data(ifo).task,options.Kg,'quiet',1,'regweight',options.reg_candidate(j), 'regweight0',options.reg0_candidate(i),'weight1',0.5,'maxitaer',options.maxitaer_mt); 
                        for ita=1:NT                            
                            [AP_va,mAP_va] = Test(mtD_tmp.bestMgen.M0+mtD_tmp.bestMgen.M{ita}, data(ifo).task(ita).xv, data(ifo).task(ita).yv, options); 
                            mtlmnn_AP_va_tmp = mtlmnn_AP_va_tmp + mean(AP_va(options.validation_measure));
                        end; 
                    end;     
                    mtlmnn_AP_va_tmp = mtlmnn_AP_va_tmp/(NF*NT);
                    all_AP_va_mtLMNN(i,j) = mtlmnn_AP_va_tmp;

                    if mtlmnn_AP_va_tmp > mtlmnn_best_AP_va
                        mtlmnn_best_AP_va = mtlmnn_AP_va_tmp; 
                        options.mtLMNN.best_reg = options.reg_candidate(j); 
                        options.mtLMNN.best_reg0 = options.reg0_candidate(i);
                    end;
                    fprintf('current_AP_va:%.3f,\tcurrent_reg0:%.3f,\tcurrent_reg:%.3f\n', mtlmnn_AP_va_tmp,options.reg0_candidate(i),options.reg_candidate(j))
                end;
            end;
            fprintf('best_AP_te:%.3f,\tbest_reg0:%.3f,\tbest_reg:%.3f\n', mtlmnn_best_AP_va,options.mtLMNN.best_reg0,options.mtLMNN.best_reg);
            best_reg = options.mtLMNN.best_reg; best_reg0 = options.mtLMNN.best_reg0;
            save(sprintf('%s/%s_mtLMNN_results_parameters.mat',results_dir,full_name), 'all_AP_va_mtLMNN','best_reg','best_reg0');
        end; 
        
        if options.mtLMNN_using_existed
            load(sprintf('%s/%s_mtLMNN_results.mat',results_dir,full_name));   
            fprintf('mtLMNN:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n',mean2(AP_tr_mtLMNN), mean2(mAP_tr_mtLMNN),mean2(AP_te_mtLMNN),mean2(mAP_te_mtLMNN));
        else
            fprintf('mtLMNN:\t');          
            for ifo = 1:NF
                fprintf('%d\t',ifo);
                
                tic
                mtD=MTlmnn_sp(data(ifo).task,options.Kg,'quiet',1,'regweight',options.mtLMNN.best_reg, 'regweight0',options.mtLMNN.best_reg0,'weight1',0.5,'maxitaer',options.maxitaer_mt); 
                time_mtLMNN(ifo) = toc;                            

                for ita=1:NT
                    mtLMNN(ifo).task(ita).M= mtD.bestMgen.M0+mtD.bestMgen.M{ita};                            
                end;                
                for ita=1:NT
                    [AP_tr,mAP_tr] = Test(mtLMNN(ifo).task(ita).M,data(ifo).task(ita).x, data(ifo).task(ita).y, options); 
                    [AP_te,mAP_te] = Test(mtLMNN(ifo).task(ita).M,data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);                     
                    AP_tr_mtLMNN(ita,:,ifo) = AP_tr;  mAP_tr_mtLMNN(ita,:,ifo) = mAP_tr; AP_te_mtLMNN(ita,:,ifo) = AP_te;   mAP_te_mtLMNN(ita,:,ifo) = mAP_te;
                end;                      
            end;
            fprintf('\n');
            fprintf('mtLMNN:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_mtLMNN),mean2(mAP_tr_mtLMNN), mean2(AP_te_mtLMNN),mean2(mAP_te_mtLMNN));
            save(sprintf('%s/%s_mtLMNN_results.mat',results_dir,full_name), 'time_mtLMNN','AP_tr_mtLMNN','mAP_tr_mtLMNN','AP_te_mtLMNN', 'mAP_te_mtLMNN','mtLMNN');
        end;
        
        %% ostRSL_global
        if options.ostRSL_global_C_search
            fprintf('\nFind best C for ostRSL_global.\n');
        	best_C = options.ostRSL_global_C;
        	best_AP_va_tmp = 0;
        	all_AP_va_ostRSL_global_C = [];
        	for idx_C=1:length(options.C_candidate)
        		options.ostRSL_global_C = options.C_candidate(idx_C);
        		AP_va_tmp = 0;
        		for ifo=1:NF
                    for ita = 1:NT                
	                    M_tmp = ostRSL(global_ins(ifo).x, global_ins(ifo).y, triGlobal(ifo).trip, options);	                        
	                    [AP_va,mAP_va] = Test(M_tmp, global_ins(ifo).xv, global_ins(ifo).yv, options);                     
	                    AP_va_tmp = AP_va_tmp + mean(AP_va(options.validation_measure));
                    end;
        		end;
        		AP_va_tmp = AP_va_tmp/(NF*NT);
        		all_AP_va_ostRSL_global_C(end+1) = AP_va_tmp;
        		if AP_va_tmp > best_AP_va_tmp
        			best_C = options.ostRSL_global_C;
        			best_AP_va_tmp = AP_va_tmp;
        		end;
        		fprintf('Current_C:%f, \t Current_AP_va:%f\n',options.ostRSL_global_C,AP_va_tmp);
        	end;
            options.ostRSL_global_C = best_C;
            fprintf('Best_C:%f, \t Best_AP_va:%f\n',best_C,best_AP_va_tmp);
        	save(sprintf('%s/%s_ostRSL_global_results_parameters.mat',results_dir,full_name),'all_AP_va_ostRSL_global_C');        	
        end;
        
        if 1
            if options.ostRSL_global_using_existed
                load(sprintf('%s/%s_ostRSL_global_results.mat',results_dir,full_name));       
                fprintf('ostRSL_global:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_ostRSL_global),mean2(mAP_tr_ostRSL_global), mean2(AP_te_ostRSL_global),mean2(mAP_te_ostRSL_global));
            else
                fprintf('ostRSL_global:\t');
                for ifo = 1:NF
                    fprintf('%d\t',ifo);
                    tic; ostRSL_global_results(ifo).M = ostRSL(global_ins(ifo).x, global_ins(ifo).y, triGlobal(ifo).trip, options);                
                    time_ostRSL_global(ifo) = toc; 
                    % [AP_tr, mAP_tr] = Test(ostRSL_global_results(ifo).M, global_ins(ifo).x, global_ins(ifo).y, options);                    
                    AP_tr = 0; mAP_tr = 0;
                    [AP_te, mAP_te] = Test(ostRSL_global_results(ifo).M, global_ins(ifo).xt, global_ins(ifo).yt, options);
                    AP_tr_ostRSL_global(:,ifo) = AP_tr;     mAP_tr_ostRSL_global(:,ifo) = mAP_tr; AP_te_ostRSL_global(:,ifo) = AP_te;     mAP_te_ostRSL_global(:,ifo) = mAP_te;  
                end;
                fprintf('\n');
                fprintf('ostRSL_global:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n',mean2(AP_tr_ostRSL_global),mean2(mAP_tr_ostRSL_global), mean2(AP_te_ostRSL_global),mean2(mAP_te_ostRSL_global));
                save(sprintf('%s/%s_ostRSL_global_results.mat',results_dir,full_name), 'time_ostRSL_global','AP_tr_ostRSL_global','mAP_tr_ostRSL_global','AP_te_ostRSL_global', 'mAP_te_ostRSL_global');
            end; 
        end;
        
        %% M = eye(d,d);
        if 0
            d = size(data(1).task(1).x,1);
            for ifo = 1:NF
                fprintf('%d\t',ifo);            
                for ita = 1:NT
                    [AP_tr, mAP_tr] = Test(eye(d,d), data(ifo).task(ita).x, data(ifo).task(ita).y, options);                    
                    [AP_te, mAP_te] = Test(eye(d,d), data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                    AP_tr_ostRSL_eye(ita,:,ifo) = AP_tr;     mAP_tr_ostRSL_eye(ita,:,ifo) = mAP_tr; AP_te_ostRSL_eye(ita,:,ifo) = AP_te;     mAP_te_ostRSL_eye(ita,:,ifo) = mAP_te;  
                end;            
            end;
            fprintf('\nostRSL_eye:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_ostRSL_eye), mean2(mAP_tr_ostRSL_eye), mean2(AP_te_ostRSL_eye),mean2(mAP_te_ostRSL_eye));
            save(sprintf('%s/%s_ostRSLE_results.mat',results_dir,full_name), 'AP_tr_ostRSL_eye', 'mAP_tr_ostRSL_eye', 'AP_te_ostRSL_eye', 'mAP_te_ostRSL_eye');
        end;
       
        %% ostRSL
        if options.ostRSL_C_search
            fprintf('\nFind best C for ostRSL.\n');
        	best_C = options.ostRSL_C;
        	best_AP_va_tmp = 0;
        	all_AP_va_ostRSL_C = [];
        	for idx_C=1:length(options.C_candidate)
        		options.ostRSL_C = options.C_candidate(idx_C);
        		AP_va_tmp = 0;
        		for ifo=1:NF
                    for ita = 1:NT                
	                    M_tmp = ostRSL([data(ifo).task(ita).x], [data(ifo).task(ita).y], triTF(ifo).task(ita).trip, options);
	                    [AP_va,mAP_va] = Test(M_tmp,data(ifo).task(ita).xv, data(ifo).task(ita).yv, options);                     
	                    AP_va_tmp = AP_va_tmp + mean(AP_va(options.validation_measure));
                    end;
        		end;
        		AP_va_tmp = AP_va_tmp/(NF*NT); all_AP_va_ostRSL_C(end+1) = AP_va_tmp;

        		if AP_va_tmp > best_AP_va_tmp
        			best_C = options.ostRSL_C; best_AP_va_tmp = AP_va_tmp;
        		end;
        		fprintf('Current_C:%f, \t Current_AP_va:%f\n',options.ostRSL_C,AP_va_tmp);
        	end;
        	options.ostRSL_C = best_C;
            fprintf('Best_C:%f, \t Best_AP_va:%f\n',best_C,best_AP_va_tmp);
            save(sprintf('%s/%s_ostRSL_results_parameters.mat',results_dir,full_name),'all_AP_va_ostRSL_C');
        end;

        if options.ostRSL.exi
            load(sprintf('%s/%s_ostRSL_results.mat',results_dir,full_name));       
            fprintf('ostRSL:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_ostRSL),mean2(mAP_tr_ostRSL), mean2(AP_te_ostRSL),mean2(mAP_te_ostRSL));
        else
            fprintf('ostRSL:\t');
            for ifo = 1:NF
                fprintf('%d\t',ifo);
                tic;                
                for ita = 1:NT                
                    ostRSL_results(ifo).task(ita).M = ostRSL( data(ifo).task(ita).x, data(ifo).task(ita).y, triTF(ifo).task(ita).trip, options);
                end;
                time_ostRSL(ifo) = toc; 

                for ita = 1:NT
                    [AP_tr, mAP_tr] = Test(ostRSL_results(ifo).task(ita).M,data(ifo).task(ita).x, data(ifo).task(ita).y, options);                    
                    [AP_te, mAP_te] = Test(ostRSL_results(ifo).task(ita).M,data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                    AP_tr_ostRSL(ita,:,ifo) = AP_tr;     mAP_tr_ostRSL(ita,:,ifo) = mAP_tr; 
                    AP_te_ostRSL(ita,:,ifo) = AP_te;     mAP_te_ostRSL(ita,:,ifo) = mAP_te;  
                end;            
            end;
            fprintf('\n');
            fprintf('ostRSL:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_ostRSL),mean2(mAP_tr_ostRSL), mean2(AP_te_ostRSL),mean2(mAP_te_ostRSL));
            save(sprintf('%s/%s_ostRSL_results.mat',results_dir,full_name), 'time_ostRSL','AP_tr_ostRSL','mAP_tr_ostRSL','AP_te_ostRSL','mAP_te_ostRSL');
        end;    
        
        %% omtRSL_cov_pa
        if options.omtRSL_CP_C_search && options.run_omtRSL_CP
            fprintf('\nFind best C for omtRSL_cov_pa.\n');
        	best_C = options.omtRSL_CP_C;
        	best_AP_va_tmp = 0;
        	all_AP_va_omtRSL_CP_C = [];
            options.rel = 'cov';
        	for idx_C=1:length(options.C_candidate)
        		options.omtRSL_CP_C = options.C_candidate(idx_C);
        		AP_va_tmp = 0;
                for ifo =1:NF                    
                    [MomtRSL_FP_tmp,que] = omtRSL(data(ifo).task,triF(ifo).trip, options);                            
                    for ita = 1:NT
                        [AP_va, mAP_va] = Test(MomtRSL_FP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, options);
                        AP_va_tmp = AP_va_tmp + mean(AP_va(options.validation_measure));                            
                    end;                            
                end;

        		AP_va_tmp = AP_va_tmp/(NF*NT);
        		all_AP_va_omtRSL_CP_C(end+1) = AP_va_tmp;

        		if AP_va_tmp > best_AP_va_tmp
        			best_C = options.omtRSL_CP_C;
        			best_AP_va_tmp = AP_va_tmp;
        		end;
        		fprintf('Current_C:%f, \t Current_AP_va:%f\n',options.omtRSL_CP_C,AP_va_tmp);
        	end;
            fprintf('Best_C:%f, \t Best_AP_va:%f\n',best_C,best_AP_va_tmp);
        	save(sprintf('%s/%s_omtRSL_CP_results_parameters.mat',results_dir,full_name),'all_AP_va_omtRSL_CP_C');
        	options.omtRSL_CP_C = best_C;
        end;

        if options.omtRSL_CP_using_existed && options.run_omtRSL_CP
            load(sprintf('%s/%s_omtRSL_CP_results.mat',results_dir,full_name));       
            fprintf('omtRSL_cov_pa:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n',mean2(AP_tr_omtRSL_CP),mean2(mAP_tr_omtRSL_CP), mean2(AP_te_omtRSL_CP),mean2(mAP_te_omtRSL_CP));
        elseif options.run_omtRSL_CP
            options.rel = 'cov';
            fprintf('omtRSL_cov_pa:\t');
            for ifo = 1:NF
                fprintf('%d\t',ifo);
                tic;
                [omtRSL_cov_pa_results(ifo).task,que_useless] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                time_omtRSL_FP(ifo) = toc; 

                for ita = 1:NT
                    [AP_tr, mAP_tr] = Test(omtRSL_cov_pa_results(ifo).task(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);   
                    [AP_te, mAP_te] = Test(omtRSL_cov_pa_results(ifo).task(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                    AP_tr_omtRSL_CP(ita,:,ifo) = AP_tr;    mAP_tr_omtRSL_CP(ita,:,ifo) = mAP_tr; AP_te_omtRSL_CP(ita,:,ifo) = AP_te;    mAP_te_omtRSL_CP(ita,:,ifo) = mAP_te;  
                end;
            end;  
            fprintf('\n');
            fprintf('omtRSL_cov_pa:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_omtRSL_CP),mean2(mAP_tr_omtRSL_CP), mean2(AP_te_omtRSL_CP),mean2(mAP_te_omtRSL_CP));
            save(sprintf('%s/%s_omtRSL_CP_results.mat',results_dir,full_name), 'time_omtRSL_CP','AP_tr_omtRSL_CP', 'mAP_tr_omtRSL_CP','AP_te_omtRSL_CP','mAP_te_omtRSL_CP');
        else
            fprintf('Ignore omtRSL_CP\n');
        end;    
        
        options.omtRSL_FP_C = options.ostRSL_C;
      
        %% omtRSL_fix_pa
        if options.omtRSL_FP_bC_search == 2     
            options.rel = 'fixed';
            fprintf('\nFind best b and C on AP_va\n');
            best_b = options.omtRSL_FP_b;
            best_C = options.omtRSL_FP_C;
            
            best_AP_va_tmp = 0;
            all_AP_va_omtRSL_FP = zeros(length(options.b_candidate),length(options.C_candidate));            
            
            for idx_b = 1:length(options.b_candidate)
                for idx_C = 1:length(options.C_candidate)
                    options.omtRSL_FP_b = options.b_candidate(idx_b);
                    options.omtRSL_FP_C = options.C_candidate(idx_C);
                    
                    AP_va_tmp = 0; 
                    for ifo =1:NF                    
                        [MomtRSL_FP_tmp,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);                            
                        for ita = 1:NT
                            [AP_va, mAP_va] = Test(MomtRSL_FP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, options);
                            AP_va_tmp = AP_va_tmp + mean(AP_va(options.validation_measure));                            
                        end;                            
                    end;
                    AP_va_tmp = AP_va_tmp/(NF*NT);                    
                    all_AP_va_omtRSL_FP(idx_b,idx_C) = AP_va_tmp;

                    if AP_va_tmp > best_AP_va_tmp
                        best_AP_va_tmp = AP_va_tmp; best_b = options.omtRSL_FP_b; best_C = options.omtRSL_FP_C;
                    end;
                    fprintf('Current_b:%f,\tCurrent_C:%f,\tCurrent_APTe:%f\n',options.omtRSL_FP_b, options.omtRSL_FP_C, AP_va_tmp);
                end;
            end;            
            options.omtRSL_FP_b = best_b; options.omtRSL_FP_C = best_C;
            fprintf('Best_b:%f,\tBest_C:%f,\tBest_AP_va:%f\n',best_b,best_C,best_AP_va_tmp);            
            save(sprintf('%s/%s_omtRSL_FP_results_parameters',results_dir,full_name),'all_AP_va_omtRSL_FP');            
            fprintf('omtRSL_FP_C:%f\t,omtRSL_FP_b:%f\n',options.omtRSL_FP_C,options.omtRSL_FP_b);
        elseif options.omtRSL_FP_bC_search == 1
            options.rel = 'fixed';
            fprintf('\nFind best b for omtRSL alone and set omtRSL_FP_C as OST_C\n');
            options.omtRSL_FP_C = options.ostRSL_C;

            best_b = options.omtRSL_FP_b;
            best_AP_va_tmp = 0;
            all_AP_va_omtRSL_FP = zeros(length(options.b_candidate),1);            
            
            for idx_b = 1:length(options.b_candidate)
                options.omtRSL_FP_b = options.b_candidate(idx_b);
                AP_va_tmp = 0; 
               
                for ifo =1:NF                    
                    [MomtRSL_FP_tmp,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);                            
                    for ita = 1:NT
                        [AP_va, mAP_va] = Test(MomtRSL_FP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, options);
                        AP_va_tmp = AP_va_tmp + mean(AP_va(options.validation_measure));                            
                    end;                            
                end;

                AP_va_tmp = AP_va_tmp/(NF*NT); all_AP_va_omtRSL_FP(idx_b) = AP_va_tmp;

                if AP_va_tmp > best_AP_va_tmp
                    best_AP_va_tmp = AP_va_tmp; best_b = options.omtRSL_FP_b;
                end;
                fprintf('Current_b:%f,\tCurrent_APTe:%f\n',options.omtRSL_FP_b, AP_va_tmp);
            end; 
             options.omtRSL_FP_b = best_b;           
            fprintf('Best_b:%f,\tBest_AP_va:%f\n',best_b, best_AP_va_tmp);            
            save(sprintf('%s/%s_omtRSL_FP_results_parameters.mat',results_dir,full_name),'all_AP_va_omtRSL_FP');            
        else
            fprintf('Using Existing omtRSL_FP_C:%f and omtRSL_FP_b:%f\n',options.omtRSL_FP_C,options.omtRSL_FP_b);
        end;

        fprintf('omtRSL_FP_C:%f\t,omtRSL_FP_b:%f\n',options.omtRSL_FP_C,options.omtRSL_FP_b);
        
        save(sprintf('%s/%s_options.mat',results_dir,full_name),'options');

        if options.omtRSL_FP_using_existed
            load(sprintf('%s/%s_omtRSL_FP_results.mat',results_dir,full_name));
            fprintf('omtRSL:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_omtRSL_FP),mean2(mAP_tr_omtRSL_FP), mean2(AP_te_omtRSL_FP),mean2(mAP_te_omtRSL_FP));
        else
            options.rel = 'fixed';
            fprintf('omtRSL_fix_pa:\t');
            for ifo = 1:NF
                fprintf('%d\t',ifo);
                tic;
                [omtRSL_results(ifo).task,que_useless] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                time_omtRSL_FP(ifo) = toc; 

                for ita = 1:NT
                    [AP_tr, mAP_tr] = Test(omtRSL_results(ifo).task(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);   
                    [AP_te, mAP_te] = Test(omtRSL_results(ifo).task(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                    AP_tr_omtRSL_FP(ita,:,ifo) = AP_tr;    mAP_tr_omtRSL_FP(ita,:,ifo) = mAP_tr; AP_te_omtRSL_FP(ita,:,ifo) = AP_te;    mAP_te_omtRSL_FP(ita,:,ifo) = mAP_te;  
                end;
            end;        
            fprintf('\n');
            fprintf('omtRSL:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_omtRSL_FP),mean2(mAP_tr_omtRSL_FP), mean2(AP_te_omtRSL_FP),mean2(mAP_te_omtRSL_FP));
            save(sprintf('%s/%s_omtRSL_FP_results.mat',results_dir,full_name), 'time_omtRSL_FP','AP_tr_omtRSL_FP', 'mAP_tr_omtRSL_FP','AP_te_omtRSL_FP', 'mAP_te_omtRSL_FP');
        end;

        %% omtRSL_fix_act and omtRSL_fix_rand
        if options.omtRSL_AR_using_existed
            load(sprintf('%s/%s_omtRSL_AR_results.mat',results_dir, full_name));
        else
            for idx_que = 1:num_que
                iq = all_iq(idx_que);
                fprintf('%d-th query\n',iq);            
                for ifo = 1:NF
                    fprintf('%d-th fold\n',ifo);

                    %% omtRSL_fix_act
                    delta_mul = options.delta_start + iq*option.delta_multiplier;
%                     delta_mul = 5;
                    options.rel = 'fixed'; options.query = 'active';  
                    tic; options.delta = 2^(delta_mul);
                    [MOMTFA,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                    time_omtRSL_FA(idx_que,ifo) = toc;  que_omtRSL_FA(idx_que,ifo) = que;

                    for ita = 1:NT
                        [AP_tr, mAP_tr] = Test(MOMTFA(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);    
                        [AP_te, mAP_te] = Test(MOMTFA(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);                
                        AP_tr_omtRSL_FA(idx_que,ita,:,ifo) = AP_tr;  mAP_tr_omtRSL_FA(idx_que,ita,:,ifo) = mAP_tr; AP_te_omtRSL_FA(idx_que,ita,:,ifo) = AP_te;  mAP_te_omtRSL_FA(idx_que,ita,:,ifo) = mAP_te;                  
                    end;
                    fprintf('omtRSL_fix_act:Que:%.3f,AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', que,mean2(AP_tr_omtRSL_FA(idx_que,:,:,ifo)), mean2(mAP_tr_omtRSL_FA(idx_que,:,:,ifo)),...
                        mean2(AP_te_omtRSL_FA(idx_que,:,:,ifo)), mean2(mAP_te_omtRSL_FA(idx_que,:,:,ifo)));

                    %% omtRSL_fix_rand
                    options.rel = 'fixed'; options.query = 'random';
                    options.delta = que_omtRSL_FA(idx_que, ifo);
                    tic; [MOMTFR,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                    time_omtRSL_FR(idx_que,ifo) = toc;  que_omtRSL_FR(idx_que,ifo) = que;

                    for ita = 1:NT
                        [AP_tr, mAP_tr] = Test(MOMTFR(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);
                        [AP_te, mAP_te] = Test(MOMTFR(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                        AP_tr_omtRSL_FR(idx_que,ita,:,ifo)  = AP_tr;  mAP_tr_omtRSL_FR(idx_que,ita,:,ifo) = mAP_tr; AP_te_omtRSL_FR(idx_que,ita,:,ifo)  = AP_te;  mAP_te_omtRSL_FR(idx_que,ita,:,ifo) = mAP_te;  
                    end;
                    fprintf('omtRSL_fix_rand:Que:%.3f,AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', que,mean2(AP_tr_omtRSL_FR(idx_que,:,:,ifo)), mean2(mAP_tr_omtRSL_FR(idx_que,:,:,ifo)),...
                        mean2(AP_te_omtRSL_FR(idx_que,:,:,ifo)), mean2(mAP_te_omtRSL_FR(idx_que,:,:,ifo)));
                    
                    
                     %% omtRSL_cov_act
                    if 0
                        options.rel = 'cov'; options.query = 'active';  
                        options.delta = 2^(-10 + iq);
                        tic; [MomtRSL_CA,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                        time_omtRSL_CA(idx_que,ifo) = toc;  que_omtRSL_CA(idx_que,ifo) = que;

                        for ita = 1:NT
                            [AP_tr, mAP_tr] = Test(MomtRSL_CA(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);    
                            [AP_te, mAP_te] = Test(MomtRSL_CA(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);                
                            AP_tr_omtRSL_CA(idx_que,ita,:,ifo) = AP_tr;  mAP_tr_omtRSL_CA(idx_que,ita,:,ifo) = mAP_tr;
                            AP_te_omtRSL_CA(idx_que,ita,:,ifo) = AP_te;  mAP_te_omtRSL_CA(idx_que,ita,:,ifo) = mAP_te;                  
                        end;
                        fprintf('omtRSL_cov_act:Que:%.3f,AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n',  que,...
                            mean2(AP_tr_omtRSL_CA(idx_que,:,:,ifo)), mean2(mAP_tr_omtRSL_CA(idx_que,:,:,ifo)),...
                            mean2(AP_te_omtRSL_CA(idx_que,:,:,ifo)), mean2(mAP_te_omtRSL_CA(idx_que,:,:,ifo)));

                        %% omtRSL_cov_rand
                        options.rel = 'cov'; options.query = 'random';
                        options.delta = queOMTCA(idx_que, ifo);
                        tic; [MOMTCR,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                        timOMTCR(idx_que,ifo) = toc;  queOMTCR(idx_que,ifo) = que;

                        for ita = 1:NT
                            [AP_tr, mAP_tr] = Test(MOMTCR(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);
                            [AP_te, mAP_te] = Test(MOMTCR(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                            AP_tr_OMTCR(idx_que,ita,:,ifo)  = AP_tr;  mAP_tr_OMTCR(idx_que,ita,:,ifo) = mAP_tr;
                            AP_te_OMTCR(idx_que,ita,:,ifo)  = AP_te;  mAP_te_OMTCR(idx_que,ita,:,ifo) = mAP_te;  
                        end;
                        fprintf('omtRSL_cov_rand:Que:%.3f,AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', que,mean2(AP_tr_OMTCR(idx_que,:,:,ifo)), mean2(mAP_tr_OMTCR(idx_que,:,:,ifo)), ...
                            mean2(AP_teOMTCR(idx_que,:,:,ifo)), mean2(mAP_teOMTCR(idx_que,:,:,ifo)));
                    end;
                end; % end of fold
            end; %end of query
             save(sprintf('%s/%s_omtRSL_AR_results.mat',results_dir, full_name),...
                'time_omtRSL_FA', 'que_omtRSL_FA', 'AP_tr_omtRSL_FA', 'mAP_tr_omtRSL_FA', 'AP_te_omtRSL_FA', 'mAP_te_omtRSL_FA',... 
                'time_omtRSL_FR', 'que_omtRSL_FR', 'AP_tr_omtRSL_FR', 'mAP_tr_omtRSL_FR', 'AP_te_omtRSL_FR', 'mAP_te_omtRSL_FR');
%                 'AP_tr_omtRSL_CA','mAP_tr_omtRSL_CA','AP_te_omtRSL_CA', 'mAP_teomtRSL_CA','timomtRSL_CA','queomtRSL_CA', 'AP_tr_OMTCR','mAP_tr_OMTCR','AP_teOMTCR', 'mAP_teOMTCR','timOMTCR','queOMTCR');
        end;   

        %% output 100$  query ratio witah different topk measure  
        for ik = 1:NTK        
            k = options.topk(ik);
            fileID = fopen(sprintf('%s/%s_fixed_top%d.txt',results_dir,full_name,k),'w');

            fprintf(fileID,'-----------------------------AP_tr-----------------------------------------\n');
            fprintf(fileID,'%s \t&\t mt-lmnn \t&\t ostRSL-rsl \t&\t ostRSL-rsl-global \t&\t omt-rsl-fixed\n', full_name);
            for ita = 1:NT
                fprintf(fileID,'%d &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f \\\\ \n',ita,...
                    mean(AP_tr_mtLMNN(ita,ik,:)), std(AP_tr_mtLMNN(ita,ik,:)),...
                    mean(AP_tr_ostRSL(ita,ik,:)),std(AP_tr_ostRSL(ita,ik,:)),...
                    mean(AP_tr_ostRSL_global(ik,:)),std(AP_tr_ostRSL_global(ik,:)),...
                    mean(AP_tr_omtRSL_FP(ita,ik,:)),std(AP_tr_omtRSL_FP(ita,ik,:)));
            end;
            fprintf(fileID,'---------------------------------------------------------------------------\n');

            fprintf(fileID,'\n-----------------------------mAP_tr-----------------------------------------\n');
            fprintf(fileID,'%s \t&\t mt-lmnn \t&\t ostRSL-rsl \t&\t ostRSL-rsl-global \t&\t omt-rsl-fixed\n', full_name);
            for ita = 1:NT
                fprintf(fileID,'%d &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f \\\\ \n',ita,...
                    mean(mAP_tr_mtLMNN(ita,ik,:)), std(mAP_tr_mtLMNN(ita,ik,:)),...
                    mean(mAP_tr_ostRSL(ita,ik,:)), std(mAP_tr_ostRSL(ita,ik,:)),...
                    mean(mAP_tr_ostRSL_global(ik,:)),std(mAP_tr_ostRSL_global(ik,:)),...
                    mean(mAP_tr_omtRSL_FP(ita,ik,:)),std(mAP_tr_omtRSL_FP(ita,ik,:)) );
            end;
            fprintf(fileID,'---------------------------------------------------------------------------\n');

            fprintf(fileID,'\n-----------------------------AP_te-----------------------------------------\n');
            fprintf(fileID,'%s \t&\t mt-lmnn \t&\t ostRSL-rsl \t&\t ostRSL-rsl-global \t&\t omt-rsl-fixed\n', full_name);
            for ita = 1:NT
                fprintf(fileID,'%d &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f \\\\ \n',ita,...
                    mean(AP_te_mtLMNN(ita,ik,:)), std(AP_te_mtLMNN(ita,ik,:)),...
                    mean(AP_te_ostRSL(ita,ik,:)), std(AP_te_ostRSL(ita,ik,:)),...
                    mean(AP_te_ostRSL_global(ik,:)), std(mAP_tr_ostRSL_global(ik,:)),...
                    mean(AP_te_omtRSL_FP(ita,ik,:)), std(AP_te_omtRSL_FP(ita,ik,:)));
            end;
            fprintf(fileID,'---------------------------------------------------------------------------\n');

            fprintf(fileID,'\n-----------------------------mAP_te-----------------------------------------\n');
            fprintf(fileID,'%s \t&\t mt-lmnn \t&\t ostRSL-rsl \t&\t ostRSL-rsl-global \t&\t omt-rsl-fixed\n', full_name);
            for ita = 1:NT
                fprintf(fileID,'%d &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f \\\\ \n',ita,...
                    mean(mAP_te_mtLMNN(ita,ik,:)), std(mAP_te_mtLMNN(ita,ik,:)),...
                    mean(mAP_te_ostRSL(ita,ik,:)), std(mAP_te_ostRSL(ita,ik,:)),...
                    mean(mAP_te_ostRSL_global(ik,:)),std(mAP_tr_ostRSL_global(ik,:)),...
                    mean(mAP_te_omtRSL_FP(ita,ik,:)),std(mAP_te_omtRSL_FP(ita,ik,:)));
            end;
            fprintf(fileID,'---------------------------------------------------------------------------\n');

            fclose(fileID);
        end;      
       
        plotTestRes(current_date,full_name)
    end;    
end; % end of datasets loop