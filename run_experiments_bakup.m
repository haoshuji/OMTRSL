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

options.tri.ratio = 1/4000; options.tri.hard= false;

options.tri.exi = false;   options.tri_time = 5;
options.mtlmnn.exi = false; options.mtlmnn.find_best_reg = false;
options.OST.exi = false;    options.find_OST_C = false;     options.OST_C = 0.01;
options.OSTG.exi = false;   options.find_OSTG_C = false;    options.OSTG_C = 1;
options.OMTFP.exi = false;  options.find_OMTFP_bC = 1;  options.OMTFP_b = 1; options.OMTFP_C =1;
options.OMTCP.exi = false;  options.find_OMTCP_C = false;   options.OMTCP_C = 10;
options.OMTAR.exi = false;
options.global.ins_exi = false;

data_names = { 'letter'};
norm_methods= {'length1' 'orignal' 'gaussian' 'scale01'}; 

options.b_candidate = [1,2,4,8,16,32,64,128,256,512];
options.C_candidate = [0.001,0.01,0.1,1,10,100,1000];
options.reg0_candidate = [1,100,300]; options.reg_candidate = [300,400,500];

%% options for omtLMNN

options.maxitaer_st = 500; options.maxitaer_mt = 500;
options.Kg = 3; options.mtlmnn.best_reg0 = 1; options.mtlmnn.best_reg = 500;

options.query = 'passive'; options.rel = 'fixed';
        
all_iq = (1:2:20);
num_que = length(all_iq);

options.delta = -1;

val_topk = (10);
%% all datasets
for idx_data=1:length(data_names)
    data_name = data_names{idx_data};
    for id = 1:length(norm_methods)
        full_name = strcat(data_name,'_',norm_methods{id});
        fprintf('Experiments on %s dataset\n',full_name);   
%         current_date = datestr(now,'mmmm-dd-yyyy-HH');
        current_date = date;

        results_dir = strcat(base_results_dir,'/',full_name,'/',current_date);
        
        if ~exist(results_dir, 'dir')  
            mkdir(results_dir);    
        end;

        load(sprintf('%s/%s_all_folds.mat',data_dir,full_name),'data');    

        NF = length(data);        NT = length(data(1).task);
        options.NT = NT;        options.data_name = full_name;        options.K = NT;
        
        %% allocating results space        
        %%time, que: NF size; others: NT * NTK * NF size
        timMTLMNN = zeros(NF,1);  APTrMTLMNN = zeros(NT,NTK,NF);  mAPTrMTLMNN = zeros(NT,NTK,NF); APTeMTLMNN  = zeros(NT,NTK,NF); mAPTeMTLMNN = zeros(NT,NTK,NF);

        timOST = zeros(NF,1);  APTrOST = zeros(NT,NTK,NF);   mAPTrOST = zeros(NT,NTK,NF); APTeOST = zeros(NT,NTK,NF); mAPTeOST = zeros(NT,NTK,NF);

        APTrOSTE  = zeros(NT,NTK,NF); mAPTrOSTE   = zeros(NT,NTK,NF); APTeOSTE  = zeros(NT,NTK,NF); mAPTeOSTE   = zeros(NT,NTK,NF);

        timOSTG = zeros(NF,1); APTrOSTG = zeros(NTK,NF); mAPTrOSTG = zeros(NTK,NF);  APTeOSTG = zeros(NTK,NF);  mAPTeOSTG = zeros(NTK,NF);

        timOMTFP = zeros(NF,1);  APTrOMTFP = zeros(NT,NTK,NF); mAPTrOMTFP = zeros(NT,NTK,NF); APTeOMTFP     = zeros(NT,NTK,NF); mAPTeOMTFP  = zeros(NT,NTK,NF); 

        timOMTCP = zeros(NF,1);  APTrOMTCP = zeros(NT,NTK,NF); mAPTrOMTCP = zeros(NT,NTK,NF); APTeOMTCP   = zeros(NT,NTK,NF);  mAPTeOMTCP   = zeros(NT,NTK,NF);  

        %%time, que: num_que *   NF, other: num_que * num_task * NTK * num_fol 
        timOMTFA      = zeros(num_que,NF);          queOMTFA  = zeros(num_que,NF);    APTrOMTFA    = zeros(num_que,NT,NTK,NF);   mAPTrOMTFA     = zeros(num_que,NT,NTK,NF);    
        APTeOMTFA    = zeros(num_que,NT,NTK,NF);   mAPTeOMTFA     = zeros(num_que,NT,NTK,NF);

        timOMTFR     = zeros(num_que,NF);          queOMTFR = zeros(num_que,NF);    APTrOMTFR   = zeros(num_que,NT,NTK,NF);   mAPTrOMTFR    = zeros(num_que,NT,NTK,NF);
        APTeOMTFR   = zeros(num_que,NT,NTK,NF);   mAPTeOMTFR    = zeros(num_que,NT,NTK,NF);

        timOMTCA      = zeros(num_que,NF);          queOMTCA  = zeros(num_que,NF);   APTrOMTCA    = zeros(num_que,NT,NTK,NF);   mAPTrOMTCA     = zeros(num_que,NT,NTK,NF);    
        APTeOMTCA    = zeros(num_que,NT,NTK,NF);   mAPTeOMTCA     = zeros(num_que,NT,NTK,NF);

        timOMTCR     = zeros(num_que,NF);          queOMTCR = zeros(num_que,NF);     APTrOMTCR   = zeros(num_que,NT,NTK,NF);   mAPTrOMTCR    = zeros(num_que,NT,NTK,NF);
        APTeOMTCR   = zeros(num_que,NT,NTK,NF);   mAPTeOMTCR    = zeros(num_que,NT,NTK,NF);
        
        %% generate global dataset for ostRSL_global
        if options.global.ins_exi
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
        if options.tri.exi
            load(sprintf('%s/%s_triplets',results_dir,full_name));
        else         
            fprintf('Generating triplets\n');     
            for ifo = 1:NF
                fprintf('%d-fold\n',ifo);                                
                tri_tmp = [];
                for ita=1:NT                    
                    triTF(ifo).task(ita).trip = GenerateTriplets(data(ifo).task(ita).x,data(ifo).task(ita).y,ita,options);                     
                    perm_ind = randperm(size(triTF(ifo).task(ita).trip,1));
                    triTF(ifo).task(ita).trip  = triTF(ifo).task(ita).trip(perm_ind,:);
                    tri_tmp =[tri_tmp;triTF(ifo).task(ita).trip];                        
                    fprintf('%d-task:%d\t',ita,size(triTF(ifo).task(ita).trip,1));
                end;
                triF(ifo).trip = tri_tmp; 	num_trip_all_tasks = size(triF(ifo).trip,1);
                fprintf('\nAll task:%d triplets\t',num_trip_all_tasks);                
                options.global_trip_num = num_trip_all_tasks;
                triGlobal(ifo).trip = GenerateTriplets(global_ins(ifo).x, global_ins(ifo).y, -1, options);
                fprintf('Global trip:%d triplets\n',size(triGlobal(ifo).trip,1)); 
            end;
            save(sprintf('%s/%s_triplets.mat',results_dir,full_name),'triTF','triF','triGlobal');
        end;
        

        %% mtLMNN    
        if options.mtlmnn.find_best_reg
            fprintf('\nFind best reg for mtlnn\n');
            mtlmnn_best_AP_va = 0;
            all_AP_va_MTLMNN = zeros(length(reg0_candidate),length(reg_candidate));
            for i=1:length(reg0_candidate)
                for j = 1:length(reg_candidate)
                    mtlmnn_AP_va_tmp = 0;
                    for ifo=1:NF               
                        mtD_tmp=MTlmnn_sp(data(ifo).task,options.Kg,'quiet',1,'regweight',reg_candidate(j), 'regweight0',reg0_candidate(i),'weight1',0.5,'maxitaer',options.maxitaer_mt); 
                        for ita=1:NT                            
                            [AP_va,mAP_va] = Test(mtD_tmp.bestMgen.M0+mtD_tmp.bestMgen.M{ita}, data(ifo).task(ita).xv, data(ifo).task(ita).yv, options); 
                            mtlmnn_AP_va_tmp = mtlmnn_AP_va_tmp + mean(AP_va(val_topk));
                        end; 
                    end;     
                    mtlmnn_AP_va_tmp = mtlmnn_AP_va_tmp/(NF*NT);
                    all_AP_va_MTLMNN(i,j) = mtlmnn_AP_va_tmp;

                    if mtlmnn_AP_va_tmp > mtlmnn_best_AP_va
                        mtlmnn_best_AP_va = mtlmnn_AP_va_tmp; options.mtlmnn.best_reg = reg_candidate(j); options.mtlmnn.best_reg0 = reg0_candidate(i);
                    end;
                    fprintf('current_AP_va:%.3f,\tcurrent_reg0:%.3f,\tcurrent_reg:%.3f\n', mtlmnn_AP_va_tmp,reg0_candidate(i),reg_candidate(j))
                end;
            end;
            fprintf('best_APTe:%.3f,\tbest_reg0:%.3f,\tbest_reg:%.3f\n', mtlmnn_best_AP_va,options.mtlmnn.best_reg0,options.mtlmnn.best_reg);
            best_reg = options.mtlmnn.best_reg; best_reg0 = options.mtlmnn.best_reg0;
            save(sprintf('%s/%s_MTLMNN_res_par.mat',results_dir,full_name), 'all_AP_va_MTLMNN');
        end; 
        
        if options.mtlmnn.exi
            load(sprintf('%s/%s_MTLMNN_res.mat',results_dir,full_name));   
            fprintf('mtLMNN:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n',mean2(APTrMTLMNN), mean2(mAPTrMTLMNN),mean2(APTeMTLMNN),mean2(mAPTeMTLMNN));
        else
            fprintf('mtLMNN:\t');          
            for ifo = 1:NF
                fprintf('%d-th fold\n',ifo);
                
                tic
                mtD=MTlmnn_sp(data(ifo).task,options.Kg,'quiet',1,'regweight',options.mtlmnn.best_reg, 'regweight0',options.mtlmnn.best_reg0,'weight1',0.5,'maxitaer',options.maxitaer_mt); 
                timMTLMNN(ifo) = toc;                            

                for ita=1:NT
                    mtLMNN(ifo).task(ita).M= mtD.bestMgen.M0+mtD.bestMgen.M{ita};                            
                end;                
                for ita=1:NT
                    [APTr,mAPTr] = Test(mtLMNN(ifo).task(ita).M,data(ifo).task(ita).x, data(ifo).task(ita).y, options); 
                    [APTe,mAPTe] = Test(mtLMNN(ifo).task(ita).M,data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);                     
                    APTrMTLMNN(ita,:,ifo) = APTr;  mAPTrMTLMNN(ita,:,ifo) = mAPTr; APTeMTLMNN(ita,:,ifo) = APTe;   mAPTeMTLMNN(ita,:,ifo) = mAPTe;
                end;                      
            end;
            fprintf('\n');
            fprintf('mtLMNN:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(APTrMTLMNN),mean2(mAPTrMTLMNN), mean2(APTeMTLMNN),mean2(mAPTeMTLMNN));
            save(sprintf('%s/%s_MTLMNN_res.mat',results_dir,full_name), 'timMTLMNN','APTrMTLMNN','mAPTrMTLMNN','APTeMTLMNN', 'mAPTeMTLMNN','mtLMNN');
        end;
        
        %% ostRSL_global
%         if options.find_OSTG_C
%             fprintf('\nFind best C for OSTG.\n');
%         	best_C = options.OSTG_C;
%         	best_AP_va_tmp = 0;
%         	all_AP_va_OSTG_C = [];
%         	for idx_C=1:length(options.C_candidate)
%         		options.OSTG_C = options.C_candidate(idx_C);
%         		AP_va_tmp = 0;
%         		for ifo=1:NF
%                     for ita = 1:NT                
% 	                    M_tmp = ostRSL(global_ins(ifo).x, global_ins(ifo).y, triGlobal(ifo).trip, options);	                        
% 	                    [AP_va,mAP_va] = Test(M_tmp, global_ins(ifo).xv, global_ins(ifo).yv, options);                     
% 	                    AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));
%                     end;
%         		end;
%         		AP_va_tmp = AP_va_tmp/(NF*NT);
%         		all_AP_va_OSTG_C(end+1) = AP_va_tmp;
% 
%         		if AP_va_tmp > best_AP_va_tmp
%         			best_C = options.OSTG_C;
%         			best_AP_va_tmp = AP_va_tmp;
%         		end;
%         		fprintf('Current_C:%f, \t Current_AP_va:%f\n',options.OSTG_C,AP_va_tmp);
%         	end;
%             options.OSTG_C = best_C;
%             fprintf('Best_C:%f, \t Best_AP_va:%f\n',best_C,best_AP_va_tmp);
%         	save(sprintf('%s/%s_OSTG_res_par.mat',results_dir,full_name),'all_AP_va_OSTG_C');        	
%         end;
%         
%         if options.OSTG.exi
%             load(sprintf('%s/%s_OSTG_res.mat',results_dir,full_name));       
%             fprintf('ostRSL_global:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(APTrOSTG),mean2(mAPTrOSTG), mean2(APTeOSTG),mean2(mAPTeOSTG));
%         else
%             fprintf('ostRSL_global:\t');
%             for ifo = 1:NF
%                 fprintf('%d\t',ifo);
%                 tic; ostRSL_global_results(ifo).M = ostRSL(global_ins(ifo).x, global_ins(ifo).y, triGlobal(ifo).trip, options);                
%                 timOSTG(ifo) = toc; 
%                 
% %                 [APTr, mAPTr] = Test(ostRSL_global_results(ifo).M, global_ins(ifo).x, global_ins(ifo).y, options);                    
%                 APTr = 0; mAPTr = 0;
%                 [APTe, mAPTe] = Test(ostRSL_global_results(ifo).M, global_ins(ifo).xt, global_ins(ifo).yt, options);
%                 APTrOSTG(:,ifo) = APTr;     mAPTrOSTG(:,ifo) = mAPTr; APTeOSTG(:,ifo) = APTe;     mAPTeOSTG(:,ifo) = mAPTe;  
% 
%             end;
%             fprintf('\n');
%             fprintf('ostRSL_global:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n',mean2(APTrOSTG),mean2(mAPTrOSTG), mean2(APTeOSTG),mean2(mAPTeOSTG));
%             save(sprintf('%s/%s_OSTG_res',results_dir,full_name), 'timOSTG','APTrOSTG','mAPTrOSTG','APTeOSTG', 'mAPTeOSTG','ostRSL_global_results');
%         end;  
        
        %% M = eye(d,d);
        d = size(data(1).task(1).x,1);
        for ifo = 1:NF
            fprintf('%d\t',ifo);            
            for ita = 1:NT
                [APTr, mAPTr] = Test(eye(d,d), data(ifo).task(ita).x, data(ifo).task(ita).y, options);                    
                [APTe, mAPTe] = Test(eye(d,d), data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                APTrOSTE(ita,:,ifo) = APTr;     mAPTrOSTE(ita,:,ifo) = mAPTr; APTeOSTE(ita,:,ifo) = APTe;     mAPTeOSTE(ita,:,ifo) = mAPTe;  
            end;            
        end;
        fprintf('\nostRSL_eye:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(APTrOSTE), mean2(mAPTrOSTE), mean2(APTeOSTE),mean2(mAPTeOSTE));
        save(sprintf('%s/%s_OSTE_res.mat',results_dir,full_name), 'APTrOSTE', 'mAPTrOSTE', 'APTeOSTE', 'mAPTeOSTE');
        
        %% ostRSL
        if options.find_OST_C
            fprintf('\nFind best C for ostRSL.\n');
        	best_C = options.OST_C;
        	best_AP_va_tmp = 0;
        	all_AP_va_OST_C = [];
        	for idx_C=1:length(options.C_candidate)
        		options.OST_C = options.C_candidate(idx_C);
        		AP_va_tmp = 0;
        		for ifo=1:NF
                    for ita = 1:NT                
	                    M_tmp = ostRSL([data(ifo).task(ita).x], [data(ifo).task(ita).y], triTF(ifo).task(ita).trip, options);
	                    [AP_va,mAP_va] = Test(M_tmp,data(ifo).task(ita).xv, data(ifo).task(ita).yv, options);                     
	                    AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));
                    end;
        		end;
        		AP_va_tmp = AP_va_tmp/(NF*NT); all_AP_va_OST_C(end+1) = AP_va_tmp;

        		if AP_va_tmp > best_AP_va_tmp
        			best_C = options.OST_C; best_AP_va_tmp = AP_va_tmp;
        		end;
        		fprintf('Current_C:%f, \t Current_AP_va:%f\n',options.OST_C,AP_va_tmp);
        	end;
        	options.OST_C = best_C;
            fprintf('Best_C:%f, \t Best_AP_va:%f\n',best_C,best_AP_va_tmp);
            save(sprintf('%s/%s_OST_res_par',results_dir,full_name),'all_AP_va_OST_C');
        end;

        if options.OST.exi
            load(sprintf('%s/%s_OST_res.mat',results_dir,full_name));       
            fprintf('ostRSL:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(APTrOST),mean2(mAPTrOST), mean2(APTeOST),mean2(mAPTeOST));
        else
            fprintf('ostRSL:\t');
            for ifo = 1:NF
                fprintf('%d\t',ifo);
                tic;                
                for ita = 1:NT                
                    ostRSL_results(ifo).task(ita).M = ostRSL( data(ifo).task(ita).x, data(ifo).task(ita).y, triTF(ifo).task(ita).trip, options);
                end;
                timOST(ifo) = toc; 

                for ita = 1:NT
                    [APTr, mAPTr] = Test(ostRSL_results(ifo).task(ita).M,data(ifo).task(ita).x, data(ifo).task(ita).y, options);                    
                    [APTe, mAPTe] = Test(ostRSL_results(ifo).task(ita).M,data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                    APTrOST(ita,:,ifo) = APTr;     mAPTrOST(ita,:,ifo) = mAPTr; APTeOST(ita,:,ifo) = APTe;     mAPTeOST(ita,:,ifo) = mAPTe;  
                end;            
            end;
            fprintf('\n');
            fprintf('ostRSL:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(APTrOST),mean2(mAPTrOST), mean2(APTeOST),mean2(mAPTeOST));
            save(sprintf('%s/%s_OST_res.mat',results_dir,full_name), 'timOST','APTrOST','mAPTrOST','APTeOST','mAPTeOST','ostRSL_results');
        end;    
        
        %% omtRSL_cov_pa
        if options.find_OMTCP_C
            fprintf('\nFind best C for omtRSL_cov_pa.\n');
        	best_C = options.OMTCP_C;
        	best_AP_va_tmp = 0;
        	all_AP_va_OMTCP_C = [];
            options.rel = 'cov';
        	for idx_C=1:length(options.C_candidate)
        		options.OMTCP_C = options.C_candidate(idx_C);
        		AP_va_tmp = 0;
                for ifo =1:NF                    
                    [MOMTFP_tmp,que] = omtRSL(data(ifo).task,triF(ifo).trip, options);                            
                    for ita = 1:NT
                        [AP_va, mAP_va] = Test(MOMTFP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, options);
                        AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));                            
                    end;                            
                end;

        		AP_va_tmp = AP_va_tmp/(NF*NT);
        		all_AP_va_OMTCP_C(end+1) = AP_va_tmp;

        		if AP_va_tmp > best_AP_va_tmp
        			best_C = options.OMTCP_C;
        			best_AP_va_tmp = AP_va_tmp;
        		end;
        		fprintf('Current_C:%f, \t Current_AP_va:%f\n',options.OMTCP_C,AP_va_tmp);
        	end;
            fprintf('Best_C:%f, \t Best_AP_va:%f\n',best_C,best_AP_va_tmp);
        	save(sprintf('%s/%s_OMTCP_res_par.mat',results_dir,full_name),'all_AP_va_OMTCP_C');
        	options.OMTCP_C = best_C;
        end;

        if options.OMTCP.exi
            load(sprintf('%s/%s_OMTCP_res.mat',results_dir,full_name));       
            fprintf('omtRSL_cov_pa:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n',mean2(APTrOMTCP),mean2(mAPTrOMTCP), mean2(APTeOMTCP),mean2(mAPTeOMTCP));
        else
            options.rel = 'cov';
            fprintf('omtRSL_cov_pa:\t');
            for ifo = 1:NF
                fprintf('%d\t',ifo);
                tic;
                [omtRSL_cov_pa_results(ifo).task,que_useless] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                timOMTFP(ifo) = toc; 

                for ita = 1:NT
                    [APTr, mAPTr] = Test(omtRSL_cov_pa_results(ifo).task(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);   
                    [APTe, mAPTe] = Test(omtRSL_cov_pa_results(ifo).task(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                    APTrOMTCP(ita,:,ifo) = APTr;    mAPTrOMTCP(ita,:,ifo) = mAPTr; APTeOMTCP(ita,:,ifo) = APTe;    mAPTeOMTCP(ita,:,ifo) = mAPTe;  
                end;
            end;  
            fprintf('\n');
            fprintf('omtRSL_cov_pa:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(APTrOMTCP),mean2(mAPTrOMTCP), mean2(APTeOMTCP),mean2(mAPTeOMTCP));
            save(sprintf('%s/%s_OMTCP_res',results_dir,full_name), 'timOMTCP','APTrOMTCP', 'mAPTrOMTCP','APTeOMTCP','mAPTeOMTCP','omtRSL_cov_pa_results');
        end;    
        
        %% omtRSL_fix_pa
        if options.find_OMTFP_bC == 2     
            options.rel = 'fixed';
            fprintf('\nFind best b and C on AP_va\n');
            best_b = options.OMTFP_b;
            best_C = options.OMTFP_C;
            
            best_AP_va_tmp = 0;
            all_AP_va_OMTFP = zeros(length(options.b_candidate),length(options.C_candidate));            
            
            for idx_b = 1:length(options.b_candidate)
                for idx_C = 1:length(options.C_candidate)
                    options.OMTFP_b = options.b_candidate(idx_b);
                    options.OMTFP_C = options.C_candidate(idx_C);
                    
                    AP_va_tmp = 0; 
                    for ifo =1:NF                    
                        [MOMTFP_tmp,que] = omtRSL(data(ifo).task,...
                            triF(ifo).trip, options);                            
                        for ita = 1:NT
                            [AP_va, mAP_va] = Test(MOMTFP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, options);
                            AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));                            
                        end;                            
                    end;
                    AP_va_tmp = AP_va_tmp/(NF*NT);                    all_AP_va_OMTFP(idx_b,idx_C) = AP_va_tmp;

                    if AP_va_tmp > best_AP_va_tmp
                        best_AP_va_tmp = AP_va_tmp; best_b = options.OMTFP_b; best_C = options.OMTFP_C;
                    end;
                    fprintf('Current_b:%f,\tCurrent_C:%f,\tCurrent_APTe:%f\n',options.OMTFP_b, options.OMTFP_C, AP_va_tmp);
                end;
            end;            
            options.OMTFP_b = best_b; options.OMTFP_C = best_C;
            fprintf('Best_b:%f,\tBest_C:%f,\tBest_AP_va:%f\n',best_b,best_C,best_AP_va_tmp);            
            save(sprintf('%s/%s_OMTFP_res_par',results_dir,full_name),'all_AP_va_OMTFP');            
            fprintf('OMTFP_C:%f\t,OMTFP_b:%f\n',options.OMTFP_C,options.OMTFP_b);
        elseif options.find_OMTFP_bC == 1
            options.rel = 'fixed';
            fprintf('\nFind best b for omtRSL alone and set OMTFP_C as OST_C\n');
            options.OMTFP_C = options.OST_C;

            best_b = options.OMTFP_b;
            best_AP_va_tmp = 0;
            all_AP_va_OMTFP = zeros(length(options.b_candidate),1);            
            
            for idx_b = 1:length(options.b_candidate)
                options.OMTFP_b = options.b_candidate(idx_b);
                AP_va_tmp = 0; 
               
                for ifo =1:NF                    
                    [MOMTFP_tmp,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);                            
                    for ita = 1:NT
                        [AP_va, mAP_va] = Test(MOMTFP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, options);
                        AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));                            
                    end;                            
                end;

                AP_va_tmp = AP_va_tmp/(NF*NT); all_AP_va_OMTFP(idx_b) = AP_va_tmp;

                if AP_va_tmp > best_AP_va_tmp
                    best_AP_va_tmp = AP_va_tmp; best_b = options.OMTFP_b;
                end;
                fprintf('Current_b:%f,\tCurrent_APTe:%f\n',options.OMTFP_b, AP_va_tmp);
            end; 
             options.OMTFP_b = best_b;           
            fprintf('Best_b:%f,\tBest_AP_va:%f\n',best_b, best_AP_va_tmp);            
            save(sprintf('%s/%s_OMTFP_res_par',results_dir,full_name),'all_AP_va_OMTFP');            
        else
            fprintf('Using Existing OMTFP_C:%f and OMTFP_b:%f\n',options.OMTFP_C,options.OMTFP_b);
        end;

        fprintf('OMTFP_C:%f\t,OMTFP_b:%f\n',options.OMTFP_C,options.OMTFP_b);
        
        save(sprintf('%s/%s_options.mat',results_dir,full_name),'options');

        if options.OMTFP.exi
            load(sprintf('%s/%s_OMTFP_res.mat',results_dir,full_name));
            fprintf('omtRSL:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(APTrOMTFP),mean2(mAPTrOMTFP), mean2(APTeOMTFP),mean2(mAPTeOMTFP));
        else
            options.rel = 'fixed';
            fprintf('omtRSL_fix_pa:\t');
            for ifo = 1:NF
                fprintf('%d\t',ifo);
                tic;
                [omtRSL_results(ifo).task,que_useless] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                timOMTFP(ifo) = toc; 

                for ita = 1:NT
                    [APTr, mAPTr] = Test(omtRSL_results(ifo).task(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);   
                    [APTe, mAPTe] = Test(omtRSL_results(ifo).task(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                    APTrOMTFP(ita,:,ifo) = APTr;    mAPTrOMTFP(ita,:,ifo) = mAPTr; APTeOMTFP(ita,:,ifo) = APTe;    mAPTeOMTFP(ita,:,ifo) = mAPTe;  
                end;
            end;        
            fprintf('\n');
            fprintf('omtRSL:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(APTrOMTFP),mean2(mAPTrOMTFP), mean2(APTeOMTFP),mean2(mAPTeOMTFP));
            save(sprintf('%s/%s_OMTFP_res.mat',results_dir,full_name), 'timOMTFP','APTrOMTFP', 'mAPTrOMTFP','APTeOMTFP', 'mAPTeOMTFP','omtRSL_results');
        end;

        %% omtRSL_fix_act and omtRSL_fix_rand
        if options.OMTAR.exi
            load(sprintf('%s/%s_OMTAR_res.mat',results_dir, full_name));
        else
            for idx_que = 1:num_que
                iq = all_iq(idx_que);
                fprintf('%d-th query\n',iq);            
                for ifo = 1:NF
                    fprintf('%d-th fold\n',ifo);

                    %% omtRSL_fix_act
                    options.rel = 'fixed'; options.query = 'active';  
                    tic; options.delta = 2^(-10 + iq);
                    [MOMTFA,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                    timOMTFA(idx_que,ifo) = toc;  queOMTFA(idx_que,ifo) = que;

                    for ita = 1:NT
                        [APTr, mAPTr] = Test(MOMTFA(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);    
                        [APTe, mAPTe] = Test(MOMTFA(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);                
                        APTrOMTFA(idx_que,ita,:,ifo) = APTr;  mAPTrOMTFA(idx_que,ita,:,ifo) = mAPTr; APTeOMTFA(idx_que,ita,:,ifo) = APTe;  mAPTeOMTFA(idx_que,ita,:,ifo) = mAPTe;                  
                    end;
                    fprintf('omtRSL_fix_act:Que:%.3f,APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', que,mean2(APTrOMTFA(idx_que,:,:,ifo)), mean2(mAPTrOMTFA(idx_que,:,:,ifo)),...
                        mean2(APTeOMTFA(idx_que,:,:,ifo)), mean2(mAPTeOMTFA(idx_que,:,:,ifo)));

                    %% omtRSL_fix_rand
                    options.rel = 'fixed'; options.query = 'random';
                    options.delta = queOMTFA(idx_que, ifo);
                    tic; [MOMTFR,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                    timOMTFR(idx_que,ifo) = toc;  queOMTFR(idx_que,ifo) = que;

                    for ita = 1:NT
                        [APTr, mAPTr] = Test(MOMTFR(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);
                        [APTe, mAPTe] = Test(MOMTFR(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                        APTrOMTFR(idx_que,ita,:,ifo)  = APTr;  mAPTrOMTFR(idx_que,ita,:,ifo) = mAPTr; APTeOMTFR(idx_que,ita,:,ifo)  = APTe;  mAPTeOMTFR(idx_que,ita,:,ifo) = mAPTe;  
                    end;
                    fprintf('omtRSL_fix_rand:Que:%.3f,APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', que,mean2(APTrOMTFR(idx_que,:,:,ifo)), mean2(mAPTrOMTFR(idx_que,:,:,ifo)),...
                        mean2(APTeOMTFR(idx_que,:,:,ifo)), mean2(mAPTeOMTFR(idx_que,:,:,ifo)));
                    
                    
                     %% omtRSL_cov_act
                    options.rel = 'cov'; options.query = 'active';  
                    options.delta = 2^(-10 + iq);
                    tic; [MOMTCA,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                    timOMTCA(idx_que,ifo) = toc;  queOMTCA(idx_que,ifo) = que;

                    for ita = 1:NT
                        [APTr, mAPTr] = Test(MOMTCA(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);    
                        [APTe, mAPTe] = Test(MOMTCA(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);                
                        APTrOMTCA(idx_que,ita,:,ifo) = APTr;  mAPTrOMTCA(idx_que,ita,:,ifo) = mAPTr;
                        APTeOMTCA(idx_que,ita,:,ifo) = APTe;  mAPTeOMTCA(idx_que,ita,:,ifo) = mAPTe;                  
                    end;
                    fprintf('omtRSL_cov_act:Que:%.3f,APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n',  que,mean2(APTrOMTCA(idx_que,:,:,ifo)),  mean2(mAPTrOMTCA(idx_que,:,:,ifo)),...
                        mean2(APTeOMTCA(idx_que,:,:,ifo)), mean2(mAPTeOMTCA(idx_que,:,:,ifo)));

                    %% omtRSL_cov_rand
                    options.rel = 'cov'; options.query = 'random';
                    options.delta = queOMTCA(idx_que, ifo);
                    tic; [MOMTCR,que] = omtRSL(data(ifo).task, triF(ifo).trip, options);
                    timOMTCR(idx_que,ifo) = toc;  queOMTCR(idx_que,ifo) = que;

                    for ita = 1:NT
                        [APTr, mAPTr] = Test(MOMTCR(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options);
                        [APTe, mAPTe] = Test(MOMTCR(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);
                        APTrOMTCR(idx_que,ita,:,ifo)  = APTr;  mAPTrOMTCR(idx_que,ita,:,ifo) = mAPTr;
                        APTeOMTCR(idx_que,ita,:,ifo)  = APTe;  mAPTeOMTCR(idx_que,ita,:,ifo) = mAPTe;  
                    end;
                    fprintf('omtRSL_cov_rand:Que:%.3f,APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', que,mean2(APTrOMTCR(idx_que,:,:,ifo)), mean2(mAPTrOMTCR(idx_que,:,:,ifo)), ...
                        mean2(APTeOMTCR(idx_que,:,:,ifo)), mean2(mAPTeOMTCR(idx_que,:,:,ifo)));
                end; % end of fold
            end; %end of query
             save(sprintf('%s/%s_OMTAR_res.mat',results_dir, full_name),...
                'APTrOMTFA','mAPTrOMTFA','APTeOMTFA', 'mAPTeOMTFA','timOMTFA','queOMTFA', 'APTrOMTFR','mAPTrOMTFR','APTeOMTFR', 'mAPTeOMTFR','timOMTFR','queOMTFR',...
                'APTrOMTCA','mAPTrOMTCA','APTeOMTCA', 'mAPTeOMTCA','timOMTCA','queOMTCA', 'APTrOMTCR','mAPTrOMTCR','APTeOMTCR', 'mAPTeOMTCR','timOMTCR','queOMTCR');
        end;   

        %% output 100$  query ratio witah different topk measure  
        for ik = 1:NTK        
            k = options.topk(ik);
            fileID = fopen(sprintf('%s/%s_fixed_top%d.txt',results_dir,full_name,k),'w');

            fprintf(fileID,'-----------------------------APTr-----------------------------------------\n');
            fprintf(fileID,'%s \t&\t mt-lmnn \t&\t ost-rsl \t&\t ost-rsl-global \t&\t omt-rsl-fixed\n', full_name);
            for ita = 1:NT
                fprintf(fileID,'%d &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f \\\\ \n',ita,...
                    mean(APTrMTLMNN(ita,ik,:)), std(APTrMTLMNN(ita,ik,:)),...
                    mean(APTrOST(ita,ik,:)),std(APTrOST(ita,ik,:)),...
                    mean(APTrOSTG(ik,:)),std(APTrOSTG(ik,:)),...
                    mean(APTrOMTFP(ita,ik,:)),std(APTrOMTFP(ita,ik,:)));
            end;
            fprintf(fileID,'---------------------------------------------------------------------------\n');

            fprintf(fileID,'\n-----------------------------mAPTr-----------------------------------------\n');
            fprintf(fileID,'%s \t&\t mt-lmnn \t&\t ost-rsl \t&\t ost-rsl-global \t&\t omt-rsl-fixed\n', full_name);
            for ita = 1:NT
                fprintf(fileID,'%d &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f \\\\ \n',ita,...
                    mean(mAPTrMTLMNN(ita,ik,:)), std(mAPTrMTLMNN(ita,ik,:)),...
                    mean(mAPTrOST(ita,ik,:)), std(mAPTrOST(ita,ik,:)),...
                    mean(mAPTrOSTG(ik,:)),std(mAPTrOSTG(ik,:)),...
                    mean(mAPTrOMTFP(ita,ik,:)),std(mAPTrOMTFP(ita,ik,:)) );
            end;
            fprintf(fileID,'---------------------------------------------------------------------------\n');

            fprintf(fileID,'\n-----------------------------APTe-----------------------------------------\n');
            fprintf(fileID,'%s \t&\t mt-lmnn \t&\t ost-rsl \t&\t ost-rsl-global \t&\t omt-rsl-fixed\n', full_name);
            for ita = 1:NT
                fprintf(fileID,'%d &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f \\\\ \n',ita,...
                    mean(APTeMTLMNN(ita,ik,:)), std(APTeMTLMNN(ita,ik,:)),...
                    mean(APTeOST(ita,ik,:)), std(APTeOST(ita,ik,:)),...
                    mean(APTeOSTG(ik,:)), std(mAPTrOSTG(ik,:)),...
                    mean(APTeOMTFP(ita,ik,:)), std(APTeOMTFP(ita,ik,:)));
            end;
            fprintf(fileID,'---------------------------------------------------------------------------\n');

            fprintf(fileID,'\n-----------------------------mAPTe-----------------------------------------\n');
            fprintf(fileID,'%s \t&\t mt-lmnn \t&\t ost-rsl \t&\t ost-rsl-global \t&\t omt-rsl-fixed\n', full_name);
            for ita = 1:NT
                fprintf(fileID,'%d &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f &\t %.3f$\\pm$%.3f \\\\ \n',ita,...
                    mean(mAPTeMTLMNN(ita,ik,:)), std(mAPTeMTLMNN(ita,ik,:)),...
                    mean(mAPTeOST(ita,ik,:)), std(mAPTeOST(ita,ik,:)),...
                    mean(mAPTeOSTG(ik,:)),std(mAPTrOSTG(ik,:)),...
                    mean(mAPTeOMTFP(ita,ik,:)),std(mAPTeOMTFP(ita,ik,:)));
            end;
            fprintf(fileID,'---------------------------------------------------------------------------\n');

            fclose(fileID);
        end;      
       
%         plotTestRes(current_date,full_name)
    end;
end; % end of datasets loop