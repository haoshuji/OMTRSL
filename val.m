clear;
% clc;
setpaths;
cur_path = pwd;
data_dir = strcat(cur_path,'/datasets');
base_results_dir = strcat(cur_path,'/results');

%% all opt
opt.use_cov = false;
opt.topk = [1,2,3,4,5,6,7,8,9,10];
ntopk = length(opt.topk);

opt.use_existing_triplets = false;   
opt.MTLMNN.use_existing = true;
opt.maxitaer_st = 500; opt.maxitaer_mt = 500;
opt.Kg = 3; opt.mtlmnn.best_reg0 = 1; opt.mtlmnn.best_reg = 500;

opt.ostRSL.use_existing = false; opt.find_ostRSL_C = true;
opt.OMTFP.use_existing = false; opt.find_omtRSL_FP_bC = 1;
opt.OMTCP.use_existing = false; opt.find_omtRSL_CP_C = true;

opt.ostRSL_C = 1;
opt.omtRSL_FP_C = 1;
opt.omtRSL_FP_b = 2;
opt.omtRSL_CP_C = 1;
opt.query = 'passive';
opt.rel = 'fixed';
opt.b_candidate = [1,2,4,8,16,32,64,128,256,512];
opt.C_candidate = [0.001,0.01,0.1,1,10,100,1000];
opt.reg0_candidate = [1,100,300]; opt.reg_candidate = [300,400,500];
%% triplets
opt.tri.ratio = 1/4000;
opt.to_generate_all_ok_triplets = false; 
opt.tri.hard = false;
opt.output_each_task_varied = false;
        
all_iq = [1:2:20];
num_que = length(all_iq);

opt.delta = -1;

ntr = 10;

val_topk = [10];
%% all datasets
data_names = {'letter'};
% norm_methods= {'orignal' 'length1' 'scale01' 'gaussian'};
norm_methods = {'orignal'};
for idx_data=1:length(data_names)
    
    data_name = data_names{idx_data};
    for id = 1:length(norm_methods)
        full_name = strcat(data_name,'_',norm_methods{id});
        fprintf('%s\n',full_name);   
%         current_date = datestr(now,'mmmm-dd-yyyy-HH');
        current_date =date;
        results_dir = strcat(base_results_dir,'/',full_name,'/',current_date);
        if ~exist(results_dir, 'dir')
            mkdir(results_dir);
        end;

        load(sprintf('%s/%s_all_folds.mat',data_dir,full_name),'data');    

        NF = length(data);
        NT = length(data(1).task);
        opt.NT = NT;
        opt.data_name = full_name;
        opt.K = NT;
        
        
        %% allocating results space        
        
        tim_MTLMNN              = zeros(NF,1);  
        AP_tr_MTLMNN           = zeros(NT,ntopk,NF);           mAP_tr_MTLMNN             = zeros(NT,ntopk,NF); 
        AP_va_MTLMNN            = zeros(NT,ntopk,NF);           mAP_va_MTLMNN             = zeros(NT,ntopk,NF);
        
        
        tim_OST              = zeros(ntr,NF,1);  
        AP_tr_OST            = zeros(ntr,NT,ntopk,NF);           mAP_tr_OST             = zeros(ntr,NT,ntopk,NF); 
        AP_va_OST            = zeros(ntr,NT,ntopk,NF);           mAP_va_OST             = zeros(ntr,NT,ntopk,NF);
        
        AP_tr_OST_eye            = zeros(NT,ntopk,NF);           mAP_tr_OST_eye             = zeros(NT,ntopk,NF); 
        AP_va_OST_eye            = zeros(NT,ntopk,NF);           mAP_va_OST_eye             = zeros(NT,ntopk,NF);
        
        tim_OMTFP       = zeros(ntr,NF,1);                         
        AP_tr_OMTFP     = zeros(ntr,NT,ntopk,NF);           mAP_tr_OMTFP      = zeros(ntr,NT,ntopk,NF);   
        AP_va_OMTFP     = zeros(ntr,NT,ntopk,NF);           mAP_va_OMTFP      = zeros(ntr,NT,ntopk,NF); 
        
        tim_OMTCP       = zeros(ntr,NF,1);                         
        AP_tr_OMTCP     = zeros(ntr,NT,ntopk,NF);           mAP_tr_OMTCP      = zeros(ntr,NT,ntopk,NF);   
        AP_va_OMTCP     = zeros(ntr,NT,ntopk,NF);           mAP_va_OMTCP      = zeros(ntr,NT,ntopk,NF);  

        %% mtlmnn
%         if opt.MTLMNN.use_existing
%             load(sprintf('%s/%s_MTLMNN_res.mat',results_dir,full_name));   
%             fprintf('mtLMNN:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n',mean2(AP_tr_MTLMNN),...
%              mean2(mAP_tr_MTLMNN),mean2(AP_va_MTLMNN),mean2(mAP_va_MTLMNN));
%         else
%             fprintf('mtLMNN:\t');          
%             for ifo = 1:NF
%                 fprintf('%d-th fold\n',ifo);
%                 
%                 tic
%                 mtD=MTlmnn_sp(data(ifo).task,opt.Kg,'quiet',1,'regweight',opt.mtlmnn.best_reg, 'regweight0',opt.mtlmnn.best_reg0,'weight1',0.5,'maxitaer',opt.maxitaer_mt); 
%                 tim_MTLMNN(ifo) = toc;                            
% 
%                 for ita=1:NT
%                     mtLMNN(ifo).task(ita).M= mtD.bestMgen.M0+mtD.bestMgen.M{ita};                            
%                 end;                
%                 for ita=1:NT
%                     [APTr,mAPTr] = Test(mtLMNN(ifo).task(ita).M,data(ifo).task(ita).x, data(ifo).task(ita).y, opt); 
%                     [APVa,mAPVa] = Test(mtLMNN(ifo).task(ita).M,data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);                     
%                     AP_tr_MTLMNN(ita,:,ifo) = APTr;  mAP_tr_MTLMNN(ita,:,ifo) = mAPTr;
%                     AP_va_MTLMNN(ita,:,ifo) = APVa;  mAP_va_MTLMNN(ita,:,ifo) = mAPVa;
%                 end;                      
%             end;
%             fprintf('\n');
%             fprintf('mtLMNN:APTr:%.3f,mAPTr:%.3f,APTe:%.3f,mAPTe:%.3f\n', mean2(AP_tr_MTLMNN),mean2(mAP_tr_MTLMNN), mean2(AP_va_MTLMNN),mean2(mAP_va_MTLMNN));
%             save(sprintf('%s/%s_MTLMNN_res.mat',results_dir,full_name), 'tim_MTLMNN','AP_tr_MTLMNN',...
%                 'mAP_tr_MTLMNN','AP_va_MTLMNN', 'mAP_va_MTLMNN');
%         end;
       
        %% M = eye(d,d);
        d = size(data(1).task(1).x,1);
        for ifo = 1:NF
            fprintf('%d\t',ifo);            
            for ita = 1:NT
                [AP_tr, mAP_tr] = Test(eye(d,d), data(ifo).task(ita).x, data(ifo).task(ita).y, opt);                    
                [AP_va, mAP_va] = Test(eye(d,d), data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);
                AP_tr_OST_eye(ita,:,ifo) = AP_tr;     mAP_tr_OST_eye(ita,:,ifo) = mAP_tr;
                AP_va_OST_eye(ita,:,ifo) = AP_va;     mAP_va_OST_eye(ita,:,ifo) = mAP_va;  
            end;            
        end;
        fprintf('\n ostRSL_eye:AP_tr:%.3f,mAP_tr:%.3f,AP_va:%.3f,mAP_va:%.3f\n', mean2(AP_tr_OST_eye), mean2(mAP_tr_OST_eye), mean2(AP_va_OST_eye),mean2(mAP_va_OST_eye));
        save(sprintf('%s/%s_OST_eye_results',results_dir,full_name), 'AP_tr_OST_eye', 'mAP_tr_OST_eye', 'AP_va_OST_eye', 'mAP_va_OST_eye');
        


        %% generate triplets
        
        num_tri_list = [];
        
        fprintf('Generating triplets\n');     
        for itr = 1:ntr
            opt.tri_time = itr;
            fprintf('%d-tri\n',itr);
            for ifo = 1:NF
                fprintf('%d-fold\n',ifo);                                
                
                all_task_trip = [];
                for ita=1:NT                    
                    tri_task(ifo).task(ita).trip = GenerateTriplets(data(ifo).task(ita).x,data(ifo).task(ita).y,ita,opt); 

                    ntr_tmp = size(tri_task(ifo).task(ita).trip,1);
                    perm_ind = randperm(ntr_tmp);
                    tri_task(ifo).task(ita).trip  = tri_task(ifo).task(ita).trip(perm_ind,:);
                    all_task_trip =[all_task_trip;tri_task(ifo).task(ita).trip];                        
                    fprintf('%d-task:%d\t',ita,size(tri_task(ifo).task(ita).trip,1));
                end;
                
                tri_fold(ifo).trip = all_task_trip;                 
                fprintf('\nAll task:%d triplets\t',size(tri_fold(ifo).trip,1));     
                
            end;       
            num_tri_list(end+1) = size(tri_fold(1).trip,1);            
            

            if opt.find_ostRSL_C && itr == 1
                fprintf('\nFind best C for ostRSL.\n');
                best_C = opt.ostRSL_C;
                best_AP_va_tmp = 0;
                all_AP_va_ostRSL_C = [];
                for idx_C=1:length(opt.C_candidate)
                    opt.ostRSL_C = opt.C_candidate(idx_C);
                    AP_va_tmp = 0;
                    for ifo=1:NF
                        for ita = 1:NT                
                            M_tmp = ostRSL([data(ifo).task(ita).x], [data(ifo).task(ita).y], tri_task(ifo).task(ita).trip, opt);
                            [AP_va,mAP_va] = Test(M_tmp,data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);                     
                            AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));
                        end;
                    end;
                    AP_va_tmp = AP_va_tmp/(NF*NT); all_AP_va_ostRSL_C(end+1) = AP_va_tmp;

                    if AP_va_tmp > best_AP_va_tmp
                        best_C = opt.ostRSL_C; best_AP_va_tmp = AP_va_tmp;
                    end;
                    fprintf('Current_C:%f, \t Current_AP_va:%f\n',opt.ostRSL_C,AP_va_tmp);
                end;
                opt.ostRSL_C = best_C;
                fprintf('Best_C:%f, \t Best_AP_va:%f\n',best_C,best_AP_va_tmp);
                save(sprintf('%s/%s_OST_res_par',results_dir,full_name),'all_AP_va_ostRSL_C','best_C');
            end;

            %% ostRSL   
            if opt.ostRSL.use_existing
                load(sprintf('%s/%s_OST_results.mat',results_dir,full_name));       
                fprintf('ostRSL:AP_tr:%.3f,mAP_tr:%.3f,AP_va:%.3f,mAP_va:%.3f\n',...
                    mean2(AP_tr_OST),mean2(mAP_tr_OST),...
                    mean2(AP_va_OST),mean2(mAP_va_OST));
            else
                fprintf('\nostRSL:\t');
                for ifo = 1:NF
                    fprintf('%d\t',ifo);

                    tic;                
                    for ita = 1:NT                
                        OST_res(ifo).task(ita).M = ostRSL(data(ifo).task(ita).x, data(ifo).task(ita).y, tri_task(ifo).task(ita).trip, opt);
                    end;
                    tim_OST(itr,ifo) = toc; 

                    for ita = 1:NT
                        [AP_tr, mAP_tr] = Test(OST_res(ifo).task(ita).M,data(ifo).task(ita).x, data(ifo).task(ita).y, opt);                    
                        [AP_va, mAP_va] = Test(OST_res(ifo).task(ita).M,data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);
                        AP_tr_OST(itr,ita,:,ifo) = AP_tr; mAP_tr_OST(itr,ita,:,ifo) = mAP_tr; 
                        AP_va_OST(itr,ita,:,ifo) = AP_va; mAP_va_OST(itr,ita,:,ifo) = mAP_va;  
                    end;            
                end;
                fprintf('\n');
                fprintf('ostRSL:AP_tr:%.3f,mAP_tr:%.3f,AP_va:%.3f,mAP_va:%.3f\n', mean2(AP_tr_OST),mean2(mAP_tr_OST), mean2(AP_va_OST),mean2(mAP_va_OST));               
            end;    

            if opt.find_omtRSL_CP_C && itr == 1
                fprintf('\nFind best C for omtRSL_cov_pa.\n');
                best_C = opt.omtRSL_CP_C;
                best_AP_va_tmp = 0;
                all_AP_va_omtRSL_CP_C = [];
                opt.rel = 'cov';
                for idx_C=1:length(opt.C_candidate)
                    opt.omtRSL_CP_C = opt.C_candidate(idx_C);
                    AP_va_tmp = 0;
                    for ifo =1:NF                    
                        [MOMTFP_tmp,que] = omtRSL(data(ifo).task,tri_fold(ifo).trip, opt);                            
                        for ita = 1:NT
                            [AP_va, mAP_va] = Test(MOMTFP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);
                            AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));                            
                        end;                            
                    end;

                    AP_va_tmp = AP_va_tmp/(NF*NT);
                    all_AP_va_omtRSL_CP_C(end+1) = AP_va_tmp;

                    if AP_va_tmp > best_AP_va_tmp
                        best_C = opt.omtRSL_CP_C;
                        best_AP_va_tmp = AP_va_tmp;
                    end;
                    fprintf('Current_C:%f, \t Current_AP_va:%f\n',opt.omtRSL_CP_C,AP_va_tmp);
                end;
                fprintf('Best_C:%f, \t Best_AP_va:%f\n',best_C,best_AP_va_tmp);
                save(sprintf('%s/%s_OMTCP_res_par.mat',results_dir,full_name),'all_AP_va_omtRSL_CP_C','best_C');
                opt.omtRSL_CP_C = best_C;
            end;

            if opt.find_omtRSL_FP_bC == 2 && itr == 1    
                opt.rel = 'fixed';
                fprintf('\nFind best b and C on AP_va\n');
                best_b = opt.omtRSL_FP_b;
                best_C = opt.omtRSL_FP_C;
                
                best_AP_va_tmp = 0;
                all_AP_va_OMTFP = zeros(length(opt.b_candidate),length(opt.C_candidate));            
                
                for idx_b = 1:length(opt.b_candidate)
                    for idx_C = 1:length(opt.C_candidate)
                        opt.omtRSL_FP_b = opt.b_candidate(idx_b);
                        opt.omtRSL_FP_C = opt.C_candidate(idx_C);
                        
                        AP_va_tmp = 0; 
                        for ifo =1:NF                    
                            [MOMTFP_tmp,que] = omtRSL(data(ifo).task, tri_fold(ifo).trip, opt);                            
                            for ita = 1:NT
                                [AP_va, mAP_va] = Test(MOMTFP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);
                                AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));                            
                            end;                            
                        end;
                        AP_va_tmp = AP_va_tmp/(NF*NT);         all_AP_va_OMTFP(idx_b,idx_C) = AP_va_tmp;

                        if AP_va_tmp > best_AP_va_tmp
                            best_AP_va_tmp = AP_va_tmp; best_b = opt.omtRSL_FP_b; best_C = opt.omtRSL_FP_C;
                        end;
                        fprintf('Current_b:%f,\tCurrent_C:%f,\tCurrent_APTe:%f\n',opt.omtRSL_FP_b, opt.omtRSL_FP_C, AP_va_tmp);
                    end;
                end;            
                opt.omtRSL_FP_b = best_b; opt.omtRSL_FP_C = best_C;
                fprintf('Best_b:%f,\tBest_C:%f,\tBest_AP_va:%f\n',best_b,best_C,best_AP_va_tmp);            
                save(sprintf('%s/%s_OMTFP_res_par',results_dir,full_name),'all_AP_va_OMTFP','best_C','best_b');            
                fprintf('omtRSL_FP_C:%f\t,omtRSL_FP_b:%f\n',opt.omtRSL_FP_C,opt.omtRSL_FP_b);
            elseif opt.find_omtRSL_FP_bC == 1 && itr == 1
                opt.rel = 'fixed';
                fprintf('\nFind best b for omtRSL alone and set omtRSL_FP_C as ostRSL_C\n');
                opt.omtRSL_FP_C = opt.ostRSL_C;

                best_b = opt.omtRSL_FP_b;
                best_AP_va_tmp = 0;
                all_AP_va_OMTFP = zeros(length(opt.b_candidate),1);            
                
                for idx_b = 1:length(opt.b_candidate)
                    opt.omtRSL_FP_b = opt.b_candidate(idx_b);
                    AP_va_tmp = 0; 
                   
                    for ifo =1:NF                    
                        [MOMTFP_tmp,que] = omtRSL(data(ifo).task, tri_fold(ifo).trip, opt);                            
                        for ita = 1:NT
                            [AP_va, mAP_va] = Test(MOMTFP_tmp(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);
                            AP_va_tmp = AP_va_tmp + mean(AP_va(val_topk));                            
                        end;                            
                    end;

                    AP_va_tmp = AP_va_tmp/(NF*NT); all_AP_va_OMTFP(idx_b) = AP_va_tmp;

                    if AP_va_tmp > best_AP_va_tmp
                        best_AP_va_tmp = AP_va_tmp; best_b = opt.omtRSL_FP_b;
                    end;
                    fprintf('Current_b:%f,\tCurrent_APTe:%f\n',opt.omtRSL_FP_b, AP_va_tmp);
                end; 
                 opt.omtRSL_FP_b = best_b;           
                fprintf('Best_b:%f,\tBest_AP_va:%f\n',best_b, best_AP_va_tmp);            
                save(sprintf('%s/%s_OMTFP_res_par',results_dir,full_name),'all_AP_va_OMTFP','best_b');            
            else
                fprintf('Using Existing omtRSL_FP_C:%f and omtRSL_FP_b:%f\n',opt.omtRSL_FP_C,opt.omtRSL_FP_b);
            end;

            %% OMTCP       
            if opt.OMTCP.use_existing
                load(sprintf('%s/%s_OMTCP_results.mat',results_dir,full_name));       
                fprintf('OMTCP:AP_tr:%.3f,mAP_tr:%.3f,AP_va:%.3f,mAP_va:%.3f\n',mean2(AP_tr_OMTCP),mean2(mAP_tr_OMTCP),mean2(AP_va_OMTCP),mean2(mAP_va_OMTCP));
            else
                opt.rel = 'cov';
                fprintf('OMTCP:\t');
                for ifo = 1:NF
                    fprintf('%d\t',ifo);
                    tic;
                    [OMTCP_results(ifo).task,que_useless] = omtRSL(data(ifo).task, tri_fold(ifo).trip, opt);
                    tim_OMTFP(itr,ifo) = toc; 

                    for ita = 1:NT
                        [AP_tr, mAP_tr] = Test(OMTCP_results(ifo).task(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, opt);   
                        [AP_va, mAP_va] = Test(OMTCP_results(ifo).task(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);
                        AP_tr_OMTCP(itr,ita,:,ifo) = AP_tr; mAP_tr_OMTCP(itr,ita,:,ifo) = mAP_tr; 
                        AP_va_OMTCP(itr,ita,:,ifo) = AP_va; mAP_va_OMTCP(itr,ita,:,ifo) = mAP_va;  
                    end;
                end;  
                fprintf('\n');
                fprintf('OMTCP:AP_tr:%.3f,mAP_tr:%.3f,AP_va:%.3f,mAP_va:%.3f\n',mean2(AP_tr_OMTCP),mean2(mAP_tr_OMTCP), mean2(AP_va_OMTCP),mean2(mAP_va_OMTCP));                
            end;    

            %% OMTFP         
            if opt.OMTFP.use_existing
                load(sprintf('%s/%s_OMTFP_results.mat',results_dir,full_name));
                fprintf('omtRSL:AP_tr:%.3f,mAP_tr:%.3f,AP_va:%.3f,mAP_va:%.3f\n',mean2(AP_tr_OMTFP),mean2(mAP_tr_OMTFP),mean2(AP_va_OMTFP),mean2(mAP_va_OMTFP));
            else
                opt.rel = 'fixed';
                fprintf('OMTFP:\t');
                for ifo = 1:NF
                    fprintf('%d\t',ifo);
                    tic;
                    [OMT_res(ifo).task,que_useless] = omtRSL(data(ifo).task, tri_fold(ifo).trip, opt);
                    tim_OMTFP(ifo) = toc; 

                    for ita = 1:NT
                        [AP_tr, mAP_tr] = Test(OMT_res(ifo).task(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, opt);   
                        [AP_va, mAP_va] = Test(OMT_res(ifo).task(ita).M, data(ifo).task(ita).xv, data(ifo).task(ita).yv, opt);
                        AP_tr_OMTFP(itr,ita,:,ifo) = AP_tr;    mAP_tr_OMTFP(itr,ita,:,ifo) = mAP_tr; 
                        AP_va_OMTFP(itr,ita,:,ifo) = AP_va;    mAP_va_OMTFP(itr,ita,:,ifo) = mAP_va;  
                    end;
                end;        
                fprintf('\n');
                fprintf('omtRSL:AP_tr:%.3f,mAP_tr:%.3f,AP_va:%.3f,mAP_va:%.3f\n', mean2(AP_tr_OMTFP),mean2(mAP_tr_OMTFP), mean2(AP_va_OMTFP),mean2(mAP_va_OMTFP));                
            end;
        
        end;
        
        save(sprintf('%s/%s_OST_results',results_dir,full_name), 'tim_OST','AP_tr_OST','mAP_tr_OST','AP_va_OST','mAP_va_OST','OST_res');
        save(sprintf('%s/%s_OMTCP_results',results_dir,full_name),'tim_OMTCP','AP_tr_OMTCP','mAP_tr_OMTCP','AP_va_OMTCP', 'mAP_va_OMTCP','OMTCP_results');
        save(sprintf('%s/%s_OMTFP_results',results_dir,full_name), 'tim_OMTFP','AP_tr_OMTFP', 'mAP_tr_OMTFP','AP_va_OMTFP', 'mAP_va_OMTFP','OMT_res');
        
        line_width=1;    font_size = 10;    font_size_legend = 10;

    %% plot
        for k=1:10
            h = figure('visible','off');

            plot(num_tri_list,repmat(mean2(mAP_tr_OST_eye(:,k,:)),[1,ntr]),'g-s');
            hold on  
            plot(num_tri_list,repmat(mean2(mAP_va_OST_eye(:,k,:)),[1,ntr]),'g-o');
            
            plot(num_tri_list,mean(mean(mAP_tr_OST(:,:,k,:),4),2),'b-X');
            plot(num_tri_list,mean(mean(mAP_va_OST(:,:,k,:),4),2),'b-<');
            
            plot(num_tri_list,mean(mean(mAP_tr_OMTCP(:,:,k,:),4),2),'m-s');
            plot(num_tri_list,mean(mean(mAP_va_OMTCP(:,:,k,:),4),2),'m-o');
                        
            plot(num_tri_list,mean(mean(mAP_tr_OMTFP(:,:,k,:),4),2),'r-s');
            plot(num_tri_list,mean(mean(mAP_va_OMTFP(:,:,k,:),4),2),'r-o');
            
            lngd = legend({'OST-eye-tr','OST-eye-va','OST-tr','OST-va','OMTCP-tr','OMTCP-va','OMTFP-tr', 'OMTFP-va'});
            set(lngd, 'Location', 'best');
            set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
            xlabel('Varied Query Ratio');
            ylabel(strcat('AP@Top',int2str(k)))
            set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
            set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
            grid
            print(h,'-depsc',strcat(results_dir,'/',data_name,'_AP_tr_va','_top_',int2str(k),'.eps'));        
            close(h);
        end;
        
        save(sprintf('%s/%s_num_tri_list.mat',results_dir,full_name),'num_tri_list');
        
        h = figure('visible','off');
%         plot(num_tri_list,repmat(mean2(mAP_tr_MTLMNN(:,:,:)),[1,ntr]),'y-s');
%         hold on
%         plot(num_tri_list,repmat(mean2(mAP_va_MTLMNN(:,:,:)),[1,ntr]),'y-o');

        plot(num_tri_list,repmat(mean2(mAP_tr_OST_eye(:,:,:)),[1,ntr]),'g-s');
        hold on
        plot(num_tri_list,repmat(mean2(mAP_va_OST_eye(:,:,:)),[1,ntr]),'g-o');

        plot(num_tri_list,mean(mean(mean(mAP_tr_OST,4),3),2),'b-s');
        plot(num_tri_list,mean(mean(mean(mAP_va_OST,4),3),2),'b-o');

        plot(num_tri_list,mean(mean(mean(mAP_tr_OMTCP,4),3),2),'m-s');
        plot(num_tri_list,mean(mean(mean(mAP_va_OMTCP,4),3),2),'m-o');

        plot(num_tri_list,mean(mean(mean(mAP_tr_OMTFP,4),3),2),'r-s');
        plot(num_tri_list,mean(mean(mean(mAP_va_OMTFP,4),3),2),'r-o');

        lngd = legend({'OST-eye-tr','OST-eye-va','OST-tr','OST-va','OMTCP-tr','OMTCP-va','OMTFP-tr', 'OMTFP-va'});
%         'MTLMNN-tr','MTLMNN-va',
        set(lngd, 'Location', 'best');
        set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        xlabel('Varied Query Ratio');
        ylabel(strcat('AP@Top',int2str(k)))
        set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        grid
        print(h,'-depsc',strcat(results_dir,'/',data_name,'_AP_tr_va','mean_top','.eps'));        
        close(h);
    end
end; % end of datasets loop