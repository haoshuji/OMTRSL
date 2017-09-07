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

options.triplets_using_existed = true; 
options.tri.hard= false;  
options.tri_time = 10;

options.mtLMNN_using_existed = false; 
options.mtLMNN_reg_search = false;
options.maxitaer_mt = 100;
options.Kg = 10; 
options.mtLMNN.best_reg0 = 1; 
options.mtLMNN.best_reg = 500;
options.Kg_candidate = [1,2,3,4,6,8,10];

data_names = { 'isolet'};%'isolet' 'letter'
norm_methods= {'length1'}; 


%% all datasets
for idx_data=1:length(data_names)
    data_name = data_names{idx_data};
    for id = 1:length(norm_methods)
        full_name = strcat(data_name,'_',norm_methods{id});
        fprintf('Experiments on %s dataset\n',full_name);   
%         current_date = datestr(now,'mmmm-dd-yyyy-HH');
        current_date =  date;%'25-Jan-2016';

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
            for iKg=1:length(options.Kg_candidate)
                time_mtLMNN = zeros(NF,1);  
                AP_tr_mtLMNN = zeros(NT,NTK,NF); mAP_tr_mtLMNN = zeros(NT,NTK,NF);
                AP_te_mtLMNN = zeros(NT,NTK,NF); mAP_te_mtLMNN = zeros(NT,NTK,NF);
                options.Kg = options.Kg_candidate(iKg);
                fprintf('%d\n',options.Kg);
                for ifo = 1:NF
                    tic
                    mtD=MTlmnn_sp(data(ifo).task,options.Kg,'quiet',1,'regweight',options.mtLMNN.best_reg, 'regweight0',options.mtLMNN.best_reg0,'weight1',0.5,'maxitaer',options.maxitaer_mt); 
                    time_mtLMNN(ifo) = toc;      
                    for ita=1:NT
                        mtLMNN_task(ita).M= mtD.bestM0{ita}+mtD.bestM{ita};                            
                    end;                
                    for ita=1:NT
                        [AP_tr,mAP_tr] = Test(mtLMNN_task(ita).M, data(ifo).task(ita).x, data(ifo).task(ita).y, options); 
                        [AP_te,mAP_te] = Test(mtLMNN_task(ita).M, data(ifo).task(ita).xt, data(ifo).task(ita).yt, options);                     
                        AP_tr_mtLMNN(ita,:,ifo) = AP_tr;  mAP_tr_mtLMNN(ita,:,ifo) = mAP_tr; AP_te_mtLMNN(ita,:,ifo) = AP_te;   mAP_te_mtLMNN(ita,:,ifo) = mAP_te;
                    end;                      
                end;
%                 fprintf('\n');
                mean(AP_te_mtLMNN(ita,:,:),3)
                fprintf('mtLMNN:AP_tr:%.3f,mAP_tr:%.3f,AP_te:%.3f,mAP_te:%.3f\n', mean2(AP_tr_mtLMNN),mean2(mAP_tr_mtLMNN), mean2(AP_te_mtLMNN),mean2(mAP_te_mtLMNN));
                save(sprintf('%s/%s_mtLMNN_results_%d.mat',results_dir,full_name,options.Kg),...
                    'time_mtLMNN','AP_tr_mtLMNN','mAP_tr_mtLMNN','AP_te_mtLMNN', 'mAP_te_mtLMNN');
            end;
        end; 
    end;    
end; % end of datasets loop