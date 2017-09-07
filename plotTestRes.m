function plotTestRes(current_date,data_name)
    % without consider ostRSL_global algorithm
    results_dir = './results/';
    %     data_names = {'landmine_scale01'};
    options.topk = [1:1:10];
    options.topk_varied = [5,10];
    ntopk = length(options.topk);
    options.b_candidate = [1,2,4,8,16,32,64,128];
    options.C_candidate = [0.001,0.01,0.1,1,10,100,1000];
    bla= [0,0,0];

    line_width=1;    font_size = 15;
    font_size_legend = 15;

    results_dir = strcat(results_dir,data_name,'/',current_date,'/');
    all_alg = {'mtLMNN','ostRSL','ostRSL_global','omtRSL_AR','omtRSL_FP'};
    for j=1:length(all_alg)
        load(sprintf('%s%s_%s_results.mat',results_dir,data_name,all_alg{j}));
    end;
    
    num_que = size(que_omtRSL_FA,1);
    
    b = options.topk;

    legend_names = {'mtLMNN'  'global'  'ostRSL' 'omtRSL'};
    legend_str = {'mtLMNN','ostRSL', 'omtRSL','omtRSL-Random','omtRSL-Active'};
  
    %% AP and mAP
    AP_y = [mean(mean(AP_te_mtLMNN(:,:,:),3),1)',mean(AP_te_ostRSL_global(:,:),2),...
        mean(mean(AP_te_ostRSL(:,:,:),3),1)', mean(mean(AP_te_omtRSL_FP(:,:,:),3),1)'];
    h = figure('visible','off');   
    min_y = min(mean(AP_te_ostRSL_global(:,:),2));
    h2 = bar(b,AP_y);
%     fH = gcf;
%     set(h2(1),'FillPattern','DiagonalLines');
%     set(h2(1),'facecolor',bla+0.9);
%     set(h2(2),'facecolor',bla+0.6);
%     set(h2(3),'facecolor','b');
%     set(h2(4),'facecolor','r');
%     applyhatch_pluscolor(fH, '\-x.', 0, [1 0 1 0], jet(4));
%     hatchfill(h2);
    ylim manual
    ylim([0.6 1]);
    xlabel('Topk');
    ylabel('Average Precision');    
    lngd = legend(h2,legend_names);
    set(lngd, 'interpreter','latex','fontsize',font_size_legend,'location','best'); 
    set(findall(gca,'-property','FontSize'),'FontSize',font_size);
    file_name = strcat(results_dir,'/',data_name,'_AP_te_bar.eps');    
    print(h,'-depsc',file_name);
    close(h);
    
    mAP_y = [mean(mean(mAP_te_mtLMNN(:,:,:),3),1)',mean(mAP_te_ostRSL_global(:,:),2),...
        mean(mean(mAP_te_ostRSL(:,:,:),3),1)', mean(mean(mAP_te_omtRSL_FP(:,:,:),3),1)'];
    h = figure('visible','off');   
    h2 = bar(b,mAP_y,'grouped');
%     set(h2(1),'facecolor','b');
%     set(h2(2),'facecolor','g');
%     set(h2(3),'facecolor','y');
%     set(h2(4),'facecolor','r');

    ylim manual
    ylim([0.6 1]);
    xlabel('Topk');
    ylabel('mean Average Precision');    
    lngd = legend(h2,legend_names);
    set(lngd, 'interpreter','latex','fontsize',font_size_legend,'location','best'); 
    set(findall(gca,'-property','FontSize'),'FontSize',font_size);
    file_name = strcat(results_dir,'/',data_name,'_mAP_te_bar.eps');    
    print(h,'-depsc',file_name);
    close(h);
    
    %% time
    
    %%bar
    time_y = [mean(time_mtLMNN), mean(time_ostRSL_global), mean(time_ostRSL), mean(time_omtRSL_FP)];
    h = figure('visible','off');   
    colors = hsv(numel(time_y));
    h2 = bar(1,time_y(1),'FaceColor','b');
    hold on;
    bar(2,time_y(2),'FaceColor','g');
    bar(3,time_y(3),'FaceColor','y');
    bar(4,time_y(4),'FaceColor','r');
    
    set(gca,'XTickLabel',{'mtLMNN', 'Global', 'ostRSL', 'omtRSL'});
    xlabel('Different algorithms');
    ylabel('Time (s)');    
    lngd = legend;
    set(lngd, 'interpreter','latex','fontsize',font_size_legend,'location','best'); 
    set(findall(gca,'-property','FontSize'),'FontSize',font_size);
    print(h,'-depsc',strcat(results_dir,'/',data_name,'_time_bar.eps'));
    close(h);
    
    %% line
    que_space = linspace(0,1,num_que);
    % draw time
    h=figure('visible','off');   
    plot(que_space,repmat(mean(time_mtLMNN),[1,num_que]),'g-s');
    hold on  
    plot(que_space,repmat(mean(time_ostRSL),[1,num_que]),'b->');
    plot(que_space,repmat(mean(time_ostRSL_global),[1,num_que]),'b-o');
    plot(que_space,repmat(mean(time_omtRSL_FP),[1,num_que]),'c-<');
    plot(mean(que_omtRSL_FR,2),mean(time_omtRSL_FR,2),'m-d');
    plot(mean(que_omtRSL_FA,2), mean(time_omtRSL_FA,2),'r-o');
    lngd = legend(legend_str);
    set(lngd, 'Location', 'best');
    set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
    xlabel('Varied Query Ratio');
    ylabel('Time(s)')
    set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
    set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
    grid
    print(h,'-depsc',strcat(results_dir,'/',data_name,'_time_line.eps'));
    close(h);
    
    %% parameter
    if 1
        
        load(sprintf('%s%s_ostRSL_results_parameters.mat',results_dir,data_name));
        load(sprintf('%s%s_omtRSL_FP_results_parameters.mat',results_dir,data_name));
    
    
        h = figure('visible','off');
        semilogx(options.C_candidate, all_AP_va_ostRSL_C,'b-o');
        lngd = legend('ostRSL');
        set(lngd, 'Location', 'best');
        set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        xlabel('Parameter C');
        ylabel(strcat('AP@Top',int2str(10)));
        ylim([0 1]);
        set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        grid
        print(h,'-depsc',strcat(results_dir,'/',data_name,'_ostRSL_C.eps'));
        close(h);
        
        h = figure('visible','off');
        plot(options.b_candidate, all_AP_va_omtRSL_FP,'r-o');
        lngd = legend('omtRSL');
        ylim([0 1]);
        set(lngd, 'Location', 'best');
        set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        xlabel('Parameter b');
        ylabel(strcat('AP@Top',int2str(10)));
        set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        grid
        print(h,'-depsc',strcat(results_dir,'/',data_name,'_omtRSL_b.eps'));
        close(h);
    end;
    
%% Active and Random
    for ik = 1:length(options.topk_varied)

        k = options.topk_varied(ik);

        % plot the results
        que_space = linspace(0,1,num_que);

        % h=figure('visible','off');   
        % plot(que_space,repmat(mean(mean(AP_tr_mtLMNN(:,k,:))),[1,num_que]),'g-s');
        % hold on  
        % plot(que_space,repmat(mean(mean(AP_tr_ostRSL(:,k,:))),[1,num_que]),'b->');
        % plot(que_space,repmat(mean(mean(AP_tr_omtRSL_FP(:,k,:))),[1,num_que]),'c-<');
        % plot(mean(queomtRSL_FR,2),mean(mean(AP_tr_omtRSL_FR(:,:,k,:),4),2),'m-d');
        % plot(mean(que_omtRSL_FA,2),mean(mean(AP_tr_omtRSL_FA(:,:,k,:),4),2),'r-o');
        % lngd = legend(legend_str);
        % set(lngd, 'Location', 'best');
        % set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        % xlabel('Varied Query Ratio');
        % ylabel(strcat('AP@Top',int2str(k)))
        % set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        % set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        % grid
        % print(h,'-depsc',strcat(results_dir,'/',data_name,'_AP_tr','_top_',int2str(k),'.eps'));
        % close(h);

        % h=figure('visible','off');   
        % plot(que_space,repmat(mean(mean(mAP_tr_mtLMNN(:,k,:))),[1,num_que]),'g-s');
        % hold on  
        % plot(que_space,repmat(mean(mean(mAP_tr_ostRSL(:,k,:))),[1,num_que]),'b->');
        % plot(que_space,repmat(mean(mean(mAP_tr_omtRSL_FP(:,k,:))),[1,num_que]),'c-<');
        % plot(mean(queomtRSL_FR,2),mean(mean(mAP_tr_omtRSL_FR(:,:,k,:),4),2),'m-d');
        % plot(mean(que_omtRSL_FA,2),mean(mean(mAP_tr_omtRSL_FA(:,:,k,:),4),2),'r-o');
        % lngd = legend(legend_str);
        % set(lngd, 'Location', 'best');
        % set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        % xlabel('Varied Query Ratio');
        % ylabel(strcat('mAP@Top',int2str(k)))
        % set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        % set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        % grid
        % print(h,'-depsc',strcat(results_dir,'/',data_name,'_mAP_tr','_top_',int2str(k),'.eps'));    
        % close(h);
        
        que_idx = [3,4,6,8,9,10,11,12,14,20];
        h=figure('visible','off');   
        plot(que_space,repmat(mean(mean(AP_te_mtLMNN(:,k,:))),[1,num_que]),'b-+');
        hold on          
%         plot(que_space,repmat(mean(AP_te_ostRSL_global(k,:)),[1,num_que]),'g');
        plot(que_space,repmat(mean(mean(AP_te_ostRSL(:,k,:))),[1,num_que]),'g-x');
        plot(que_space,repmat(mean(mean(AP_te_omtRSL_FP(:,k,:))),[1,num_que]),'k-*');
        que_omtRSL_FR_tmp = mean(que_omtRSL_FR,2); AP_tmp = mean(mean(AP_te_omtRSL_FR(:,:,k,:),4),2);
        plot(que_omtRSL_FR_tmp(que_idx),AP_tmp(que_idx),'m-d');
        que_omtRSL_FA_tmp = mean(que_omtRSL_FA,2);AP_tmp = mean(mean(AP_te_omtRSL_FA(:,:,k,:),4),2);
        plot(que_omtRSL_FA_tmp(que_idx),AP_tmp(que_idx),'r-o');
        lngd = legend(legend_str);
        set(lngd, 'Location', 'best');
        set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        xlabel('Varied Query Ratio');
        ylabel(strcat('AP@Top',int2str(k)))
        set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        grid
        print(h,'-depsc',strcat(results_dir,'/',data_name,'_AP_te','_top_',int2str(k),'.eps'));
        close(h);
        

        que_idx = [3,4,6,8,9,10,11,12,14,20];
        h=figure('visible','off');   
        plot(que_space,repmat(mean(mean(mAP_te_mtLMNN(:,k,:))),[1,num_que]),'b-+');
        hold on          
%         plot(que_space,repmat(mean(mAP_te_ostRSL_global(k,:)),[1,num_que]),'g');
        plot(que_space,repmat(mean(mean(mAP_te_ostRSL(:,k,:))),[1,num_que]),'g-x');
        plot(que_space,repmat(mean(mean(mAP_te_omtRSL_FP(:,k,:))),[1,num_que]),'k-*');
        que_omtRSL_FR_tmp = mean(que_omtRSL_FR,2); mAP_tmp = mean(mean(mAP_te_omtRSL_FR(:,:,k,:),4),2);
        plot(que_omtRSL_FR_tmp(que_idx),mAP_tmp(que_idx),'m-d');
        que_omtRSL_FA_tmp = mean(que_omtRSL_FA,2);mAP_tmp = mean(mean(mAP_te_omtRSL_FA(:,:,k,:),4),2);
        plot(que_omtRSL_FA_tmp(que_idx),mAP_tmp(que_idx),'r-o');
        lngd = legend(legend_str);
        set(lngd, 'Location', 'best');
        set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        xlabel('Varied Query Ratio');
        ylabel(strcat('mAP@Top',int2str(k)))
        set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        grid
        print(h,'-depsc',strcat(results_dir,'/',data_name,'_mAP_te','_top_',int2str(k),'.eps'));
        close(h);
    end;
