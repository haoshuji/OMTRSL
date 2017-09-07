function plot_results4(current_date,data_name)
    % without consider OST_global algorithm
    results_dir = './results/';
    %     data_names = {'landmine_scale01'};
    options.topk = [1:1:10];
    options.topk_varied = [1,5,10];
    ntopk = length(options.topk);


    line_width=1;    font_size = 20;
    font_size_legend = 25;

    results_dir = strcat(results_dir,data_name,'/',current_date,'/');
    all_alg = {'MTLMNN','OST','OSTE','OMTAR','OMTFP','OMTCP'};
    for j=1:length(all_alg)
        load(sprintf('%s%s_%s_res.mat',results_dir,data_name,all_alg{j}));
    end;
    num_que = size(queOMTFA,1);
    
    b = options.topk;
    
%% AP_tr
    % h=figure('visible','off');   
    % plot(options.topk,mean(mean(APTrMTLMNN(:,:,:),3),1),'g-s');
    % hold on  
    % plot(options.topk,mean(mean(APTrOST(:,:,:),3),1),'b->');
    % plot(options.topk,mean(mean(APTrOMTFP(:,:,:),3),1),'r-s');    
    % lngd = legend('MTLMNN','OST','omtRSL', 'location','best');
    % set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
    % xlabel('Topk');
    % ylabel('Average Precision')
    % set(findall(gca,'-property','FontSize'),'FontSize',font_size);
    % set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
    % grid
    % print(h,'-depsc',strcat(results_dir,'/',data_name,'_AP_tr','.eps'));   
    % close(h);
    
    % his_y = [mean(mean(APTrMTLMNN(:,:,:),3),1)',mean(mean(APTrOST(:,:,:),3),1)',...
    %    mean(mean(APTrOMTFP(:,:,:),3),1)'];
    % h = figure('visible','off'); 
    % h2 = bar(b,his_y,'grouped');
    % xlabel('Topk');
    % ylabel('Average Precision');    
    % l = cell(1,3);
    % l{1}='MTLMNN'; l{2}='OST'; l{3}='omtRSL';
    % lngd = legend(h2,l);
    % set(lngd, 'interpreter','latex','fontsize',font_size_legend,'location','best'); 
    % set(findall(gca,'-property','FontSize'),'FontSize',font_size);
    % print(h,'-depsc',strcat(results_dir,'/',data_name,'_APTr','.eps'));
    % close(h);

% %% mAP_tr
%     h=figure('visible','off');   
%     plot(options.topk,mean(mean(mAPTrMTLMNN(:,:,:),1),3),'g-s');
%     hold on  
%     plot(options.topk,mean(mean(mAPTrOST(:,:,:),1),3),'b->');
%     plot(options.topk,mean(mean(mAPTrOMTFP(:,:,:),1),3),'r-s');    
%     lngd = legend({'MTLMNN','OST','omtRSL'});
%     set(lngd, 'Location', 'best');
%     set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
%     xlabel('Topk');
%     ylabel('Mean Average Precision')
%     set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
%     set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
%     grid
%     print(h,'-depsc',strcat(results_dir,'/',data_name,'_mAP_tr','.eps'));
%     close(h);
    
    % his_y = [mean(mean(mAPTrMTLMNN(:,:,:),3),1)',mean(mean(mAPTrOST(:,:,:),3),1)', mean(mean(mAPTrOMTFP(:,:,:),3),1)'];
    % h = figure('visible','off');  
    % h2 = bar(b,his_y,'grouped');
    % xlabel('Topk');
    % ylabel('Average Precision');    
    % l = cell(1,3);
    % l{1}='MTLMNN'; l{2}='OST'; l{3}='omtRSL';
    % lngd = legend(h2,l);
    % set(lngd, 'interpreter','latex','fontsize',font_size_legend,'location','best'); 
    % set(findall(gca,'-property','FontSize'),'FontSize',font_size);
    % print(h,'-depsc',strcat(results_dir,'/',data_name,'_mAPTrhis','.eps'));
    % close(h);
    
%% AP_te    
    % h=figure('visible','off');   
    % plot(options.topk,mean(mean(APTeMTLMNN(:,:,:),3),1),'g-s');
    % hold on  
    % plot(options.topk,mean(mean(APTeOST(:,:,:),3),1),'b->');
    % plot(options.topk,mean(mean(APTeOMTFP(:,:,:),3),1),'r-s');    
    % lngd = legend({'MTLMNN','OST','omtRSL'});
    % set(lngd, 'Location', 'best');
    % set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
    % xlabel('Topk');
    % ylabel('Average Precision')
    % set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
    % set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
    % grid
    % print(h,'-depsc',strcat(results_dir,'/',data_name,'_AP_te','.eps'));
    % close(h);
    
    his_y = [mean(mean(APTeMTLMNN(:,:,:),3),1)',mean(mean(APTeOST_eye(:,:,:),3),1)',...
        mean(mean(APTeOST(:,:,:),3),1)',...
        mean(mean(APTeOMTFP(:,:,:),3),1)'...
        mean(mean(APTeomtRSL_cov_pa(:,:,:),3),1)'];
    h = figure('visible','off');   
    h2 = bar(b,his_y,'grouped');
    xlabel('Topk');
    ylabel('Average Precision');    
    l = cell(1,5);
    l{1}='MTLMNN';l{2}='OST-eye'; l{3}='OST'; l{4}='omtRSL-Fix',l{5}='omtRSL-Cov';
    lngd = legend(h2,l);
    set(lngd, 'interpreter','latex','fontsize',font_size_legend,'location','best'); 
    set(findall(gca,'-property','FontSize'),'FontSize',font_size);
    print(h,'-depsc',strcat(results_dir,'/',data_name,'_APTehis','.eps'));
    saveas(h,strcat(results_dir,'/',data_name,'_APTe','.fig'));
    close(h);
    
%% mAP_te    
    % h=figure('visible','off');   
    % plot(options.topk,mean(mean(mAPTeMTLMNN(:,:,:),3),1),'g-s');
    % hold on  
    % plot(options.topk,mean(mean(mAPTeOST(:,:,:),3),1),'b->');
    % plot(options.topk,mean(mean(mAPTeOMTFP(:,:,:),3),1),'r-s');    
    % lngd = legend({'MTLMNN','OST','omtRSL'});
    % set(lngd, 'Location', 'best');
    
    % set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
    % xlabel('Topk');
    % ylabel('Mean Average Precision')
    % set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
    % set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
    % grid
    % print(h,'-depsc',strcat(results_dir,'/',data_name,'_mAP_te','.eps'));
    % % print(h,'-dpdf',strcat(results_dir,'/',data_name,'_mAP_te','.pdf'));
    % saveas(h,strcat(results_dir,'/',data_name,'_mAP_te','.fig'));
    % close(h);
    
    his_y = [mean(mean(mAPTeMTLMNN(:,:,:),3),1)',mean(mean(mAPTeOST(:,:,:),3),1)', mean(mean(mAPTeOMTFP(:,:,:),3),1)'];
    h = figure('visible','off');   
    h2 = bar(b,his_y,'grouped');
    xlabel('Topk');
    ylabel('Average Precision');    
    l = cell(1,3);
    l{1}='MTLMNN'; l{2}='OST'; l{3}='omtRSL';
    lngd = legend(h2,l);
    set(lngd, 'interpreter','latex','fontsize',font_size_legend,'location','best'); 
    set(findall(gca,'-property','FontSize'),'FontSize',font_size);
    print(h,'-depsc',strcat(results_dir,'/',data_name,'_mAPTehis','.eps'));
    saveas(h,strcat(results_dir,'/',data_name,'_mAPTehis','.fig'));
    close(h);
    
    
    que_space = linspace(0,1,num_que);
    % draw time
    h=figure('visible','off');   
    loglog(que_space,repmat(mean(timMTLMNN),[1,num_que]),'g-s');
    hold on  
    loglog(que_space,repmat(mean(timOST),[1,num_que]),'b->');
    loglog(que_space,repmat(mean(timOMTFP),[1,num_que]),'c-<');
    loglog(mean(queOMTFR,2),mean(timOMTFR,2),'m-d');
    loglog(mean(queOMTFA,2), mean(timOMTFA,2),'r-o');
    lngd = legend(legend_st);%'OOSL', 'OROSL', 'OAOSL', 
    set(lngd, 'Location', 'best');
    set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
    xlabel('Varied Query Ratio');
    ylabel('Time(s)')
    set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
    set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
    % ylim([7 100])
    grid
    print(h,'-depsc',strcat(results_dir,'/',data_name,'_time','.eps'));
    % print(h,'-dpdf',strcat(results_dir,'/',data_name,'_time','.pdf'));
    saveas(h,strcat(results_dir,'/',data_name,'_time','.fig'));

    close(h);

    legend_str = {'MTLMNN','ostRSL','omtRSL-Passive','omtRSL-Random','omtRSL-Active'};
    legend_str2 = {'MTLMNN','ostRSL','OMTFP','OMTFR','OMTFA','OMTCR','OMTCA'};
%% Active and Random
    for ik = 1:length(options.topk_varied)

        k = options.topk_varied(ik);

        % plot the results
        que_space = linspace(0,1,num_que);

        h=figure('visible','off');   
        plot(que_space,repmat(mean(mean(APTrMTLMNN(:,k,:))),[1,num_que]),'g-s');
        hold on  
        plot(que_space,repmat(mean(mean(APTrOST(:,k,:))),[1,num_que]),'b->');
        plot(que_space,repmat(mean(mean(APTrOMTFP(:,k,:))),[1,num_que]),'c-<');
        plot(mean(queOMTFR,2),mean(mean(APTrOMTFR(:,:,k,:),4),2),'m-d');
        plot(mean(queOMTFA,2),mean(mean(APTrOMTFA(:,:,k,:),4),2),'r-o');
        lngd = legend(legend_str);
        set(lngd, 'Location', 'best');
        set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        xlabel('Varied Query Ratio');
        ylabel(strcat('AP@Top',int2str(k)))
        set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        grid
        print(h,'-depsc',strcat(results_dir,'/',data_name,'_AP_tr','_top_',int2str(k),'.eps'));
        close(h);

        h=figure('visible','off');   
        plot(que_space,repmat(mean(mean(mAPTrMTLMNN(:,k,:))),[1,num_que]),'g-s');
        hold on  
        plot(que_space,repmat(mean(mean(mAPTrOST(:,k,:))),[1,num_que]),'b->');
        plot(que_space,repmat(mean(mean(mAPTrOMTFP(:,k,:))),[1,num_que]),'c-<');
        plot(mean(queOMTFR,2),mean(mean(mAPTrOMTFR(:,:,k,:),4),2),'m-d');
        plot(mean(queOMTFA,2),mean(mean(mAPTrOMTFA(:,:,k,:),4),2),'r-o');
        lngd = legend(legend_str);
        set(lngd, 'Location', 'best');
        set(lngd, 'interpreter','latex','fontsize',font_size_legend); 
        xlabel('Varied Query Ratio');
        ylabel(strcat('mAP@Top',int2str(k)))
        set(findall(gca,'-property','FontSize'),'FontSize',font_size); 
        set(findall(gca, '-property', 'linewidth'),'linewidth',line_width);
        grid
        print(h,'-depsc',strcat(results_dir,'/',data_name,'_mAP_tr','_top_',int2str(k),'.eps'));    
        close(h);

        h=figure('visible','off');   
        plot(que_space,repmat(mean(mean(APTeMTLMNN(:,k,:))),[1,num_que]),'g-s');
        hold on  
        plot(que_space,repmat(mean(mean(APTeOST(:,k,:))),[1,num_que]),'b->');
        plot(que_space,repmat(mean(mean(APTeOMTFP(:,k,:))),[1,num_que]),'c-<');
        % plot(mean(queOMTFR,2),mean(mean(APTeOMTFR(:,:,k,:),4),2),'m-x');
        % plot(mean(queOMTFA,2),mean(mean(APTeOMTFA(:,:,k,:),4),2),'m-o');
        plot(mean(queOMTCR,2),mean(mean(APTeOMTCR(:,:,k,:),4),2),'m-o');
        plot(mean(queOMTCA,2),mean(mean(APTeOMTCA(:,:,k,:),4),2),'r-o');
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
        

        h=figure('visible','off');   
        plot(que_space,repmat(mean(mean(mAPTeMTLMNN(:,k,:))),[1,num_que]),'g-s');
        hold on  
        plot(que_space,repmat(mean(mean(mAPTeOST(:,k,:))),[1,num_que]),'b->');
        plot(que_space,repmat(mean(mean(mAPTeOMTFP(:,k,:))),[1,num_que]),'c-<');
        % plot(mean(queOMTFR,2),mean(mean(mAPTeOMTFR(:,:,k,:),4),2),'m-x');
        % plot(mean(queOMTFA,2), mean(mean(mAPTeOMTFA(:,:,k,:),4),2),'m-o');
        plot(mean(queOMTCR,2),mean(mean(mAPTeOMTCR(:,:,k,:),4),2),'m-d');
        plot(mean(queOMTCA,2),mean(mean(mAPTeOMTCA(:,:,k,:),4),2),'r-o');
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
