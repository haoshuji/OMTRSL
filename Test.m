function [AP, mAP] = Test(M, X, Y,options)
    AP = [];
    mAP = [];
	topk_list = options.topk;
    [d,n] = size(X);
    M = reshape(M,[d,d]); 
    X_sparce = sparse(X);
    tic;
    sim = X_sparce'*(M*X_sparce);    
%     fprintf('%f time cost\n',toc);
%     tic;
%     sim = zeros(n,n);
%     for i=1:n
%         for j=i:n
%             sim(i,j) = X_sparce(:,i)'*M*X_sparce(:,j);
%             sim(j,i) = sim(i,j);
%         end;
%     end;
%     fprintf('%f time cost\n',toc);
%     sim = X'*(M*X);
    tic;
    sim = full(sim);
    [sorted_sim,sorted_ind] = sort(sim,2,'descend'); %sort each row in descend, ,2,'descend'
%     rem_div = round(n/10);    
% %     tic;
%     sorted_ind = zeros(n,10);
%     max_values = maxk(sim,10,2,'sorting',false);
%     for i=1:n
%         uni_top_k_element = sort(unique(maxk(sim(i,:),10, 'sorting', false)),'descend');        
%         for j=1:length(uni_top_k_element)
%             find_ind = find(sim(i,:)==uni_top_k_element(j));
%             sorted_ind(i,end+1:end+length(find_ind)) = find_ind;
%         end;
%         if ~rem(i,rem_div)
%             fprintf('.');
%         end;
%     end;
%     fprintf('%f time cost\n',toc);
    
    for itopk=1:length(topk_list)
        k = topk_list(itopk);
		sum_AP = 0;
		sum_mAP = 0;
		for i=1:length(Y)
			num_each_query = 0;
			avg_pre_each = 0;
			y_temp = Y(sorted_ind(i,1:k)) == Y(i);
			for j=1:k
				if y_temp(j) == 1
					num_each_query = num_each_query + 1;	
					avg_pre_each = avg_pre_each + num_each_query/j;
				end;
			end;
			sum_AP = sum_AP + num_each_query/k;
			sum_mAP = sum_mAP + avg_pre_each/k;
		end;
		AP(itopk) = sum_AP / length(Y);
		mAP(itopk) = sum_mAP / length(Y);
	end;
