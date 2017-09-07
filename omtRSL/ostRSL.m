function M = ostRSL(X, Y, triplets, options)
%online single task Relative Similarity Learning
if triplets(1,5) == -1
    C = options.ostRSL_global_C;
else
    C = options.ostRSL_C;
end;

[d,n] = size(X);
M = eye(d);
num_triplets = size(triplets,1);
% fprintf('%d triplets\n',num_triplets);
rem_div = round(num_triplets/10);
num_update = 0;
for ind_tri=1:num_triplets
    trip = triplets(ind_tri,:);
    trip_cell = num2cell(trip);
	[x_ind, x_1_ind, x_2_ind, y, it_useless] = trip_cell{:};
	X_t = X(:,x_ind)*(X(:,x_1_ind) - X(:,x_2_ind))';
	
%     X_t = X_t/norm(X_t,'fro');
    
    p = trace(M*X_t');

	ell = max(0, 1-y*p);
	if ell > 0.0
        num_update = num_update + 1;
		tau = min(C, ell/norm(X_t,'fro')^2);
		M = M + tau*y*X_t;
	end;
%     if ~rem(ind_tri,rem_div)
%         fprintf('.');
%     end;
end;
% fprintf('num_update/num_tri = %f\n',num_update/num_triplets);