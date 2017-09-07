dim = size(task(1).x,1);
% Regular LMNN
for i=1:length(task)

    % Validation set: final L
    task(i).errI_1=knnclassify(eye(dim),task(i).x,task(i).y,task(i).xv,task(i).yv,Kg).*100;
    task(i).errL_1=knnclassify(D(i).minL,task(i).x,task(i).y,task(i).xv,task(i).yv,Kg).*100;

    % Validation set: best L
    if isfield(D(i),'bestL')
     task(i).errL_2=knnclassify(D(i).bestL,task(i).x,task(i).y,task(i).xv,task(i).yv,Kg).*100;
    end

    % Test set: final L
    task(i).errI_1_te=knnclassify(eye(dim),[task(i).x task(i).xv],[task(i).y task(i).yv],task(i).xTe,task(i).yTe,Kg).*100;
    task(i).errL_1_te=knnclassify(D(i).minL,[task(i).x task(i).xv],[task(i).y task(i).yv],task(i).xTe,task(i).yTe,Kg).*100;

    % Test set: best L
    if isfield(D(i),'bestL')
     task(i).errL_2_te=knnclassify(D(i).bestL,[task(i).x task(i).xv],[task(i).y task(i).yv],task(i).xTe,task(i).yTe,Kg).*100;
    end;
    
end;

errs_1=[[task.errI_1];[task.errL_1]];
errs_2=[[task.errI_1];[task.errL_2]];

errs_1_te=[[task.errI_1_te];[task.errL_1_te]];
errs_2_te=[[task.errI_1_te];[task.errL_2_te]];

for i=1:length(reg0)
    for j=1:length(reg)
        % Checking if I got the right Dets
        if mtD{i,j}.pars.regweight ~= reg(j) || mtD{i,j}.pars.regweight0 ~= reg0(i)
            error('Wrong mtD.')
        end

        % run mtLMNN
        for st=1:length(task)

            % Validation set: specific M0,Mi
            task(st).errMTL_1=knnclassify( decompose( mtD{i,j}.bestM0{st} + mtD{i,j}.bestM{st}), task(st).x,task(st).y,task(st).xv,task(st).yv,Kg).*100;
            task(st).errMTL0_1=knnclassify( decompose( mtD{i,j}.M0 + mtD{i,j}.M{st}), task(st).x,task(st).y,task(st).xv,task(st).yv,Kg).*100;
            
            % Validation set: general M0 + Mi
            task(st).errMTL_2=knnclassify(decompose( mtD{i,j}.bestMgen.M0 + mtD{i,j}.bestMgen.M{st}),task(st).x,task(st).y,task(st).xv,task(st).yv,Kg).*100;
            

            % Test set: specific M0,Mi
            task(st).errMTL_1_te=knnclassify(decompose( mtD{i,j}.bestM0{st}+ mtD{i,j}.bestM{st}),[task(st).x task(st).xv],[task(st).y task(st).yv],task(st).xTe,task(st).yTe,Kg).*100;
            task(st).errMTL0_1_te=knnclassify(decompose( mtD{i,j}.M0 + mtD{i,j}.M{st}) ,[task(st).x task(st).xv],[task(st).y task(st).yv],task(st).xTe,task(st).yTe,Kg).*100;

            % Test set: general M0 + Mi
            task(st).errMTL_2_te=knnclassify(decompose( mtD{i,j}.bestMgen.M0 + mtD{i,j}.bestMgen.M{st}),[task(st).x task(st).xv],[task(st).y task(st).yv],task(st).xTe,task(st).yTe,Kg).*100;
            
        end;

        % Validation set Errors
        errs_1=[errs_1;[task.errMTL0_1];[task.errMTL_1]];
        errs_2=[errs_2;[task.errMTL0_1];[task.errMTL_2]];

        % Test set errors
        errs_1_te=[errs_1_te;[task.errMTL0_1_te];[task.errMTL_1_te]];
        errs_2_te=[errs_2_te;[task.errMTL0_1_te];[task.errMTL_2_te]];
  end;
end;
