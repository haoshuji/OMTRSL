function [Det]=mtLMNN(task,varargin)
%
% function [L,Det]=mtLMNN(x,y,Kg,'Parameter1',Value1,'Parameter2',Value2,...);
%
% Input:
%
% x = input matrix (each column is an input vector) 
% y = labels 
% (*optional*)  L = initial transformation matrix (e.g eye(size(x,1)))
% (*optional*) Kg = attract Kg nearest similar labeled vectos
% 
% Parameters:
% stepsize = (default 1e-09)
% save = (def 0) save the initial computation 
% skip = (def 0) loads the initial computation instead of
%        recomputing (only works if previous run was on exactly the same data)
% correction = (def 15) how many steps between each update 
%              The number of impostors are fixed for until next "correction"
% factor = (def 1.1) multiplicative factor by which the
%         "correction" gab increases
% obj = (def 1) if 1, solver solves in L, if 0, solver solves in L'*L
% thresho = (def 1e-9) cut off for change in objective function (if
%           improvement is less, stop)
% thresha = (def 1e-22) cut off for stepsize, if stepsize is
%           smaller stop
% validation = (def 0) fraction of training data to be used as
%               validation set (best output is stored in Det.bestL)
% validationstep = (def 50) every "valcount" steps do validation
% maxiter = maximum number of iterations (default: 10000)
% scale = (def. 0) if 1, all data gets re-scaled s.t. average
%         distance to closest neighbor is 1
% quiet = {0,1} surpress output (default=0)  
%
%
% Output:
%
% L = linear transformation xnew=L*x
%    
% Det.obj = objective function over time
% Det.nimp = number of impostors over time
% Det.pars = all parameters used in run
% Det.time = time needed for computation
% Det.iter = number of iterations
% Det.verify = verify (results of validation - if used)
%  
% Version 1.0
% copyright by Kilian Q. Weinbergerr (2005)
% University of Pennsylvania
% contact kilianw@seas.upenn.edu
%

% print help screen if no parameters are given                                                                        

global L0 M0 M0old doingcorrection bestMgen;        
L0=[];      
M0=[];            
M0old=M0;          
bestMgen.M0 = M0;

tempnum=num2str(round(rand(1).*1000));
tempname=sprintf('temp_%s.mat',tempnum);
% fprintf('Filename: %s \n',tempname);

% fprintf('LMNN stable\n');
if(nargin==0)
 help lmnn;
 return;
end;
                            
% backwards compatibility 
if ~isstruct(task)
    x=task;
    task=[];
    task(1).x=x;
    task(1).y=varargin{1};    
    varargin(1)=[];       
end;
                           
% read in neighborhood size Kg and initial matrix L0
while(~isempty(varargin) && isnumeric(varargin{1}))
 % pars neighborhood size Kg
 if prod(size(varargin{1}))==1,
    Kg=varargin{1};
    varargin(1)=[];
    continue;
 end;     

 % pars initial matrix L0
 if prod(size(varargin{1}))>1,
    L0=varargin{1};
    varargin(1)=[];
 end;    
end;
                                    
% If not specified set to default values
if(isempty(L0))
 % fprintf(['Initial starting point not specified.\nStarting with identity matrix.\n']);  
 L0=eye(size(task(1).x,1));
end;      

% compute actual Mahalanobis matrix
M0=L0'*L0;

if(exist('Kg','var')~=1)                                                
 Kg=3;
 % fprintf(['Neighborhood size not specified.\nSetting it to Kg=%i.\n'],Kg);  
end;
tic

% Do a dimensionality check
for i=1:length(task)                       
    % Initializationip
    % fprintf('Task %i has %i dimensions and %i input vectors\n',i,size(task(i).x));
    if(size(task(i).x,1)~=size(L0,2)) error('x and L must have matching dimensions!\n');end;
end;



% set parameters
pars.stepsize=1e-07;
pars.minstepsize=0;
pars.maxiter=1000;
pars.factor=1.1;
pars.correction=15;
pars.thresho=1e-7;
pars.thresha=1e-22;
pars.scale=0;
pars.obj=0;
pars.quiet=0;
pars.classsplit=0;
pars.validation=0;
pars.validationstep=10;
pars.earlystopping=0;
pars.valrand=1;

pars.stepgrowth=1.01;
pars.weight1=0.5;

% by shuji
% pars.maximp=100000;
pars.maximp=1000;
pars.maximp0=2000;
% pars.maximp0=1000000;

pars.treesize=20;
pars.checkup=2; %0=notree 1=tree 2=choose
pars.targetlabels=[];

pars.regweight=1.0;
pars.regweight0=1.0;

%pars.mtweight = 1.0; 

% keyboard;

pars=extractpars(varargin,pars);

doingcorrection=pars.correction;

% verification dataset
if(pars.validation<0 | pars.validation>1)
    error('validation parameter should be >0 and <1. Thanks.');
end;

if(~isempty(pars.targetlabels))
    pars.targetlabels=pars.targetlabels(itr);
end;


obj=zeros(1,pars.maxiter);
nimp=zeros(1,pars.maxiter);
verify=zeros(length(task),0);
bestMeanErr = Inf;
bestiter = 0;

% lazy workaround 
oldval = pars.validation;
pars.validation = 0;

for i=1:length(task) 
 if isfield(task(i),'L'),L=task.L; else L=L0; end; %zeros(size(L0));end;
 if isfield(task(i),'L_for_nn'), pars.L_for_nn = task(i).L_for_nn; end;
 if isfield(task(i),'Kg'),k=task.Kg; else k=Kg;end;
 state(i)=initstate(L,task(i).x,task(i).y,k,pars);
 state(i).xv = task(i).xv;
 state(i).yv = task(i).yv;
 bestMgen.M{i} = state(i).M;
 state(i).obj = obj;
end;

% lazy workaround ending
pars.validation = oldval;
        
stepsize0=1e-07;
                                               
clear('Kg')      

% Main Loop
for iter=1:pars.maxiter 

     % save old position
     oldstate=state;
     M0old=M0;

     % take gradient step for M0        
     if iter>1,                                   
       M0=M0-stepsize0.*pars.regweight0.*2.*(M0-eye(size(M0)));
       for st=1:length(state),
           % global metric
           M0=M0-stepsize0.*mat(state(st).Gradient);
       end;                   
       M0=makepsd(M0,pars); % project onto PSD cone

       for st=1:length(state),       
           % local metric
           % state(st).stepsize=stepsize0;                               
           state(st).Gradient=state(st).Gradient+2*pars.regweight.*vec(state(st).M);
    %       state(st).Gradient=state(st).Gradient+pars.regweight.*ones(size(state(st).Gradient));
           state(st)=step(state(st),pars,stepsize0);
        end;
     end;                       


     % validation
     for st=1:length(state)
       % disp(iter)
       [state(st),dobreak_1(st),new_verify(st,:)]=validate(iter,state(st),verify(st,:),pars,st);
     end
     verify = new_verify;
     clear new_verify;
     if(pars.validation>0 && (mod(iter,pars.validationstep)==0 || iter==1) && mean(verify(:,end)) <= bestMeanErr)
         % disp(verify)
         bestiter = iter;
         bestMeanErr = mean(verify(:,end));
         bestMgen.M0 = M0;
         for st=1:length(state)
           bestMgen.M{st} = state(st).M;
         end
     end
    if all(dobreak_1),break;end;


     for st=1:length(state)
      % compute gradient
      state(st)=computegradient(state(st),pars);

      % compute objective and # of impostors
      o=(state(st).Gradient)'*vec( state(st).M(:)+ M0(:))+state(st).totalactive.*(1-pars.weight1);
      o=o+pars.regweight.*sum(sum((state(st).M).^2));
      state(st).obj(iter) = o;
      obj(iter)=obj(iter)+o;
      % fprintf(' %2.4f ',o);
      % add regularizer
      nimp(iter)=nimp(iter)+state(st).totalactive;

      % print status message
      % printstatus(iter,obj,state(st),pars);
     end;                            

     obj(iter)=obj(iter)+pars.regweight0.*sum(sum((M0-eye(size(M0))).^2));
     printstatus(iter,obj(iter),nimp(iter),0,stepsize0,pars)     

     % if iter >= 245, keyboard; end

     % increase stepsize or (decrease stepsize and restore state)
     [obj,state,stepsize0]=updatestepsize(iter,obj,state,oldstate,stepsize0,pars);

     % update state stepsizes: shibin
     for st = 1:length(state)
       %state(st).correction = doingcorrection;
       %[state(st).obj,state(st)]=updatestatestepsize(iter,state(st).obj,state(st),oldstate(st),pars);   
       state(st).stepsize = stepsize0;
     end

     % check if we should terminate
    for st = 1:length(state)
     [state(st),pars,dobreak_2(st)]=terminationcheck(iter,obj,state(st),pars);
    end 
    if all(dobreak_2),break;end;

      % update impostors 
      doingcorrection=doingcorrection-1;
      if doingcorrection==0,
          for st=1:length(state)               
              [state(st),added(st),pars]=findimpostors(state(st),pars);
           end;
           if sum(added <= 0)
              pars.correction=min(pars.correction*2+2,300);
              doingcorrection=pars.correction - 1;
           else
              pars.correction=round(pars.correction*pars.factor); 
              doingcorrection=pars.correction;
           end     
      end;
      % save file edited by Shuji_Hao
      % if(mod(iter,10)==0) 
      %   save(tempname,'state','iter','obj','pars','verify');
      % end;
end;
%keyboard;
% Output
Det.obj=obj(1:iter);
Det.nimp=nimp(1:iter);
Det.pars=pars;
Det.time=toc;
Det.iter=iter;
Det.verify=verify;

Det.bestMgen = bestMgen;
Det.bestitergen = bestiter;   
Det.M0=M0;
Det.M={};
for i=1:length(task)
    Det.M{i}=state(i).M;
    Det.bestM{i} = state(i).bestM;
    Det.bestM0{i} = state(i).bestM0;
    Det.bestiter(i) = state(i).bestiter;
end;    
            
% delete(tempname)
% fprintf('\n');     
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     
     
     
     
     
                              
function state=initstate(L,x,y,Kg,pars)
global M0
    
% generate state 
%[itr,ite]=makesplits(y,1-pars.validation,1,pars.classsplit,Kg+1,pars.valrand);
%state.xv=x(:,ite);
%state.yv=y(:,ite);
%state.x=x(:,itr);
%state.y=y(itr);
state.x = x;
state.y = y;
state.xv = [];
state.yv = [];
state.obj = [];
state.correction = pars.correction;
if isfield(pars,'L_for_nn')
%    fprintf('here \n')
  [gen,NN]=getGenLS(pars.L_for_nn*state.x,state.y,Kg,pars);
else
  [gen,NN]=getGenLS(state.x,state.y,Kg,pars);
end
state.NN=NN;
state.dfG=vec(SOD(state.x,gen(1,:),gen(2,:)));
state.Kg=Kg;
state.L=L;              
state.M=state.L'*state.L;
state.bestM=state.M;
state.bestM0=M0;
state.stepsize=pars.stepsize;
state=flushstate(state,pars);
state.totalactive=0;                                      
state.besterr=inf;
state.earlycounter=0;
state.bestiter=0; 

state.imp=checkup(decompose( state.M + M0 ),state.x,state.y,state.NN(end,:),pars,1);       
if(size(state.imp,2)>pars.maximp0)
 ip=randperm(size(state.imp,2));
 ip=ip(1:pars.maximp0);
 state.imp=state.imp(:,ip);    
end;
if(pars.scale)
 Lx=state.L*x;
 sc=sqrt(mean(sum( ((Lx-Lx(:,state.NN(end,:)))).^2,1)));
 state.L=state.L./(sc+2);
end;



   

function [state,dobreak,verify]=validate(iter,state,verify,pars,st)                 
global M0;                  
dobreak=false;
if pars.validation==0,state.bestM=state.M;return;end; %state.bestL=state.L;return;end;              
if(pars.validation>0 && (mod(iter,pars.validationstep)==0 || iter==1))
 %% if iter == 1 disp('here'); keyboard; end
 temp = knnclassify(decompose(M0+state.M),state.x,state.y,state.xv,state.yv,state.Kg).*100;
 verify=[verify temp(2)];
 %fprintf('kNN validation error: %2.2f ',verify(end)*100);
 
 if(verify(end)<=state.besterr) 
     %fprintf('<= %2.2f   :-) \n',state.besterr*100);
     state.besterr=verify(end);     
     state.bestM = state.M;
     state.bestM0 = M0;
     state.bestiter=iter;
     state.earlycounter=0;
 else
     %fprintf('> %2.2f   :-( %i/%i\n',state.besterr*100,state.earlycounter,pars.earlystopping);
     state.earlycounter=state.earlycounter+1;
 end;                                         
 if(pars.earlystopping>0 && state.earlycounter>pars.earlystopping)
       % fprintf('\nTask %d: Validation error is no longer improving!\n',st);
       dobreak=true;
 end;
end;


    
            
function printstatus(iter,obj,nimp,timp,stepsize,pars)    
    if(~pars.quiet)
        fprintf(['%i.  Obj:%2.2f Nimp:%i/%i  Stepsize:%e           \n   '],iter,obj,nimp,timp,stepsize);
     end;
            
function [state,added,pars]=findimpostors(state,pars)
global M0 doingcorrection;
    

   if(~pars.quiet)fprintf('');end;
   Vio=setdiff(checkup(decompose(state.M+M0),state.x,state.y,state.NN(state.Kg,:),pars)',state.imp','rows')';
   if(pars.maximp<inf)
     i=randperm(size(Vio,2));
     Vio=Vio(:,i(1:min(pars.maximp,size(Vio,2))));
   end;

   ol=size(state.imp,2);
   [imp i1 i2]=unique([state.imp Vio].','rows');
   state.imp=imp.';
   if(size(state.imp,2)~=ol)
      for nnid=1:state.Kg;
	   state.a1{nnid}=i2(state.a1{nnid});
	   state.a2{nnid}=i2(state.a2{nnid});
      end;
   end;                                                                                               
   added=size(state.imp,2)-ol;
   % put next correction a little more into the future if no new impostors were added
   if(added<=0) 
      %pars.correction=min(pars.correction*2+2,300);
      % doingcorrection=pars.correction-1;
   else
      %pars.correction=round(pars.correction*pars.factor); %shibin
      % doingcorrection=pars.correction;
   end;
if ~pars.quiet,fprintf('\nAdded %i new impostors\n',added);end;
                   
function [state,pars,dobreak]=terminationcheck(iter,obj,state,pars)
% check if we have terminated       
global doingcorrection;
dobreak=false;
if(iter<=10) return;end;
 if (max(abs(diff(obj(iter-10:iter))))<pars.thresho*obj(iter)  | state.stepsize<pars.thresha)
  if(pars.correction-doingcorrection>=5) 
     doingcorrection=1;
  else
    switch(pars.obj)
     case 0
      if(~pars.quiet)fprintf('Stepsize too small. No more progress!\n');end;
      dobreak=true;
     case 1
      pars.obj=0;
      pars.correction=15;
      pars.stepsize=1e-9;
      doingcorrection=0;
      state=flushstate(state,pars);
      if(~pars.quiet) 
        fprintf('\nVerifying solution! %i\n',obj(iter)); 
      end;
     case 3
      if(~pars.quiet)fprintf('Stepsize too small. No more progress!\n');end;
      dobreak=true;
     end;
  end;
 end;
                                  
function state=flushstate(state,pars)
    for nnid=1:state.Kg; 
        state.a1{nnid}=[];
        state.a2{nnid}=[];
    end;
    state.df=zeros(size(state.dfG));       
    state.Gradient=state.dfG .* pars.weight1;
  
function [obj,state]=updatestatestepsize(iter,obj,state,oldstate,pars)
    % either increases or decreases stepsize by a multiplicative factor
global doingcorrection;
                                   
delta=obj(iter)-obj(max(iter-1,1));            
if(iter>1 & delta>0 & state.correction~=pars.correction) 
        % correct stepsize and undo last step
        state.stepsize=state.stepsize*0.5;
        % fprintf('***correcting stepsize***\n');
        if(state.stepsize<pars.minstepsize) state.stepsize=pars.minstepsize;end;
        state=oldstate;
        obj(iter)=obj(iter-1);
else 
        if(state.correction~=pars.correction) state.stepsize=state.stepsize*pars.stepgrowth;end;
end;
                              
function [obj,state,stepsize0]=updatestepsize(iter,obj,state,oldstate,stepsize0,pars)
    % either increases or decreases stepsize by a multiplicative factor
global doingcorrection M0 M0old;
% if iter >= 16, keyboard; end
                                   
delta=obj(iter)-obj(max(iter-1,1));            
if(iter>1 & delta>0 & doingcorrection~=pars.correction) 
        % correct stepsize and undo last step
        % fprintf('***correcting stepsize***\n');
        if(stepsize0<pars.minstepsize) stepsize0=pars.minstepsize;end;
        state=oldstate;
        obj(iter)=obj(iter-1);
        M0=M0old;
        stepsize0=stepsize0*0.5;
else 
        if(doingcorrection~=pars.correction) stepsize0=stepsize0*pars.stepgrowth;end;
end;
                                                                       

    
    

function state=computegradient(state,pars)
global M0;
                
    L=decompose( state.M + M0 );
    N=length(state.y);        
    Lx=L*state.x;
    state.totalactive=0;
                                                 
    g0=cdist(Lx,state.imp(1,:),state.imp(2,:));

    kk=1;
    for nnid=kk:state.Kg
     Ni(nnid,1:N)=(sum((Lx-Lx(:,state.NN(nnid,:))).^2,1)+1);
    end;

     g1=Ni(:,state.imp(1,:)); 
     g2=Ni(:,state.imp(2,:)); 
     act1=[];act2=[];

    for nnid=state.Kg:-1:kk
      act1=find(g0<g1(nnid,:)); 
      act2=find(g0<g2(nnid,:)); 
      active=[act1 act2];
     if(~isempty(state.a1{nnid}) | ~isempty(state.a2{nnid}))
      [plus1,minus1]=sd(act1(:)',state.a1{nnid}(:)');
      [plus2,minus2]=sd(act2(:)',state.a2{nnid}(:)');
     else
      plus1=act1;plus2=act2;
      minus1=[];minus2=[];
     end;

     MINUS1a=[state.imp(1,minus1) state.imp(2,minus2)]; MINUS1b=[state.imp(1,[plus1 plus2])];
     MINUS2a=[state.NN(nnid,state.imp(1,minus1)) state.NN(nnid,state.imp(2,minus2))]; MINUS2b=[state.imp(2,[plus1 plus2])];

     [isplus2,i]= sort(state.imp(2,plus2));plus2=plus2(i);
     PLUS1a=[state.imp(1,plus1) isplus2]; PLUS1b=[state.imp(1,[minus1 minus2])];
     PLUS2a=[state.NN(nnid,state.imp(1,plus1)) state.NN(nnid,isplus2)]; PLUS2b=[state.imp(2,[minus1 minus2])];

     loss1=max(g1(nnid,:)-g0,0);
     loss2=max(g2(nnid,:)-g0,0);

     [PLUS ,pweight]=count([PLUS1a;PLUS2a]);
     [MINUS,mweight]=count([MINUS1a;MINUS2a]);


     df2=SODW(state.x,PLUS(1,:),PLUS(2,:),pweight)-SODW(state.x,MINUS(1,:),MINUS(2,:),mweight);
     df4=SOD(state.x,PLUS1b,PLUS2b)-SOD(state.x,MINUS1b,MINUS2b);
     state.df=state.df+vec(df2+df4);

     state.a1{nnid}=act1;state.a2{nnid}=act2;
     state.totalactive=state.totalactive+length(active);

    end;
state.Gradient=state.dfG.*pars.weight1+state.df.*(1-pars.weight1);



function state=step(state,pars,stepsize0);
global M0;                    
                    
G=mat(state.Gradient);
% do step in gradient direction
switch(pars.obj)
  case 0    % updating Q
     Q=state.M-state.stepsize.*G;
   case 1   % updating L
     G=2.*(state.L*G);
     state.L=state.L-state.stepsize.*G;     
     state.M=state.L'*state.L;
     return;
  case 3
     Q=state.L'*state.L;
	 Q=Q-stepsize.*G;
	 Q=diag(Q);
 	 state.L=diag(sqrt(max(Q,0)));
     state.M=state.L'*state.L;
     return;
  otherwise
   error('Objective function has to be 0,1,2\n');
end;                     
% Make sure Q+M0 is psd
[M_psd]=makepsd( Q + M0 , pars);

state.M = M_psd - M0;

state.L=decompose(state.M);

        

function [M,L]=makepsd(Q,pars);
    % decompose Q
    [L,dd]=eig(Q);
    dd=real(diag(dd));
    L=real(L);
    % reassemble Q (ignore negative eigenvalues)
    j=find(dd<1e-10);
    if(~isempty(j)) 
        if(~pars.quiet)fprintf('[%i]',length(j));end;
    end;
    dd(j)=0;
    [temp,ii]=sort(-dd);
    L=L(:,ii);
    dd=dd(ii);
    
    L=(L*diag(sqrt(dd)))';
    M=L'*L;






function [gen,NN]=getGenLS(x,y,Kg,pars);
if(~pars.quiet);fprintf('Computing nearest neighbors ...\n');end; %#ok<SEPEX>
[D,N]=size(x);

if(~isempty(pars.targetlabels))
    y=pars.targetlabels;
end;

un=unique(y);
Gnn=zeros(Kg,N);
for c=un
%  fprintf('%i nearest genuine neighbors for class %i:',Kg,c);
 i=find(y==c);
 nn=LSKnn(x(:,i),x(:,i),2:Kg+1);
 Gnn(:,i)=i(nn);
%  fprintf('\r');
end;

% fprintf('\n');

NN=Gnn;
gen1=vec(Gnn(1:Kg,:)')';
gen2=vec(repmat(1:N,Kg,1)')';

gen=[gen1;gen2];



function imp=checkup(L,x,y,NN,pars,first)
persistent treetime notreetime;
if(nargin==6)
    treetime=-1;
    notreetime=-1;
else
    first=0;
end;
t1=toc;
if(treetime<notreetime)
  imp=checkupmtree(L,x,y,NN,pars,first);treetime=toc-t1;
else
  imp=checkupnotree(L,x,y,NN,pars,first);notreetime=toc-t1;
end;
 

function imp=checkupmtree(L,x,y,NN,pars,first)
if(~pars.quiet);fprintf('Computing nearest neighbors ...\n');end;
[D,N]=size(x);

mL=max(L');
L=L(find(mL),:);
Lx=L*x;
Ni=sum((Lx-Lx(:,NN)).^2,1)+1;
un=unique(y);

% build up ball trees
for c=1:length(un)
 classindex{c}=find(y==un(c));
 forest{c}.tree=buildmtreemex(Lx(:,classindex{c}),pars.treesize);
end;
imp=[];
for c=1:length(un)-1
if(~pars.quiet)fprintf('All impostors for class %i    \r',c);end;
 for c2=c+1:length(un)
      limps=findNimex(forest{c2}.tree,Lx(:,classindex{c2}),Lx(:,classindex{c}),Ni(classindex{c2}),Ni(classindex{c}));         
      
%    keyboard;
 	if(size(limps,2)>pars.maximp)
  	 ip=randperm(size(limps,2));
  	 ip=ip(1:pars.maximp);
  	 limps=limps(:,ip);
  	end;
	limps=[classindex{c}(limps(1,:));classindex{c2}(limps(2,:))];
    imp=[imp limps];
   end;
end;
try
 imp=unique(sort(imp)','rows')';
catch
 fprintf('Sorry, probably ran out of memory!');
 keyboard;  
end;



function imp=checkupnotree(L,x,y,NN,pars,first)
if(~pars.quiet) fprintf('Computing nearest neighbors ...\n');end;
[D,N]=size(x);

Lx=L*x;
Ni=sum((Lx-Lx(:,NN)).^2,1)+1;

un=unique(y);
imp=[];
index=1:N;
for c=un(1:end-1)
 if(~pars.quiet)fprintf('All nearest impostor neighbors for class %i :',c);end;
 i=index(find(y(index)==c));
 index=index(find(y(index)~=c));
 limps=LSImps2(Lx(:,index),Lx(:,i),Ni(index),Ni(i),pars);
 if(size(limps,2)>pars.maximp)
  ip=randperm(size(limps,2));
  ip=ip(1:pars.maximp);
  limps=limps(:,ip);
 end;
 imp=[imp [i(limps(2,:));index(limps(1,:))]];

 if(~pars.quiet)fprintf('\r');end;
end;
try
 imp=unique(sort(imp)','rows')';
catch
 fprintf('Sorry, probably ran out of memory!');
 keyboard;  
end;



function limps=LSImps2(X1,X2,Thresh1,Thresh2,pars);
B=750;
[D,N2]=size(X2);
N1=size(X1,2);
limps=[];
for i=1:B:N2
  BB=min(B,N2-i);
  try
  newlimps=findimps3Dac(X1,X2(:,i:i+BB), Thresh1,Thresh2(i:i+BB));
%  t1=toc;D=distance(X1,X2(:,i:i+BB));newlimps=findimps3Ddist(D,Thresh1,Thresh2(i:i+BB));toc-t1
  if(~isempty(newlimps) & newlimps(end)==0)    
    [minv,endpoint]=min(min(newlimps));
    newlimps=newlimps(:,1:endpoint-1);
  end;
  newlimps=unique(newlimps','rows')';
  catch
    keyboard;
    end;
  newlimps(2,:)=newlimps(2,:)+i-1;
  limps=[limps newlimps];
  if(~pars.quiet)fprintf('(%i%%) ',round((i+BB)/N2*100)); end;
end;
if(~pars.quiet)fprintf(' [%i] ',size(limps,2));end;




function NN=LSKnn(X1,X2,ks,pars);
B=750;
[D,N]=size(X2);
NN=zeros(length(ks),N);
DD=zeros(length(ks),N);

for i=1:B:N
  BB=min(B,N-i);
%   fprintf('.');
  Dist=distance(X1,X2(:,i:i+BB));
%   fprintf('.');
  [dist,nn]=mink(Dist,max(ks));
  clear('Dist');
%   fprintf('.'); 
  NN(:,i:i+BB)=nn(ks,:);
  clear('nn','dist');
%   fprintf('(%i%%) ',round((i+BB)/N*100)); 
end;
  

function [L,dd]=decompose(M);
% function [L,dd]=decompose(M);
%
% decomposees the (positive semidefinite) matrix M into L and dd such that M=L'*L
% where the ROWS of L are the eigenvectors of M
% sorted with decreasing eigenvalues (which are stored in dd)
% 
% copyright by Kilian Q. Weinberger, 2006
%

[L,dd]=eig(M);
L=real(transpose(L*sqrt(max(dd,0))));
dd=diag(dd);
[temp,ii]=sort(-dd);
L=L(ii,:);
dd=dd(ii);


