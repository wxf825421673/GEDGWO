%**************************************************************************
% GEDGWO                                                                  %
% The output function is:                                                 %
% Global_score: the value of best fitness                                 %
% Global_pos:   the best solution for problem                             %
% Convergence_curve: the best fitness value of each iteration             %
%**************************************************************************

function  [Global_score,Global_pos,Convergence_curve]=GEDGWO(SearchAgents_no,FEsmax,lb,ub,dim,fobj,fnum)
%======================Initialization======================================
lb=lb.*ones(1,dim);ub=ub.*ones(1,dim);
Positions=repmat(lb,SearchAgents_no,1)+rand(SearchAgents_no,dim).*(repmat(ub-lb,SearchAgents_no,1));
FEs = 0;% function evaluation
SEL=floor(SearchAgents_no/2); 
rank=zeros(1,SearchAgents_no);
Positions_new=zeros(SearchAgents_no,dim);
Convergence_curve=[];
% =======================Calculating the fitness=========================== 
    fitness  = fobj(Positions',fnum);FEs  =  FEs+SearchAgents_no;       
    [~, index]   = sort(fitness);
    % weights Eq.(10)
    weights = log(SEL+1/2)-log(1:SEL)';
    weights = weights/sum(weights);
    flag = 1;  % if stagnates, flag=0, otherwise flag=1   
    mean_Pop_former=inf;
% ===============================Main loop=================================   
while FEs<FEsmax
    
    i = 1:SearchAgents_no;
    rank(index(i)) = (SearchAgents_no-i+1)/SearchAgents_no ; % Eq.(13): rank the population according to their fitness
    %select alpha belta delta wolf
    Alpha_pos =Positions(index(1),:); 
    Beta_pos = Positions(index(2),:); 
    Delta_pos = Positions(index(3),:);
   
    if flag~=0
    % calculate Covariance matrix Eq.(11)
    Xsel = Positions(index(1:SEL), :);
    xmean = weights'*Xsel;
    C = 1/(SEL)*(Xsel - xmean(ones(SEL,1), :))'*(Xsel - xmean(ones(SEL,1), :));
    C = triu(C) + transpose(triu(C,1)); % enforce symmetry
    [B,D] = eig(C);
    D  = sqrt(diag(D+0.0000001*eye(dim)));   
    end
    
    a=2-2*(FEs/FEsmax); % Eq.(5)  
    z=randn(SearchAgents_no,dim);
for i=1:SearchAgents_no 
    e=floor(3*rand+1);
    leader=index(e); %select leader from alpha,belta,delta wolf randomly
    
    if flag==0 % if stagnates, 
        % =========================Gaussian random walk====================
            Positions_new(i,:) = Positions(leader,:)+z(i,:)*a/2.*abs(Positions(i,:)-Positions(leader,:))+rand*Positions(leader,:)-rand*Positions(i,:); % Eq.(21), r5=rand,r6=rand
    else  % flag=1
      if rank(i)<0.5  %  inferior solutions       
         if rand<a/2  % r4=rand
         % ===============================basic GWO========================   
            A1=2*a*rand(1,dim)-a; % Eqs.(3-8)
            C1=2*rand(1,dim); 
            D_alpha=abs(C1.*Alpha_pos-Positions(i,:)); 
            X1=Alpha_pos-A1.*D_alpha;
            A2=2*a*rand(1,dim)-a; 
            C2=2*rand(1,dim); 
            D_beta=abs(C2.*Beta_pos-Positions(i,:)); 
            X2=Beta_pos-A2.*D_beta;     
            A3=2*a*rand(1,dim)-a; 
            C3=2*rand(1,dim);
            D_delta=abs(C3.*Delta_pos-Positions(i,:));
            X3=Delta_pos-A3.*D_delta;             
            Positions_new(i,:)=(X1+X2+X3)/3;  
         else 
         % ===============================ISR strategy=====================
            meanISR= Positions(leader,:)+rand*(Positions(leader,:)-Positions(i,:)); % Eqs.(19-20), r3=rand 
            Positions_new(i,:)=meanISR+randn*a/2*abs(Positions(i,:)-Positions(leader,:)); % Eq.(17) 
         end  
      else  %  superior solutions
      % =======================GED with shifted mean=======================            
             means=(xmean+Positions(leader,:)+Positions(i,:))/3; % Eq.(15)
             Positions_new(i,:)=means+(B * (D .*z(i,:)'))'; % Eq.(16) 
      end
    end
    % ==========================boundary control===========================
    for j=1:dim
      if Positions_new(i,j)<lb(j)|| ~(isreal(Positions_new(i,j))),   Positions_new(i,j)=rand*(ub(j)-lb(j))+lb(j);  end
      if Positions_new(i,j)>ub(j)|| ~(isreal(Positions_new(i,j))),   Positions_new(i,j)=rand*(ub(j)-lb(j))+lb(j);  end
    end        
end
             fitness_new = fobj(Positions_new',fnum);FEs  =  FEs+SearchAgents_no;
                % ==========================elitism mechanism==============
                elit=fitness_new<fitness ;    % Eq.(23)
                fitness(elit)=fitness_new(elit);
                Positions(elit,:) = Positions_new(elit,:);

         
    % update best solution 
    [~, index]  = sort(fitness);  
     Global_pos  = Positions(index(1),:);
     Global_score= fitness(index(1));
    % stagnation check
    mean_Pop=mean(fitness(index(1:SEL)));
    if mean_Pop==mean_Pop_former 
    flag=0;
    else
    flag=1;
    end
    mean_Pop_former=mean_Pop;
    Convergence_curve= [Convergence_curve,Global_score];  %record best fitness score
end

 