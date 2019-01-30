clear 
close all
clc

mex cec14_func.cpp -DWINDOWS
SearchAgents_no = 500;       % Population size (you are advised to enlarge the population size for high-dimensional test)
func_num        = 1;         % Number of the benchmarks that can be selected from 1 to 30 
dim             = 30;        % fixed dimensionality of the benchmarks 30D test
FEsmax          = dim*10000; % Maximum function evaluations
lowB            = -100;      % lower boundary
upB             = 100;       % Upper boundary
fobj            = str2func('cec14_func');
optimum         = 100*func_num;
tic;
[Best_score,Best_pos,GEDGWO_cg_curve]=GEDGWO(SearchAgents_no,FEsmax,lowB,upB,dim,fobj,func_num);%Run GEDGWO
toc;
Best_score=Best_score-optimum;if Best_score<1e-08;Best_score=0;end
GEDGWO_cg_curve_error=GEDGWO_cg_curve-optimum;
GEDGWO_cg_curve_error(GEDGWO_cg_curve_error<1e-08)=0;
display(['The function number is : ', num2str(func_num)]);
display(['The best solution obtained by GEDGWO is : ', num2str(Best_pos)]);
display(['The best optimal value of the problem obtained by GEDGWO is : ', num2str(Best_score)]);

%Plot the error value obtained by GEDGWO
[~,n] = size(GEDGWO_cg_curve_error);
x=1:FEsmax/n:FEsmax;
semilogy(x,GEDGWO_cg_curve_error,'r','LineWidth',3);    
xlabel('FEs','FontSize',16 );
ylabel('log(f(x)-f(x^*))','FontSize',16);
set(gca,'FontSize',16);


