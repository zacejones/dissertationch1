clear;
clc;

global c_gamma c_rho I_N IS_bar L_N L_S c_sigma gL c_lambda Phi_IS


factual_N_import_share=0.2270857*100; %Calculated in R
factual_S_import_share=0.1388399*100; %Calculated in R
factual_omega=131047.3/31416.3; %Data from World Bank USGDPWorker/CHNGDPWorker

loop_T_N=[0.1:0.1:1.5];
loop_T_S=[0.1:0.1:1.5];
loop_alpha=[2.1:0.1:4];

diff_matrix=0;
N_SHARE_MATRIX=0;
S_SHARE_MATRIX=0;

for i = 1:length(loop_T_N)
    for j = 1:length(loop_T_S)
        for k = 1:length(loop_alpha)
            
PARAMETERS_SEARCH_3;



t_N = loop_T_N(i);
t_S = loop_T_S(j);
T_N = 1 + t_N;
T_S = 1 + t_S;

c_alpha = loop_alpha(k);

%myfminunc_handle=@(x) myconstraints(x,m,N);
%[x]=fminunc(myfminunc_handle,wage_0,options);

tracker=[t_N,t_S,c_alpha]
tol=1e-12;
options=optimset('display','off','TolFun',tol,'MaxFunEvals',inf,'MaxIter',5000);%'plotfcns',{@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize, @optimplotfirstorderopt});
steadystate_handle = @(x) steady_state(x,T_S,T_N,c_alpha);
x0 = [3,0.1,3,3,1,0.2,0.2];
ub=[10,1,20,10,10,10,10];
lb=[0,0,0,0,0,0,0];
[x] = fmincon(steadystate_handle,x0,[],[],[],[],lb,ub,@constraint,options);

x_N = x(1);
I_F = x(2);
c_omega = x(3);
c_N = x(4);
c_S = x(5);
Psi_N = x(6);
Psi_S = x(7);

test1 = c_alpha/c_gamma - (((c_rho+I_N)/(c_rho+I_N+IS_bar))*((c_N*L_N/Psi_N)*T_S^(-c_sigma)+(c_S*L_S/Psi_S))/((c_N*L_N/Psi_N)+(c_S*L_S/Psi_S)*T_N^(-c_sigma))*c_omega^(c_sigma) - c_omega);
test2 = c_N - (c_omega+(c_rho-gL)*c_omega*c_gamma*x_N*((c_lambda*I_N/(c_lambda*I_N+I_F))+((c_lambda*I_N/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)))));
test3 = c_S - (1 + (c_rho-gL)*c_alpha*(c_lambda*I_N/(c_lambda*I_N+IS_bar))*(1/((c_sigma-1)*(c_alpha+c_gamma*c_omega)*(c_rho+I_N+IS_bar)*(Phi_IS)+c_alpha*c_lambda*I_N)));
test4 = Psi_N - ((c_sigma*c_omega/(c_sigma-1))^(1-c_sigma)*c_lambda*I_N/(c_lambda*I_N+I_F)+T_N^(1-c_sigma)*(c_sigma/(c_sigma-1))^(1-c_sigma)*(c_lambda*I_N/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)) + T_N^(1-c_sigma)*(IS_bar/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)));
test5 = Psi_S - ((c_sigma*c_omega*T_S/(c_sigma-1))^(1-c_sigma)*c_lambda*I_N/(c_lambda*I_N+I_F)+(c_sigma/(c_sigma-1))^(1-c_sigma)*(c_lambda*I_N/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)) + (IS_bar/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)));
test6 = 1 - (c_gamma*x_N*((c_sigma-1)*(c_rho+I_N)*(c_lambda*I_N/(c_lambda*I_N+I_F))+I_N));
test7 = 1 - (x_N*(L_N/L_S)*(I_F/(c_lambda*I_N + I_F))*((c_sigma-1)*(c_alpha+c_gamma*c_omega)*(c_rho+I_N+IS_bar)*Phi_IS +c_alpha*c_lambda*I_N));

testsolution=[
test1,
test2,
test3,
test4,
test5,
test6,
test7];

testsolution=double(testsolution);
sum(testsolution.^2);

nN = I_N/(I_N+I_F);
nF = I_F/(I_N+I_F)*I_N/(I_N+IS_bar);
nS = I_F/(I_N+I_F)*IS_bar/(I_N+IS_bar);

S_import_share = 1/Psi_S *(T_S)^(1-c_sigma)*(c_omega*c_sigma/(c_sigma-1))^(1-c_sigma)*(c_lambda*I_N/(c_lambda*I_N + I_F))*100;
N_import_share = 1/Psi_N *(T_N)^(1-c_sigma)* (((c_sigma/(c_sigma-1))^(1-c_sigma)*(I_F/(c_lambda*I_N+I_F))*(c_lambda*I_N/(c_lambda*I_N+IS_bar)))+I_F/(c_lambda*I_N +I_F)*IS_bar/(c_lambda*I_N+IS_bar))*100;

diff_N = (S_import_share-factual_S_import_share)^2;
diff_S = (N_import_share-factual_N_import_share)^2;
diff_omega = (factual_omega - c_omega)^2;

N_SHARE_MATRIX(i,j,k)=N_import_share;
S_SHARE_MATRIX(i,j,k)=S_import_share;
alpha_MATRIX(i,j,k)=c_omega;

diff_MATRIX(i,j,k)=diff_N + diff_S + diff_omega; %See if moments are matched
if sum(testsolution.^2)>1e-8 %Eliminate minimums where the objective function isn't satisfied
    diff_MATRIX(i,j,k)=10;
end
solution_MATRIX(i,j,k)=sum(testsolution.^2); %See if objective is solved.


        end
    end
end


