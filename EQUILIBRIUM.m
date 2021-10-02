clear;
clc;

global c_alpha c_gamma c_rho I_N IS_bar L_N L_S c_sigma gL c_lambda Phi_IS

PARAMETERS;

factual_N_import_share=0.2270857;
factual_S_import_share=0.1388399;
factual_omega=131047.3/31416.3; %Data from World Bank USGDPWorker/CHNGDPWorker


tol=1e-12;
options=optimset('display','off','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000,'plotfcns',{@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize, @optimplotfirstorderopt});
steadystate_handle = @(x) steady_state(x,T_S,T_N,c_alpha);
x0 = [3,0.1,3,3,1,0.2,0.2];
ub=[10,1,20,10,10,10,10];
lb=[0,0,0,0,0,0,0];
[x,fval] = fmincon(steadystate_handle,x0,[],[],[],[],lb,ub,@constraint,options);

x_N = x(1);
I_F = x(2);
c_omega = x(3);
c_N = x(4);
c_S = x(5);
Psi_N = x(6);
Psi_S = x(7);

nN = I_N/(I_N+I_F);
nF = I_F/(I_N+I_F)*I_N/(I_N+IS_bar);
nS = I_F/(I_N+I_F)*IS_bar/(I_N+IS_bar);

S_import_share = 1/Psi_S *(T_S)^(1-c_sigma)*(c_omega*c_sigma/(c_sigma-1))^(1-c_sigma)*(c_lambda*I_N/(c_lambda*I_N + I_F));
N_import_share = 1/Psi_N *(T_N)^(1-c_sigma)* (((c_sigma/(c_sigma-1))^(1-c_sigma)*(I_F/(c_lambda*I_N+I_F))*(c_lambda*I_N/(c_lambda*I_N+IS_bar)))+I_F/(c_lambda*I_N +I_F)*IS_bar/(c_lambda*I_N+IS_bar));

diff_N = (S_import_share-factual_S_import_share)^2;
diff_S = (N_import_share-factual_N_import_share)^2;
diff_omega=(factual_omega-c_omega)^2;
diff_total= diff_N + diff_S+diff_omega;

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

testsolution=double(testsolution)
sum(testsolution.^2)

solution=[x_N,
I_F,
c_omega,
c_N,
c_S,
Psi_N,
Psi_S]

diff_N = (S_import_share-factual_S_import_share)^2;
diff_S = (N_import_share-factual_N_import_share)^2;

    