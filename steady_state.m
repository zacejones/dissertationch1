function y = steady_state(x,T_S,T_N,c_alpha)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

x_N=x(1);
I_F=x(2);
c_omega=x(3);
c_N=x(4);
c_S=x(5);
Psi_N=x(6);
Psi_S=x(7);

global c_gamma c_rho I_N IS_bar L_N L_S c_sigma gL c_lambda Phi_IS

y(1) = c_alpha/c_gamma - (((c_rho+I_N)/(c_rho+I_N+IS_bar))*((c_N*L_N/Psi_N)*T_S^(-c_sigma)+(c_S*L_S/Psi_S))/((c_N*L_N/Psi_N)+(c_S*L_S/Psi_S)*T_N^(-c_sigma))*c_omega^(c_sigma) - c_omega);

y(2) = c_N - (c_omega+(c_rho-gL)*c_omega*c_gamma*x_N*((c_lambda*I_N/(c_lambda*I_N+I_F))+((c_lambda*I_N/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)))));

y(3) = c_S - (1 + (c_rho-gL)*c_alpha*(c_lambda*I_N/(c_lambda*I_N+IS_bar))*(1/((c_sigma-1)*(c_alpha+c_gamma*c_omega)*(c_rho+I_N+IS_bar)*(Phi_IS)+c_alpha*c_lambda*I_N)));

y(4) = Psi_N - ((c_sigma*c_omega/(c_sigma-1))^(1-c_sigma)*c_lambda*I_N/(c_lambda*I_N+I_F)+T_N^(1-c_sigma)*(c_sigma/(c_sigma-1))^(1-c_sigma)*(c_lambda*I_N/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)) + T_N^(1-c_sigma)*(IS_bar/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)));

y(5) = Psi_S - ((c_sigma*c_omega*T_S/(c_sigma-1))^(1-c_sigma)*c_lambda*I_N/(c_lambda*I_N+I_F)+(c_sigma/(c_sigma-1))^(1-c_sigma)*(c_lambda*I_N/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)) + (IS_bar/(c_lambda*I_N+IS_bar))*(I_F/(c_lambda*I_N+I_F)));

y(6) = 1 - (c_gamma*x_N*((c_sigma-1)*(c_rho+I_N)*(c_lambda*I_N/(c_lambda*I_N+I_F))+I_N));

y(7) = 1 - (x_N*(L_N/L_S)*(I_F/(c_lambda*I_N + I_F))*((c_sigma-1)*(c_alpha+c_gamma*c_omega)*(c_rho+I_N+IS_bar)*Phi_IS +c_alpha*c_lambda*I_N));


y=sum(y.^2);
end

