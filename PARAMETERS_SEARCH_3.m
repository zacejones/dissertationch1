%%%%%%%%%%%%%%%%%
%%% PARAMETERS %%
%%%%%%%%%%%%%%%%%
%DS is Dinopoulos and Segerstrom 2007 (unpublished)

global c_gamma c_rho I_N IS_bar L_N L_S c_sigma gL c_lambda Phi_IS

%Population Growth Rate (Sener = 0.01) (DS = 0.014)
gL = 0.014;

% Innovation size (Sener=1.25) (DS=1.7)
c_delta = 1.25;

% Discount Rate (Sener = 0.07) (DS = 0.07)
c_rho = 0.07;

% CES (Dinopoulous, Heins, Unel = 2) (DS=1.5) (3.42 Ossa) (4 Broda,
% Weinstein) (5 ~Soddenberry) (5.26 Soddenberry Table 4)
c_sigma = 5.26;

% South Copying Probability (D&S 2007 Footnote 12)
IS_bar = 0.05;

% Innovative R+D Productivity (Sener=1.25, to calibrate growth = 0.5%) (DS = 1)
c_gamma = 1;

% Product Quality
c_lambda = c_delta^(c_sigma-1);

% Innovative R+D Rate
I_N = gL/(c_lambda-1);

%Labor (DS=1,2 Respectively)
L_N = 1;
L_S = 2;

%Phi_IS, an endogenous parameter
Phi_IS=(c_lambda*I_N/(c_lambda*I_N+IS_bar))+(c_sigma/(c_sigma-1))*(IS_bar/(c_lambda*I_N+IS_bar));





