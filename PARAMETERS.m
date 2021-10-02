%%%%%%%%%%%%%%%%%
%%% PARAMETERS %%
%%%%%%%%%%%%%%%%%
%DS is Dinopoulos and Segerstrom 2007 (unpublished)

global c_alpha c_gamma c_rho I_N IS_bar L_N L_S c_sigma gL c_lambda Phi_IS

%Population Growth Rate (Sener = 0.01) (DS = 0.014)
gL = 0.014;

% Innovation size (Sener=1.25) (DS=1.7)
c_delta = 1.25;

% Discount Rate (Sener = 0.07) (DS = 0.07)
c_rho = 0.07;

% CES (Dinopoulous, Heins, Unel = 2) (DS=1.5)(3.42 Ossa), (4 Broda,
% Weinstein) (5 ~Soddenberry) (5.26 Soddenberry Table 4)
c_sigma = 5.26;


% South Copying Probability (D&S 2007 Footnote 12)
IS_bar = 0.05;

% Innovative R+D Productivity (Sener=1.25, to calibrate growth = 0.5%) (DS = 1)
c_gamma = 1;

% Adaptive R+D Productivity (for now just copying innovative productivity,
% but this should really receive its own calibration (DS = 3.5) (ZJ=10.4
% calibration for sigma=3.42) (7.4 for sigma=4 ZJ) (3.1 for sigma=5.26 ZJ_)
c_alpha = 3.1;

% Product Quality
c_lambda = c_delta^(c_sigma-1);

% Innovative R+D Rate
I_N = gL/(c_lambda-1);

%Labor (DS=1,2 Respectively)
L_N = 1;
L_S = 2;

%Phi_IS, an endogenous parameter
Phi_IS=(c_lambda*I_N/(c_lambda*I_N+IS_bar))+(c_sigma/(c_sigma-1))*(IS_bar/(c_lambda*I_N+IS_bar));

%t_N=1.2 ZJ Calibration s=3.42, 1 ZJ Calibration s=4 0.7 for z=5.26 ZJ 
%t_S=0.6 ZJ Calibration s=3.42, .04 ZJ Calibration s=4, .2 for s=5.26 ZJ
t_N = 0.7*1.25;
t_S = 0.2;
T_N = 1 + t_N;
T_S = 1 + t_S;
% Tariffs, to be added later



