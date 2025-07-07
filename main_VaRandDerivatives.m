clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VaR and Derivatives
%%%%%%%%%%%%%%%%%%%%%%%%%%
imgDir = 'Images/'; % Directory for saving figures
txtDir = 'Results/'; % Directory for saving results
txtFilename = fullfile(txtDir, 'VaR_Derivatives.txt'); % Output file for results

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%
P0 = 305;        % Underlying stock price
K = 300;        % Option strike price
TTM = 4/12;        % Option time to maturity (in years)
sg = 0.25;       % Annualized implied volatility
rf = 0.08;      % Risk-free rate
q = 0.03;          % Dividend yield
alpha = 0.95;   % Confidence level
days = 80;     % VaR horizon in days
dt = days / 250; % Time step for VaR horizon
M = 1000000;    % Number of simulations
C0 = blsprice(P0, K, rf, TTM, sg, q); % Price the option using Black-Scholes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHOD 1: VaR Exact
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
worst_r = norminv(1 - alpha) * sg * sqrt(dt); % Worst-case return
worst_P = P0 * exp(worst_r); % Worst-case stock price
worst_C = blsprice(worst_P, K, rf, TTM - dt, sg, q); % Worst-case option price
VaR_C = C0 - worst_C; % Exact VaR
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHOD 2: Monte Carlo + Full Revaluation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% Simulate the risk factor (returns)
r = sg * sqrt(dt) * randn(M, 1);
% Simulate future stock prices
PT = P0 * exp(r);
% Reprice the option at future prices
CT = blsprice(PT, K, rf, TTM - dt, sg, q);
% Compute option profit and loss
C_PL = CT - C0;
% Compute Expected Shortfall (ES) using a helper function
[VaR_MC, ES_MC] = get_riskmeasures('NP', C_PL, alpha);
toc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METHOD 3: Delta-Gamma Monte Carlo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
% Step 1: Compute the Greeks
Theta = blstheta(P0, K, rf, TTM, sg, q);
Delta = blsdelta(P0, K, rf, TTM, sg, q);
Gamma = blsgamma(P0, K, rf, TTM, sg, q);

% Step 2: Simulate returns
r = sg * sqrt(dt) * randn(M, 1);

% Step 3: Simulate future stock prices
PT = P0 * exp(r);

% Step 4: Simulate option P&L using Greeks
% 4.1 Delta-only method
CT_d = C0 + Theta * dt + Delta * (PT - P0);

% 4.2 Delta-Gamma method
CT_dg = C0 + Theta * dt + Delta * (PT - P0) + 0.5 * Gamma * (PT - P0).^2;

% Step 5: Compute VaR and ES using sample quantities
% 5.1 VaR-Delta MC and Exact
C_PL_d = CT_d - C0;
[VaR_d_MC, ES_d_MC] = get_riskmeasures('NP', C_PL_d, alpha);
P_worst = quantile(PT, 1 - alpha);
VaR_d_ex = -(Theta * dt + Delta * (P_worst - P0));

% 5.2 VaR-Delta-Gamma MC and Exact
C_PL_dg = CT_dg - C0;
[VaR_dg_MC, ES_dg_MC] = get_riskmeasures('NP', C_PL_dg, alpha);
VaR_dg_ex = -(Theta * dt + Delta * (P_worst - P0) + 0.5 * Gamma * (P_worst - P0).^2);
toc;

% Compare the results
Output_VaR = table(days, VaR_C, VaR_MC, VaR_d_ex, VaR_d_MC, VaR_dg_ex, VaR_dg_MC)
Output_ES = table(days, ES_MC, ES_d_MC, ES_dg_MC)

%Write to a text file
log_to_file(Output_VaR, txtFilename);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting the distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minPL = min(C_PL_d);
perPL = prctile(C_PL_d, 30);
maxPL = max(C_PL_d);

h1 = figure('Color', [1 1 1]);

% Full Monte Carlo
subplot(3, 2, 1);
histogram(C_PL, 100, 'Normalization', 'pdf');
xlabel('Option P\&L', 'Interpreter', 'latex');
title('Full Monte Carlo', 'Interpreter', 'latex');
hold on;
plot(-VaR_MC, 0, 'r*');
plot(-VaR_C, 0, 'g*');
fplot(@(x) normpdf(x, mean(C_PL), std(C_PL)), [minPL, maxPL])
xlim([minPL, maxPL+0.1]);

% Left tail (Full Monte Carlo)
subplot(3, 2, 2);
histogram(C_PL, 100, 'Normalization', 'pdf');
xlim([minPL, perPL]);
hold on;
plot(-VaR_MC, 0, 'r*');
plot(-VaR_C, 0, 'g*');
xlabel('Option P\&L', 'Interpreter', 'latex');
title('Left Tail', 'Interpreter', 'latex');

% Delta Monte Carlo
subplot(3, 2, 3);
histogram(CT_d - C0, 100, 'Normalization', 'pdf');
hold on;
plot(-VaR_d_ex, 0, 'r*');
plot(-VaR_C, 0, 'g*');
xlabel('Option P\&L', 'Interpreter', 'latex');
title('Delta Monte Carlo', 'Interpreter', 'latex');
xlim([minPL, maxPL]);

% Left tail (Delta Monte Carlo)
subplot(3, 2, 4);
histogram(CT_d - C0, 100, 'Normalization', 'pdf');
hold on;
plot(-VaR_d_MC, 0, 'r*');
plot(-VaR_C, 0, 'g*');
xlim([minPL, perPL]);
xlabel('Option P\&L', 'Interpreter', 'latex');
title('Left Tail', 'Interpreter', 'latex');

% Delta-Gamma Monte Carlo
subplot(3, 2, 5);
histogram(CT_dg - C0, 100, 'Normalization', 'pdf');
hold on;
plot(-VaR_dg_ex, 0, 'r*');
plot(-VaR_C, 0, 'g*');
xlabel('Option P\&L', 'Interpreter', 'latex');
title('Delta-Gamma Monte Carlo', 'Interpreter', 'latex');
xlim([minPL, maxPL]);

% Left tail (Delta-Gamma Monte Carlo)
subplot(3, 2, 6);
histogram(CT_dg - C0, 100, 'Normalization', 'pdf');
hold on;
plot(-VaR_dg_MC, 0, 'r*');
plot(-VaR_C, 0, 'g*');
xlim([minPL, perPL]);
xlabel('Option P\&L', 'Interpreter', 'latex');
title('Left Tail', 'Interpreter', 'latex');
sgtitle(['P\&L distribution (Horizon: ' num2str(days) ' days)'], 'interpreter', 'latex');
saveas(h1, fullfile(imgDir, ['RA_fig_Comparing_VaR_horizon_' num2str(days) '.png']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Approximating the Greeks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 1e-5; % Small perturbation
C0 = blsprice(P0, K, rf, TTM, sg, q);
C0p = blsprice(P0 + h, K, rf, TTM, sg, q);
C0m = blsprice(P0 - h, K, rf, TTM, sg, q);

% Finite difference approximations
Fwd_diff = (C0p - C0) / h;
Bwd_diff = (C0 - C0m) / h;
Cent_diff = (C0p - C0m) / (2 * h);

% Compare with analytical Delta
check = table(Delta, Fwd_diff, Bwd_diff, Cent_diff);
disp(check);

SecDer = (C0p - 2 * C0 + C0m)/h^2
Gamma