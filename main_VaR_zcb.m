
% Clear workspace, close all figures, and clear command window
clear all;
close all;
clc;

% Initial parameters
y0 = 0.1;               % Initial yield
TTM = 1;                % Time to maturity (in years)
BondPrice = @(y, TTM) exp(-y * TTM); % Bond price function
alpha = 0.99;           % Confidence level for VaR
M = 100000;             % Number of Monte Carlo simulations

% Time step and bond price at t=0
dt = 10 / 250;          % Time step (10 days, assuming 250 trading days/year)
B0 = BondPrice(y0, TTM); % Initial bond price

% Yield change parameters
mu_dy = 0.01 * dt;      % Mean yield change
vol_dy = 0.03 * sqrt(dt); % Volatility of yield change

% Calculate worst-case yield and corresponding bond price
best_dy = mu_dy + norminv(alpha) * vol_dy; % Worst-case yield change
best_y = y0 + best_dy;                     % Worst-case yield
worst_B = BondPrice(best_y, TTM - dt);     % Worst-case bond price
VaR_B = B0 - worst_B;                       % Value at Risk (VaR)

%% Monte Carlo Simulation with Delta-Gamma Approximation
% Step 1: Compute Greeks using finite differences
dt_fd = 1 / 250;        % Time step for finite differences
dy_fd = 1e-5;           % Small change in yield for finite differences

% Bond prices for finite differences
Bp = BondPrice(y0 + dy_fd, TTM); % Bond price at y0 + dy_fd
Bm = BondPrice(y0 - dy_fd, TTM); % Bond price at y0 - dy_fd

% Compute Greeks
Theta = (BondPrice(y0, TTM - dt_fd) - B0) / dt_fd; % Theta
Delta = (Bp - Bm) / (2 * dy_fd);                   % Delta
Gamma = (Bp - 2 * B0 + Bm) / (dy_fd)^2;            % Gamma

% Step 2: Simulate yield changes
yT = y0 + mu_dy + vol_dy * randn(M, 1); % Simulated yields

% Step 3: Reprice the bond using Delta-Gamma approximation
BT = B0 + Theta * dt + Delta * (yT - y0) + 0.5 * Gamma * (yT - y0).^2;
PL = BT - B0; % Profit/Loss

% Compute VaR and Expected Shortfall (ES)
VaR_MC_DG = -prctile(PL, (1 - alpha) * 100); % Monte Carlo VaR
ExcessLoss = find(PL < -VaR_MC_DG);          % Losses exceeding VaR
ES_MC_DG = -mean(PL(ExcessLoss));            % Expected Shortfall (ES)

% Display results
check = [VaR_B, VaR_MC_DG]; % Compare analytical and MC VaR
disp('VaR (Analytical vs Monte Carlo):');
disp(check);

% Plot histogram of P&L distribution
h = figure('Color', [1 1 1]);
histogram(PL, 100, 'Normalization', 'pdf'); % Histogram of P&L
hold on;
plot(-VaR_MC_DG, 0, 'r*', 'MarkerSize', 10); % Mark VaR
plot(-ES_MC_DG, 0, 'g*', 'MarkerSize', 10);  % Mark ES
xlabel('Profit/Loss');
ylabel('Density');
title('P&L Distribution with VaR and ES');
legend('P&L Distribution', 'VaR', 'ES');
hold off;
