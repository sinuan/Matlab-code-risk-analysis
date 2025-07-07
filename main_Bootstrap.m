clear; close all; clc;

%**************************************************%
%%%%%%%% BOTTOM-UP APPROACH TO PORTFOLIO RISK %%%%%%
%%%%%%%% BOOTSTRAP SIMULATION %%%%%%%%%%%%%%%%%%%%%
%**************************************************%

%% ============================
%  Load Data and Define Market
% ============================
filename = 'DOW30_merge.csv'; % Input file name
marketName = 'DOW30'; % Market name
imgDir = 'Images/'; % Directory for saving figures
txtDir = 'Results/'; % Directory for saving results
txtFilename = fullfile(txtDir, 'TopDown_RiskAnalysis.txt'); % Output file for results

% Ensure directories exist
if ~exist(imgDir, 'dir'), mkdir(imgDir); end
if ~exist(txtDir, 'dir'), mkdir(txtDir); end

% Load dataset
dataset = readtable(filename, 'MissingRule', 'omitrow'); % Read data
colLabels = dataset.Properties.VariableNames; % Column labels
tickers = colLabels(2:end); % Extract tickers (asset names)
histPrices = dataset{:, 2:end}; % Historical prices
histDates = dataset{:, 1}; % Historical dates

[NObs, NAsset] = size(histPrices); % Number of observations and assets

% Compute Asset Log-Returns
logRet = log(histPrices(2:end, :) ./ histPrices(1:end-1, :)); % Log returns
T = size(logRet, 1); % Number of time periods

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Bootstrap Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
alpha = 0.90; % Confidence level
Nb = 500; % Number of bootstrap samples
w = ones(NAsset, 1) / NAsset; % Equal-weighted portfolio

% Preallocate arrays for bootstrap results
VaRb = zeros(Nb, 1); % Bootstrap VaR values
ESb = zeros(Nb, 1); % Bootstrap ES values

% Run the bootstrap simulation
for i = 1:Nb
    % Randomly sample returns with replacement
    U = randi(T, T, 1); % Random indices
    simRet = logRet(U, :); % Bootstrap sample of returns
    
    % Portfolio returns for the bootstrap sample
    simRetPortfolio = simRet * w;
    
    % Compute VaR and ES for the bootstrap sample
    [VaRb(i), ESb(i)] = get_riskmeasures('NP', simRetPortfolio, alpha);
end

% Bootstrap VaR Estimate & Confidence Intervals
VaRboot = mean(VaRb); % Bootstrap VaR estimate
VaRbootCI = prctile(VaRb, [5, 95]); % 90% confidence interval for VaR

% Bootstrap ES Estimate & Confidence Intervals
ESboot = mean(ESb); % Bootstrap ES estimate
ESbootCI = prctile(ESb, [5, 95]); % 90% confidence interval for ES

% Display results
VaRBootstrap = [VaRboot; VaRbootCI(:)]; % VaR results
ESBootstrap = [ESboot; ESbootCI(:)]; % ES results
Quantity = {'Bootstrap Estimate'; 'Lower CI'; 'Upper CI'}; % Labels
Synthesis = table(Quantity, VaRBootstrap, ESBootstrap); % Create table
disp(Synthesis); % Display table

% Save results to a text file
writetable(Synthesis, txtFilename, 'Delimiter', 'tab');

%% ============================
%  Plot Bootstrap Distributions
% ============================
% Plot VaR distribution
h1 = figure('Color', [1 1 1]);
histogram(VaRb, 30, 'Normalization', 'pdf'); % Histogram of VaR
hold on;
plot(VaRboot, 0, 'b*', 'MarkerSize', 10); % Bootstrap VaR estimate
plot(VaRbootCI, [0, 0], 'r*', 'MarkerSize', 10); % Confidence intervals
title('Bootstrap Distribution of VaR', 'Interpreter', 'latex');
xlabel('VaR');
ylabel('Density');
saveas(h1, fullfile(imgDir, 'Port_HS_VaRDist.png')); % Save figure

% Plot ES distribution
h2 = figure('Color', [1 1 1]);
histogram(ESb, 30, 'Normalization', 'pdf'); % Histogram of ES
hold on;
plot(ESboot, 0, 'b*', 'MarkerSize', 10); % Bootstrap ES estimate
plot(ESbootCI, [0, 0], 'r*', 'MarkerSize', 10); % Confidence intervals
title('Bootstrap Distribution of ES', 'Interpreter', 'latex');
xlabel('ES');
ylabel('Density');
saveas(h2, fullfile(imgDir, 'Port_HS_ESDist.png')); % Save figure