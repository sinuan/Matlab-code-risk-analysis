clear; close all; clc;

%**************************************************%
%%%%%%%% BOTTOM-UP APPROACH TO PORTFOLIO RISK %%%%%%
%%%%%%%% BOOTSTRAP SIMULATION N-PERIODS     %%%%%%%%%
%**************************************************%

%% ============================
%  Load Data and Define Market
% ============================
filename = 'DOW30_merge.csv'; % Input file name
marketName = 'DOW30'; % Market name
imgDir = 'Images/'; % Directory for saving figures
txtDir = 'Results/'; % Directory for saving results
txtFilename = fullfile(txtDir, 'Bootstrap_Multiperiod.txt'); % Output file for results

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

% Compute Asset Log-Returns for a selected asset
PickAsset = 2; % Index of the selected asset
logRet = log(histPrices(2:end, PickAsset) ./ histPrices(1:end-1, PickAsset)); % Log returns
T = size(logRet, 1); % Number of time periods

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Bootstrap Estimates for N-Days VaR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndays = 10; % VaR horizon (in days)
Nb = 1000; % Number of bootstrap samples
alpha = 0.99; % Confidence level

% Preallocate arrays for bootstrap results
VaRNdays = zeros(Nb, Ndays); % Bootstrap VaR values for each horizon
ESNdays = zeros(Nb, Ndays); % Bootstrap ES values for each horizon

% Repeat Bootstrap simulations
for i = 1:Nb
    % Simulate T x Ndays cumulative returns
    U = randi(T, T, Ndays); % Random indices for bootstrapping
    simLogRetT = cumsum(logRet(U), 2); % Cumulative returns for each horizon
    
    % Compute VaR and ES for each horizon
    for nday = 1:Ndays
        [VaRNdays(i, nday), ESNdays(i, nday)] = get_riskmeasures('NP', simLogRetT(:, nday), alpha);
    end
end

% Bootstrap estimates and confidence intervals
VaRBNdays = mean(VaRNdays); % Mean VaR for each horizon
VaRBNdaysCI = prctile(VaRNdays, [5, 95]); % 90% confidence intervals for VaR

ESBNdays = mean(ESNdays); % Mean ES for each horizon
ESBNdaysCI = prctile(ESNdays, [5, 95]); % 90% confidence intervals for ES

% Create a table for results
Bootstrap = table((1:Ndays)', VaRBNdays', ESBNdays', ...
    'VariableNames', {'Horizon', 'VaR', 'ExpShortfall'});
disp(Bootstrap); % Display table

% Save results to a text file
writetable(Bootstrap, txtFilename, 'Delimiter', 'tab');

%% ============================
%  Plot Bootstrap Results
% ============================
% Plot VaR and ES over the horizon
h1 = figure('Color', [1 1 1]);
plot(1:Ndays, VaRBNdays, 'g*', 'LineWidth', 1.5); % VaR estimates
hold on;
plot(1:Ndays, VaRBNdaysCI, 'g', 'LineWidth', 1.5); % VaR confidence intervals
plot(1:Ndays, ESBNdays, 'r*', 'LineWidth', 1.5); % ES estimates
plot(1:Ndays, ESBNdaysCI, 'r', 'LineWidth', 1.5); % ES confidence intervals
xlabel('Horizon (days)', 'Interpreter', 'latex');
title('Bootstrap Estimates of n-days VaR and ES', 'Interpreter', 'latex');
legend('VaR', 'VaR CI (5%)', 'VaR CI (95%)', 'ES', 'ES CI (5%)', 'ES CI (95%)', ...
    'Location', 'best', 'Interpreter', 'latex');
saveas(h1, fullfile(imgDir, [tickers{PickAsset} '_HS_ndays_VaR_ES.png']));

% Plot distributions of 10-day VaR and ES
h2 = figure('Color', [1 1 1]);
subplot(1, 2, 1);
histogram(VaRNdays(:, end), 30, 'Normalization', 'pdf'); % VaR distribution
hold on;
plot(VaRBNdays(end), 0, 'g*', 'MarkerSize', 10); % VaR estimate
plot(VaRBNdaysCI(:, end), [0, 0], 'r*', 'MarkerSize', 10); % VaR confidence intervals
title('Bootstrap Distribution of 10-day VaR', 'Interpreter', 'latex');
xlabel('VaR', 'Interpreter', 'latex');
xlim([min(min(VaRNdays(:, end), ESNdays(:, end))), max(max(VaRNdays(:, end), ESNdays(:, end)))]);

subplot(1, 2, 2);
histogram(ESNdays(:, end), 30, 'Normalization', 'pdf'); % ES distribution
hold on;
plot(ESBNdays(end), 0, 'g*', 'MarkerSize', 10); % ES estimate
plot(ESBNdaysCI(:, end), [0, 0], 'r*', 'MarkerSize', 10); % ES confidence intervals
title('Bootstrap Distribution of 10-day ES', 'Interpreter', 'latex');
xlabel('ES', 'Interpreter', 'latex');
xlim([min(min(VaRNdays(:, end), ESNdays(:, end))), max(max(VaRNdays(:, end), ESNdays(:, end)))]);

saveas(h2, fullfile(imgDir, [tickers{PickAsset} '_HS_ndays_VaR_ES_distribution.png']));