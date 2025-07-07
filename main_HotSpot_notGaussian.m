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

% Compute Asset Log-Returns
logRet = log(histPrices(2:end, :) ./ histPrices(1:end-1, :)); % Log returns
T = size(logRet, 1); % Number of time periods

% Assign Confidence level
alpha = 0.99; % Confidence level for VaR and ES

% Fix the horizon
horizon = 1; % Time horizon for risk measures

% Assign weights (equal-weighted portfolio)
w = ones(NAsset, 1) / NAsset; % Portfolio weights

% Compute portfolio returns
logRetPortfolio = logRet * w; % Portfolio log returns

% Estimate the mean vector using the sample mean
mu = mean(logRet, 1)' * horizon; % Mean returns scaled by horizon

% Estimate the covariance matrix using the sample covariance matrix
Sigma = cov(logRet) * horizon; % Covariance matrix scaled by horizon

% Fix the number of simulations
M = 1000000; % Number of Monte Carlo simulations

% Simulate asset returns
simReturns = mvnrnd(mu, Sigma, M); % Simulated returns

% Simulate portfolio returns
simPortfolioReturns = simReturns * w; % Simulated portfolio returns

% Compute VaR and ES of the portfolio
[VaR, ES] = get_riskmeasures('NP', simPortfolioReturns, alpha); % Non-parametric VaR and ES

% Computing Marginal VaR (MVaR)
eps = 0.01 * VaR; % Error margin for VaR interval
lowerBound = -VaR - eps; % Lower bound for VaR interval
upperBound = -VaR + eps; % Upper bound for VaR interval
pos = find((simPortfolioReturns >= lowerBound) & (simPortfolioReturns <= upperBound)); % Positions within VaR interval

% Extract simulations that satisfy the condition
condReturns = simReturns(pos, :); % Conditional returns

% Compute MVaR: take the conditional mean
MVaR = -mean(condReturns)'; % Marginal VaR

% Compute Component VaR (CVaR)
CVaR = w .* MVaR; % Component VaR
CVaRp = CVaR / sum(CVaR); % Normalized Component VaR
checkVaR = [sum(CVaR), VaR]; % Check consistency

% Compute Marginal Expected Shortfall (MES)
posES = find(simPortfolioReturns <= -VaR); % Positions for ES
condReturnsES = simReturns(posES, :); % Conditional returns for ES

% Compute MES: take the conditional mean
MES = -mean(condReturnsES)'; % Marginal Expected Shortfall

% Compute Component ES (CES)
CES = w .* MES; % Component Expected Shortfall
CESp = CES / sum(CES); % Normalized Component ES
checkES = [sum(CES), ES]; % Check consistency

% Plot CVaR and CES contributions
threshold = 0.05; % Threshold for identifying hot spots
X = categorical(tickers); % Ticker names for plotting

h1 = figure('Color', [1 1 1]);
subplot(2, 1, 1);
bar(X, CVaRp); % Bar plot of CVaR contributions
hold on;
bar(X, CVaRp .* (CVaRp > threshold), 'r'); % Highlight contributions above threshold
yline(threshold, 'r', 'LineWidth', 1.5); % Threshold line
xlabel('Stock', 'Interpreter', 'latex');
ylabel('CVaR\%', 'Interpreter', 'latex');
title('Identifying Hot Spots (CVaR)', 'Interpreter', 'latex');

subplot(2, 1, 2);
bar(X, CESp); % Bar plot of CES contributions
hold on;
bar(X, CESp .* (CESp > threshold), 'r'); % Highlight contributions above threshold
yline(threshold, 'r', 'LineWidth', 1.5); % Threshold line
xlabel('Stock', 'Interpreter', 'latex');
ylabel('CES\%', 'Interpreter', 'latex');
title('Identifying Hot Spots (CES)', 'Interpreter', 'latex');

% Save figure
saveas(h1, fullfile(imgDir, 'HotSpots_CVaR_CES.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Historical risk decomposition: HORIZON 1 PERIOD
[HS_VaR, HS_ES] = get_riskmeasures('NP', logRetPortfolio, alpha); % Historical VaR and ES

% Computing Historical MVaR
eps = 0.1 * HS_VaR; % Error margin for VaR interval
lowerBound = -HS_VaR - eps; % Lower bound for VaR interval
upperBound = -HS_VaR + eps; % Upper bound for VaR interval
pos = find((logRetPortfolio >= lowerBound) & (logRetPortfolio <= upperBound)); % Positions within VaR interval

% Extract simulations that satisfy the condition
condReturnsHS = logRet(pos, :); % Conditional historical returns

% Compute Historical MVaR: take the conditional mean
HS_MVaR = -mean(condReturnsHS)'; % Historical Marginal VaR

% Compute Historical CVaR
HS_CVaR = w .* HS_MVaR; % Historical Component VaR
HS_CVaRp = HS_CVaR / sum(HS_CVaR); % Normalized Historical Component VaR
checkVaR = [sum(HS_CVaR), HS_VaR]; % Check consistency

% Compute Historical MES
posES = find(logRetPortfolio <= -HS_VaR); % Positions for ES
HS_MES = -mean(logRet(posES, :))'; % Historical Marginal Expected Shortfall

% Compute Historical CES
HS_CES = w .* HS_MES; % Historical Component Expected Shortfall
HS_CESp = HS_CES / sum(HS_CES); % Normalized Historical Component ES
checkES = [sum(HS_CES), HS_ES]; % Check consistency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare Gaussian and Non-Parametric Risk Contributions
z = norminv(1 - alpha); % Z-score for Gaussian VaR
volPortfolio = sqrt(w' * Sigma * w); % Portfolio volatility
MVaRg = -(mu + z * Sigma * w / volPortfolio); % Gaussian Marginal VaR
CVaRg = w .* MVaRg; % Gaussian Component VaR
CVaRpg = CVaRg / sum(CVaRg); % Normalized Gaussian Component VaR

% Plot CVaR contributions
h2 = figure('Color', [1 1 1]);
bar(X, [CVaRp, CVaRpg, HS_CVaRp]); % Bar plot of CVaR contributions
yline(threshold, 'r', 'LineWidth', 1.5); % Threshold line
xlabel('Stock', 'Interpreter', 'latex');
ylabel('CVaR\%', 'Interpreter', 'latex');
legend('Monte Carlo', 'Gaussian', 'Historical', 'Interpreter', 'latex');
title('Identifying Hot Spots (CVaR Contributions)', 'Interpreter', 'latex');

% Save figure
saveas(h2, fullfile(imgDir, 'CVaR_Contributions_Comparison.png'));