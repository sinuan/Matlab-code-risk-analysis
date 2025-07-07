clear; close all; clc;

%**************************************************%
%%%%%%%% BOTTOM-UP APPROACH TO PORTFOLIO RISK %%%%%%
%%%%%%%% PARAMETRIC APPROACHES %%%%%%%%%%%%%%%%%%%%%
%**************************************************%

%% ============================
%  Load Data and Define Market
% ============================
filename = 'DOW30_merge.csv';
marketName = 'DOW30';
imgDir = 'Images/'; % Directory for saving figures (must exist)
txtDir = 'Results/'; % Directory for saving results
txtFilename = fullfile(txtDir, 'TopDown_RiskAnalysis.txt');

% Ensure directories exist
if ~exist(imgDir, 'dir'), mkdir(imgDir); end
if ~exist(txtDir, 'dir'), mkdir(txtDir); end

% Load dataset
dataset = readtable(filename, 'MissingRule', 'omitrow');
colLabels = dataset.Properties.VariableNames;
tickers = colLabels(2:end); % Extract tickers
histPrices = dataset{:, 2:end}; % Historical prices
histDates = dataset{:, 1}; % Historical dates

[NObs, NAsset] = size(histPrices);

% Compute Asset Log-Returns
logRet = log(histPrices(2:end, :) ./ histPrices(1:end-1, :));

% Assign Portfolio Composition: Equally Weighted
w_eq = ones(NAsset, 1) / NAsset;

% Assign Confidence Level
alpha = 0.99;

% Estimate Mean Vector and Covariance Matrix
MeanV = mean(logRet)';
Sigma = cov(logRet);

% Compute portfolio variance
sg2p = w_eq' * Sigma * w_eq;

% Compute portfolio VaR (assume zero mean)
z = norminv(1 - alpha, 0, 1);
VaR_w = -z * sqrt(sg2p);

% Compute Marginal VaR
MVaR_w = -z * Sigma * w_eq / sqrt(sg2p);

% Compute Component VaR
CVaR_w = w_eq .* MVaR_w;
chk = [sum(CVaR_w), VaR_w];

% Compute Component VaR in percentage
CVaR_w_p = CVaR_w / sum(CVaR_w);

% Create table for Component VaR results
T_Component_g = table(w_eq, MVaR_w, CVaR_w, CVaR_w_p);
T_Component_g.Properties.VariableNames = {'Weights', 'MVaR', 'CVaR', 'CVaR%'};
T_Component_g.Properties.RowNames = tickers;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Handle functions %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sg2p = @(x, Sigma) x' * Sigma * x;
MVaR = @(x, Sigma) Sigma * x / sqrt(sg2p(x, Sigma));
CVaR = @(x, Sigma) x .* MVaR(x, Sigma);
Conv = @(x, Sigma) sqrt(x.^2 .* diag(Sigma) - (x .* CVaR(x, Sigma)).^2);

% Build the minimum variance portfolio
x0 = ones(NAsset, 1) / NAsset;
w_mv = fmincon(@(x) sg2p(x, Sigma), x0, [], [], ones(1, NAsset), 1, zeros(NAsset, 1), ones(NAsset, 1));
sg2mv = sg2p(w_mv, Sigma);
MVaR_mv = MVaR(w_mv, Sigma);
CVaR_mv = CVaR(w_mv, Sigma);
CVaR_mv_p = CVaR_mv / sum(CVaR_mv);
Conv_mv = Conv(w_mv, Sigma);

% Build the risk-parity portfolio
x0 = w_eq;
w_rp = fmincon(@(x) std(CVaR(x, Sigma)), x0, [], [], ones(1, NAsset), 1, zeros(NAsset, 1), ones(NAsset, 1));
sg2rp = sg2p(w_rp, Sigma);
MVaR_rp = MVaR(w_rp, Sigma);
CVaR_rp = CVaR(w_rp, Sigma);
CVaR_rp_p = CVaR_rp / sum(CVaR_rp);
Conv_rp = Conv(w_rp, Sigma);

% Build the conviction-parity portfolio
x0 = w_rp;
w_cp = fmincon(@(x) std(Conv(x, Sigma)), x0, [], [], ones(1, NAsset), 1, zeros(NAsset, 1), ones(NAsset, 1));
sg2cp = sg2p(w_cp, Sigma);
MVaR_cp = MVaR(w_cp, Sigma);
CVaR_cp = CVaR(w_cp, Sigma);
CVaR_cp_p = CVaR_cp / sum(CVaR_cp);
Conv_cp = Conv(w_cp, Sigma);

% Prepare tables for the report
T_w = table(w_eq, w_mv, w_rp, w_cp);
T_w.Properties.VariableNames = {'EQ', 'GMV', 'RP', 'Conv'};
T_w.Properties.RowNames = tickers;

T_MRisk = table(MVaR(w_eq, Sigma), MVaR_mv, MVaR_rp, MVaR_cp);
T_MRisk.Properties.VariableNames = {'EQ', 'GMV', 'RP', 'Conv'};
T_MRisk.Properties.RowNames = tickers;

T_CRisk = table(CVaR_w_p, CVaR_mv_p, CVaR_rp_p, CVaR_cp_p);
T_CRisk.Properties.VariableNames = {'EQ', 'GMV', 'RP', 'Conv'};
T_CRisk.Properties.RowNames = tickers;

T_Conv = table(Conv(w_eq, Sigma), Conv_mv, Conv_rp, Conv_cp);
T_Conv.Properties.VariableNames = {'EQ', 'GMV', 'RP', 'Conv'};
T_Conv.Properties.RowNames = tickers;

% Plot portfolio weights
h = figure('Color', [1 1 1]);
StockNames = categorical(tickers);
bar(StockNames, [w_eq, w_mv, w_rp, w_cp]);
legend('E.W.', 'M.V.', 'R.P.', 'C.P.', 'interpreter', 'latex', 'location', 'northwest');
title('Portfolio weights', 'interpreter', 'latex');
print(h, '-dpng', fullfile(imgDir, 'Fig_HotSpot_EW_MV_RP'));

% Plot Marginal VaR
h = figure('Color', [1 1 1]);
bar(StockNames, [MVaR(w_eq, Sigma), MVaR_mv, MVaR_rp, MVaR_cp]);
legend('E.W.', 'M.V.', 'R.P.', 'interpreter', 'latex', 'location', 'northwest');
title('Marginal VaR', 'interpreter', 'latex');
print(h, '-dpng', fullfile(imgDir, 'Fig_HotSpot_MVaR_EW_MV_RP'));

% Plot Component VaR%
h = figure('Color', [1 1 1]);
bar(StockNames, [CVaR_w_p, CVaR_mv_p, CVaR_rp_p, CVaR_cp_p]);
legend('E.W.', 'M.V.', 'R.P.', 'interpreter', 'latex', 'location', 'northwest');
title('Component VaR\%', 'interpreter', 'latex');
print(h, '-dpng', fullfile(imgDir, 'Fig_HotSpot_CVaR_EW_MV_RP'));

% Plot Conviction
h = figure('Color', [1 1 1]);
bar(StockNames, [Conv(w_eq, Sigma), Conv_mv, Conv_rp, Conv_cp]);
legend('E.W.', 'M.V.', 'R.P.', 'C.P.', 'interpreter', 'latex', 'location', 'northwest');
title('Conviction', 'interpreter', 'latex');
print(h, '-dpng', fullfile(imgDir, 'Fig_HotSpot_Conv_EW_MV_RP'));

% Assess the strategies in terms of incremental VaR
IVaR_mb = MVaR_w' * (w_mv - w_eq);
IVaR_rp = MVaR_rp' * (w_rp - w_eq);
IVaR_conv = MVaR_cp' * (w_rp - w_eq);

% Plot Incremental VaR
Strategy = categorical({'Min. Var.', 'Risk Parity', 'Conv. Parity'});
h = figure('Color', [1 1 1]);
bar(Strategy, [IVaR_mb, IVaR_rp, IVaR_conv]);
title('Incremental VaR', 'interpreter', 'latex');
print(h, '-dpng', fullfile(imgDir, 'Fig_IVaR_EW_MV_RP'));