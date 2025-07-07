    %% ===================================================
%  COMPARING COVARIANCE ESTIMATES
%  ===================================================
%% ============================
%  Load Data and Define Market
% ============================
filename = 'DOW30_merge.csv';
fileindex = 'SPIndex_merge.csv';
marketName = 'DOW30';
imgDir = 'Images/'; % Directory for saving figures (must exist)
txtDir = 'Results/'; % Directory for saving results
txtFilename = fullfile(txtDir, 'TopDown_RiskAnalysis.txt');

% Ensure directories exist
if ~exist(imgDir, 'dir'), mkdir(imgDir); end
if ~exist(txtDir, 'dir'), mkdir(txtDir); end

% Load stock dataset 
dataset = readtable(filename, 'MissingRule', 'omitrow');
colLabels = dataset.Properties.VariableNames;
tickers = colLabels(2:end); % Extract tickers
histPrices = dataset{:, 2:end}; % Historical prices
histDates = dataset{:, 1}; % Historical dates

[NObs, NAsset] = size(histPrices);

% Load index dataset 
dataset = readtable(fileindex, 'MissingRule', 'omitrow');
histSPPrices = dataset{:, 2:end}; % Historical prices
histSPDates = dataset{:, 1}; % Historical dates

% Compute Asset Log-Returns
LogRet = log(histPrices(2:end, :) ./ histPrices(1:end-1, :));

% Estimating the covariance matrix
% Sample covariance
Sigma_Sample = cov(LogRet);
[Sigma_LW, shrinkage] = get_LedoitWolfCov(LogRet(1:j, :));

% Initialize variables for rolling window analysis
k = 1;
CN_sample = []; % Condition number for sample covariance
CN_LW = []; % Condition number for Ledoit-Wolf covariance

for j = 250:250:NObs
    % Sample covariance for rolling window
    Sigma_Sample_r = cov(LogRet(1:j, :));

    % Ledoit-Wolf covariance for rolling window
    [Sigma_LW_r, shrinkage(k)] = get_LedoitWolfCov(LogRet(1:j, :));

    % Compute condition numbers
    CN_sample(k) = max(eig(Sigma_Sample_r)) / min(eig(Sigma_Sample_r));
    CN_LW(k) = max(eig(Sigma_LW_r)) / min(eig(Sigma_LW_r));

    k = k + 1;
end

% Display condition numbers for comparison
disp([(250:250:NObs)', CN_sample', CN_LW']);

% Compute S&P Index Log-Returns
SPLogRet = log(histSPPrices(2:end, :) ./ histSPPrices(1:end-1, :));

% Find common dates between stock and index data
[common_dates, iSP, iStock] = intersect(histSPDates(2:end), histDates(2:end));

StockRet = LogRet(iStock, :);
SPRet = SPLogRet(iSP, :);

% Estimate market model (CAPM) for each stock
beta_0 = zeros(NAsset, 1); % Intercept
beta_1 = zeros(NAsset, 1); % Slope coefficient
res_var = zeros(NAsset, 1); % Residual variance
R2 = zeros(NAsset, 1); % R-squared

for i = 1:NAsset
    % Linear regression model: r_i = a_i + b_i * rm + eps_i
    model = fitlm(StockRet(:, i), SPRet);
    beta = model.Coefficients.Estimate;
    beta_0(i, 1) = beta(1); % Intercept
    beta_1(i, 1) = beta(2); % Slope coefficient
    res_var(i, 1) = model.MSE; % Residual variance
    R2(i, 1) = model.Rsquared.Ordinary; % R-squared
end

% Estimate covariance matrix using the market model
Sigma_1F = beta_1 * beta_1' * var(SPRet) + diag(res_var);

% Download Fama-French 3 Factors data
url = 'https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_daily_CSV.zip';
outputFile = 'F-F_Research_Data_Factors_daily_CSV.zip';
websave(outputFile, url); % Download the ZIP file
unzip(outputFile); % Extract the contents

% Import Fama-French factors data
FF_factors_daily = readtable('F-F_Research_Data_Factors_daily.csv');
disp(head(FF_factors_daily)); % Display the first few rows

% Convert dates to numeric format (YYYYMMDD)
Dates_datetime = datetime(histDates, 'InputFormat', 'MM/dd/yyyy');
Dates_numeric = year(Dates_datetime) * 10000 + month(Dates_datetime) * 100 + day(Dates_datetime);

% Find common dates between stock and Fama-French data
[common_dates, iFactors, iStock] = intersect(FF_factors_daily.Var1, Dates_numeric(2:end));

StockRet = LogRet(iStock, :);
FactRet = log(1 + FF_factors_daily{iFactors, 2:4} / 100); % Factor returns
Rf = log(1 + FF_factors_daily{iFactors, 5} / 100); % Risk-free rate

% Estimate 3-factor model for each stock
beta3F_0 = zeros(NAsset, 1); % Intercept
beta3F_1 = zeros(NAsset, 3); % Factor loadings (Mkt, SMB, HML)
res3F_var = zeros(NAsset, 1); % Residual variance
R23F = zeros(NAsset, 1); % R-squared

for i = 1:NAsset
    % Prepare data for regression
    T = array2table([StockRet(:, i), FactRet], ...
        'VariableNames', {'StockRet', 'Mkt_RF', 'SMB', 'HML'});
    
    % Linear regression model: r_i = a_i + b1_i * Mkt + b2_i * SMB + b3_i * HML + eps_i
    model = fitlm(T, 'StockRet ~ Mkt_RF + SMB + HML');
    beta = model.Coefficients.Estimate;
    beta3F_0(i, 1) = beta(1); % Intercept
    beta3F_1(i, :) = beta(2:4); % Factor loadings
    res3F_var(i, 1) = model.MSE; % Residual variance
    R23F(i, 1) = model.Rsquared.Ordinary; % R-squared
end

% Estimate covariance matrix using the 3-factor model
Sigma_3F = beta3F_1 * cov(FactRet(:, 1:3)) * beta3F_1' + diag(res3F_var);

% Compute condition numbers for covariance matrices
CN1F = max(eig(Sigma_1F)) / min(eig(Sigma_1F));
CN3F = max(eig(Sigma_3F)) / min(eig(Sigma_3F));
CN3Sample = max(eig(Sigma_Sample)) / min(eig(Sigma_Sample));
CN3LW = max(eig(Sigma_LW)) / min(eig(Sigma_LW));

% Define a function for weighted covariance matrix
MySigma = @(a) a * Sigma_3F + (1 - a) * Sigma_Sample;
CNMy = @(a) max(eig(MySigma(a))) / min(eig(MySigma(a)));

% Display results
disp('Condition Numbers:');
disp(['1-Factor Model: ', num2str(CN1F)]);
disp(['3-Factor Model: ', num2str(CN3F)]);
disp(['Sample Covariance: ', num2str(CN3Sample)]);
disp(['Ledoit-Wolf Covariance: ', num2str(CN3LW)]);


% %Global Min Var Portfolio
% wg_sample  = inv(Sigma_Sample)*ones(NAsset,1)/sum(inv(Sigma_Sample)*ones(NAsset,1))
% wg_LW  = inv(Sigma_LW)*ones(NAsset,1)/sum(inv(Sigma_LW)*ones(NAsset,1))
% wg_1F  = inv(Sigma_1F)*ones(NAsset,1)/sum(inv(Sigma_1F)*ones(NAsset,1))
% wg_3F  = inv(Sigma_3F)*ones(NAsset,1)/sum(inv(Sigma_3F)*ones(NAsset,1))
% 
% bar([wg_sample, wg_LW, wg_1F, wg_3F])