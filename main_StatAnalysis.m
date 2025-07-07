clear; close all; clc;

%% ============================
%  Load Data and Define Market
% ============================
filename = 'SPIndex_merge.csv'; 
myMarket = 'SP500';  
img_directory = 'Images/'; % Directory for saving figures (must exist)
txt_directory  = "Results/"; % Directory for saving results 

%Check if imgdirectory and txt_directory  exist, otherwise create them
get_check_directory(img_directory)
get_check_directory(txt_directory)

% Load dataset
dataset = readtable(filename, 'MissingRule', 'omitrow');
ColLabels = dataset.Properties.VariableNames;
Tickers = ColLabels(2:end); % Extract tickers
HistPrices = dataset{:, 2:end}; % Historical prices
HistDates = dataset{:, 1}; % Historical dates

[NObs, NAsset] = size(HistPrices);

% Select Asset
PickAsset = randi(NAsset); % to make random choice 
symbol = Tickers{PickAsset};
disp(['Selected Asset: ', symbol]);


% Compute Log Returns
LogRet = log(HistPrices(2:end, PickAsset) ./ HistPrices(1:end-1, PickAsset)); 
Dates = HistDates(2:end);

%% ============================
%  Plot Price and LogRet series
% ============================
h = figure('Color', [1 1 1]);
subplot(1,2,1)
plot(HistDates, HistPrices(:, PickAsset))
xlabel('Time', 'interpreter','latex')
ylabel('SP500 Index', 'interpreter','latex')

subplot(1,2,2)
plot(Dates, LogRet)
xlabel('Time', 'interpreter','latex')
ylabel('SP500 returns', 'interpreter','latex')
print(h, [img_directory, 'Ticker_', symbol, '_price_and_returns'], '-dpng');

%% ============================
%  Compute Basic Statistics
% ============================
nobs = length(LogRet);
mean_r = mean(LogRet);
var_r = var(LogRet);
sd_r = sqrt(var_r);
min_r = min(LogRet);
max_r = max(LogRet);
sk_r = skewness(LogRet);
k_r = kurtosis(LogRet);

% Create Summary Table
Synthesis = table(nobs, mean_r, var_r, sd_r, min_r, max_r, sk_r, k_r);
Synthesis.Properties.VariableNames = {'N_Obs', 'Mean', 'Var', 'St_Dev', 'Min', 'Max', 'Skew', 'Kurt'};
Synthesis.Properties.RowNames = {symbol};



%% ============================
%  Plot Histogram vs Gaussian PDF
% ============================
IQR = prctile(LogRet, 75) - prctile(LogRet, 25);
binsize = 2*IQR/nobs^(1/3);
nbins = round((max(LogRet) - min(LogRet)) / binsize );

h = figure('Color', [1 1 1]);
subplot(2,2,[1:2])
histogram(LogRet, nbins, 'Normalization', 'pdf')
hold on
fplot(@(x) normpdf(x, mean_r, sd_r), [min_r max_r])
ksdensity(LogRet)
xlabel('Log Returns', 'interpreter', 'latex')
title([symbol ': Empirical vs Gaussian Density'], 'interpreter', 'latex')
xlim([mean_r-6*sd_r mean_r+6*sd_r])
legend('Empirical', 'Normal', 'Kernel Density', 'interpreter', 'latex')

subplot(2,2,3)
histogram(LogRet, nbins, 'Normalization','pdf')
hold on
fplot(@(x) normpdf(x, mean_r, sd_r), [prctile(LogRet,1) prctile(LogRet,99)])
ksdensity(LogRet)
xlim([prctile(LogRet, 1) prctile(LogRet, 50)])
xlabel('Log Returns', 'interpreter', 'latex')
title([symbol ': Left Tail'], 'interpreter', 'latex')
legend('Empirical', 'Normal', 'Kernel Density', 'interpreter', 'latex')

subplot(2,2,4)
histogram(LogRet, nbins, 'Normalization','pdf')
hold on
fplot(@(x) normpdf(x, mean_r, sd_r), [prctile(LogRet,60) prctile(LogRet,99)])
ksdensity(LogRet)
xlim([prctile(LogRet,50) prctile(LogRet,99)])
xlabel('Log Returns', 'interpreter', 'latex')
legend('Empirical', 'Normal', 'Kernel Density', 'interpreter', 'latex')
title([symbol ': Right Tail'], 'interpreter', 'latex')

print(h, [img_directory, 'Ticker_', symbol, '_histogram_vs_gaussian'], '-dpng');

%% ============================
%  Theoretical vs Empirical Frequencies
% ============================
Z = (LogRet - mean(LogRet)) / std(LogRet); % Standardized Returns

% Define intervals
Range = [0 0.25;
         0.25 0.5;
         0.5 1;
         1 1.5;
         1.5 2;
         2 3;
         3 100];

for j = 1:size(Range, 1)  
    TheorPr(j,1) = 2 * (normcdf(Range(j,2)) - normcdf(Range(j,1)));
    EmpFr(j,1) = sum((abs(Z) < Range(j,2)) .* (abs(Z) > Range(j,1))) / nobs;
end

% Display Table
FreqTable = table(Range(:,1), Range(:,2), round(EmpFr*100,4), round(TheorPr*100,4));
FreqTable.Properties.VariableNames = {'Lower Bound', 'Upper Bound', 'Empirical Freq (%)', 'Theoretical Freq (%)'};

Xrange=[1:nobs];
h=figure('Color',[1 1 1]);
X = categorical(Range);
bar(X(:,1), [EmpFr*100 TheorPr*100]) 
ylabel('Frequencies (\%)', 'interpreter','latex')
legend('Empirical Frequency','Theoretical Freq.','interpreter','latex')
% Set the centered title for the entire figure
sgtitle([symbol ': Theoretical vs Empirical Frequencies'], 'interpreter', 'latex')

print(h,[img_directory, 'Ticker_', symbol '_emp_vs_theor'],'-dpng')

%% ============================
%  Q-Q Plot
% ============================
Z = (LogRet - mean_r) / sd_r;
thquantile = icdf('normal', [0.5:nobs-0.5] / nobs, 0, 1)';
emquantile = sort(Z);

h = figure('Color', [1 1 1]);
plot(thquantile, [thquantile, emquantile])
axis equal;
grid on;
title([symbol ': QQPlot'], 'interpreter', 'latex')
xlabel('Theoretical Gaussian Quantiles', 'interpreter', 'latex');
ylabel('Sample Quantiles', 'interpreter', 'latex');
print(h, [img_directory, 'Ticker_', symbol, '_qq_plot'], '-dpng');

%% ============================
%  Perform Jarque-Bera Normality Test
% ============================
JB = (nobs/6) * sk_r^2 + (nobs/24) * (k_r - 3)^2;
CV = icdf('chi', 0.95, 2);
pvalue = 1 - cdf('chi', JB, 2);

%Use Matlab built-in function
[H0, pvalue, jbstat, critval] = jbtest(Z);

% Create Summary Table
JBSynthesis = table(nobs, sk_r, k_r, JB, CV, pvalue);
JBSynthesis.Properties.VariableNames = {'N. Obs', 'Skewness', 'Kurtosis', ...
                    'JB', 'CV', 'PValue'};
JBSynthesis.Properties.RowNames = {symbol};



%% ============================
%  Clustering of Extreme Losses/Gains
% ============================
ExtremeG = zeros(nobs,1);
ExtremeL = zeros(nobs,1);
ExtremeG(Z > 3) = Z(Z > 3);
ExtremeL(Z < -3) = Z(Z < -3);

h = figure('Color', [1, 1, 1]);
subplot(1,2,1)
plot(datenum(Dates), Z,'.')
hold on
plot(datenum(Dates(Z>3)), Z(Z>3),'r.')
plot(datenum(Dates(Z<-3)), Z(Z<-3),'r.')
yline(3)
yline(-3)
xlabel('Time', 'interpreter','latex')
ylabel('Standardized returns', 'interpreter','latex')
xlim([datenum(Dates(2)) datenum(Dates(end))])
dateaxis('x', 12, Dates(2)) 

subplot(1,2,2)
plot(Dates, ExtremeG, 'g');
hold on
plot(Dates, abs(ExtremeL), 'r');
xlabel('Time', 'Interpreter', 'latex');
ylabel('Gains and Losses larger than 3 st. dev.', 'Interpreter', 'latex');
legend('Gains','Losses', 'Interpreter', 'latex')
% Set the centered title for the entire figure
sgtitle([symbol ': Clustering of Extreme Losses/Gains'], 'Interpreter', 'latex');
grid on;

print(h, [img_directory, 'Ticker_', symbol '_extreme_losses_gains'], '-dpng');

%% ============================
%  Test for Autocorrelation
% ============================
maxlags = 20;
[acfValues,  ~, ~] = autocorr(Z, 'NumLags', maxlags );
[acfSqValues,  ~, ~] = autocorr(Z.^2, 'NumLags', maxlags );
[acfAbsValue,  ~, ~] = autocorr(abs(Z), 'NumLags', maxlags );

% Create Summary Table
ACFSynthesis = table([1:maxlags]', round(acfValues(2:end),4), ...
    round(acfSqValues(2:end),4), round(acfAbsValue(2:end),4));
ACFSynthesis.Properties.VariableNames = {'Lag', 'ACF Ret', 'ACF SQ. Ret', 'ACF Abs. Ret'};


h = figure('Color',[1 1 1]);
autocorr(Z, 'NumLags', maxlags);
xlim([0.5, maxlags+0.5])
xlabel('Lag','interpreter','latex')
ylabel('Sample autocorrelation','interpreter','latex')
title([symbol ': Autocorrelation of Log-Returns'], 'interpreter', 'latex');
print(h, [img_directory, 'Ticker_', symbol, '_acf'], '-dpng');

h = figure('Color',[1 1 1]);
autocorr(Z.^2, 'NumLags', maxlags);
xlim([0.5, maxlags+0.5])
xlabel('Lag','interpreter','latex')
ylabel('Sample autocorrelation','interpreter','latex')
title([symbol ': Autocorrelation of Squared Log-Returns'], 'interpreter', 'latex');
print(h, [img_directory, 'Ticker_', symbol, '_acf2'], '-dpng');

h = figure('Color',[1 1 1]);
autocorr(abs(Z), 'NumLags', maxlags);
xlim([0.5, maxlags+0.5])
xlabel('Lag','interpreter','latex')
ylabel('Sample autocorrelation','interpreter','latex')
title([symbol ': Autocorrelation of Absolute Log-Returns'], 'interpreter', 'latex');
print(h, [img_directory, 'Ticker_', symbol, '_acfabs'], '-dpng');

%% ============================
%  Ljung-Box and Box-Pierce Test
% ============================

% Compute Test Statistics
BoxPierce = sum(acfValues(2:end).^2) * nobs;
LjungBox = sum(acfValues(2:end).^2 ./ (nobs - (1:maxlags)')) * nobs * (nobs + 2);
pvalBP = 1 - cdf('chi', BoxPierce, maxlags)
pvalLB = 1 - cdf('chi', LjungBox, maxlags)

BoxPierce2 = sum(acfSqValues(2:end).^2) * nobs;
LjungBox2 = sum(acfSqValues(2:end).^2 ./ (nobs - (1:maxlags)')) * nobs * (nobs + 2);
pvalBP2 = 1 - cdf('chi', BoxPierce2, maxlags)
pvalLB2 = 1 - cdf('chi', LjungBox2, maxlags)

BoxPierceAbs = sum(acfAbsValue(2:end).^2) * nobs;
LjungBoxAbs = sum(acfAbsValue(2:end).^2 ./ (nobs - (1:maxlags)')) * nobs * (nobs + 2);
pvalBPA = 1 - cdf('chi', BoxPierceAbs, maxlags)
pvalLBA = 1 - cdf('chi', LjungBoxAbs, maxlags)

CV = icdf('chisquare', 0.95, maxlags);

% Create Summary Table
PB1Synthesis = table(nobs, maxlags, BoxPierce, CV, pvalBP);
PB1Synthesis.Properties.VariableNames = {'N. Obs', 'Lags', 'Box-Pierce', 'CV', 'PValue'};
PB1Synthesis.Properties.RowNames = {symbol};
LB1Synthesis = table(nobs, maxlags, LjungBox, CV, pvalLB);
LB1Synthesis.Properties.VariableNames = {'N. Obs', 'Lags', 'Ljung-Box', 'CV', 'PValue'};
LB1Synthesis.Properties.RowNames = {symbol};


PB2Synthesis = table(nobs, maxlags, BoxPierce2, CV, pvalBP2);
PB2Synthesis.Properties.VariableNames = {'N. Obs', 'Lags', 'Box-Pierce', 'CV', 'PValue'};
PB2Synthesis.Properties.RowNames = {symbol};
LB2Synthesis = table(nobs, maxlags, LjungBox2, CV, pvalLB2);
LB2Synthesis.Properties.VariableNames = {'N. Obs', 'Lags', 'Ljung-Box', 'CV', 'PValue'};
LB2Synthesis.Properties.RowNames = {symbol};

PB3Synthesis = table(nobs, maxlags, BoxPierceAbs, CV, pvalBPA);
PB3Synthesis.Properties.VariableNames = {'N. Obs', 'Lags', 'Box-Pierce', 'CV', 'PValue'};
PB3Synthesis.Properties.RowNames = {symbol};

LB3Synthesis = table(nobs, maxlags, LjungBoxAbs, CV, pvalLBA);
LB3Synthesis.Properties.VariableNames = {'N. Obs', 'Lags', 'Ljung-Box', 'CV', 'PValue'};
LB3Synthesis.Properties.RowNames = {symbol};

% ============================
% Get the dates of the 5 largest losses and gains
% ============================

% Sort indices to find the top 5 gains and losses
[~, sortedIdx] = sort(LogRet); % Sort log returns in ascending order

% Get indices of 5 largest gains (last 5 elements)
top_gains_idx = sortedIdx(end-4:end); 

% Get indices of 5 largest losses (first 5 elements)
top_losses_idx = sortedIdx(1:5); 

% Extract corresponding dates and log return values
gains_returns = num2str(LogRet(top_gains_idx));
gains_dates = datestr(Dates(top_gains_idx), 'dd/mm/yyyy');
losses_returns = num2str(LogRet(top_losses_idx));
losses_dates = datestr(Dates(top_losses_idx), 'dd/mm/yyyy');

top_losses = table(losses_dates, losses_returns, ...
                        'VariableNames', {'Date', 'Log_Return'});
top_gains = table(gains_dates, gains_returns, ...
                        'VariableNames', {'Date', 'Log_Return'});


%Update filename according to the chosen ticker
txtfilename = txt_directory + "StatAnal_" + symbol + '_.txt'
log_to_file("# ========================================================", txtfilename)
log_to_file(strjoin(["Statistical analysis of ", Tickers{PickAsset}, " log-return series"], ''), txtfilename);
log_to_file("# ========================================================", txtfilename)

% Display the results
log_to_file(strjoin(["Starting date: ",  datestr(Dates(1))]), txtfilename)
log_to_file(strjoin(["Ending date: " datestr(Dates(end))]), txtfilename)
log_to_file("Log-return Summary Statistics", txtfilename)
log_to_file(Synthesis, txtfilename)

% Empirical and Theoretical Frequencies
log_to_file("--------------------------", txtfilename)
log_to_file('Top 5 Largest Losses', txtfilename)
log_to_file(top_losses, txtfilename, 1)

log_to_file("--------------------------", txtfilename)
log_to_file("Top 5 Largest Gains", txtfilename)
log_to_file(top_gains, txtfilename, 1)

% Empirical and Theoretical Frequencies
log_to_file("--------------------------", txtfilename)
log_to_file("Empirical and Theoretical Frequencies", txtfilename)
log_to_file(FreqTable, txtfilename, 1)

%Jarque-Bera Normality Test
log_to_file("--------------------------", txtfilename)
log_to_file("Jarque-Bera Normality Test", txtfilename)
log_to_file(JBSynthesis, txtfilename)

if JB > CV
    log_to_file('Reject null of normality', txtfilename);
else
    log_to_file('Accept null of normality', txtfilename);
end

log_to_file("--------------------------", txtfilename)
log_to_file("Autocorrelation Analysis", txtfilename)

log_to_file(ACFSynthesis, txtfilename, 1)

log_to_file("--------------------------", txtfilename)
log_to_file("Autocorrelation Tests of returns", txtfilename)
log_to_file('Box Pierce Test', txtfilename)
log_to_file(PB1Synthesis, txtfilename)

log_to_file("Ljung Box Test", txtfilename)

log_to_file(LB1Synthesis, txtfilename)
if LjungBox > CV
    log_to_file('There is evidence of serial correlation in log-returns', txtfilename);
else
    log_to_file('There is no evidence of serial correlation in log-returns', txtfilename);
end

log_to_file("--------------------------", txtfilename)
log_to_file("Autocorrelation of squared returns", txtfilename)

log_to_file("Box Pierce Test", txtfilename)
log_to_file(PB2Synthesis, txtfilename)

log_to_file("Ljung Box Test", txtfilename)
log_to_file(LB2Synthesis, txtfilename)
if LjungBox2 > CV
    log_to_file('There is evidence of serial correlation in squared log-returns', txtfilename);
else
    log_to_file('There is no evidence of serial correlation in squared log-returns', txtfilename);
end

log_to_file("--------------------------", txtfilename)
log_to_file("Autocorrelation of absolute returns", txtfilename)

log_to_file("Box Pierce Test", txtfilename)
log_to_file(PB3Synthesis, txtfilename)
log_to_file("Ljung Box Test", txtfilename)
log_to_file(LB3Synthesis, txtfilename)
if LjungBoxAbs > CV
    log_to_file('There is evidence of serial correlation in absolute log-returns', txtfilename);
else
    log_to_file('There is no evidence of serial correlation in absolute log-returns', txtfilename);
end

log_to_file("# Analysis Completed", txtfilename)
log_to_file("# ========================================================", txtfilename)
