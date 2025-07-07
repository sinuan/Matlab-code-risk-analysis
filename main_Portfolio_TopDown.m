clear; close all; clc;

%**************************************************%
%%%%%%%%TOP DOWN APPROACH TO PORTFOLIO RISK%%%%%%%%%
%COMPARING PARAMETRIC AND NOT PARAMETRIC APPROACHES%
%**************************************************%

%% ============================
%  Load Data and Define Market
% ============================
filename = 'DOW30_merge.csv'; 
myMarket = 'DOW30';  
img_directory = 'Images/'; % Directory for saving figures (must exist)
txt_directory  = "Results/"; % Directory for saving results 
txtfilename = txt_directory + 'TopDown_RiskAnalysis.txt';

% Ensure directories exist
if ~exist(img_directory, 'dir')
    mkdir(img_directory);
end
if ~exist(txt_directory, 'dir')
    mkdir(txt_directory);
end

% Load dataset
dataset = readtable(filename, 'MissingRule', 'omitrow');
ColLabels = dataset.Properties.VariableNames;
Tickers = ColLabels(2:end); % Extract tickers
HistPrices = dataset{:, 2:end}; % Historical prices
HistDates = dataset{:, 1}; % Historical dates

[NObs, NAsset] = size(HistPrices);

% Assign Portfolio quantities
Q = ones(NAsset,1);

% Determine the portfolio value over time
Port = HistPrices * Q;

% Log returns are computed using the natural logarithm of price ratios
LogRet_p = log(Port(2:end) ./ Port(1:end-1));
Dates = HistDates(2:end);
Dates = datetime(Dates, 'InputFormat', 'dd/MM/yyyy');

% Plot Portfolio log-return
figure('Color',[1 1 1]);
plot(Dates, LogRet_p)
title('Portfolio log-return','interpreter','latex')
xlim([Dates(1) Dates(end)])
dateaxis('x', 12)
xlabel('Time','interpreter','latex')
print(gcf, [img_directory, 'TopDown_PortfolioReturns'], '-dpng')
% Define functions for parametric and non-parametric VaR
%Gaussian case
VaR_g = @(ConfLevel, LogRet_p) - norminv(1-ConfLevel, mean(LogRet_p), std(LogRet_p, 1));
ES_g = @(ConfLevel, LogRet_p) -( mean(LogRet_p) - std(LogRet_p, 1)*normpdf(norminv(1-ConfLevel, 0, 1), 0, 1 )/(1-ConfLevel));
%Non parametric case
VaR_np = @(ConfLevel, LogRet_p) -prctile(LogRet_p, (1-ConfLevel)*100); 
ES_np  = @(ConfLevel, LogRet_p)  -mean(LogRet_p(LogRet_p<=-VaR_np(ConfLevel, LogRet_p)))

jj = 1;
ConfLevelRange = [0.90:0.01:0.99];
for aj= ConfLevelRange
   Risk(jj,1) = VaR_g(aj, LogRet_p);
   Risk(jj,2) = VaR_np(aj, LogRet_p);
   Risk(jj,3) = ES_g(aj, LogRet_p);
   Risk(jj,4) = ES_np(aj, LogRet_p);
   jj = jj + 1;
end

% Plot VaR and Expected Shortfall
figure('Color',[1 1 1]);
plot(ConfLevelRange,Risk)
legend('Gaussian VaR','Non-Param. VaR', 'Gaussian ES','Non-Param. ES', ...
    'interpreter','latex','Location','northwest')
xlabel('$\alpha$ (Conf. Level)','interpreter','latex')
title('Portfolio Risk','interpreter','latex')
print(gcf, [img_directory, 'TopDown_Risk_ConfLevel'], '-dpng')

% Plot VaR
% figure('Color',[1 1 1]);
% fplot(@(a) [VaR_g(a, LogRet_p) VaR_np(a, LogRet_p)], [0.9, 0.99])
% hold on
% fplot(@(a) [ES_g(a, LogRet_p) ES_np(a, LogRet_p)], [0.9, 0.99])
% legend('Gaussian VaR','Non-Param. VaR', 'Gaussian ES','Non-Param. ES','interpreter','latex')
% xlabel('$\alpha$ (Conf. Level)','interpreter','latex')
% title('Portfolio Value at Risk','interpreter','latex')
% print(gcf, [img_directory, 'TopDown_Portfolio_VaR'], '-dpng')

% Set confidence level
ConfLevel = 0.90;

% VaR and ES Calculation Horizon n periods

HoldingPeriod = 1:100;
mu_d = mean(LogRet_p); std_d = std(LogRet_p);
skew_d = skewness(LogRet_p); kurt_d = kurtosis(LogRet_p);

for j = 1:length(HoldingPeriod)
    range = sort([NObs:-HoldingPeriod(j):1]);
    Port_hp = Port(range, :);
    LogRet_hp = log(Port_hp(2:end,:) ./ Port_hp(1:end-1,:));
    fitstat(j, :) = [length(LogRet_hp), mean(LogRet_hp), var(LogRet_hp), skewness(LogRet_hp), kurtosis(LogRet_hp)];
    NObs_hp(j, 1) = length(LogRet_hp);
    VaR_gaussian_hp(j, 1) = VaR_g(ConfLevel, LogRet_hp);
    VaR_gaussian_sr(j, 1) = -norminv(1-ConfLevel, mu_d*HoldingPeriod(j), std_d*sqrt(HoldingPeriod(j)));
    VaR_nonparam_hp(j, 1) = VaR_np(ConfLevel, LogRet_hp);
    ES_gaussian_hp(j, 1) = ES_g(ConfLevel, LogRet_hp);
    ES_nonparam_hp(j, 1) = ES_np(ConfLevel, LogRet_hp);
end

xxx = length(HoldingPeriod);
% Create Summary Table
StatSynthesis = table(fitstat(:,1), fitstat(:,2), fitstat(:,3), fitstat(:,4), fitstat(:,5));
StatSynthesis .Properties.VariableNames = {'N. Obs', 'Exp. Return', 'Variance', 'Skewness', ...
                    'Kurtosis'};

% Create Summary Table
VaRSynthesis = table((1:xxx)', NObs_hp(:,1), VaR_gaussian_hp(:,1), VaR_gaussian_sr(:,1), VaR_nonparam_hp(:, 1), ES_gaussian_hp(:,1), ES_nonparam_hp(:,1));
VaRSynthesis .Properties.VariableNames = {'Holding Period', 'N. Obs', 'Gaussian VaR', 'Gaussian VaR (Sq.Root Rule)', 'Non Parametric VaR', 'Gaussian ES', 'Non Parametric ES'};

% Plot Moments Scaling
figure('Color', [1 1 1])
subplot(2, 2, 1)
plot(1:xxx, [fitstat(:,2) mu_d*(1:xxx)'])
legend('Expected Return', 'Scaling using $n$','interpreter','latex','Location','best')
xlabel('Days', 'Interpreter','latex')
ylabel('Expected Return', 'Interpreter','latex')

subplot(2, 2, 2)
plot(1:xxx, [fitstat(:,3) std_d^2*(1:xxx)'])
legend('Sample Variance', 'Scaling using $n$','interpreter','latex','Location','best')
xlabel('Days', 'Interpreter','latex')
ylabel('Variance', 'Interpreter','latex')

subplot(2, 2, 3)
plot(1:xxx, [fitstat(:, 4) skew_d./sqrt((1:xxx)')])
xlabel('Days', 'Interpreter','latex')
ylabel('Skewness', 'Interpreter','latex')
legend('Sample Skewness', 'Scaling using $1/\sqrt{n}$','interpreter','latex','Location','best')

subplot(2, 2, 4)
plot(1:xxx, [fitstat(:, 5) kurt_d./(1:xxx)'])
xlabel('Days', 'Interpreter','latex')
ylabel('Kurtosis', 'Interpreter','latex')
legend('Sample Kurtosis', 'Scaling using $1/n$','interpreter','latex','Location','best')
print(gcf, [img_directory, 'TopDown_Scaling_nperiods'], '-dpng')

% Plot VaR & ES for different holding periods
figure('Color', [1 1 1])
plot(HoldingPeriod, [VaR_gaussian_hp VaR_gaussian_sr VaR_nonparam_hp ES_gaussian_hp ES_nonparam_hp])
legend('Gaussian VaR', 'Gaussian VaR (Sq. Root Rule)','Non-Parametric VaR','Gaussian ES','Non-Parametric ES','interpreter','latex','Location','best')
xlabel('Holding Period $n$','interpreter','latex')
title(['Portfolio Risk $(\alpha:$ ' num2str(ConfLevel*100) '\%)'],'interpreter','latex')
print(gcf, [img_directory, 'TopDown_Portfolio_VaR_nperiods'], '-dpng')

%%%%%%%SAVE RESULTS TO FILE%%%%%%%%%%%%%
log_to_file("# ========================================================", txtfilename)
log_to_file("TOP DOWN APPROACH TO RISK ", txtfilename);
log_to_file("# ========================================================", txtfilename)
log_to_file(strjoin(["Starting date: ",  datestr(Dates(1))]), txtfilename)
log_to_file(strjoin(["Ending date: " datestr(Dates(end))]), txtfilename)

log_to_file(strjoin(["Number of Assets: " num2str(NAsset)]), txtfilename, 1)

log_to_file('', txtfilename, 1)
log_to_file('-------------------------', txtfilename, 1)
log_to_file('Risk measures for different confidence levels', txtfilename, 1)
log_to_file('Conf.Level VaRG VaRNP ESG ESNP', txtfilename, 1)
log_to_file([ConfLevelRange' Risk], txtfilename, 1)
log_to_file('-------------------------', txtfilename, 1)
log_to_file('Stat measures for different holding periods', txtfilename, 1)
log_to_file(StatSynthesis, txtfilename, 1)
log_to_file('', txtfilename, 1)
log_to_file('-------------------------', txtfilename, 1)
log_to_file('Risk measures for different holding periods', txtfilename, 1)
log_to_file(VaRSynthesis, txtfilename, 1)
log_to_file('', txtfilename, 1)
log_to_file("# Analysis Completed", txtfilename)
log_to_file("# ========================================================", txtfilename)


