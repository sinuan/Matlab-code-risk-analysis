close all
clear all
clc

%%%%%INPUT
%Threshold in percentage
Threshold = 0.95 
%Fix time horizon
Horizon = [1:5000];
%Number of simulations
M = 1000

%Choose dataset and directory where to store figures
filename = 'DOW30merge.csv';
myMarket = 'DOW30'
imgdirectory = 'Images/'%this dir must be created

%Load the the data set and obtain relevant statitistical information 
Dataset = readtable(filename, MissingRule = 'omitrow');
ColLabels = Dataset.Properties.VariableNames;
StockTicker = ColLabels(2:end);

%Load historical prices
StockPrices = Dataset(:,2:end).Variables;

%Determine n. of obs. (NObs) and n. of assets (NAsset)
[NObs NAsset] = size(StockPrices);

%Choose at random an Asset
PickAsset = randi(NAsset); %1
StockName = StockTicker{PickAsset}
display(['The chosen stock is: ' StockName])

% Compute daily log-returns and estimate mean, and volatility
LogRet = log(StockPrices (2:end,PickAsset)./StockPrices (1:end-1,PickAsset))
LogRetMean = mean(LogRet);
LogRetStd = std(LogRet);

%Find current price
P_0 = StockPrices(end)

%Find threshold price
P_Threshod = Threshold*P_0;
Loss = P_0*(1-Threshold)

for n = Horizon
    Mean_n = LogRetMean*n;
    Vol_n = LogRetStd*sqrt(n);%
    Z = (log(Threshold) - Mean_n)/Vol_n ;
    ThProb(n, 1) = normcdf(Z);
    r_n = Mean_n + Vol_n *randn(M,1);
    P_n = P_0*exp(r_n);
    MCProb(n, 1) = mean(P_n <= P_Threshod); 
end

%Horizon with maximum probability
nhat = -log(Threshold)/(LogRetMean);
Mean_nhat = LogRetMean * nhat;
Vol_nhat = LogRetStd*nhat^0.5;%
Z_nhat = (log(Threshold) - Mean_nhat )/Vol_nhat ;
MaXProb = normcdf(Z_nhat);

figtitle = [StockName ': Probability of losing more than ' num2str(Loss) ' USD at different horizons']
h = figure('Color', [1 1 1])
plot(Horizon, [ThProb MCProb])
hold on
plot(nhat, MaXProb, 'b*')
xlim([0 500])
xlabel('Horizon (days)', 'Interpreter', 'latex')
ylabel('Probability', 'Interpreter', 'latex')
legend('Exact', ['MC (M=' num2str(M) ')'],  'Interpreter', 'latex')
title(figtitle,  'Interpreter', 'latex')
figname = ['Images/ProbGaining_Sim_500_M' num2str(M)]
print(h, figname, '-dpng')

h = figure('Color', [1 1 1])
plot(Horizon, [ThProb MCProb])
hold on
plot(nhat, MaXProb, 'b*')
xlabel('Horizon (days)', 'Interpreter', 'latex')
ylabel('Probability', 'Interpreter', 'latex')
legend('Exact', ['MC (M=' num2str(M) ')'],  'Interpreter', 'latex')
title(figtitle,  'Interpreter', 'latex')
figname = ['Images/ProbGaining_Sim_5000_M' num2str(M)]
print(h, figname, '-dpng')

%Results
format long
Range = [10 50 100 200 250 500 1000 2000 3000 4000 5000]
Results = [Horizon(Range)' ThProb(Range) MCProb(Range)]



