%% ===================================================
%  STUDENT'S T DISTRIBUTION FITTING AND RISK ANALYSIS
%  ===================================================
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

% Randomly Select an Asset
PickAsset = randi(NAsset, 1);
disp(['Asset: ', Tickers{PickAsset}]);

% Compute Log-Returns
LogRet = log(HistPrices(2:end, PickAsset) ./ HistPrices(1:end-1, PickAsset)); 


%% ===================================================
%  PLOT LEFT TAIL OF STUDENT'S T DISTRIBUTION
%  ===================================================
h = figure('Color',[1 1 1]);
hold on;
range = [-5 -2];

% Plot PDFs for different degrees of freedom
v_values = [4, 10, 20, 100]; 
for v = v_values
    fplot(@(x) pdf('tLocationScale', x, 0, 1, v), range);
end
fplot(@(x) pdf('norm', x, 0, 1), range); % Normal distribution

legend('$\nu=4$', '$\nu=10$', '$\nu=20$', '$\nu=100$', '$\nu=\infty$', ...
       'Location', 'best', 'Interpreter', 'latex');
xlabel('Returns', 'Interpreter', 'latex');
ylabel('Student t PDF', 'Interpreter', 'latex');
print(h, '-dpng', fullfile(imgdirectory, 'Tut_StudentT_Tailpdf'));

%% ===================================================
%  FIT STUDENT'S T DISTRIBUTION USING METHOD OF MOMENTS
%  ===================================================
nobs = length(LogRet);
mu_mm = mean(LogRet); 
nu_mm = 4 + 6 / (kurtosis(LogRet) - 3); % Degrees of freedom estimate
sg_mm = sqrt(((nu_mm - 2) / nu_mm) * var(LogRet));

% Display Results
disp(table(nobs, mu_mm, sg_mm, nu_mm, 'VariableNames', {'nobs', 'mu', 'sigma', 'nu'}, 'RowNames', {'MM'}));

%% ===================================================
%  FIT STUDENT'S T DISTRIBUTION USING MAXIMUM LIKELIHOOD
%  ===================================================
phat = mle(LogRet, 'distribution', 'tLocationScale');
mu_ml = phat(1);
sg_ml = phat(2);
nu_ml = phat(3);

% Display Results
disp(table(nobs, mu_ml, sg_ml, nu_ml, 'VariableNames', {'nobs', 'mu', 'sigma', 'nu'}, 'RowNames', {'ML'}));

%% ===================================================
%  COMPARE EMPIRICAL, GAUSSIAN, AND STUDENT T PDFS
%  ===================================================
h = figure('Color',[1 1 1]);
hold on;

% Plot Student's T (MM and ML)
fplot(@(x) pdf('tLocationScale', x, mu_mm, sg_mm, nu_mm), [min(LogRet), max(LogRet)], 'r');
fplot(@(x) pdf('tLocationScale', x, mu_ml, sg_ml, nu_ml), [min(LogRet), max(LogRet)], 'b');
fplot(@(x) normpdf(x, mean(LogRet), std(LogRet)), [min(LogRet), max(LogRet)]);

% Histogram of Empirical Returns
histogram(LogRet, round(sqrt(nobs)), 'Normalization', 'pdf');
legend('Student''s T (MM)', 'Student''s T (MLE)', 'Gaussian', 'Empirical', 'Location', 'best', 'Interpreter', 'latex');
xlabel('Log-return', 'Interpreter', 'latex');
ylabel('PDF', 'Interpreter', 'latex');
print(h, '-dpng', fullfile(imgdirectory, 'Tut_StudentT_pdf.png'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QQPlot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prc = linspace(1,100,100)/100;
prc_t= -(mu_ml+sg_ml*tinv(1-prc, nu_ml));
prc_sample = prctile(LogRet, prc*100);

h=figure('Color',[1 1 1])
plot(prc_t,prc_t,'b')
hold on
plot(prc_t,prc_sample, 'r.')
xlabel('Theoretical percentiles','interpreter','latex')
ylabel('Empirical percentiles','interpreter','latex')
print(h,'Tut_StudentT_qqplot.png','-dpng')

%% ===================================================
%  CREATE HANDLE FUNCTIONS
%  ===================================================

VaR_g = @(a, mu, sg) - (mu + sg * icdf('norm', 1 - a, 0, 1));
ES_g = @(a, mu, sg) - (mu - pdf('norm', icdf('norm', 1 - a, 0, 1), 0, 1) * sg / (1 - a));

VaR_t = @(a, mu, sg, nu) - (mu + sg * icdf('T', 1 - a, nu));
ES_t = @(a , mu, sg, nu) - quadgk(@(x) x .* pdf('tLocationScale', x, mu, sg, nu), -inf, -VaR_t(a, mu, sg, nu)) / (1 - a);

VaR_np = @(alpha, LogRet)  -prctile(LogRet, (1-alpha)*100);
ES_np = @(alpha, LogRet) -mean(LogRet.*(LogRet<-VaR_np(alpha, LogRet)))/(1-alpha) + ...
        ((1-alpha)-mean((LogRet<-VaR_np(alpha, LogRet)))).*VaR_np(alpha, LogRet)/(1-alpha)

%% ===================================================
%  COMPUTE AND COMPARE VALUE-AT-RISK (VaR) AND EXPECTED SHORTFALL (ES)
%  ===================================================
alpha = [0.9, 0.95, 0.99];

for j = 1:length(alpha)
    VaR(j, :) = [1 - alpha(j), VaR_g(alpha(j), mu_mm, sg_mm), ...
                VaR_t(alpha(j), mu_mm, sg_mm, nu_mm) , ...
                VaR_t(alpha(j), mu_ml, sg_ml, nu_ml), ...
                VaR_np(alpha(j), LogRet)];
    
    ES(j, :) = [1 - alpha(j), ES_g(alpha(j), mu_mm, sg_mm), ...
                ES_t(alpha(j), mu_mm, sg_mm, nu_mm), ...
                ES_t(alpha(j), mu_ml, sg_ml, nu_ml), ...
                ES_np(alpha(j), LogRet)];
end

disp(array2table(VaR, 'VariableNames', {'ConfidenceLevel', 'VaR G', 'VaR T (MM)' , 'VaR T (ML)', 'VaR NP'}, 'RowNames', cellstr(num2str(alpha(:)))));

disp(array2table(ES, 'VariableNames', {'ConfidenceLevel', 'ES_Gaussian', 'ES T (MM)', 'ES T (ML)', 'ES NP'}, 'RowNames', cellstr(num2str(alpha(:)))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulating a Student t at different horizons 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NSim=100000; %number of simulations
Horizon=10; %number of periods
Tret=mu_mm+sg_mm*trnd(nu_mm,NSim,Horizon); %1 period simulated returns
CumRet=[zeros(NSim,1), cumsum(Tret,2)]; %simulated cumulative returns

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute VaR and ES using simulated returns at different horizons 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.99;
for j=0:Horizon
    VaR_g_n(j+1) = VaR_g(alpha, mean(LogRet)*j, std(LogRet)*j^0.5); 
    ES_g_n(j+1) = ES_g(alpha, mean(LogRet)*j, std(LogRet)*j^0.5);
    VaR_t_n(j+1) = VaR_np(alpha, CumRet(:,j+1));
    ES_t_n(j+1) = ES_np(alpha, CumRet(:,j+1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Term structure of VaR and ES using Gaussian and Student's T distribution
h=figure('Color',[1 1 1])
plot([0:Horizon], [VaR_g_n', ES_g_n', VaR_t_n', ES_t_n'])
xlabel('Horizon (days)','interpreter','latex')
title('Risk Measures','interpreter','latex')
legend('VaR G','ES G','VaR T', 'ES T', 'interpreter','latex','location','best')
print(h, '-dpng', fullfile(imgDir, 'Tut_StudentT_VaR_horizon.png'))

%Simulated distribution at the horizon
h=figure('Color',[1 1 1])
histogram(CumRet(:,end), round(sqrt(NSim)),'normalization','pdf')
hold on
xmin=mean(LogRet)*Horizon-5*std(LogRet)*Horizon^0.5;
xmax=mean(LogRet)*Horizon+5*std(LogRet)*Horizon^0.5;
fplot(@(x) normpdf(x,mean(LogRet)*Horizon,std(LogRet)*Horizon^0.5), ...
    [xmin xmax])
xlim([xmin xmax])
xlabel('Cumulative return', 'interpreter','latex')
title(['Simulated PDF (Horizon: ', num2str(Horizon) ' day)'], 'interpreter','latex')
legend('Simulated Student''s T','Gaussian', 'interpreter','latex','location','northwest')
print(h, '-dpng', fullfile(imgDir, 'Tut_StudentT_MC.png'))

%Simulated returns at different horizons and 
h=figure('Color',[1 1 1])
plot([0:Horizon],CumRet')
hold on
plot([0:Horizon],prctile(CumRet, [alpha*100, (1-alpha)*100]),'g*')
hold on
plot([0:Horizon],[min(CumRet)' max(CumRet)'],'r*')
xlabel('Horizon (days)','interpreter','latex')
title('Simulated Cumulative Returns','interpreter','latex')
print(h, '-dpng', fullfile(imgDir, 'Tut_StudentT_ScenarioMC.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kurtosis at different horizons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure('Color',[1 1 1])
plot([0:Horizon],kurtosis(CumRet))
xlabel('Horizon (days)','interpreter','latex')
title('Kurtosis of cum. returns','interpreter','latex')
print(h, '-dpng', fullfile(imgDir, 'Tut_KurtosisT_MC.png'));

