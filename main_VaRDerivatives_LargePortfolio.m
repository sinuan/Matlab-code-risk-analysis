clear all 
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Stock parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
P0 = 300; %underlying stock price
TTM = 1; %option time to maturity
sg = 0.2; %annualized implied volatility
rf = 0.05; %risk-free rate
q = 0; %dividend yield

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% VaR/ES parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = 0.99; %confidence level
days  = 50 ; 
VaRHorizon = days/250; %VaR Horizon: 10, 30, 90, 245

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MC parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 10000; %Number of simulations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FixedOrRandom = 2;
if FixedOrRandom  == 1
    cp = [1 -1]%calls and puts
    units = [1 1]; %number of options
    longshort = units./abs(units);
    Strike = [310 290];
else 
    NumRnd = 30;
    u = rand(1,NumRnd);
    cp = (u<0.5) - 1*(u>0.5);%calls and puts
    units = randi([-10, 10], 1, NumRnd) + 10^(-6); %number of options; avoid exactly zero positions
    longshort = units./abs(units);
    Strike = randi([0.5*P0, 2*P0], 1, NumRnd);
end
NumOpts = length(Strike);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option characteristics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
payoff = @(cp,S,K) max(cp*(S-K),0)
SRange = [100:500];
for j = 1:NumOpts 
    OptPayoff(j,:) = payoff(cp(j), SRange, Strike(j));
end
PortPayoff = units*OptPayoff; 

h = figure('Color',[1 1 1])
plot(SRange, PortPayoff,'.')
xlabel('$S(T)$','interpreter','latex')
ylabel('Portfolio payoff','interpreter','latex')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CALCULATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Price Options and compute portfolio value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:NumOpts 
    optprice(j,1) = get_optionprice(cp(j), P0, Strike(j), rf, TTM, sg);
end
Port0 = units*optprice;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate returns
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r  = randn(M,1)*sg*VaRHorizon^0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate prices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PT = P0*exp(r); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Full Reevaluation
% Reprice options value and compute their P&L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1:NumOpts 
    optpriceHor(j,:) = get_optionprice(cp(j), PT, Strike(j), rf, TTM-VaRHorizon, sg);
    PL_opt(j,:) = longshort(j)*(optpriceHor(j,:) - optprice(j));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulate portfolio value and P&L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PortHor = units*optpriceHor;%(1xNumOpts) (NumOptsxM)=(1XM)
PL_port = PortHor-Port0;


h = figure('Color',[1 1 1])
plotmatrix([PT, PortHor' ])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute risk measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[VaR_opt  ES_opt] = get_risk_sample(PL_opt', alpha);
[VaR_port  ES_port] = get_risk_sample(PL_port', alpha)
VaRanalysis = table(VaR_opt  , VaR_port)
ESanalysis = table(ES_opt, ES_port)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulated distributions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if FixedOrRandom  == 1
    MinAll = min([min(PL_opt(1,:)), min(PL_opt(2,:)), min(PL_port)]);
    MaxAll = max([max(PL_opt(1,:)), max(PL_opt(2,:)), max(PL_port)]);
    h = figure('Color',[1 1 1])
    subplot(3,1,1)
    histogram(PL_opt(1,:), 'normalization','pdf')
    xlabel('P\&L (contract 1)','interpreter','latex')
    ylabel('PDF','interpreter','latex')
    xlim([MinAll MaxAll]) 
    title('Simulated PDF','interpreter','latex')
    subplot(3,1,2)
    histogram(PL_opt(2,:), 'normalization','pdf')
    xlabel('P\&L (contract 2)','interpreter','latex')
    ylabel('PDF','interpreter','latex')
    xlim([MinAll MaxAll]) 
    subplot(3,1,3)
    histogram(PL_port, 'normalization','pdf')
    xlabel('P\&L (portfolio)','interpreter','latex')
    ylabel('PDF','interpreter','latex')
    xlim([MinAll MaxAll]) 
else
    MinAll = min(PL_port);
    MaxAll = max(PL_port);
    h = figure('Color',[1 1 1])
    histogram(PL_port, 'normalization','pdf')
    hold on
    fplot(@(X) normpdf(X, mean(PL_port), std(PL_port)), [MinAll ,MaxAll])
    xlabel('P\&L (portfolio)','interpreter','latex')
    ylabel('PDF','interpreter','latex')
    xlim([MinAll MaxAll]) 
    legend('Simulated','Gaussian','interpreter','latex')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compare VaR of individual options with Portfolio VaR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[units' VaR_opt'  ES_opt']
[VaR_port sum(abs(units'.*VaR_opt')) ES_port sum(abs(units'.*ES_opt'))]

%stop

h = figure('Color',[1 1 1])
X = categorical([0:NumOpts] )
bar(X, [VaR_port abs(units'.*VaR_opt')'])
hold on
bar([VaR_port],'r')
xlabel('Components','interpreter','latex')
ylabel('VaR','interpreter','latex')

h = figure('Color',[1 1 1])
subplot(2,1,1)
X = categorical({'Portfolio', 'Sum of VaRs'} )
bar(X, [VaR_port sum(abs(units'.*VaR_opt'))])
ylabel('VaR','interpreter','latex')

subplot(2,1,2)
X = categorical({'Portfolio', 'Sum of ESs'} )
bar(X, [ES_port sum(abs(units'.*ES_opt'))])
ylabel('Expected Shortfall','interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Risk contribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute Marginal VaR
eps = VaR_port*0.01;
position = find((PL_port<-VaR_port+eps).*(PL_port>-VaR_port-eps));
for j = 1:NumOpts 
    MVaR(j)  = -mean(PL_opt(j,position));
    CVaR(j)  = MVaR(j)*units(j);
end
%check
[sum(CVaR) VaR_port]
CVaRp = CVaR/sum(CVaR)

%Compute Marginal Expected Shortfall
position = find((PL_port<-VaR_port));
for j = 1:length(Strike)
    MES(j)  = -mean(PL_opt(j,position));
    CES(j)  = MES(j)*units(j);
end
[sum(CES) ES_port]
CESp = CES/sum(CES);



if FixedOrRandom  == 1
    h = figure('Color',[1 1 1])
    subplot(2,1,1)
    X = categorical({'Contract 1','Contract 2'});
    bar(X, CVaRp,0.4 )
    subplot(2,1,2)
    X = categorical({'Contract 1','Contract 2'});
    bar(X, CESp,0.4 )
    title('Risk Contribution','interpreter','latex')
else
    X = [1:NumRnd];
    h = figure('Color',[1 1 1])
    subplot(2,1,1)
    bar(X, CVaRp,0.4 )
    subplot(2,1,2)
    bar(X, CESp,0.4 )
    title('Risk Contribution','interpreter','latex')

end