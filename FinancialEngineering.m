close all; clear all; clc;

stock='SPY';

lookback = 252;
%cfrom = 'Jan 1 2009';
%cto = 'Aug 27 2013';
cfrom = datestr(today-lookback);
cto = datestr(today);

date_from=datenum(cfrom);
date_to=datenum(cto);
  
adjClose = fetch(yahoo,stock,'adj close',date_from,date_to);
div = fetch(yahoo,stock,date_from,date_to,'v')
returns= (adjClose(2:end,2)./adjClose(1:end-1,2)-1);
logreturns=log(adjClose(2:end,2)./adjClose(1:end-1,2)-1);

volat_d = std(returns);             % Daily volatility
volat = volat_d * sqrt(lookback)    % Annualized volatility

% plot adjusted Close price of  and mark days when dividends
% have been announced
plot(adjClose(:,1),adjClose(:,2),'color',[0.6 0.6 0.6])
hold on;
plot(div(:,1),min(adjClose(:,2))+10,'ob');
ylabel('SPY (US$)');
clable=strcat(cfrom,' to ',cto);
xlabel(clable);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation inputs
%side = 'call';                             % Option side
%style = 'american';                        % Option style
%price = adjClose(length(adjClose),2);    % Current instrument price (147.31, as of 2012/11/05)
%strike = 140;			% Strike price
%riskfree = .0007;		% Risk-free rate, Yield on 3m US Treasury Yields, as of 2012/11/05
%divyield = .0199;		% Dividend yield on S&P 500 (IVV), as of 2012/11/05
%tte = 30; %(datetime(2012,12,22) - datetime(2012,11, 6)).days  	# Time to expiration in days
side = 'call';  style = 'american';  price = 142.410;
strike = 140.0; riskfree = .0001; divyield = .0199; tte = 46;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear all; close all; clc;
%c = yahoo;
%AdjPrice = fetch(c,'IBM','Adj Close','08/01/10','08/25/15')
%plot(AdjPrice(2))
%tday = datestr(AdjPrice(1:end,1),'yyyymmdd')
%tday=datestr(datenum(tday, 'mm/dd/yyyy'), 'yyyymmdd'); % convert the format into yyyymmdd.

%% Financial Engineering in Matlab

%% Numerical Methods (Option Pricing Model)
% Pricing European options with the Binomial Tree 
clear
s.M = s.N + 1 % Number of terminal nodes of tree
s.u = 1 + s.pu % Expected value in the up state
s.d = 1 - self.pd % Expected value in the down state
s.qu = (math.exp((self.r-self.div)*self.dt) - self.d) / (self.u-self.d)
s.qd = 1-self.qu
% Pricing American options with the Binomial Tree
% Overview of Numerical Methods 
% Finite-difference Methods for One-factor Models 
% Further Finite-difference Methods for One-factor Models
% Monte Carlo Simulation
% https://github.com/alapala/finite_difference_methods
%%%%%%%%%%%%%%% Problem and method parameters %%%%%%%%%%%%%%%%% 
s = 0; k = 0; r = 0; t = 0; sd = 0; n = 0;
u = exp(sd*sqrt(t/n));
d = 1/u;
p = (exp(r*t/n) - d) / (u - d);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over each node and calculate the Cox Ross Rubinstein underlying price tree
priceTree = nan(steps+1,steps+1);
priceTree(1,1) = S0;
for idx = 2:steps+1
    priceTree(1:idx-1,idx) = priceTree(1:idx-1,idx-1)*u;
    priceTree(idx,idx) = priceTree(idx-1,idx-1)*d;
end

% Calculate the value at expiry
valueTree = nan(size(priceTree));
switch oType
    case 'PUT'
        valueTree(:,end) = max(X-priceTree(:,end),0);
    case 'CALL'
        valueTree(:,end) = max(priceTree(:,end)-X,0);
end

% Loop backwards to get values at the earlier times
steps = size(priceTree,2)-1;
for idx = steps:-1:1
    valueTree(1:idx,idx) = ...
        exp(-r*dt)*(p*valueTree(1:idx,idx+1) ...
        + (1-p)*valueTree(2:idx+1,idx+1));
    if earlyExercise
        switch oType
            case 'PUT'
                valueTree(1:idx,idx) = ...
                    max(X-priceTree(1:idx,idx),valueTree(1:idx,idx));
            case 'CALL'
                valueTree(1:idx,idx) = ...
                    max(priceTree(1:idx,idx)-X,valueTree(1:idx,idx));
        end
    end
end

% Output the option price
oPrice = valueTree(1)

% Explicit finite difference

% Implicit finite difference

% Crank-Nicolson finite difference


%% Time Series Analysis
%% Preductive Model

%% The value of time

% 3.1 Present value
% An investment project has an investment cost of 100 today, 
% and produces cash flows of 75 each of the next two years.
% What is the Net Present Value of the project?
C=[-100 75 75];
t=[0 1 2];
r = 0.1;
d=(1/(1+r)).^t
NPV = C*d'

% 3.2.1 Internal rate of return.
% We are considering an investment with the following cash flows at dates 0, 1 and 2:
% C0 = -100; C1 = 10; C2 = 110
% 1. The current interest rate (with discrete, annual compounding) is 5%. Determine 
%    the present value of the cash flows.
% 2. Find the internal rate of return of this sequence of cash flows.
C=[-100 10 110];
t=[0 1 2];
r = 0.05;
d=(1/(1+r)).^t;
NPV=C*d';
IRR = irr(C)

% 3.3 Continously compounded interest
% 1. Given a 15% interest rate with monthly compounding, 
% calculate the equivalent interest rate with continuous compounding.
% 2. Given a 12% interest rate with continuous compounding, find the 
% equivalent interest rate with quarterly compounding.
% Carrying out the calculations:
r = 12 * log( 1+0.15/12)
r4 = 4 * ( exp(0.12/4)-1)
