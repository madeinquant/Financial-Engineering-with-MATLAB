%% Financial Engineering with MATLAB(R)

%% Load in some data
% Again we'll import Bund data sampled minutely
close all; clear; clc;

stock='SPY';

lookback = 252;
%cfrom = 'Jan 1 2009';
%cto = 'Aug 27 2013';
cfrom = datestr(today-lookback);
cto = datestr(today);

date_from=datenum(cfrom);
date_to=datenum(cto);
  
adjClose = fetch(yahoo,stock,'adj close',date_from,date_to);
div = fetch(yahoo,stock,date_from,date_to,'v');
returns= (adjClose(2:end,2)./adjClose(1:end-1,2)-1);
logreturns=log(adjClose(2:end,2)./adjClose(1:end-1,2)-1);

volat_d = std(returns);             % Daily volatility
volat = volat_d * sqrt(lookback);    % Annualized volatility

% plot adjusted Close price of  and mark days when dividends
% have been announced
plot(adjClose(:,1),adjClose(:,2),'color',[0.6 0.6 0.6])
hold on;
plot(div(:,1),min(adjClose(:,2))+10,'ob');
ylabel('SPY (US$)');
clable=strcat(cfrom,' to ',cto);
xlabel(clable);

%clear all; close all; clc;
%c = yahoo;
%AdjPrice = fetch(c,'IBM','Adj Close','08/01/10','08/25/15')
%plot(AdjPrice(2))
%tday = datestr(AdjPrice(1:end,1),'yyyymmdd')
%tday=datestr(datenum(tday, 'mm/dd/yyyy'), 'yyyymmdd'); % convert the format into yyyymmdd.

%% Cox-Ross-Rubinstein
% Pricing american options with Cox-Ross-Rubinstein method 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation inputs
%side = 'call';                             % Option side
%style = 'american';                        % Option style
%price = adjClose(length(adjClose),2);    % Current instrument price (147.31, as of 2012/11/05)
%strike = 140;			% Strike price
%riskfree = .0007;		% Risk-free rate, Yield on 3m US Treasury Yields, as of 2012/11/05
%divyield = .0199;		% Dividend yield on S&P 500 (IVV), as of 2012/11/05
%tte = 30; %(datetime(2012,12,22) - datetime(2012,11, 6)).days  	# Time to expiration in days
close all; clear; clc;
side = 'call';  style = 'american';  price = 142.410;
strike = 140.0; riskfree = .0001; divyield = .0199; tte = 46;
volat = 0.182;
%Pre-processing of inputs and calculation of per-step figures
n = 8;							    % Depth of binomial tree (levels are numbered from 0 to n)
tdelta = tte / (n * 365);		    % Time delta per one step (as fraction of year)
u = exp(volat * sqrt(tdelta));		% Up movement per step
d = 1/u;							% Down movement per step
rf = exp(riskfree * tdelta) - 1;	% Risk-free rate per step
dy = exp(divyield * tdelta) - 1;    % Dividend yield per step
pu = (1+rf-dy-d) / (u-d);		    % Probability of up movement
pd = 1 - pu;			            % Probability of down movement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate terminal nodes of binomial tree
level = {};
for i = 0:n
    pr = (price * (d^i) * (u^(n-i))); 	
	% Option value at the node (depending on side)
    if strcmp(side,'call')
        ov = max(0.0, pr-strike);
    else
        ov = max(0.0, strike-pr);
    end;
    level{i+1} = {pr, ov};
    %level = [level, ov];%level.append((pr, ov));
	fprintf('Node Price %.3f, Option Value %.3f\n', pr, ov)
end;

levels = {'None','None','None'}; % Remember levels 0,1,2 for the greeks

% reduce binomial tree
fprintf('Reduce Binomial Tree\n')
for i = n-1:-1:0
    % fprintf('i: %d\n',i)
    levelNext = {};
    fprintf('Tree level %i\n', i)
    for j = 0:i     % Iterate through nodes from highest to lowest price
        % fprintf('j: %d\n',j)
        node_u = level{j+1};
        node_d = level{j+2};
        % Instrument's price at the node
        pr = node_d{1} / d;  
        % Option value at the node (depending on side)
        ov = (node_d{2} * pd + node_u{2} * pu) / (1 + rf);	
        if strcmp(style,'american') % American options can be exercised anytime
            if strcmp(side,'call')
                ov = max(ov, pr-strike);
            else 
                ov = max(ov, strike-pr);
            end;
        end;
        levelNext{j+1} = {pr, ov};
        fprintf('Node Price %.3f, Option Value %.3f\n', pr, ov)
    end;
    level = levelNext;
    if j<=2
        if j>0
            levels{j}=level; % save level 0,1,2 of the tree
        end;
    end;
end;

%% Black-Scholes-Merton
%% Monte Carlo Simulation
%% Explicit Finite Difference Method/American option
close all; clear; clc;

S0 = 50; K = 50; r = 0.1; T = 5/12;
sigma = 0.4; Smax = 100; M = 100; N = 1000; 
is_call = false;

dS = Smax / M; dt = T /  N;
i_values = 0:1:M-1; 
j_values = 0:1:N-1; 
grid = zeros(M+1, N+1);
boundary_conds = linspace(0,Smax,M+1);

% setup_boundary_conditions
if is_call
    grid(:, end) = max(boundary_conds - K, 0);
    grid(end, 1:end-1) = (Smax - K) * exp(-r * dt * (N-j_values));
else
    grid(:, end) = max(K-boundary_conds, 0);
    grid(1, 1:end-1) = (K - Smax) * exp(-r * dt * (N-j_values));
   
end;

% setup_coefficients
a = 0.5*dt*((sigma^2) * (i_values.^2) - r*i_values);
b = 1 - dt*((sigma^2) * (i_values.^2) + r);
c = 0.5*dt*((sigma^2) * (i_values.^2) + r*i_values);

% traverse_grid
for j = length(j_values):-1:1
    for i = 3:M
        grid(i,j) = a(i)*grid(i-1,j+1) + b(i)*grid(i,j+1) + c(i)*grid(i+1,j+1);
    end;
end;

% interpolate
% """ Use piecewise linear interpolation on the initial grid column to get the closest price at S0.
%
interp1q(boundary_conds, grid(:, 1), S0)

%% Implicit Finite Difference Method/American option
close all; clear; clc;

S0 = 50; K = 50; r = 0.1; T = 5/12;
sigma = 0.4; Smax = 100; M = 100; N = 1000; 
is_call = false;

dS = Smax / M; dt = T /  N;
i_values = 0:1:M-1; 
j_values = 0:1:N-1; 
grid = zeros(M+1, N+1);
boundary_conds = linspace(0,Smax,M+1);

% setup_boundary_conditions
if is_call
    grid(:, end) = max(boundary_conds - K, 0);
    grid(end, 1:end-1) = (Smax - K) * exp(-r * dt * (N-j_values));
else
    grid(:, end) = max(K-boundary_conds, 0);
    grid(1, 1:end-1) = (K - Smax) * exp(-r * dt * (N-j_values));
end;

% setup_coefficients
a = 0.5*(r*dt*i_values - (sigma^2)*dt*(i_values.^2));
b = 1 + (sigma^2)*dt*(i_values.^2) + r*dt;
c = -0.5*(r * dt*i_values + (sigma^2)*dt*(i_values.^2));

coeffs = diag(a(3:M), -1) + diag(b(2:M)) + diag(c(2:M-1), 1);

% traverse_grid
% """ Solve using linear systems of equations """
[P, L, U]= lu(coeffs);
aux = zeros(M-1);
for j=N:-1:1
    aux(1) = dot(-a(2), grid(1, j));
    x1 = L \ (grid(2:M, j+1)+aux(1));
    x2 = U \ x1;
    grid(2:M, j) = x2;
end;

% interpolate
% """ Use piecewise linear interpolation on the initial grid column to get the closest price at S0.
%
interp1q(boundary_conds, grid(:, 1), S0)

%% Crank-Nicolson/American option
% setup_coefficients
alpha = 0.25*dt*( (sigma^2)*(i_values.^2) - r*i_values);
beta = -dt*0.5*( (sigma^2)*(i_values.^2) + r);
gamma = 0.25*dt*( (sigma^2)*(i_values.^2) + r*i_values);

M1 = -diag(alpha(2:M), -1) + diag(1-beta(1:M)) - diag(gamma(1:M-1), 1);
M2 = diag(alpha(2:M), -1) + diag(1+beta(1:M)) + diag(gamma(1:M-1), 1);

% traverse_grid
% """ Solve using linear systems of equations """
[P, L, U]= lu(coeffs);
aux = zeros(M-1);
for j=N:-1:1
    aux(0) = dot(-a(1), grid(0, j));
    x1 = solve(L, grid(1:M, j+1)+aux);
    x2 = solve(U, x1);
    grid(1:M, j) = x2;
end;
    
% interpolate
% """ Use piecewise linear interpolation on the initial grid column to get the closest price at S0.
%
interp1q(boundary_conds, grid(:, 1), S0)

%% Non-path-dependent interest rate product
%% Path-dependent interest rate product

%% Two-factor explicit
%% Two-factor implicit

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
