%% Pre-initialize
    clear;
    clc;
    close all;

%% Parameters

    lambda = 50*rand(1);                % Regularaization factor
    n = 1000;                            % Number of measurements
    lb = -2*lambda;                     % Lower bound
    ub = 2*lambda;                      % Upper bound

    y = (ub-lb).*rand(n,1) + lb;        % Measurements
    v = ones(n-1,1);        
    D = diag(v,1);
    D = D(1:n-1,:);
    D(logical(eye(size(D)))) = -1*ones(n-1,1);  % Differencing matrix
    tol = 1e-5;                         % Tolerance
    
%% Solve total variation problem
    l1 = TV_BB(D, y, lambda, tol);
    l2 = TV_l2(D, y, lambda, tol);

%% Plot
    stairs(y,':','LineWidth',1.5);
    hold on;
    stairs(l1.x,'LineWidth',1.5);
    stairs(l2.x,'k');
    legend('Original signal','L1 signal','L2 signal');
