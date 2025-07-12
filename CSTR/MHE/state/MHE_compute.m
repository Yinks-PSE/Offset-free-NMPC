function [x_MHE, dist_state, dist_output] = MHE_compute(ys_meas, u_first_moves, N_MHE, y_meas)

%global Nsteps
%This is the Moving Horizon Estimator Function.
%The length of the horizon is 5

hstep = 0.005;
%Nsteps = 3;
Ts = 0.05;

global x_apriori x_aposteriori P Rx Rdx Ry dist_states;

% Initialization of the aposteriori state, covariance, and weighting matrices
if isempty(x_apriori)
    x_aposteriori = [0.877; 324.5];  % 4 states
    P = diag([1e-6;1e-6]);  % Update for 4 states
    dist_states = [0; 0];  % 4 state disturbances
    
    Rx = diag([1e-4; 1e-6]);  % Update for 4 states
    Ry = diag(1);  % 2 outputs
    Rdx = diag([1; 1]);  % Update for 4 states

    P = covariance_update(P, Rx, Ry, x_aposteriori, u_first_moves(end, :), Ts);
end

persistent d_prev
if isempty(d_prev)
    d_prev = zeros(2 * N_MHE, 1);  % Update for 4 states
end

% Predict the apriori state based on the aposteriori state
x_apriori = state_predictor(x_aposteriori, u_first_moves(end, :), 1, 1, hstep, dist_states)';

% Initial decision variables for the optimizer
dec_var_initial = [x_apriori; zeros(2 * N_MHE, 1)];
dec_var_LB = [0; 300; -0.01 * ones(2 * N_MHE, 1)];
dec_var_UB = [1; 360; 0.01 * ones(2 * N_MHE, 1)];

% Optimizer call using fmincon
[dec_var, cost_MHE] = fmincon(@(dec_var) MHE_cost_fun(dec_var, x_apriori, d_prev, ys_meas, u_first_moves, P, Rdx, Ry, hstep, N_MHE), ...
    dec_var_initial, [], [], [], [], dec_var_LB, dec_var_UB);

disp(cost_MHE);
d_prev = dec_var(3:end);

% Predict the state throughout the horizon
x_MHE = ones(N_MHE, 2);  % 4 states
x_past_init = dec_var(1:2);  % Initialize with first 4 states

j = 3;
for i = 1:N_MHE
    x_MHE(i, :) = state_predictor(x_past_init, u_first_moves(i, :), 1, 1, hstep, dec_var(j:j+1));
    j = j + 2;
    x_past_init = x_MHE(i, :)';
end

x_MHE = x_MHE(end, :)';

% Update aposteriori state
x_aposteriori = dec_var(1:2);

% Retrieve state disturbances from the decision vector
dist_states = dec_var((2 * N_MHE + 1):(2 * N_MHE + 2));
dist_state = dist_states;

%disp(size(dist_state))

% Calculate output disturbance
dist_output = ys_meas(N_MHE, :)' - x_MHE(2);

% Update covariance matrix
P = covariance_update(P, Rx, Ry, x_MHE, u_first_moves(end, :), Ts);

end

function cost = MHE_cost_fun(dec_var, x_apriori, d_prev, ys_meas, u_first_moves, P, Rdx, Ry, hstep, N_MHE)

j = 3;
x_pasts = ones(N_MHE, 2);  % 4 states
x_past_init = dec_var(1:2);

%disp(size(x_past_init))

for i = 1:N_MHE
    x_pasts(i, :) = state_predictor(x_past_init, u_first_moves(i, :), 1, 1, hstep, dec_var(j:j+1));
    j = j + 2;
    x_past_init = x_pasts(i, :)';
end

y1_past = x_pasts(:, 1); 
y1_meas = ys_meas(:, 1);

% Stage cost based on the output prediction
PIy = Ry^(-1);
PIdx = Rdx^(-1);

stage_cost_1 = (y1_meas - y1_past)' * PIy(1, 1) * (y1_meas - y1_past);

% Disturbance cost
dec_var_x1_dist = (dec_var(3:2:end) - d_prev(1:2:end));
dec_var_x2_dist = (dec_var(4:2:end) - d_prev(2:2:end));

stage_cost_2 = PIdx(1, 1) * (dec_var_x1_dist' * dec_var_x1_dist) + ...
               PIdx(2, 2) * (dec_var_x2_dist' * dec_var_x2_dist);

% Total stage cost
stage_cost = stage_cost_1 + stage_cost_2;

% Arrival cost
arrival_cost = (x_apriori - dec_var(1:2))' * (P^(-1)) * (x_apriori - dec_var(1:2));

% Calculate total cost
if isreal(dec_var) && isreal(x_pasts)
    cost = stage_cost + arrival_cost;
else
    cost = inf;
end

end

function P_matrix = covariance_update(P_matrix, Rx, Ry, x_MHE, u, Ts)

%values of parameters
 F = 100;
    V = 100;
    k0 = 7.2e10;
    E_R = 8750;
    DeltaH = -5e4;
    rho = 1000;
    Cp = 0.239;
    UA = 5e4;
    CAF = 1.0;
    T_f = 350;

    %%%%%%INPUT PARAMETER
Tc= u(1);

Ca = x_MHE(1);
T = x_MHE(2);

%Ts=1 ;



% calculate the Jacobians for state and measurement equations

a11 = 1 + Ts*(-F/V - k0 * exp(-E_R / T));
a12 = Ts*(Ca * k0 * (E_R / T^2) * exp(-E_R / T));

a21 = Ts*(-DeltaH * k0 * exp(-E_R / T)) / (rho * Cp);
a22 = 1 + Ts*(-F/V - (DeltaH * Ca * k0 * E_R / (rho * Cp * T^2)) * exp(-E_R / T) - UA / (V * rho * Cp));

A = [a11,a12;
    a21,a22];

    C = [0, 1];  % Measurement matrix, assuming direct state observation

    G=[1, 0;
        0, 1]  ;
    


P_matrix = G*Rx*G' + A*P_matrix*A' - A*P_matrix*C'*(Ry + C*P_matrix*C')^(-1)*C*P_matrix*A';

end


function x = state_predictor(x_estim,u,Nu,Np,hstep,dist_state)

%This function predicts values of states throughout the prediction horizon. i.e.
%From k+1 to Np+1 using Runge-Kutta discretization.
%This function requires x_estim to be a column vector

%This function predicts values of states throughout the prediction horizon. i.e.
%From k+1 to Np+1 using Runge-Kutta discretization.
%This function requires x_estim to be a column vector.
%h = 0.5;

k1 = hstep*rate_of_state_change(x_estim',u(1));
k2 = hstep*rate_of_state_change(x_estim'+(k1/2),u(1));
k3 = hstep*rate_of_state_change(x_estim'+(k2/2),u(1));
k4 = hstep*rate_of_state_change(x_estim'+k3,u(1));
x(1,:) = x_estim' + (1/6)*(k1 + (2*k2) + (2*k3) + k4)' + dist_state';

for i = 2:Nu 
    
k1 = hstep*rate_of_state_change(x(i-1,:),u(i));
k2 = hstep*rate_of_state_change(x(i-1,:)+(k1/2),u(i));
k3 = hstep*rate_of_state_change(x(i-1,:)+(k2/2),u(i));
k4 = hstep*rate_of_state_change(x(i-1,:)+k3,u(i));
x(i,:) = x(i-1,:) + (1/6)*(k1 + (2*k2) + (2*k3) + k4) + dist_state';
end

for i = Nu+1:Np
k1 = hstep*rate_of_state_change(x(i-1,:),u(Nu));
k2 = hstep*rate_of_state_change(x(i-1,:)+(k1/2),u(Nu));
k3 = hstep*rate_of_state_change(x(i-1,:)+(k2/2),u(Nu));
k4 = hstep*rate_of_state_change(x(i-1,:)+k3,u(Nu));
x(i,:) = x(i-1,:) + (1/6)*(k1 + (2*k2) + (2*k3) + k4) + dist_state';    
end

end


function dxdt = rate_of_state_change(x_estim, u)

%values of parameters
F = 100;     % Flow rate [L/min]
    V = 100;     % Volume [L]
    k0 = 7.2e10; % Pre-exponential factor [min^-1]
    E_R = 8750;  % Activation energy divided by R [K]
    delta_H = -5e4; % Heat of reaction [J/mol]
    rho = 1000;  % Density [g/L]
    Cp = 0.239;  % Heat capacity [J/g/K]
    U_A = 5e4;     % Heat transfer coefficient [J/min/K]
    Caf = 1.0;   % Feed concentration [mol/L]
    Tf = 350;    % Feed temperature [K]

%%% OUTPUT PARAMETERS
Ca = x_estim(1);
    T = x_estim(2);


%%% input parameters


Tc= u(1);

k = k0 * exp(-E_R / T);

    dCAdt = F/V * (Caf - Ca) - k * Ca;

    dTdt = F/V * (Tf - T) + (-delta_H / (rho * Cp)) * k * Ca + (U_A / (rho * Cp * V)) * (Tc - T);

    dxdt = [dCAdt; dTdt];
end
