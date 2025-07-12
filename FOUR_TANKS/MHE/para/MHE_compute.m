function [x_MHE, dist_state, dist_output] = MHE_compute(ys_meas, u_first_moves, N_MHE, y_meas)

global Nsteps
%This is the Moving Horizon Estimator Function.
%The length of the horizon is 5

hstep = 1;
Nsteps = 3;
Ts = 3;

global x_apriori x_aposteriori P Rx Rdx Ry dist_states;

% Initialization of the aposteriori state, covariance, and weighting matrices
if isempty(x_apriori)
    x_aposteriori = [4; 5; 21; 19];  % 4 states
    P = 1e-3*ones(4, 4);  % Update for 4 states
    dist_states = [0; 0; 0; 0];  % 4 state disturbances
    
    Rx = diag([1e-6; 1e-6; 1e-6; 1e-6]);  % Update for 4 states
    Ry = diag([1e-3; 1e-3]);  % 2 outputs
    Rdx = diag([1e-8; 1e-8; 1e-8; 1e-8]);  % Update for 4 states

    P = covariance_update(P, Rx, Ry, x_aposteriori, u_first_moves(end, :), Ts);
end

persistent d_prev
if isempty(d_prev)
    d_prev = zeros(4 * N_MHE, 1);  % Update for 4 states
end

% Predict the apriori state based on the aposteriori state
x_apriori = state_predictor(x_aposteriori, u_first_moves(end, :), 1, 1, hstep, dist_states)';

% Initial decision variables for the optimizer
dec_var_initial = [x_apriori; zeros(4 * N_MHE, 1)];
dec_var_LB = [0; 0; 0; 0; -10 * ones(4 * N_MHE, 1)];
dec_var_UB = [30; 30; 30; 30; 10 * ones(4 * N_MHE, 1)];

% Optimizer call using fmincon
[dec_var, cost_MHE] = fmincon(@(dec_var) MHE_cost_fun(dec_var, x_apriori, d_prev, ys_meas, u_first_moves, P, Rdx, Ry, hstep, N_MHE), ...
    dec_var_initial, [], [], [], [], dec_var_LB, dec_var_UB);

disp(cost_MHE);
d_prev = dec_var(5:end);

% Predict the state throughout the horizon
x_MHE = ones(N_MHE, 4);  % 4 states
x_past_init = dec_var(1:4);  % Initialize with first 4 states

j = 5;
for i = 1:N_MHE
    x_MHE(i, :) = state_predictor(x_past_init, u_first_moves(i, :), 1, 1, hstep, dec_var(j:j+3));
    j = j + 4;
    x_past_init = x_MHE(i, :)';
end

x_MHE = x_MHE(end, :)';

% Update aposteriori state
x_aposteriori = dec_var(1:4);

% Retrieve state disturbances from the decision vector
dist_states = dec_var((4 * N_MHE + 1):(4 * N_MHE + 4));
dist_state = dist_states;

% Calculate output disturbance
dist_output = ys_meas(N_MHE, :)' - x_MHE(1:2);

% Update covariance matrix
P = covariance_update(P, Rx, Ry, x_MHE, u_first_moves(end, :), Ts);

end

function cost = MHE_cost_fun(dec_var, x_apriori, d_prev, ys_meas, u_first_moves, P, Rdx, Ry, hstep, N_MHE)

j = 5;
x_pasts = ones(N_MHE, 4);  % 4 states
x_past_init = dec_var(1:4);

for i = 1:N_MHE
    x_pasts(i, :) = state_predictor(x_past_init, u_first_moves(i, :), 1, 1, hstep, dec_var(j:j+3));
    j = j + 4;
    x_past_init = x_pasts(i, :)';
end

y1_past = x_pasts(:, 1); 
y1_meas = ys_meas(:, 1);
y2_past = x_pasts(:, 2); 
y2_meas = ys_meas(:, 2);

% Stage cost based on the output prediction
PIy = Ry^(-1);
PIdx = Rdx^(-1);

stage_cost_1 = (y1_meas - y1_past)' * PIy(1, 1) * (y1_meas - y1_past) + ...
               (y2_meas - y2_past)' * PIy(2, 2) * (y2_meas - y2_past);

% Disturbance cost
dec_var_x1_dist = (dec_var(5:4:end) - d_prev(1:4:end));
dec_var_x2_dist = (dec_var(6:4:end) - d_prev(2:4:end));
dec_var_x3_dist = (dec_var(7:4:end) - d_prev(3:4:end));
dec_var_x4_dist = (dec_var(8:4:end) - d_prev(4:4:end));

stage_cost_2 = PIdx(1, 1) * (dec_var_x1_dist' * dec_var_x1_dist) + ...
               PIdx(2, 2) * (dec_var_x2_dist' * dec_var_x2_dist) + ...
               PIdx(3, 3) * (dec_var_x3_dist' * dec_var_x3_dist) + ...
               PIdx(4, 4) * (dec_var_x4_dist' * dec_var_x4_dist);

% Total stage cost
stage_cost = stage_cost_1 + stage_cost_2;

% Arrival cost
arrival_cost = (x_apriori - dec_var(1:4))' * (P^(-1)) * (x_apriori - dec_var(1:4));

% Calculate total cost
if isreal(dec_var) && isreal(x_pasts)
    cost = stage_cost + arrival_cost;
else
    cost = inf;
end

end

function P = covariance_update(P, Rx, Ry, x_MHE, u, Ts)

g = 981;  % Gravitational constant (m/s^2)
A1 = 35;  
A2 = 39;  
A3 = 35;  
A4 = 39;

a1 = 0.97;  
a2 = 0.74;  
a3 = 0.42;  
a4 = 0.426;

alpha1 = 0.8542;
alpha2 = 0.7722;
alpha3 = 0.7126;
alpha4 = 0.5149;

% States
h1 = x_MHE(1);
h2 = x_MHE(2);
h3 = x_MHE(3);
h4 = x_MHE(4);

% Inputs
q1 = u(1);
q2 = u(2);

% Jacobian matrix A (linearized state transition matrix)
A = [1 + Ts * (a1 * alpha1 * g / A1) * (2 * g * h1)^-0.5, 0, Ts * (a3 * alpha3 * g / A1) * (2 * g * h3)^-0.5, 0;
     0, 1 + Ts * (a2 * alpha2 * g / A2) * (2 * g * h2)^-0.5, 0, Ts * (a4 * alpha4 * g / A2) * (2 * g * h4)^-0.5;
     0, 0, 1 + Ts * (a3 * alpha3 * g / A3) * (2 * g * h3)^-0.5, 0;
     0, 0, 0, 1 + Ts * (a4 * alpha4 * g / A4) * (2 * g * h4)^-0.5];

C = [1, 0, 0, 0;
     0, 1, 0, 0];

P = Rx + A * P * A' - A * P * C' * (Ry + C * P * C')^(-1) * C * P * A';

end

function x = state_predictor(x_estim, u, Nu, Np, hstep, dist_state)

%disp(size(dist_state))

x(1, :) = RK4(x_estim, [u(1); u(Nu+1)], hstep) + dist_state';

for i = 2:Nu 
    x(i, :) = RK4(x(i-1, :), [u(i); u(Nu+i)], hstep) + dist_state';
end

for i = Nu+1:Np
    x(i, :) = RK4(x(i-1, :), [u(i); u(2*Nu)], hstep) + dist_state';
end

end

function [xpred] = RK4(xo, uo, hstep)

global Nsteps

h = hstep;
xk1 = xo;

for k = 1:Nsteps
    k1 = h * rate_of_state_change(xk1, uo);
    k2 = h * rate_of_state_change(xk1 + k1/2, uo);
    k3 = h * rate_of_state_change(xk1 + k2/2, uo);
    k4 = h * rate_of_state_change(xk1 + k3, uo);
    xk1 = xk1 + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
end

xpred = xk1';

end

function dhdt = rate_of_state_change(x_estim, u)

g = 981;  % Gravitational constant (m/s^2)
A1 = 35;  
A2 = 39;  
A3 = 35;  
A4 = 39;

a1 = 0.97;  
a2 = 0.74;  
a3 = 0.42;  
a4 = 0.426;

alpha1 = 0.8542*0.95;
alpha2 = 0.7722*0.95;
alpha3 = 0.7126*0.95;
alpha4 = 0.5149*0.95;

h1 = x_estim(1);
h2 = x_estim(2);
h3 = x_estim(3);
h4 = x_estim(4);

q1 = u(1);  
q2 = u(2);  

dh1dt = -(a1 * alpha1 / A1) * sqrt(2 * g * h1) + (a3 * alpha3 / A1) * sqrt(2 * g * h3) + (q1 / A1);
dh2dt = -(a2 * alpha2 / A2) * sqrt(2 * g * h2) + (a4 * alpha4 / A2) * sqrt(2 * g * h4) + (q2 / A2);
dh3dt = -(a3 * alpha3 / A3) * sqrt(2 * g * h3) + ((1 - alpha2) * q2) / A3;
dh4dt = -(a4 * alpha4 / A4) * sqrt(2 * g * h4) + ((1 - alpha1) * q1) / A4;

dhdt = [dh1dt; dh2dt; dh3dt; dh4dt];

end
