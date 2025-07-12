function [x_MHE, dist_state, dist_output] = MHE_compute2(ys_meas, u_first_moves, N_MHE, y_meas)

global Nsteps
%This is the Moving Horizon Estimator Function.
%The length of the horizon is 5

hstep = 0.5;
Nsteps = 3;
Ts = 0.1;

global x_apriori x_aposteriori P Rx Rdx Ry dist_states;

% Initialization of the aposteriori state, covariance, and weighting matrices
if isempty(x_apriori)
    x_aposteriori = [6; 5; 19.14];  % 4 states
    P = diag([1e-7;1e-7;1e-7]);  % Update for 4 states
    dist_states = [0; 0; 0];  % 4 state disturbances
    
    Rx = diag([1e-3; 1e-3; 5e-1]);  % Update for 4 states
    Ry = diag(1);  % 2 outputs
    Rdx = diag([1e-2; 1e-2; 1e-2]);  % Update for 4 states

    P = covariance_update(P, Rx, Ry, x_aposteriori, u_first_moves(end, :), Ts);
end

persistent d_prev
if isempty(d_prev)
    d_prev = zeros(3 * N_MHE, 1);  % Update for 4 states
end

% Predict the apriori state based on the aposteriori state
x_apriori = state_predictor(x_aposteriori, u_first_moves(end, :), 1, 1, hstep, dist_states)';

% Initial decision variables for the optimizer
dec_var_initial = [x_apriori; zeros(3 * N_MHE, 1)];
dec_var_LB = [0; 0; 0; -3 * ones(3 * N_MHE, 1)];
dec_var_UB = [20; 60; 60; 3 * ones(3 * N_MHE, 1)];

% Optimizer call using fmincon
[dec_var, cost_MHE] = fmincon(@(dec_var) MHE_cost_fun(dec_var, x_apriori, d_prev, ys_meas, u_first_moves, P, Rdx, Ry, hstep, N_MHE), ...
    dec_var_initial, [], [], [], [], dec_var_LB, dec_var_UB);

disp(cost_MHE);
d_prev = dec_var(4:end);

% Predict the state throughout the horizon
x_MHE = ones(N_MHE, 3);  % 4 states
x_past_init = dec_var(1:3);  % Initialize with first 4 states

j = 4;
for i = 1:N_MHE
    x_MHE(i, :) = state_predictor(x_past_init, u_first_moves(i, :), 1, 1, hstep, dec_var(j:j+2));
    j = j + 3;
    x_past_init = x_MHE(i, :)';
end

x_MHE = x_MHE(end, :)';

% Update aposteriori state
x_aposteriori = dec_var(1:3);

% Retrieve state disturbances from the decision vector
dist_states = dec_var((3 * N_MHE + 1):(3 * N_MHE + 3));
dist_state = dist_states;

% Calculate output disturbance
dist_output = ys_meas(N_MHE, :)' - x_MHE(1);

% Update covariance matrix
P = covariance_update(P, Rx, Ry, x_MHE, u_first_moves(end, :), Ts);

end

function cost = MHE_cost_fun(dec_var, x_apriori, d_prev, ys_meas, u_first_moves, P, Rdx, Ry, hstep, N_MHE)

j = 4;
x_pasts = ones(N_MHE, 3);  % 4 states
x_past_init = dec_var(1:3);

for i = 1:N_MHE
    x_pasts(i, :) = state_predictor(x_past_init, u_first_moves(i, :), 1, 1, hstep, dec_var(j:j+2));
    j = j + 3;
    x_past_init = x_pasts(i, :)';
end

y1_past = x_pasts(:, 1); 
y1_meas = ys_meas(:, 1);

% Stage cost based on the output prediction
PIy = Ry^(-1);
PIdx = Rdx^(-1);

stage_cost_1 = (y1_meas - y1_past)' * PIy(1, 1) * (y1_meas - y1_past);

% Disturbance cost
dec_var_x1_dist = (dec_var(4:3:end) - d_prev(1:3:end));
dec_var_x2_dist = (dec_var(5:3:end) - d_prev(2:3:end));
dec_var_x3_dist = (dec_var(6:3:end) - d_prev(3:3:end));

stage_cost_2 = PIdx(1, 1) * (dec_var_x1_dist' * dec_var_x1_dist) + ...
               PIdx(2, 2) * (dec_var_x2_dist' * dec_var_x2_dist) + ...
               PIdx(3, 3) * (dec_var_x3_dist' * dec_var_x3_dist) ;

% Total stage cost
stage_cost = stage_cost_1 + stage_cost_2;

% Arrival cost
arrival_cost = (x_apriori - dec_var(1:3))' * (P^(-1)) * (x_apriori - dec_var(1:3));

% Calculate total cost
if isreal(dec_var) && isreal(x_pasts)
    cost = stage_cost + arrival_cost;
else
    cost = inf;
end

end

function P_matrix = covariance_update(P_matrix, Rx, Ry, x_MHE, u, Ts)

%definig the parameters
 Ki = 22; %Substrate inhibition constatnt
 Km= 1.2;     % substrate saturation constant
    Pm = 50;     % product saturation constant
    Sf = 20; % substrate concentration in the feed
    Yxs = 0.4; %cell-mass yeild
    alpha = 2.2; %kinematic parameter
    beta = 0.2; %kinemtic prameter
    miu_m = 0.48; %the maximum growth rate

% States
X = x_MHE(1);
S = x_MHE(2);
P = x_MHE(3);

% Inputs
D = u(1);

% Jacobian matrix A (linearized state transition matrix)

a11 = 1 - Ts*(D - (miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki)));
a12 = 0;
a13 = 0;

a21=-(Ts*(miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki)))/Yxs;
a22=1 - D*Ts;
a23=0;

a31= -Ts*(beta + alpha*(miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki)));
a32=0;
a33 = 1 - D*Ts;

A = [a11,a12,a13;a21,a22,a23;a31,a32,a33];


C = [1, 0, 0];

G = [1,0,0;
    0,1,0;
    0,0,1];

P_matrix = G*Rx*G' + A*P_matrix*A' - A*P_matrix*C'*(Ry + C*P_matrix*C')^(-1)*C*P_matrix*A';

end

function x = state_predictor(x_estim, u, Nu, Np, hstep, dist_state)

%disp(size(dist_state))

x(1, :) = RK4(x_estim, u(1), hstep) + dist_state';

for i = 2:Nu 
    x(i, :) = RK4(x(i-1, :), u(i), hstep) + dist_state';
end

for i = Nu+1:Np
    x(i, :) = RK4(x(i-1, :), u(i), hstep) + dist_state';
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

function dxdt = rate_of_state_change(x_estim, u)

%definig the parameters
 Ki = 22; %Substrate inhibition constatnt
 Km= 1.2;     % substrate saturation constant
    Pm = 50;     % product saturation constant
    Sf = 20; % substrate concentration in the feed
    Yxs = 0.4; %cell-mass yeild
    alpha = 2.2; %kinematic parameter
    beta = 0.2; %kinemtic prameter
    miu_m = 0.48; %the maximum growth rate


X = x_estim(1);
S = x_estim(2);
P = x_estim(3);

D = u(1);   


miu = (miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki))  ;

dx1dt = -D * X + miu * X   ;

dx2dt = D * (Sf - S) - (1/Yxs) * miu * X ;

dx3dt = -D * P + (alpha * miu + beta) * X  ;

dxdt = [dx1dt; dx2dt; dx3dt];

end
