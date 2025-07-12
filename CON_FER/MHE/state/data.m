signals_data = out.X.Data;
X = out.x.Data;
S = out.S.Data;
P = out.P.Data;


state = out.state.Data;

input = out.D.Data;

input = squeeze(input);


time = out.t_save.time;

X_pred = X(:, 1); 
X_est = X(:, 2);

S_pred = S(:, 1);
S_est = S(:, 2);

P_pred = P(:, 1);
P_est = P(:, 2);

setpoint = signals_data(:, 1);
predicted = signals_data(:, 2);

state1 = state(:, 1);
state2 = state(:, 2);
state3 = state(:, 3);

save NMPC_EKF_SD.mat setpoint predicted input X_est X_pred S_est S_pred P_est P_pred state1 state2 state3 time



