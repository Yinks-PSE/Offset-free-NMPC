signals_data = out.X.Data;
X = out.x.Data;
S = out.S.Data;
P = out.P.Data;


Yxs = out.Yxs.Data;

input = out.D.Data;


time = out.t_save.time;

X_pred = X(:, 1); 
X_est = X(:, 2);

S_pred = S(:, 1);
S_est = S(:, 2);

P_pred = P(:, 1);
P_est = P(:, 2);

setpoint = signals_data(:, 1);
predicted = signals_data(:, 2);

actual = Yxs(:, 1);
ekf = Yxs(:, 2);

save NMPC_EKF_SD.mat setpoint predicted input X_est X_pred S_est S_pred P_est P_pred actual ekf time



