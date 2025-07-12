signals_data = out.T.Data;
state_est = out.est.Data;
para_est = out.state.Data;

input = out.Tc.Data;


time = out.t_save.time;

setpoint = signals_data(:, 1);
predicted = signals_data(:, 2);

para_actual = para_est(:, 1);
para_mhe = para_est(:, 2);

actual = state_est(:, 1);
mhe = state_est(:, 2);

save NMPC_EKF_SD.mat setpoint predicted input actual mhe para_actual para_mhe time



