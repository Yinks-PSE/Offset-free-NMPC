signals_data = out.T.Data;
state_est = out.est.Data;
state_dist = out.state.Data;

input = out.Tc.Data;

ref = out.r_save.Data(:,1);
time = out.t_save.time;

setpoint = signals_data(:, 1);
predicted = signals_data(:, 2);

state1 = state_dist(:, 1);
state2 = state_dist(:, 2);

actual = state_est(:, 1);
ekf = state_est(:, 2);

save NMPC_EKF_SD.mat setpoint predicted input actual ekf state1 state2 time


%nom stands for nominal
%predicted_nom = out.T_nom.Data(:, 2);
%output_nom = out.Tc_nom.Data;


