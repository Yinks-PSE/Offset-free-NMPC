signals_data = out.T.Data;
state_est = out.est.Data;
E_R = out.E_R.Data;

input = out.Tc.Data;

%ref = out.r_save.Data(:,1);
time = out.t_save.time;

setpoint = signals_data(:, 1);
predicted = signals_data(:, 2);

actual = state_est(:, 1);
ekf = state_est(:, 2);

E_R_sp = E_R(:, 1);
E_R_pred = E_R(:, 2);

save NMPC_EKF_SD.mat setpoint predicted input actual ekf E_R_sp E_R_pred time


%nom stands for nominal
%predicted_nom = out.T_nom.Data(:, 2);
%output_nom = out.Tc_nom.Data;


