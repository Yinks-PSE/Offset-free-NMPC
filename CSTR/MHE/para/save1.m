signals_data = out.T_nom.Data;
input_nmpc = out.Tc_nom.Data;

time = out.t_save.time;

predicted_nmpc = signals_data(:, 2);

save NMPC.mat predicted_nmpc input_nmpc time