H1 = out.H1.Data;
H2 = out.H2.Data;
h1 = out.h1.Data;
h2 = out.h2.Data;
h3 = out.h3.Data;
h4 = out.h4.Data;


state = out.state.Data;
state = squeeze(state);

input = out.V.Data;

input = squeeze(input);


time = out.t_save.time;

h1_pred = h1(:, 1); 
h1_est = h1(:, 2);

h2_pred = h2(:, 1);
h2_est = h2(:, 2);

h3_pred = h3(:, 1);
h3_est = h3(:, 2);

h4_pred = h4(:, 1);
h4_est = h4(:, 2);

H1_sp = H1(:, 1);
H1_pred = H1(:, 2);

H2_sp = H2(:, 1);
H2_pred = H2(:, 2);

V1 = input(1, :);
V2 = input(2, :);

state1 = state(1, :);
state2 = state(2, :);
state3 = state(3, :);
state4 = state(4, :);

save NMPC_EKF_SD.mat H1_sp H1_pred H2_sp H2_pred V1 V2 h1_pred h1_est h2_pred h2_est h3_pred h3_est h4_pred h4_est state1 state2 state3 state4 time



