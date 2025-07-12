H1 = out.H1_nom.Data;
H2 = out.H2_nom.Data;
input_nmpc = out.V_nom.Data;

time = out.t_save.time;


H1_nom_pred = H1(:, 2);
H2_nom_pred = H2(:, 2);

V1_nom = input_nmpc(:, 1);
V2_nom = input_nmpc(:, 2);

save NMPC.mat V1_nom V2_nom H1_nom_pred H2_nom_pred time