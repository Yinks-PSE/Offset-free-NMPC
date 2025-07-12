clc
load('NMPC_EKF_SD')
load('NMPC')

figure(1)
subplot(2,2,1)
plot(time, H1_sp, 'r', time, H1_nom_pred, 'b-', time,H1_pred,'g-',  'LineWidth',2);
ylabel('Height H1, [cm]','FontSize' , 8)
xlim([0,1200])
xlabel('Time [Seconds]','FontSize',8)
legend('Set point','NMPC','NMPC+EKF+SD')

subplot(2,2,2)
plot(time, H2_sp, 'r', time, H2_nom_pred, 'b-', time,H2_pred,'g-',  'LineWidth',2);
ylabel('Height H2, [cm]','FontSize' , 8)
xlim([0,1200])
xlabel('Time [Seconds]','FontSize',8)
legend('Set point','NMPC','NMPC+EKF+SD')

subplot(2,2,3)
plot(time, V1, 'b', time, V1_nom,'k-',  'LineWidth',2);
ylabel('Flow rate q1, [cm^3/s]','FontSize' , 8)
xlim([0,1200])
xlabel('Time [Seconds]','FontSize',8)
legend('NMPC','NMPC+EKF+SD')

subplot(2,2,4)
plot(time, V2, 'b', time, V2_nom,'k-',  'LineWidth',2);
ylabel('Flow rate q2, [cm^3/s]','FontSize' , 8)
xlim([0,1200])
xlabel('Time [Seconds]','FontSize',8)
legend('NMPC','NMPC+EKF+SD')

figure(2)
subplot(4,1,1)
plot(time, h1_pred, 'r', time, h1_est, 'k:',  'LineWidth',2);
title('Height, H1 [cm]')
xlabel('Time [Seconds]','FontSize',8)
legend('Actual state','EstimatedNMPC+EKF+SD')

subplot(4,1,2)
plot(time, h2_pred, 'r', time, h2_est, 'k:',  'LineWidth',2);
title('Height, H2 [cm]')
xlabel('Time [Seconds]','FontSize',8)
legend('Actual state','EstimatedNMPC+EKF+SD')

subplot(4,1,3)
plot(time, h3_pred, 'r', time, h3_est, 'k:',  'LineWidth',2);
title('Height, H3 [cm]')
xlabel('Time [Seconds]','FontSize',8)
legend('Actual state','EstimatedNMPC+EKF+SD')

subplot(4,1,4)
plot(time, h4_pred, 'r', time, h4_est, 'k:',  'LineWidth',2);
title('Height, H4 [cm]')
xlabel('Time [Seconds]','FontSize',8)
legend('Actual state','EstimatedNMPC+EKF+SD')

figure(3)
subplot(1,1,1)
plot(time, state1, 'g', time, state2, 'r', time, state3,'y', time, state4, 'b', 'LineWidth',2);
title('State disturbance')
xlabel('Time [Seconds]','FontSize',8)
legend('h1','h2', 'h3', 'h4')
