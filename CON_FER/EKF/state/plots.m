clc
load('NMPC_EKF_SD')
load('NMPC')

figure(1)
subplot(2,1,1)
plot(time, setpoint, 'r', time, predicted_nmpc, 'b-', time,predicted,'g-',  'LineWidth',2);
ylabel('Cell concentration, g/L-1','FontSize' , 8)
xlabel('Time [Seconds]','FontSize',8)
legend('Set point','NMPC','NMPC+EKF+SD')

subplot(2,1,2)
plot(time, input, 'b', time,input_nmpc,'k-',  'LineWidth',2);
ylabel('Dilution rate, h-1','FontSize' , 8)
xlabel('Time [Seconds]','FontSize',8)
legend('NMPC','NMPC+EKF+SD')

figure(2)
subplot(3,1,1)
plot(time, X_pred, 'r', time, X_est, 'k:',  'LineWidth',2);
title('Cell concentration [X]')
xlabel('Time [Seconds]','FontSize',8)
legend('Actual state','EstimatedNMPC+EKF+SD')

subplot(3,1,2)
plot(time, S_pred, 'r', time, S_est, 'k:',  'LineWidth',2);
title('Substrate con. [S]')
xlabel('Time [Seconds]','FontSize',8)
legend('Actual state','EstimatedNMPC+EKF+SD')

subplot(3,1,3)
plot(time, P_pred, 'r', time, P_est, 'k:',  'LineWidth',2);
title('Product conc. [P]')
xlabel('Time [Seconds]','FontSize',8)
legend('Actual state','EstimatedNMPC+EKF+SD')

figure(3)
subplot(1,1,1)
plot(time, state1, 'g', time, state2, 'r', time, state3, 'p', 'LineWidth',2);
title('State disturbance')
xlabel('Time [Seconds]','FontSize',8)
legend('state X','state S', 'State P')
