clc
load('NMPC_EKF_SD')
load('NMPC')

figure(1)
subplot(2,1,1)
plot(time, setpoint, 'r', time, predicted_nmpc, 'b-', time,predicted,'g-',  'LineWidth',2);
ylabel('Temperature in the reactor, T(K)','FontSize' , 8)
xlabel('Time [Seconds]','FontSize',8)
legend('Set point','NMPC','NMPC+MHE+SD')

subplot(2,1,2)
plot(time, input, 'b', time,input_nmpc,'k-',  'LineWidth',2);
ylabel('Temperature of the coolant stream, T(K)','FontSize' , 8)
xlabel('Time [Seconds]','FontSize',8)
legend('NMPC','NMPC+MHE+SD')

figure(2)
subplot(1,1,1)
plot(time, actual, 'r', time, mhe, 'k:',  'LineWidth',2);
title('Estimated state')
xlabel('Time [Seconds]','FontSize',8)
legend('Actual state','EstimatedNMPC+MHE+SD')

figure(3)
subplot(1,1,1)
plot(time, para_actual, 'g', time, para_mhe, 'r',  'LineWidth',2);
title('Estimated parameter')
xlabel('Time [Seconds]','FontSize',8)
legend('state Ca','state T')
