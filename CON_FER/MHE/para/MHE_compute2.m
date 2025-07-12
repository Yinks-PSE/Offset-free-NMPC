function [x_MHE, param_MHE, dist_output] = MHE_compute2(ys_meas, u_first_moves, N_MHE, y_meas)

%global Nsteps
%This is the Moving Horizon Estimator Function.
%The length of the horizon is 5

%Q = Rx
%R = Ry


%hstep = 0.5;
%Nsteps = 5;
Ts = 0.5;
dist_states = [0;0;0];

global x_apriori x_aposteriori params_MHE P Q R;

%Initialization of the aposteriori state, covariance and weighting matrices.

if isempty(x_apriori)
    x_aposteriori = [6; 5 ;19.14];
    params_MHE = 0.4;
    P = diag([1e-3;1e-3;1e-3]);
    Q = diag([1e2;1e2;1e2]);
    R = 1e10;
    Pw = covariance_update(P,Q,R,x_aposteriori,u_first_moves,params_MHE,Ts);
end

% Ensure Yxs is initialized
  %  if isempty(params_MHE)
        %params_MHE = 0.4; % Default value
    %end

%disp(eye(3,3)/P);

%This creates the apriori state using predictions based on the past
%aposteriori state.

x_apriori = state_predictor(x_aposteriori,u_first_moves(end,:),1,1,params_MHE)';

%Preparation of the optimizer, starting point, lower and upper bounds are
%defined.
dec_var_initial = [x_apriori; 0.4*ones(N_MHE,1)];
%dec_var_initial = [x_apriori; repmat(params_MHE,N_MHE,1)];
dec_var_LB = [0;0;0;0];
dec_var_UB = [20;60;60;1];
%dec_var_LB = [0;0;0;0*ones(N_MHE,1)];
%dec_var_UB = [20;60;60;1*ones(N_MHE,1)];

options = optimoptions('fmincon','Display','final','Algorithm','sqp');

%rng default % For reproducibility
%x0 = dec_var_initial;
%gs = GlobalSearch('NumTrialPoints',500);
%cost = @(dec_var)MHE_cost_fun(dec_var,x_apriori,ys_meas,u_first_moves,P,R,N_MHE);
%problem = createOptimProblem('fmincon','x0',x0,'objective',cost,'lb',dec_var_LB,'ub',dec_var_UB,'nonlcon',@(dec_var)nlcf(dec_var,N_MHE));
%[dec_var,cost_MHE] = run(gs,problem);
%options = optimoptions('fmincon','Display','final','Algorithm','interior-point');

%Note that the decision vector consists of the aposteriori state
%variables and the state disturbances throughout the prediction horizon.

[dec_var,cost_MHE] = fmincon(@(dec_var)MHE_cost_fun(dec_var,x_apriori,ys_meas,u_first_moves,P,R,N_MHE),dec_var_initial,[],[],[],[],dec_var_LB,dec_var_UB,[],options);

%Displays the value of the objective function in the command window.

disp(cost_MHE);
%disp('u_first_moves')
%disp(u_first_moves);

%Uses the aposteriori state gotten from optimization and the other decision
%variables to predict the current state.

%j = 4;
x_MHE = ones(N_MHE,3);
x_past_init = dec_var(1:3);
disp(x_past_init)
for i = 1:N_MHE
    x_MHE(i,:) = state_predictor(x_past_init,u_first_moves(i,:),1,1,dec_var(4));
    %j = j+1;
    x_past_init = x_MHE(i,:)';
end

x_MHE = x_MHE(end,:)';

%Aposteriori state is updated.

x_aposteriori = dec_var(1:3);

%Estimated parameters are updated

params_MHE = dec_var(end);
param_MHE = params_MHE;

disp('params_MHE');
disp(dec_var(4:end));

%Disturbance in outputs is calculated.

dist_output = ys_meas(N_MHE,:)' - x_MHE(1);

%Covariance and weighting matrices are updated.

Pw = covariance_update(P,Q,R,x_MHE,u_first_moves,params_MHE,Ts);

end

function cost = MHE_cost_fun(dec_var,x_apriori,ys_meas,u_first_moves,P,R,N_MHE)

%The objective function utilized by the moving horizon estimator is a sum of the
%stage cost and the arrival cost. The stage cost penalizes the difference between 
%the measured state and outputs and the estimated ones. 
%The arrival cost penalizes the difference between the apriori state and
%the aposteriori state.
%Note that the weights used for each penalty term are either the error noise covariance
% of the outputs and states or are derived from them in some cases.


%The next few lines calculates the state and output vector throughout the
% horizon of the estimator using the aposteriori state
%and the state disturbances. i.e. The decision vector.

%j = 4;
x_pasts = ones(N_MHE,3);
x_past_init = dec_var(1:3);
for i = 1:N_MHE
    x_pasts(i,:) = state_predictor(x_past_init,u_first_moves(i,:),1,1,dec_var(4));
    %j = j+1;
    x_past_init = x_pasts(i,:)';
end
 
y1_past = x_pasts(:,1); y1_meas = ys_meas(:,1);
%y2_past = x_pasts(:,2); y2_meas = ys_meas(:,2);

%First stage cost is a weighted sum of the differences between the measured
%output and those calculated using the decision vector.

%disp(size(R))

stage_cost_1 = (y1_meas - y1_past)'*R(1,1)*(y1_meas - y1_past);% + (y2_meas - y2_past)'*R(2,2)*(y2_meas - y2_past);

%Second stage cost is a weighted sum of the estimated parameter. However it
%is assigned a weight of zero.
 
dec_var_param_1 = ones(N_MHE,1);
%dec_var_param_2 = ones(N_MHE,1);

j = 1;
for i = 4:1:(N_MHE+3)
    dec_var_param_1(j,:) = dec_var(i);
   % dec_var_param_2(j,:) = dec_var(i+1);
    j = j+1;
end

stage_cost_2 = ((dec_var_param_1'*dec_var_param_1)); %+ (0*(dec_var_param_2'*dec_var_param_2));
                
stage_cost = stage_cost_1 + stage_cost_2;

arrival_cost = (x_apriori - dec_var(1:3))'*(eye(3,3)/P)*(x_apriori - dec_var(1:3));


%if eq==0
    cost = (stage_cost + arrival_cost);
%else
    %cost = inf;

%if isinteger(arrival_cost) && isinteger(stage_cost)
    
%else
%endt = 1e15;

end

function P_matrix = covariance_update(P_matrix,Q,R,x_MHE,u,param_MHE,Ts)

alpha=2.2;
beta=0.2;
miu_m=0.48;
Pm=50;
Km=1.2;
Ki=22;

D = u(1);
Yxs= param_MHE;

%X = x_MHE(1);
S = x_MHE(2);
P = x_MHE(3);


%%% assigning values to the state variables


miu=(miu_m*(1-P/Pm))*S/(Km+S+(S^2/Ki)) ;

%disp(Yxs);

%COMPUTATION OF JACOBI;ANS

%The Jacobians of the vector functions  f,b and h, F = J(f)
%and H = J(h) and B = J(b) are computed.

%Note that f is a vector function mapping the previous value of the
%state vector x_k-1 to its present value x_k.
%h is a vector function mapping the present value of the state vector x_k to
%the output (measured variable used for state estimation) y_k.
%b is a vector function mapping the inputs to the outputs.

a11 = 1 - Ts*(D - miu);
a12 = 0;
a13 = 0;

a21=-(Ts*miu)/Yxs;
a22=1 - D*Ts;
a23=0;

a31= -Ts*(beta + alpha*miu);
a32=0;
a33 = 1 - D*Ts;

%disp(a11);
%disp(a21);
%disp(a22);
%disp(a31);
%disp(a33);

  A = [a11,a12,a13;a21,a22,a23;a31,a32,a33];
  
 
h11 = 1;    h21 = 0;   
h12 = 0;    h22 = 1;
h13 = 0;    h23 = 0;


H = [h11,h12,h13; h21,h22,h23];

%B =[[        (Ts*g1)/A1,                 0]
    %[                 0,        (Ts*g2)/A2]
    %[                 0, -(Ts*(g2 - 1))/A3]
   %[ -(Ts*(g1 - 1))/A4,                 0]];

%P is updated according to the equations below.
V = R + H*Q*H';
I = eye(3,3);

P_matrix = A*P_matrix*A' + I*Q*I + A*Q*H'*(V^-1)*H*Q*A;
end


function x = state_predictor(x_estim,u,Nu,Np,Yxs)

%This function predicts values of states throughout the prediction horizon. i.e.
%From k+1 to Np+1 using Runge-Kutta discretization.
%This function requires x_estim to be a column vector.
h = 1;

k1 = h*rate_of_state_change(x_estim',u(1),Yxs);
k2 = h*rate_of_state_change(x_estim'+(k1/2),u(1),Yxs);
k3 = h*rate_of_state_change(x_estim'+(k2/2),u(1),Yxs);
k4 = h*rate_of_state_change(x_estim'+k3,u(1),Yxs);
x(1,:) = x_estim' + (1/6)*(k1 + (2*k2) + (2*k3) + k4);

for i = 2:Nu 
    
k1 = h*rate_of_state_change(x(i-1,:),u(i),Yxs);
k2 = h*rate_of_state_change(x(i-1,:)+(k1/2),u(i),Yxs);
k3 = h*rate_of_state_change(x(i-1,:)+(k2/2),u(i),Yxs);
k4 = h*rate_of_state_change(x(i-1,:)+k3,u(i),Yxs);
x(i,:) = x(i-1,:) + (1/6)*(k1 + (2*k2) + (2*k3) + k4);
end

for i = Nu+1:Np
k1 = h*rate_of_state_change(x(i-1,:),u(Nu),Yxs);
k2 = h*rate_of_state_change(x(i-1,:)+(k1/2),u(Nu),Yxs);
k3 = h*rate_of_state_change(x(i-1,:)+(k2/2),u(Nu),Yxs);
k4 = h*rate_of_state_change(x(i-1,:)+k3,u(Nu),Yxs);
x(i,:) = x(i-1,:) + (1/6)*(k1 + (2*k2) + (2*k3) + k4);    
end

end

function dxdt=rate_of_state_change(x_estim,u,Yxs)
%%%%%%%%%%%%%%%%%%%

 % Debugging
    

%values of parameters
%Yxs=0.4 ;
alpha=2.2  ;
beta=0.2  ;
miu_m=0.48   ;
Pm=50   ;
Km=1.2  ;
Ki=22  ;

%%% assigning values to the state variables
X = x_estim(1);
S = x_estim(2);
P = x_estim(3);

%%% input parameters
D  = u(1);
Sf = 20;

%%disp("Dimensions:");
%disp(["X:", size(X)]);
%disp(["S:", size(S)]);
%disp(["P:", size(P)]);
%disp(["D:", size(D)]);
%disp(["Sf:", size(Sf)]);
%disp(["Yxs:", size(Yxs)]);


miu=(miu_m*(1-(P/Pm)))*S/(Km+S+(S^2/Ki))  ;
dx1dt=-D*X + miu*X     ;
dx2dt=D*(Sf-S)- (1/Yxs)*miu*X   ;
dx3dt=-D*P+(alpha*miu+beta)*X   ;

dxdt=[dx1dt dx2dt dx3dt];
end
