clc, clear all

%% Simulink Model
% To run this example, Simulink(R) is required.
if ~mpcchecktoolboxinstalled('simulink')
    disp('Simulink is required to run this example.')
    return
end

%%initialize NMPC object and set properties
nx = 3;
ny = 1;


%'UD', [5,6,7,8,9,10]

nlobj = nlmpc(nx, ny, "MV", 1, 'UD', [2]);

%The prediction model sample time is the same as the controller sample time
Ts = 0.05; %Ts = 2
nlobj.Ts = Ts;

%Prediction and control horizon
nlobj.PredictionHorizon = 5; % Increased for better long-term prediction
nlobj.ControlHorizon = 1;

%weights on output and manipulated variables
nlobj.Weights.OutputVariables = 50; % Balanced weight for faster tracking %0.1
nlobj.Weights.ManipulatedVariables = 0.1; % Increased to reduce control effort penalty %0.05;
nlobj.Weights.ManipulatedVariablesRate = 0.01; % Reduced to improve responsiveness


%constraints on manipulated variables
nlobj.MV(1).RateMin = -10;

nlobj.MV(1).RateMax = 10;

nlobj.MV(1).Min = 280;

nlobj.MV(1).Max = 350;

%constraint on output variables to keep them within reasonable bounds
nlobj.OV(1).Min = 300;

nlobj.OV(1).Max = 360;


%% specify the nonlinear state and output functions
nlobj.Model.StateFcn = 'cstrStateFcnCT';
nlobj.Model.OutputFcn = 'cstrOutputFcn';


%% validate the nonlinear MPC object with initial conditions
x0 = [0.877;324.5;8750];

y0 = [0.877;324.5];  %[0.0033;920.86; 0.0104; 0.0085; 0.0118; 798];
u0 = 300;
validateFcns(nlobj, x0, u0);
