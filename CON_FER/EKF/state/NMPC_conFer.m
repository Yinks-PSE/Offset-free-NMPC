clc, clear all

%% Simulink Model
% To run this example, Simulink(R) is required.
if ~mpcchecktoolboxinstalled('simulink')
    disp('Simulink is required to run this example.')
    return
end

%%initialize NMPC object and set properties
nx = 6;
ny = 1;


%'UD', [5,6,7,8,9,10]

nlobj = nlmpc(nx, ny, "MV", 1, 'UD', [2,3,4]);

%The prediction model sample time is the same as the controller sample time
Ts = 0.1; %Ts = 2
nlobj.Ts = Ts;

%Prediction and control horizon
nlobj.PredictionHorizon = 20; % Increased for better long-term prediction
nlobj.ControlHorizon = 2;

%weights on output and manipulated variables
nlobj.Weights.OutputVariables = 10; % Balanced weight for faster tracking %0.1
nlobj.Weights.ManipulatedVariables = 0; % Increased to reduce control effort penalty %0.05;
nlobj.Weights.ManipulatedVariablesRate = 0.1; % Reduced to improve responsiveness


%constraints on manipulated variables
nlobj.MV(1).RateMin = -0.1;

nlobj.MV(1).RateMax = 0.1;

nlobj.MV(1).Min = 0;

nlobj.MV(1).Max = 5;

%constraint on output variables to keep them within reasonable bounds
nlobj.OV(1).Min = 0;

nlobj.OV(1).Max = 20;


%% specify the nonlinear state and output functions
nlobj.Model.StateFcn = 'conFerStateFcnCT';
nlobj.Model.OutputFcn = 'conFerOutputFcn';


%% validate the nonlinear MPC object with initial conditions
x0 = [6.0;5.0;19.14;0;0;0];

y0 = 6.0;  %[0.0033;920.86; 0.0104; 0.0085; 0.0118; 798];
u0 = 0.202;
validateFcns(nlobj, x0, u0);
