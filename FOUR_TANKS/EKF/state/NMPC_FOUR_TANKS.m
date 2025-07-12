clc, clear all

%% Simulink Model
% To run this example, Simulink(R) is required.
if ~mpcchecktoolboxinstalled('simulink')
    disp('Simulink is required to run this example.')
    return
end


%%initialize NMPC object and set properties
nx = 8;
ny = 2;
%nu = 2;

%'UD', [5,6,7,8,9,10]

nlobj = nlmpc(nx, ny, 'MV', [1,2], 'UD', [3,4,5,6]);

%The prediction model sample time is the same as the controller sample time
Ts = 3; %Ts = 2
nlobj.Ts = Ts;

%Prediction and control horizon
nlobj.PredictionHorizon = 10;
nlobj.ControlHorizon = 2;

%weights on output and manipulated variables
nlobj.Weights.OutputVariables =  [1,1]; %0.1
nlobj.Weights.ManipulatedVariables = [0,0]; %0.05;
nlobj.Weights.ManipulatedVariablesRate = [0.1,0.1];

%constraints on manipulated variables
nlobj.MV(1).RateMin = -10;
nlobj.MV(2).RateMin = -10;

nlobj.MV(1).RateMax = 10;
nlobj.MV(2).RateMax = 10;

nlobj.MV(1).Min = 15;
nlobj.MV(2).Min = 10;

nlobj.MV(1).Max = 150;
nlobj.MV(2).Max = 98;

%constraint on output variables to keep them within reasonable bounds
nlobj.OV(1).Min = 0;
nlobj.OV(2).Min = 0;
%nlobj.OV(3).Min = 2.5;
%nlobj.OV(4).Min = 2.5;


nlobj.OV(1).Max = 32;
nlobj.OV(2).Max = 32;
%nlobj.OV(3).Max = 8;
%nlobj.OV(4).Max = 8;


%% specify the nonlinear state and output functions
nlobj.Model.StateFcn = 'four_tankStateFcnCT';
nlobj.Model.OutputFcn = 'four_tankOutputFcn';


%SolverOptions.Algorithm = 'sqp';


%% validate the nonlinear MPC object with initial conditions
h0 = [4;5;21;19;0;0;0;0];  %[0.0033;920.86; 0.0104; 0.0085; 0.0118; 798];
u0 = [55;75];
validateFcns(nlobj, h0, u0);