%% 🚀 How to Use

%% Step 1: Run the NMPC Initialization
%Open MATLAB and set the working directory to:
   % 1. CON_FER/MHE/state/ OR CON_FER/EKF/state/ does not matter which on to open first just maintain directory
 
%2. Open `NMPC_conFer.m` and run the script.

%3. You should see the following confirmation in the command window:
   
   %Model.StateFcn is OK.
   %Model.OutputFcn is OK.
   

%% Step 2: Run the Simulink Model
%1. Still in the same directory, open:
   
  % SIM_continuousFermenter.slx
   
   %This will launch Simulink.

%2. Change the **Stop Time** (top toolbar) to: 500
  

%3. Click **Run** to start the simulation.

%% Step 3: Check Results
%Once simulation completes:
%- Return to the MATLAB window.
 % - A **response graph** showing system behavior under disturbances
%- You can also view saved images: `estimated_state.jpg`, `state_dist.jpg`, `response.jpg`.
