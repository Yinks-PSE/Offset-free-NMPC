function dxdt = cstrStateFcnCT(x, u)


%definig the parameters
 F = 100;     % Flow rate [L/min]
    V = 100;     % Volume [L]
    k0 = 7.2e10; % Pre-exponential factor [min^-1]
    %E_R = 8750;  % Activation energy divided by R [K]
    delta_H = -5e4; % Heat of reaction [J/mol]
    rho = 1000;  % Density [g/L]
    Cp = 0.239;  % Heat capacity [J/g/K]
    U_A = 5e4;     % Heat transfer coefficient [J/min/K]
    Caf = 1.0;   % Feed concentration [mol/L]
    Tf = 350;    % Feed temperature [K]
    E_R = x(3);


%%% state parameters
Ca = x(1);
T = x(2);

%%% input parameters
Tc = u(1);


dxdt = zeros(3,1);

% Rate constant (Arrhenius equation)
        k = k0 * exp(-E_R / T);
        
        % Mass balance for Ca
        dxdt(1) = F/V * (Caf - Ca) - k * Ca;
        
        % Energy balance for T
        dxdt(2) = F/V * (Tf - T) + (-delta_H / (rho * Cp)) * k * Ca + (U_A / (rho * Cp * V)) * (Tc - T);

        dxdt(3) = u(2);


 
