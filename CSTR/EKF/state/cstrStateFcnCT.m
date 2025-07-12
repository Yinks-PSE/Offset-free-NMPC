function dydt = cstrStateFcnCT(y, u)


%definig the parameters
 F = 100;     % Flow rate [L/min]
    V = 100;     % Volume [L]
    k0 = 7.2e10; % Pre-exponential factor [min^-1]
    E_R = 8750;  % Activation energy divided by R [K]
    delta_H = 5e4; % Heat of reaction [J/mol]
    rho = 1000;  % Density [g/L]
    Cp = 0.239;  % Heat capacity [J/g/K]
    U_A = 5e4;     % Heat transfer coefficient [J/min/K]
    Caf = 1.0;   % Feed concentration [mol/L]
    Tf = 350;    % Feed temperature [K]


%%% state parameters
Ca = y(1);
T = y(2);

%%% input parameters
Tc = u(1);

dy1 = y(3);
dy2 = y(4);

dydt = zeros(4,1);

% Rate constant (Arrhenius equation)
        k = k0 * exp(-E_R / T);
        
        % Mass balance for Ca
        dydt(1) = F/V * (Caf - Ca) - k * Ca + dy1;
        
        % Energy balance for T
        dydt(2) = F/V * (Tf - T) + (delta_H / (rho * Cp)) * k * Ca + (U_A / (rho * Cp * V)) * (Tc - T) + dy2;

        dydt(3) = u(2);

        dydt(4) = u(3);


 
