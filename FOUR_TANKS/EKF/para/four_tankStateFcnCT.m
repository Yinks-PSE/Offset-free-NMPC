function dhdt = four_tankStateFcnCT(h, u)
    % Define system parameters

    %disp(h); disp(u)
    g = 981;  % gravitational constant (m/s^2)
    
    % Tank areas (in cm^2)
    A1 = 35;  
    A2 = 39;  
    A3 = 35;  
    A4 = 39;  
    
    % Outlet areas (in cm^2)
    a1 = 0.97;  
    a3 = 0.42;  
    a2 = 0.74;  
    a4 = 0.426;  

    alpha1 = 0.8542 * 0.95;
    alpha3 = 0.7126 * 0.95;
    alpha2 = 0.7722 * 0.95;
    alpha4 = 0.5149 * 0.95;

    
    % Flow distribution coefficients
    gamma1 = 0.23;  
    gamma2 = 0.19;  
    
    % State variables (water levels in each tank)
    h1 = h(1);  
    h2 = h(2);  
    h3 = h(3);  
    h4 = h(4);  
    
    % Inputs (pump flow rates)
    q1 = u(1);  
    q2 = u(2);  
    %alpha1 = h(5);

    dh1 = h(5);
    dh2 = h(6);
    dh3 = h(7);
    dh4 = h(8);

    dhdt = zeros(8,1);
    % Equations from the journal
    dhdt(1) = -(a1*alpha1/A1) * sqrt(2 * g * h1) + (a3*alpha3/A1) * sqrt(2 * g * h3) + (gamma1 * q1) / A1 +dh1;
    dhdt(2) = -(a2*alpha2/A2) * sqrt(2 * g * h2) + (a4*alpha4/A2) * sqrt(2 * g * h4) + (gamma2 * q2) / A2 +dh2;
    dhdt(3) = -(a3*alpha3/A3) * sqrt(2 * g * h3) + ((1 - gamma2) * q2) / A3 +dh3;
    dhdt(4) = -(a4*alpha4/A4) * sqrt(2 * g * h4) + ((1 - gamma1) * q1) / A4 + dh4;
    
    %dhdt(5) = u(3);
    dhdt(5) = u(3);
    dhdt(6) = u(4);
    dhdt(7) = u(5);
    dhdt(8) = u(6);

