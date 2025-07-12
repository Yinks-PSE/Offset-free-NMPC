function dxdt = conFerStateFcnCT(x, u)


%definig the parameters
 Ki = 22; %Substrate inhibition constatnt
 Km= 1.2;     % substrate saturation constant
    Pm = 50;     % product saturation constant
    Sf = 20; % substrate concentration in the feed
    %Yxs = 0.4; %cell-mass yeild
    alpha = 2.2; %kinematic parameter
    beta = 0.2; %kinemtic prameter
    miu_m = 0.48; %the maximum growth rate
    Yxs = x(4);

%%% state parameters
X = x(1);
S = x(2);
P = x(3);

%%% input parameters
D = u(1);



dxdt = zeros(4,1);


miu = (miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki))  ;

dxdt(1) = -D * X + miu * X  ;

dxdt(2) = D * (Sf - S) - (1/Yxs) * miu * X ;

dxdt(3) = -D * P + (alpha * miu + beta) * X;

dxdt(4) = u(2);

 
