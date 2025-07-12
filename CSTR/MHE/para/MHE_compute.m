function [Xobs,dist_output,ER_est,fcost,Ecount,exitflag]=MHE_compute(uo,um,ymeas...
                                        ,yplant,Ts,hstep,Nsteps,Rp,Rx,Ry...
                                        ,d_initial, d_MIN_colvector,d_MAX_colvector)

% DEFINING PERSISTENT VARIABLES

persistent PIo
if isempty(PIo)

PIo=1e-4*eye(2) ; 

end

persistent xest_prev
if isempty(xest_prev)
xest_prev = [0.877 ;324.5];
end

persistent p_bar
if isempty(p_bar)
p_bar = 8750 ;
end

% % OPTIONS = optimset('Algorithm','sqp','Display','final');
% OPTIONS = optimset('Algorithm','active-set');
% % OPTIONS = optimset('Algorithm','Interior-Point');
% 
% tic

OPTIONS = optimset('Algorithm','sqp','Display','final','TolFun', 1e-6);
%OPTIONS = optimset('Algorithm','active-set');
%OPTIONS = optimset('Algorithm','Interior-Point');

tic

% tic
% OPT1 = optimset('Algorithm','sqp');

%[dopt,fvalp] =fmincon(@(p) obj_est(um,hstep,Nsteps,p,ym,xest_prev,p_prev,PIo,Ry,Rp),[0.12466;0.74068;300],[],[],[],[],[0;0;100],[1;1;600],[],OPT1) 
try
    [dopt,fval,exitflag] = fmincon(@(d) obj_est(um,hstep,Nsteps,d,ymeas,xest_prev,p_bar,PIo,Rp,Ry), ...
        d_initial, [], [], [], [], d_MIN_colvector, d_MAX_colvector, [], OPTIONS);
catch ME
    disp('Error in fmincon:');
    disp(ME.message);
    disp('Stack trace:');
    disp(ME.stack);
    return;
end   ;

disp('dopt:'), disp(dopt);
if any(isnan(dopt))
    error('fmincon returned NaN in dopt.');
end

%disp(dopt)
Ecount = toc ;
fcost=fval ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
ER_est=[dopt(3)]  ;
E_R=ER_est;
p_bar = E_R ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE Xobs from estimated state
[Xobs]=state_cal(um,hstep,Nsteps,dopt)  ;
%xo=Xobs  ; 
N_MHE = 5 ;

yprediction = Xobs(2:2);


%[yprediction]=Euler(yplant,uo,E_R,hstep,Nsteps);

% function [yprediction]=ymodel(xo,uo,ko,hstep,Nsteps)
% 
% %[xk1]=Euler(xo,uo,ko,hstep,Nsteps);
% 
% xpred = Euler(xo,uo,ko,hstep,Nsteps)  ;
%     
% yprediction=xpred  ;

% [yprediction]=RK4(yplant,uo,ko,hstep,Nsteps);
dist_output = yplant - yprediction  ;

%PIo=co_var(xo,uo,hstep,PIo,ko,Rx,Ry)  ;
PIo=co_var(Xobs,uo,hstep,PIo,E_R,Rx,Ry)

%PIo=co_var(Xobs,uo,ko,Ts,PIo,Rx,Ry)     
%disp(PIo)

[Xprevious]=apriori_state(um,hstep,Nsteps,dopt)    ;
xest_prev=Xprevious  ;

end
%%% MHE-COMPUTE ENDS HERE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STARTING FROM HERE ARE THE INGREDIENTS FOR COMPUTING MHE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Xobs]=state_cal(um,hstep,Nsteps,dopt);

%xo=[dopt(1);dopt(2)]  ;
xo = [dopt(1);dopt(2)]  ;

disp('size xo in state_cal');disp(size(xo))
E_R = dopt(3)  ;

if ~isvector(xo) || length(xo) ~= 2
    error('Error: xo in state_cal is not a 2-element vector.');
end


%STEP 1
    
    xpred1=Euler(xo,um(:,1),E_R,hstep,Nsteps)   ;
    
 %STEP 2
    
    xpred2=Euler(xpred1,um(:,2),E_R,hstep,Nsteps)   ;
    %xpred2=xpred1+1/6*(k1_2+2*k2_2+2*k3_2+k4_2) ;
    
  %STEP 3
    
    xpred3=Euler(xpred2,um(:,3),E_R,hstep,Nsteps)   ;
    
    %STEP 4

    xpred4=Euler(xpred3,um(:,4),E_R,hstep,Nsteps)   ;
    
   %STEP 5
    
    xpred5=Euler(xpred4,um(:,5),E_R,hstep,Nsteps)   ;
    
    
    %STEP 6
    
    xpred6=Euler(xpred5,um(:,6),E_R,hstep,Nsteps)   ;
    
    
[Xobs]=xpred6 ; 
%disp(Xobs)


end
function PIk=co_var(ymeas,input,Ts,PIo,E_R,Rx,Ry)

%values of parameters
 F = 100;
    V = 100;
    k0 = 7.2e10;
    %E_R = 8750;
    DeltaH = -5e4;
    rho = 1000;
    Cp = 0.239;
    UA = 5e4;
    CAF = 1.0;
    T_f = 350;

    %%%%%%INPUT PARAMETER
Tc= input(1);
xhat=ymeas  ;

%Ts=1 ;



% calculate the Jacobians for state and measurement equations

a11 = 1 + Ts*(-F/V - k0 * exp(-E_R / xhat(2)));
a12 = Ts*(xhat(1) * k0 * (E_R / xhat(2)^2) * exp(-E_R / xhat(2)));

a21 = Ts*(-DeltaH * k0 * exp(-E_R / xhat(2))) / (rho * Cp);
a22 = 1 + Ts*(-F/V - (DeltaH * xhat(1) * k0 * E_R / (rho * Cp * xhat(2)^2)) * exp(-E_R / xhat(2)) - UA / (V * rho * Cp));

A = [a11,a12;
    a21,a22];

    C = [0, 1];  % Measurement matrix, assuming direct state observation

    G=[1, 0;
        0, 1]  ;
    

    %y + C*PIo*C'
%A*PIo*C'*(Ry + C*PIo*C')^(-1)*C*PIo*A' 
%A*PIo*A'
%G*Rx*G' 

PIk=G*Rx*G' + A*PIo*A' - A*PIo*C'*(Ry + C*PIo*C')^(-1)*C*PIo*A' ;

end

function [Xprevious]=apriori_state(um,hstep,Nsteps,dopt);


xo=[dopt(1);dopt(2)]  ;
E_R=dopt(3)  ;

%STEP 1
    
    xp=Euler(xo,um(:,1),E_R,hstep,Nsteps);
    xpred1= xp ;
    %xpred1=xo+1/6*(k1_1+2*k2_1+2*k3_1+k4_1) ;
    
    
Xprevious=xpred1  ;
end

% INGREDIENTS STOPS HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STARTING FROM HERE ARE FILES FOR FMINCON FOR MHE ESTIMATION

function er_fun=obj_est(um,hstep,Nsteps,d,ymeas,xest_prev,p_bar,PIo,Rp,Ry) 

[ypred]=est_fun(um,hstep,Nsteps,d);


ypredr=ypred(end,:)'  ;
ymr=ymeas'  ;

piy=Ry^(-1)  ;


Wx1=piy(1,1)  ;

Wp=Rp^(-1)  ;  


xbar=[(d(1)-xest_prev(1)) ;  (d(2)-xest_prev(2)) ];

er_fun=(ymr-ypredr)'*Wx1*(ymr-ypredr)  +  Wp*(d(3)- p_bar)^2  + xbar'*PIo^(-1)*xbar   ;
%er_fun=(ym1-y1)'*Wx1*(ym1-y1) + (ym2-y2)'*Wx2*(ym2-y2)  +  Wp*(p(3)- p_prev)^2  + xbar'*PIo^(-1)*xbar   ;



%er_fun=(ym1-y1)'*Wx1*(ym1-y1) + (ym2-y2)'*Wx2*(ym2-y2)  +  Wp*p(3)^2  + xbar'*PI*xbar ;
%er_fun=(ym1-y1)'*Wx1*(ym1-y1) + (ym2-y2)'*Wx2*(ym2-y2)  +  Wx*xd 
%er_fun=(ym1-y1)'*Wx1*(ym1-y1) + (ym2-y2)'*Wx2*(ym2-y2)   ;







end

function [ypred]=est_fun(um,hstep,Nsteps,d);

xo=[d(1);d(2)]  ;

disp('size xo in est_fun');disp(size(xo))

E_R=d(3)   ;
if ~isvector(xo) || length(xo) ~= 2
    error('Error: xo in est_fun is not a 2-element vector.');
end

%STEP 1
    
    xpred1=Euler(xo,um(:,1),E_R,hstep,Nsteps)   ;
    
 %STEP 2
    
    xpred2=Euler(xpred1,um(:,2),E_R,hstep,Nsteps)   ;
    %xpred2=xpred1+1/6*(k1_2+2*k2_2+2*k3_2+k4_2) ;
    
  %STEP 3
    
    xpred3=Euler(xpred2,um(:,3),E_R,hstep,Nsteps)   ;
    
    %STEP 4

    xpred4=Euler(xpred3,um(:,4),E_R,hstep,Nsteps)   ;
    
   %STEP 5
   
    xpred5=Euler(xpred4,um(:,5),E_R,hstep,Nsteps)   ;
    
    
    %STEP 6
    
    xpred6=Euler(xpred5,um(:,6),E_R,hstep,Nsteps)   ;
    
  

%ypred=[xpred1 +  [d(1);d(2)] ,...
%       xpred2 +  [d(3);d(4)] ,...
 %      xpred3 +  [d(5);d(6)] ,...
 %      xpred4 +  [d(7);d(8)] ,...
 %      xpred5 +  [d(9);d(10)] ,...
  %     xpred6 +  [d(11);d(12)] ]  ;
   
    
    
    
%ypred=[xpred1(1);xpred2(1);xpred3(1);xpred4(1);xpred5(1);xpred6(1)]  
%ypred=[xpred1, xpred2, xpred3, xpred4, xpred5, xpred6] ; 
ypred=[xpred1, xpred2, xpred3, xpred4,xpred5,xpred6] ;
end

function [xk1]=Euler(xo,uo,E_R,hstep,Nsteps);

disp('size xo in Euler');disp(size(xo));
disp('length xo in Euler');disp(length(xo));

if length(xo) ~= 2
    error('Error: xo in Euler is not a 2-element vector.');
end


%xk=RK4(xo,uo,du,h,Np,Nu,ko,dx1,dx2);
%disp(xo)
%delta=Ts/Nsteps  ;
xk1=xo ;
    
for k=1:Nsteps
    
   % xk1=xk1 + hstep*odefun(xk1,uo,ko);  %u(:,k)
    k1=hstep*odefun(xk1,uo,E_R);
    k2=hstep*odefun(xk1+k1/2,uo,E_R);
    k3=hstep*odefun(xk1+k2/2,uo,E_R);
    k4=hstep*odefun(xk1+k3,uo,E_R);
    
    xk1 = xk1 + 1/6*(k1+2*k2+2*k3+k4);
   
   
   
end




end


 


function xdot=odefun(x,u,E_R);



%%%%%%%%%%%%%%%%%%%
%disp(x)
%values of parameters
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

%%% OUTPUT PARAMETERS
Ca = x(1);
    T = x(2);


%%% input parameters


Tc= u(1);

k = k0 * exp(-E_R / T);

    dCAdt = F/V * (Caf - Ca) - k * Ca;

    dTdt = F/V * (Tf - T) + (-delta_H / (rho * Cp)) * k * Ca + (U_A / (rho * Cp * V)) * (Tc - T);

    xdot = [dCAdt; dTdt];





end
