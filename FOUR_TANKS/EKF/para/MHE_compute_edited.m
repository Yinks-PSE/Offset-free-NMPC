function [Xobs,fcost,Ecount,exitflag]=MHE_compute_edited(uo,um,ymeas,yplant,Ts, hstep,Nsteps,Rp,Rx,Ry...
    ,d_initial,d_MIN_colvector,d_MAX_colvector)

%defining persistent variables
persistent PIo
if isempty(PIo)

    PIo = 1e-5*eye(6);
end

persistent xest_prev
if isempty(xest_prev)
    xest_prev = [0.002156;427.2;0.2136;0.03598;0.002156;449.3];
end

%disp(d_initial)
%disp(um);disp(ymeas);disp(Rp);disp(Rx);disp(Ry);disp(d_MIN_colvector);disp(d_MAX_colvector)


OPTIONS = optimset('Algorithm', 'sqp','Display','final');

[dopt,fval,exitflag] = fmincon(@(d) obj_est(um,hstep,Nsteps,d,ymeas,xest_prev,PIo,Rp,Ry)...
    ,d_initial,[],[],[],[],d_MIN_colvector,d_MAX_colvector,[],OPTIONS);

disp(dopt); %display the value for variable dopt
Ecount = toc;
fcost = fval;


[Xobs] = state_cal(um,hstep,Nsteps,dopt);

disp(Xobs)

N_MHE = 6;

[yprediction] = RK4(yplant,uo,hstep,Nsteps);

%ESTIMATE THE CO-VARIANCE MATRIX
PIo = co_var(Xobs,uo,Ts,PIo,Rx,Ry);
disp(PIo)

%ESTIMATE THE APRIORI STATE
[Xprevious] = apriori_state(um,hstep,Nsteps,dopt);
xest_prev=Xprevious;

end


function [Xobs] = state_cal(um,hstep,Nsteps,dopt)

xo = [dopt(1);dopt(2);dopt(3);dopt(4);dopt(5);dopt(6)];

%STEP1
xpred1=RK4(xo,um(:,1),hstep,Nsteps);

%STEP2
xpred2=RK4(xpred1,um(:,2),hstep,Nsteps);

%STEP3
xpred3=RK4(xpred2,um(:,3),hstep,Nsteps);

%STEP4
xpred4=RK4(xpred3,um(:,4),hstep,Nsteps);

%STEP5
xpred5=RK4(xpred4,um(:,5),hstep,Nsteps);

%STEP6
xpred6=RK4(xpred5,um(:,6),hstep,Nsteps);

%STEP7
xpred7=RK4(xpred6,um(:,7),hstep,Nsteps);

%STEP8
xpred8=RK4(xpred7,um(:,8),hstep,Nsteps);

%STEP9
xpred9=RK4(xpred8,um(:,9),hstep,Nsteps);

%STEP10
xpred10=RK4(xpred9,um(:,10),hstep,Nsteps);

[Xobs] = xpred10;
disp(Xobs)

end


function PIk = co_var(Ts,PIo,Rx,Ry)


A = eye(6);

%A=[             1 - (Ts*(2*a1*c1*g*sp*(2*g*h1)^(a1 - 1) + 2*a4*c4*g*sp*(2*g*h1)^(a4 - 1)))/(a*w),                                                                                                                                                                                                                          0,                                                                                                                                                                                                                                          0;
 %(Ts*(2*a1*c1*g*sp*(2*g*h1)^(a1 - 1) + 2*a4*c4*g*sp*(2*g*h1)^(a4 - 1)))/(c*w + (b*h2*w)/Hmax), 1 - (Ts*b*w*(c1*sp*(2*g*h1)^a1 - c2*sp*(2*g*h2)^a5 + c4*sp*(2*g*h1)^a4 - sp*u2*(2*g*h2)^a2))/(Hmax*(c*w + (b*h2*w)/Hmax)^2) - (Ts*(2*a5*c2*g*sp*(2*g*h2)^(a5 - 1) + 2*a2*g*sp*u2*(2*g*h2)^(a2 - 1)))/(c*w + (b*h2*w)/Hmax),                                                                                                                                                                                                                                          0;
  %                                                                                          0,                                                                                                                        (Ts*(2*a5*c2*g*sp*(2*g*h2)^(a5 - 1) + 2*a2*g*sp*u2*(2*g*h2)^(a2 - 1)))/(w*(R^2 - (R - h3)^2)^(1/2)), 1 - (Ts*(2*R - 2*h3)*(c2*sp*(2*g*h2)^a5 - c3*sp*(2*g*h3)^a6 + sp*u2*(2*g*h2)^a2 - sp*u3*(2*g*h3)^a3))/(2*w*(R^2 - (R - h3)^2)^(3/2)) - (Ts*(2*a6*c3*g*sp*(2*g*h3)^(a6 - 1) + 2*a3*g*sp*u3*(2*g*h3)^(a3 - 1)))/(w*(R^2 - (R - h3)^2)^(1/2))]   ;      

C= eye(6) ;
  

%G=[1, 0, 0;
  % 0, 1, 0;
   %0, 0, 1]  ;

%  old arrival cost
%PIk=G*Q*G' + A*PIo*A' - A*PIo*C'*(Rw + C*PIo*C')^(-1)*C*PIo*A'   ;

%PIk=G*Rx*G' + A*PIo*A' - A*PIo*C'*(Ry + C*PIo*C')^(-1)*C*PIo*A'   ;
%PIk=real(PIk) ;


%% new arrival cost update
PIk = Rx + A*PIo*A' - A*PIo*C'*(Ry + C*PIo*C')^(-1)*C*PIo*A'   ;


end

function [Xprevious] = apriori_state(um,hstep,Nsteps,dopt)

xo = [dopt(1);dopt(2);dopt(3);dopt(4);dopt(5);dopt(6)];

disp('xo:'); disp(xo)

xp = RK4(xo,um(:,1),hstep,Nsteps);
Xprevious = xp;

end

function er_fun = obj_est(um,hstep,Nsteps,d,ymeas,xest_prev,PIo,Rp,Ry)

[ypred] = est_fun(um,d,hstep,Nsteps);

ypredr = ypred';
ymr=ymeas';

y1=ypredr(:,1); ym1=ymr(:,1);
y2=ypredr(:,2); ym2=ymr(:,2);
y3=ypredr(:,3); ym3=ymr(:,3);
y4=ypredr(:,4); ym4=ymr(:,4);

%disp('ypredr:'); disp(ypredr);
%disp('ymr:'); disp(ymr)

PI_p = Rp^-1;
PI_y = Ry^-1;

c_arrival = [(d(1)-xest_prev(1));(d(2)-xest_prev(2));(d(3)-xest_prev(3))...
    ;(d(4)-xest_prev(4));(d(5)-xest_prev(5));(d(6)-xest_prev(6))];

%disp('d:'); disp(d);
    %disp('c_arrival:'); disp(c_arrival);
    %disp('ym1-y1:'); disp(ym1-y1);
    %disp('ym2-y2:'); disp(ym2-y2);
    %disp('ym3-y3:'); disp(ym3-y3);
    %disp('ym4-y4:'); disp(ym4-y4);

er_fun=(ym1-y1)'*PI_y(1,1)*(ym1-y1) + (ym2-y2)'*PI_y(2,2)*(ym2-y2) + ...
    (ym3-y3)'*PI_y(3,3)*(ym3-y3) + (ym4-y4)'*PI_y(4,4)*(ym4-y4) + ...
      c_arrival'*PIo^(-1)*c_arrival;

%disp(err_fun)

end


function [ypred] = est_fun(um,d,hstep,Nsteps)
xo = [d(1);d(2);d(3);d(4);d(5);d(6)];

%disp('xo:'); disp(xo)
%disp('um:'); disp(um)

%STEP 1
    
    xpred1=RK4(xo,um(:,1),hstep,Nsteps)   ;  

    %disp('xpred1:'); disp(xpred1)
    
 %STEP 2
    
    xpred2=RK4(xpred1,um(:,2),hstep,Nsteps)   ;
   
    
  %STEP 3
    
    xpred3=RK4(xpred2,um(:,3),hstep,Nsteps)   ;
    
    
      %STEP 4
    
    xpred4=RK4(xpred3,um(:,4),hstep,Nsteps)   ;
    
      %STEP 5
    
    xpred5=RK4(xpred4,um(:,5),hstep,Nsteps)   ;
    
      %STEP 6
    
    xpred6=RK4(xpred5,um(:,6),hstep,Nsteps)   ;
    
          %STEP 7
    
    xpred7=RK4(xpred6,um(:,7),hstep,Nsteps)   ;
    
          %STEP 8
    
    xpred8=RK4(xpred7,um(:,8),hstep,Nsteps)   ;

          %STEP 9
    
    xpred9=RK4(xpred8,um(:,9),hstep,Nsteps)   ;

          %STEP 10
    
    xpred10=RK4(xpred9,um(:,10),hstep,Nsteps)   ;

    ypred=[xpred1,xpred2,xpred3,xpred4,xpred5,xpred6,xpred7,xpred8,xpred9,xpred10];

end


function [xpred]=RK4(xo,uo,hstep,Nsteps)
%disp('uo:'); disp(uo);
%disp('hstep:'); disp(hstep);
%disp('Nsteps:'); disp(Nsteps);
%disp('xo:'); disp(xo)

h=hstep;
xk1=xo;

for k=1:Nsteps
    k1=h*odefun(xk1,uo);
    k2=h*odefun(xk1+k1/2,uo);
    k3=h*odefun(xk1+k2/2,uo);
    k4=h*odefun(xk1+k3,uo);
    
   xk1 = xk1 + 1/6*(k1+2*k2+2*k3+k4);

  % disp('xk1:'); disp(xk1);
end

xpred = xk1;

end

function dydt = odefun(x,u)

%definig the parameters
deltaHfv = 498;
deltaHcr = 506;
deltaHrg = -394200;
cpo = 3.1335;
cpa = 1.074;
cps = 1.005;
sf = 14.57;
sc = 787.9601;
Ecc = 41.79;
Ecr = 101.5;
Ecb = 158.6;
Din = 915.1;
R = 0.0083143;
%kfo = 0.8;
kcc = 0.0189;
kcr = 19620;
Oin = 0.2136;
kcb = 187940;
ksig = 0.1197;
%alpha = 0.12;
Vc = 175738;
Va = 20;
Pris = 100000;
Vris = 2500;
%Frc = 13028;
%Fin = 1720;
%Fa = 2080;
%Toil = 503;
Ta = 298;
%yf0 = 0.8;
%tc = 1.453;
%phi = 0.724;

%%% state parameters
Crc = x(1);
Treg = x(2);
Ofg = x(3);
Ccat = x(4);
Csc = x(5);
Tris = x(6);

%%% input parameters
Fa = u(1);
Fin = u(2);
Toil = u(3);
Frc = u(4);


if any(x <= 0)
        disp('Warning: Non-positive state in odefun:'); disp('x:'); disp(x);
end
    
    if any(isnan(x)) || any(isnan(u))
        disp('NaN detected in odefun input:'); disp(u);
    end

Kr = kcr*exp(-Ecr/R*(Tris - 60))/(Ccat*Crc^0.15);


dydt(1) = ((Frc*(Csc - Crc))/Vc) - kcb*Ofg*Crc*exp(-Ecb/(R*(Treg+510)));

dydt(2) = ((Frc*(Tris - Treg))/Vc) + (((cpa/cps)*Fa*(Ta - Treg)))/(Vc) +  ((-deltaHrg*kcb*Ofg*Crc*(exp(-Ecb/(R*Treg)))))/(sf);

dydt(3) = ((0.0313*Fa*(Oin - Ofg)))/(Va) - (ksig*kcb*Ofg*Crc*(Vc/Va)*exp(-Ecb/(R*(Treg+510))));

%disp('Ccat:');disp(dydt(2))

dydt(4) = (-Frc*Ccat/Vris) + (kcc*Pris*exp(-Ecc/(R*(Tris- 60)))/((Ccat*Crc^0.06)));

dydt(5) = Frc*(Crc - Csc)/Vris + kcc*Pris*exp(-Ecc/(R*(Tris- 60))/(Ccat*Crc^0.06));

dydt(6) = (Frc*(Treg - Tris))/(Vris) + (cpo/(cps*Vris))*Fin*(Toil - Tris) + (0.875*(-deltaHfv*Fin))/(sc*Vris) + (-deltaHcr*Pris*Din*Fin*Kr)/(2*sc*(Fin + 0.0313*Vris*Pris*Din*Kr));


 %if any(isnan(dydt)) || any(~isreal(dydt))
 %       disp('NaN or complex detected in odefun output:'); disp(dydt);
 %end

dydt = [dydt(1);dydt(2);dydt(3);dydt(4);dydt(5);dydt(6)];

end