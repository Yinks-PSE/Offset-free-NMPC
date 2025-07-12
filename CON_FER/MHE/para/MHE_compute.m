function [Xobs,dist_output,Yxs_est,fcost,Ecount,exitflag]=MHE_compute(uo,um,ymeas...
                                        ,yplant,Ts,hstep,Nsteps,Rp,Rx,Ry...
                                        ,d_initial, d_MIN_colvector,d_MAX_colvector)

% DEFINING PERSISTENT VARIABLES

persistent PIo
if isempty(PIo)

PIo=1e-4*eye(3) ; 

end

persistent xest_prev
if isempty(xest_prev)
xest_prev = [6.0;5.0;19.14];
end

persistent p_bar
if isempty(p_bar)
p_bar = 0.4 ;
end

% % OPTIONS = optimset('Algorithm','sqp','Display','final');
% OPTIONS = optimset('Algorithm','active-set');
% % OPTIONS = optimset('Algorithm','Interior-Point');
% 
% tic

OPTIONS = optimset('Algorithm','sqp','Display','final');
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
%disp(dopt)
Ecount = toc ;
fcost=fval ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
Yxs_est=[dopt(4)]  ;
Yxs=Yxs_est;
p_bar = Yxs ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE Xobs from estimated state
[Xobs]=state_cal(um,hstep,Nsteps,dopt)  ;
%xo=Xobs  ; 
N_MHE = 5 ;

yprediction = Xobs(1:1);


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
PIo=co_var(Xobs,uo,hstep,PIo,Yxs,Rx,Ry)

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
xo = [dopt(1);dopt(2);dopt(3)]  ;

%disp('size xo in state_cal');disp(size(xo))
Yxs = dopt(4)  ;


%STEP 1
    
    xpred1=Euler(xo,um(:,1),Yxs,hstep,Nsteps)   ;
    
 %STEP 2
    
    xpred2=Euler(xpred1,um(:,2),Yxs,hstep,Nsteps)   ;
    %xpred2=xpred1+1/6*(k1_2+2*k2_2+2*k3_2+k4_2) ;
    
  %STEP 3
    
    xpred3=Euler(xpred2,um(:,3),Yxs,hstep,Nsteps)   ;
    
    %STEP 4

    xpred4=Euler(xpred3,um(:,4),Yxs,hstep,Nsteps)   ;
    
   %STEP 5
    
    xpred5=Euler(xpred4,um(:,5),Yxs,hstep,Nsteps)   ;
    
    
    %STEP 6
    
    xpred6=Euler(xpred5,um(:,6),Yxs,hstep,Nsteps)   ;

    xpred7=Euler(xpred6,um(:,7),Yxs,hstep,Nsteps);
    xpred8=Euler(xpred7,um(:,8),Yxs,hstep,Nsteps);
    xpred9=Euler(xpred8,um(:,9),Yxs,hstep,Nsteps);
    xpred10=Euler(xpred9,um(:,10),Yxs,hstep,Nsteps);
    
    
[Xobs]=xpred10 ; 
%disp(Xobs)


end
function PIk=co_var(ymeas,input,Ts,PIo,Yxs,Rx,Ry)

%definig the parameters
 Ki = 22; %Substrate inhibition constatnt
 Km= 1.2;     % substrate saturation constant
    Pm = 50;     % product saturation constant
    Sf = 20; % substrate concentration in the feed
    %Yxs = 0.4; %cell-mass yeild
    alpha = 2.2; %kinematic parameter
    beta = 0.2; %kinemtic prameter
    miu_m = 0.48; %the maximum growth rate

    %%%%%%INPUT PARAMETER
D = input(1);
xhat=ymeas  ;

%Ts=1 ;

X = xhat(1);
S = xhat(2); 
P = xhat(3);
% Jacobian matrix A (linearized state transition matrix)

a11 = 1 - Ts*(D - (miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki)));
a12 = 0;
a13 = 0;

a21=-(Ts*(miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki)))/Yxs;
a22=1 - D*Ts;
a23=0;

a31= -Ts*(beta + alpha*(miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki)));
a32=0;
a33 = 1 - D*Ts;

A = [a11,a12,a13;a21,a22,a23;a31,a32,a33];


C = [1, 0, 0];

G = [1,0,0;
    0,1,0;
    0,0,1];
    

    %y + C*PIo*C'
%A*PIo*C'*(Ry + C*PIo*C')^(-1)*C*PIo*A' 
%A*PIo*A'
%G*Rx*G' 

PIk=G*Rx*G' + A*PIo*A' - A*PIo*C'*(Ry + C*PIo*C')^(-1)*C*PIo*A' ;

end

function [Xprevious]=apriori_state(um,hstep,Nsteps,dopt);


xo=[dopt(1);dopt(2);dopt(3)]  ;
Yxs=dopt(4)  ;

%STEP 1
    
    xp=Euler(xo,um(:,1),Yxs,hstep,Nsteps);
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


xbar=[(d(1)-xest_prev(1)) ;  (d(2)-xest_prev(2)); (d(3)-xest_prev(3)) ];

er_fun=(ymr-ypredr)'*Wx1*(ymr-ypredr)  +  Wp*(d(4)- p_bar)^2  + xbar'*PIo^(-1)*xbar   ;
%er_fun=(ym1-y1)'*Wx1*(ym1-y1) + (ym2-y2)'*Wx2*(ym2-y2)  +  Wp*(p(3)- p_prev)^2  + xbar'*PIo^(-1)*xbar   ;



%er_fun=(ym1-y1)'*Wx1*(ym1-y1) + (ym2-y2)'*Wx2*(ym2-y2)  +  Wp*p(3)^2  + xbar'*PI*xbar ;
%er_fun=(ym1-y1)'*Wx1*(ym1-y1) + (ym2-y2)'*Wx2*(ym2-y2)  +  Wx*xd 
%er_fun=(ym1-y1)'*Wx1*(ym1-y1) + (ym2-y2)'*Wx2*(ym2-y2)   ;







end

function [ypred]=est_fun(um,hstep,Nsteps,d);

xo=[d(1);d(2);d(3)]  ;

%disp('size xo in est_fun');disp(size(xo))

Yxs=d(4)   ;


%STEP 1
    
    xpred1=Euler(xo,um(:,1),Yxs,hstep,Nsteps)   ;
    
 %STEP 2
    
    xpred2=Euler(xpred1,um(:,2),Yxs,hstep,Nsteps)   ;
    %xpred2=xpred1+1/6*(k1_2+2*k2_2+2*k3_2+k4_2) ;
    
  %STEP 3
    
    xpred3=Euler(xpred2,um(:,3),Yxs,hstep,Nsteps)   ;
    
    %STEP 4

    xpred4=Euler(xpred3,um(:,4),Yxs,hstep,Nsteps)   ;
    
   %STEP 5
   
    xpred5=Euler(xpred4,um(:,5),Yxs,hstep,Nsteps)   ;
    
    
    %STEP 6
    
    xpred6=Euler(xpred5,um(:,6),Yxs,hstep,Nsteps)   ;

    xpred7=Euler(xpred6,um(:,7),Yxs,hstep,Nsteps);
    xpred8=Euler(xpred7,um(:,8),Yxs,hstep,Nsteps);
    xpred9=Euler(xpred8,um(:,9),Yxs,hstep,Nsteps);
    xpred10=Euler(xpred9,um(:,10),Yxs,hstep,Nsteps);

    
  

%ypred=[xpred1 +  [d(1);d(2)] ,...
%       xpred2 +  [d(3);d(4)] ,...
 %      xpred3 +  [d(5);d(6)] ,...
 %      xpred4 +  [d(7);d(8)] ,...
 %      xpred5 +  [d(9);d(10)] ,...
  %     xpred6 +  [d(11);d(12)] ]  ;
   
    
    
    
%ypred=[xpred1(1);xpred2(1);xpred3(1);xpred4(1);xpred5(1);xpred6(1)]  
%ypred=[xpred1, xpred2, xpred3, xpred4, xpred5, xpred6] ; 
ypred=[xpred1, xpred2, xpred3, xpred4,xpred5,xpred6,xpred7,xpred8,xpred9,xpred10] ;
end

function [xpred]=Euler(xo,uo,Yxs,hstep,Nsteps);


%xk=RK4(xo,uo,du,h,Np,Nu,ko,dx1,dx2);
%disp(xo)
%delta=Ts/Nsteps  ;
xk1=xo ;
    
for k=1:Nsteps
    
   % xk1=xk1 + hstep*odefun(xk1,uo,ko);  %u(:,k)
    k1=hstep*odefun(xk1,uo,Yxs);
    k2=hstep*odefun(xk1+k1/2,uo,Yxs);
    k3=hstep*odefun(xk1+k2/2,uo,Yxs);
    k4=hstep*odefun(xk1+k3,uo,Yxs);
    
    xk1 = xk1 + 1/6*(k1+2*k2+2*k3+k4);
   
   
   
end

xpred = xk1;


end


 


function xdot=odefun(x,u,Yxs);

%definig the parameters
 Ki = 22; %Substrate inhibition constatnt
 Km= 1.2;     % substrate saturation constant
    Pm = 50;     % product saturation constant
    Sf = 20; % substrate concentration in the feed
    %Yxs = 0.4; %cell-mass yeild
    alpha = 2.2; %kinematic parameter
    beta = 0.2; %kinemtic prameter
    miu_m = 0.48; %the maximum growth rate


X = x(1);
S = x(2);
P = x(3);

D = u(1);   


miu = (miu_m * (1 - (P/Pm))) * S / (Km + S + (S^2 / Ki))  ;

dx1dt = -D * X + miu * X   ;

dx2dt = D * (Sf - S) - (1/Yxs) * miu * X ;

dx3dt = -D * P + (alpha * miu + beta) * X  ;

xdot = [dx1dt; dx2dt; dx3dt];

end

