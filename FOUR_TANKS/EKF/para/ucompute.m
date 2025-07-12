function [u_first_move,fcost,Ecount,exitflag,ypred] = ucompute(r1,r2,r3,r4,u1set,u2set,u3set,u4set,Wy1,Wy2,Wy3,Wy4,Wu1,Wu2,Wu3,Wu4...
    ,Wdu1,Wdu2,Wdu3,Wdu4,xo,uo,hstep,Nsteps,Np,Nu,P,alpha,x1min,x2min,x3min,x4min,x1max,x2max,x3max,x4max,du1min,du2min...
    ,du3min,du4min,du1max,du2max,du3max,du4max,u1min,u2min,u3min,u4min,u1max,u2max,u3max,u4max,U_initail,UMIN_colvector,UMAX_colvector)


OPTIONS = optimset('Algorithm', 'sqp', 'Display', 'final');


% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
%fun = objective function
%A = linear inequality constraints
%b = linear inequality constriants
%Aeq = linear equaality constraints
%beq = linear equality constraints
%lb = lower bounds
%ub = upper bounds
%nonlcon = Nonlinear constraints - nlfun
%options = optimization options

[uopt,fval,exitflag] = fmincon(@(u) obj(r1,r2,r3,r4,u1set,u2set,u3set,u4set,Wy1,Wy2,Wy3,Wy4,Wu1,Wu2,Wu3,Wu4...
    ,Wdu1,Wdu2,Wdu3,Wdu4,u,xo,uo,hstep,Nsteps,Np,Nu,P),U_initial,[],[],[],[],UMIN_colvector,UMAX_colvector...
    , @(u) nlfun(uo,xo,u,hstep,Nsteps,Np,Nu,x1min,x2min,x3min,x4min,x1max,x2max,x3max,x4max,du1min,du2min,du3min,du4min...
    ,du1max,du2max,du3max,du4max,P,alpha,r1,r2,r3,r4), OPTIONS);

Ecount = toc;

u1_first_move = uopt(1, :); %all the element in the first row are selected and assign to u_first_move
u2_first_move = uopt(Nu+1, :);
u3_first_move = uopt(2*Nu+1, :);
u4_first_move = uopt(3*Nu+1, :);

u_first_move = [u1_first_move;u2_first_move;u3_first_move;u4_first_move];
fcost = fval;

ypred = RK4(xo,uo,hstep,Nsteps);

[xPx, Tconst] = terminal_const(uo,xo,uopt,hstep,Nsteps,Np,Nu,P,alpha,r1,r2,r3,r4);



end

function cost = obj(r1,r2,r3,r4,u1set,u2set,u3set,u4set,Wy1,Wy2,Wy3,Wy4,Wu1,Wu2,Wu3,Wu4,Wdu1,Wdu2,Wdu3...
    ,Wdu4,u,xo,uo,hstep,Nsteps,Np,Nu,P)


yk = output_prediction(xo,uo,u,hstep,Nsteps,Np,Nu);  %call the 'output_orediction' function to predict the future outputs


%initialization of rate of change of control inputs(du1,du2,du3)and the
%control inputs(u1,u2,u3) to zero vectors of length Nu.

du1 = zeros(Nu,1);
du2 = zeros(Nu,1);
du3 = zeros(Nu,1);
du4 = zeros(Nu,1);

u1 = zeros(1,Nu); 
u2 = zeros(1,Nu);  
u3 = zeros(1,Nu);
u4 = zeros(1,Nu);


%calculate the changes in the control inputs over time(du1,du2,du3)
du1(1) = u1(1)-uo(1); 
du2(1) = u2(1)-uo(2); 
du3(1) = u3(1)-uo(3);
du4(1) = u4(1)-uo(4);

for k=2:Nu-1
    du1(k) = u1(k)-u1(k-1);
    du2(k) = u2(k)-u2(k-1);
    du3(k) = u3(k)-u3(k-1);
    du4(k) = u4(k)-u4(k-1);
end


%OUTPUT COST CALCULATION
%calculate the differences between the setpoints and the predicted outputs
r1set = r1;
r2set = r2;
r3set = r3;
r4set = r4;


fy1 = (r1set*ones(Np-1, 1)-yk(1:Np-1, 1));
fy2 = (r2set*ones(Np-1, 1)-yk(1:Np-1, 2));
fy3 = (r3set*ones(Np-1, 1)-yk(1:Np-1, 3));
fy4 = (r4set*ones(Np-1, 1)-yk(1:Np-1, 4));



%INPUT COST CALCULATIONS
%calculate the difference between the setpoints and the control inputs
fu1 = (u1set*ones(Nu, 1)-u1(1:Nu));
fu2 = (u2set*ones(Nu, 1)-u2(1:Nu));
fu3 = (u3set*ones(Nu, 1)-u3(1:Nu));
fu4 = (u4set*ones(Nu, 1)-u4(1:Nu));

%INPUT RATE COST CALCULATION
%assign the rate changes of the control inputs to fdu1 fdu2 fdu3
fdu1 = du1;
fdu2 = du2;
fdu3 = du3;
fdu4 = du4;


%TOTAL COST CALCULATION
cost1 = fy1'*Wy1*fy1 + fy2'*Wy2*fy2 + fy3'*Wy3*fy3 + fy4'*Wy4*fy4;
cost2 = fu1'*Wu1*fu2 + fu2'*Wu2*fu2 + fu3'*Wu3*fu3 + fu4'*Wu4*fu4;
cost3 = fdu1'*Wdu1*fdu1 + fdu2'*Wdu2*fdu2 + fdu3'*Wdu3*fdu3 + fdu4'*Wdu4*fdu4;


%TERMINAL COST CALCULATION
cost_terminal = (yk(Np,:)-[r1,r2,r3,r4])*P*(yk(Np,:)-[r1,r2,r3,r4])';

cost = cost1 + cost2 + cost3 + cost_terminal;
end

function [c, ceq] = nlfun(uo,xo,u,hstep,Nsteps,Np,Nu,x1min,x2min,x3min,x4min,x1max,x2max,x3max,x4max,du1min,du2min,du3min,du4min...
    ,du1max,du2max,du3max,du4max,P,alpha,r1,r2,r3,r4)

yk = output_prediction(xo,uo,u,hstep,Nsteps,Np,Nu);


du1 = zeros(Nu,1);
du2 = zeros(Nu,1);
du3 = zeros(Nu,1);
du4 = zeros(Nu,1);

u1 = zeros(1,Nu); 
u2 = zeros(1,Nu);  
u3 = zeros(1,Nu);
u4 = zeros(1,Nu);


%calculate the changes in the control inputs over time(du1,du2,du3)
du1(1) = u1(1)-uo(1); 
du2(1) = u2(1)-uo(2); 
du3(1) = u3(1)-uo(3);
du4(1) = u4(1)-uo(4);

for k=2:Nu
    du1(k) = u1(k)-u1(k-1);
    du2(k) = u2(k)-u2(k-1);
    du3(k) = u3(k)-u3(k-1);
    du4(k) = u4(k)-u4(k-1);
end

y1k = yk(:,1);
y2k = yk(:,2);
y3k = yk(:,3);
y4k = yk(:,4);

y1_min_colvector = x1min*ones(Np,1);
y1_max_colvector = x1max*ones(Np,1);
du1_min_colvector = du1min*ones(Nu,1);
du1_max_colvector = du1max*ones(Nu,1);

y2_min_colvector = x2min*ones(Np,1);
y2_max_colvector = x2max*ones(Np,1);
du2_min_colvector = du2min*ones(Nu,1);
du2_max_colvector = du2max*ones(Nu,1);

y3_min_colvector = x3min*ones(Np,1);
y3_max_colvector = x3max*ones(Np,1);
du3_min_colvector = du3min*ones(Nu,1);
du3_max_colvector = du3max*ones(Nu,1);

y4_min_colvector = x4min*ones(Np,1);
y4_max_colvector = x4max*ones(Np,1);
du4_min_colvector = du4min*ones(Nu,1);
du4_max_colvector = du4max*ones(Nu,1);

%The inequality constraint functions are defined
%ineq_du_min is the inequality constraint function on du_min
%ineq_du_max is the inequality constraint function on du_min
%ineq_y_min is the inequality constraint function on y_min
%ineq_y_max is the inequality constraint function on y_max




ineq_y1_min = y1_min_colvector-y1k;
ineq_y1_max = y1k-y1_max_colvector;
ineq_du1_min = du1_min_colvector-du1;
ineq_du1_max = du1-du1_max_colvector;

ineq_y2_min = y2_min_colvector-y2k;
ineq_y2_max = y2k-y2_max_colvector;
ineq_du2_min = du2_min_colvector-du2;
ineq_du2_max = du2-du2_max_colvector;

ineq_y3_min = y3_min_colvector-y3k;
ineq_y3_max = y3k-y3_max_colvector;
ineq_du3_min = du3_min_colvector-du3;
ineq_du3_max = du3-du3_max_colvector;

ineq_y4_min = y4_min_colvector-y4k;
ineq_y4_max = y4k-y4_max_colvector;
ineq_du4_min = du4_min_colvector-du4;
ineq_du4_max = du4-du4_max_colvector;


Terminal_constraint = (yk(Np,:)-[r1,r2,r3,r4])*P*(yk(Np,:)-[r1,r2,r3,r4])' - alpha;

c = [ineq_y1_min; ineq_y1_max; ineq_du1_min; ineq_du1_max;...
     ineq_y2_min; ineq_y2_max; ineq_du2_min; ineq_du2_max;...
     ineq_y3_min; ineq_y3_max; ineq_du3_min; ineq_du3_max;...
     ineq_y4_min; ineq_y4_max; ineq_du4_min; ineq_du4_max;...
     Terminal_constraint];

ceq = []; %ceq = [] means no equality constraints are defined for this problem
end


function [xPx, Tconst] = terminal_const(uo,xo,u,hstep,Nsteps,Np,Nu,P,alpha,r1,r2,r3,r4)

yk = output_prediction(xo,uo,u,hstep,Nsteps,Np,Nu);
Np = Np-1;

xPx = [yk(Np,:)-[r1,r2,r3,r4]]*P*[yk(Np,:)-[r1,r2,r3,r4]]' ;
Tconst = [yk(Np,:)-[r1,r2,r3,r4]]*P*[yk(Np,:)-[r1,r2,r3,r4]]' - alpha;

end

function yk = output_prediction(xo,uo,u,hstep,Nsteps,Np,Nu)
xk = state_prediction(xo,uo,u,hstep,Nsteps,Np,Nu);
yk = xk + repmat(Np,1);

%repmat- repeat copes of array

end

function xk = state_prediction(xo,uo,u,hstep,Nsteps,Np,Nu)

xo=xo(1:3);

du1 = zeros(Nu,1);
du2 = zeros(Nu,1);
du3 = zeros(Nu,1);
du4 = zeros(Nu,1);

u1 = zeros(Nu,1); 
u2 = zeros(Nu,1);  
u3 = zeros(Nu,1);
u4 = zeros(Nu,1);

u = [u1'; u2'; u3'; u4'];


du1(1) = u1(1)-uo(1);
du2(1) = u2(1)-uo(2);
du3(1) = u3(1)-uo(3);
du4(1) = u4(1)-uo(4);

for k=2:Nu
    du1(k) = u1(k)-u1(k-1);
    du2(k) = u2(k)-u2(k-1);
    du3(k) = u3(k)-u3(k-1);
    du4(k) = u4(k)-u4(k-1);
end

%STEP 1
%k=0;
    x(:,1) = RK4(xo,uo,hstep,Nsteps)  ;
    
    
for k=1:Nu-1
    
     x(:,k+1) = RK4(x(:,k),u(:,k),hstep,Nsteps)  ;
   % x(:,k+1)=x(:,k)+h*odefun(x(:,k),u(:,k));  %u(:,k)
end


for k=Nu:Np-1

    
    x(:,k+1) = RK4(x(:,k),u(:,Nu),hstep,Nsteps)  ;
    
end





xk=x';

end

function [xpred] = RK4(xo,uo,hstep,Nsteps)
h=hstep;
xk1=xo;

for k=1:Nsteps

    k1=h*odefun(xk1,uo);
    k2=h*odefun(xk1+k1/2,uo);
    k3=h*odefun(xk1+k2/2,uo);
    k4=h*odefun(xk1+k3,uo);

    xk1 = xk1 + 1/6*(k1+2*k2+2*k3+k4);
end
xpred = xk1;
    

end


function dydt = odefun(x, u)



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
alpha = 0.12;
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

Kr = kcr*exp(-Ecr/R*(Tris - 60))/(Ccat*Crc^0.15);


dydt(1) = ((Frc*(Csc - Crc))/Vc) - kcb*Ofg*Crc*exp(-Ecb/(R*(Treg+510)));

dydt(2) = ((Frc*(Tris - Treg))/Vc) + (((cpa/cps)*Fa*(Ta - Treg)))/(Vc) +  ((-deltaHrg*kcb*Ofg*Crc*(exp(-Ecb/(R*Treg)))))/(sf);

dydt(3) = ((0.0313*Fa*(Oin - Ofg)))/(Va) - (ksig*kcb*Ofg*Crc*(Vc/Va)*exp(-Ecb/(R*(Treg+510))));

%disp('Ccat:');disp(dydt(2))

dydt(4) = (-Frc*Ccat/Vris) + (kcc*Pris*exp(-Ecc/(R*(Tris- 60)))/((Ccat*Crc^0.06)));

dydt(5) = Frc*(Crc - Csc)/Vris + kcc*Pris*exp(-Ecc/(R*(Tris- 60))/(Ccat*Crc^0.06));

dydt(6) = (Frc*(Treg - Tris))/(Vris) + (cpo/(cps*Vris))*Fin*(Toil - Tris) + (0.875*(-deltaHfv*Fin))/(sc*Vris) + (-deltaHcr*Pris*Din*Fin*Kr)/(2*sc*(Fin + 0.0313*Vris*Pris*Din*Kr));


dydt = [dydt(1);dydt(2);dydt(3);dydt(4);dydt(5);dydt(6)];

end
 




