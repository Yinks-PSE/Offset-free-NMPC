function [u_first_move,fcost,Ecount,exitflag,ypred]=ucompute(r1,r2,u1set,u2set,Wy1,Wy2...
         ,Wu1,Wu2,Wdu1,Wdu2,xo,uo,hstep,Nsteps,Np,d1,d2,d3,d4,cof,dx...
         ,x1min,x2min,x1max,x2max,u1min,u2min,u1max,u2max...
          ,du1min, du2min, du1max, du2max...
         ,U_initial,DUMIN_colvector,DUMAX_colvector)

    
OPTIONS = optimset('Algorithm','sqp','Display','final');
%OPTIONS = optimset('Algorithm','active-set');

tic

[duopt,fval,exitflag] =fmincon(@(du) obj(r1,r2,u1set,u2set,Wy1,Wy2...
    ,Wu1,Wu2,Wdu1,Wdu2,du,xo,uo,hstep,Nsteps,Np,d1,d2,d3,d4,cof,dx)...
    ,U_initial,[],[],[],[],DUMIN_colvector,DUMAX_colvector,@(du) nlfun(uo,xo,du...
    ,hstep,Nsteps,Np,x1min,x2min...
     ,x1max,x2max,u1min,u2min,u1max,u2max...
      ,du1min, du2min, du1max, du2max...
      ,cof,dx),OPTIONS)   ;     



 Ecount = toc ;
 
u1_first_move=duopt(1,:)      ;
u2_first_move=duopt(Np+1,:)   ;



u_first_move=[u1_first_move ;u2_first_move ]   ;
fcost=fval ;

ypred = RK4(xo,uo,cof,dx,hstep,Nsteps)  ;
ypred = ypred(1:2);

end


function cost=obj(r1,r2,u1set,u2set,Wy1,Wy2...
    ,Wu1,Wu2,Wdu1,Wdu2,du,xo,uo,hstep,Nsteps,Np,d1,d2,d3,d4,cof,dx);


xk=yprediction(xo,uo,du,hstep,Nsteps,Np,cof,dx)   ;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% du1 = zeros(Nu,1); 
% du2 = zeros(Nu,1);  
% du3 = zeros(Nu,1);

u1 = zeros(1,Np); 
u2 = zeros(1,Np);  

du1 = zeros(Np,1); 
du2 = zeros(Np,1);  


for i = 1:Np
du1(i,:) = du(i);
end

for i = 1:Np
du2(i,:) = du(Np+i);
end



%%%%%%%%%%%%%%%%%%%%%%%%%
    u1(1)=du1(1)+uo(1)  ;   u2(1)=du2(1)+uo(2)  ;  u(:,1)=[u1(1),u2(1)] ;
for k=2:Np
    u1(k)=du1(k)+u1(k-1)  ;  u2(k)=du2(k)+u2(k-1)  ; u(:,k)=[u1(k),u2(k)] ;
end



% OUTPUT COST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1set=r1-d1  ;  r2set=r2-d2   ; 
%xk(:,3)
%xk(1:Np-1,3)

%fy1=(r1set*ones(Np-1,1)-xk(:,1));
%fy2=(r2set*ones(Np-1,1)-xk(:,2));
%fy3=(r3set*ones(Np-1,1)-xk(:,3));

fy1=(r1set*ones(Np-1,1)-xk(1:Np-1,1));
fy2=(r2set*ones(Np-1,1)-xk(1:Np-1,2));

% INPUT COST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fu1=(u1set*ones(Np,1)-u1');
fu2=(u2set*ones(Np,1)-u2');


% INPUT RATE COST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fdu1=du1;    fdu2=du2;     


%THE COSTS
cost1=fy1'*Wy1*fy1     + fy2'*Wy2*fy2        ;
cost2=fu1'*Wu1*fu1     + fu2'*Wu2*fu2          ;
cost3=fdu1'*Wdu1*fdu1  + fdu2'*Wdu2*fdu2       ;
%cost3=n1'*Wdu1*n1  + n2'*Wdu2*n2      ;

%  The terminal cost


%cost_terminal=[xk(Np,:)-[r1,r2,r3,r4]]*P*[xk(Np,:)-[r1,r2,r3,r4]]'     ;
%cost_terminal=[xk(Np,:)]*P*[xk(Np,:)]'     ;

cost=cost1  +   cost2 +  cost3    ;      


end


function [c,ceq] = nlfun(uo,xo,du...
    ,hstep,Nsteps,Np,x1min,x2min...
     ,x1max,x2max,u1min,u2min,u1max,u2max...
     ,du1min, du2min, du1max, du2max...
      ,cof,dx);

xk=yprediction(xo,uo,du,hstep,Nsteps,Np,cof,dx)   ;

du1 = zeros(Np,1); 
du2 = zeros(Np,1);  


for i = 1:Np
du1(i,:) = du(i);
end

for i = 1:Np
du2(i,:) = du(Np+i);
end



%%%%%%%%%%%%%%%%%%%%%%%%%
    u1(1)=du1(1)+uo(1)  ;   u2(1)=du2(1)+uo(2)  ;  u(:,1)=[u1(1),u2(1)] ;
for k=2:Np
    u1(k)=du1(k)+u1(k-1)  ;  u2(k)=du2(k)+u2(k-1)  ; u(:,k)=[u1(k),u2(k)] ;
end


y1k=xk(:,1)   ; %output 1
y2k=xk(:,2)    ; %output 2

%%%%%%%%%%%%%%

y1_min_colvector = x1min*ones(Np,1);
y1_max_colvector = x1max*ones(Np,1);
u1_min_colvector = u1min*ones(Np,1);
u1_max_colvector = u1max*ones(Np,1);
du1_min_colvector = du1min*ones(Np,1);
du1_max_colvector = du1max*ones(Np,1);


y2_min_colvector = x2min*ones(Np,1);
y2_max_colvector = x2max*ones(Np,1);
u2_min_colvector = u2min*ones(Np,1);
u2_max_colvector = u2max*ones(Np,1);
du2_min_colvector = du2min*ones(Np,1);
du2_max_colvector = du2max*ones(Np,1);


%The inequality constraint functions are defined
%ineq_du_min is the inequality constraint function on du_min
%ineq_du_max is the inequality constraint function on du_min
%ineq_y_min is the inequality constraint function on y_min
%ineq_y_max is the inequality constraint function on y_max

%size(y1_min_colvector)
%size(y1k)


ineq_du1_min = du1_min_colvector - du1 ;
ineq_du1_max = du1 - du1_max_colvector;

ineq_y1_min = y1_min_colvector - y1k   ;
ineq_y1_max = y1k - y1_max_colvector  ;
%disp(u1_min_colvector)
%disp(u1)
ineq_u1_min = u1_min_colvector - u1';
ineq_u1_max = u1' - u1_max_colvector;

ineq_du2_min = du2_min_colvector - du2 ;
ineq_du2_max = du2 - du2_max_colvector;
ineq_y2_min = y2_min_colvector - y2k;
ineq_y2_max = y2k - y2_max_colvector;
ineq_u2_min = u2_min_colvector - u2';
ineq_u2_max = u2' - u2_max_colvector;

%c is a vector function of all the inequality cocnstraints
%ceq = [] means no equality constraints are defined for this problem



%Terminal_constraint=[xk(Np,:)-[r1,r2,r3,r4]]*P*[xk(Np,:)-[r1,r2,r3,r4]]' -  alpha    ;
%cterm=X*P*X' - 0.7  ;
%cost_terminal=[xk(Np,:)-[r1,r2,r3]]*P*[xk(Np,:)-[r1,r2,r3]]'     ;
%cost_terminal=[xk(Np,:)]*P*[xk(Np,:)]'     ;


c = [ineq_du1_min; ineq_du1_max; ineq_y1_min; ineq_y1_max; ineq_u1_min; ineq_u1_max;...
    ineq_du2_min;  ineq_du2_max; ineq_y2_min; ineq_y2_max; ineq_u2_min; ineq_u2_max]      ;
 
 
% c = [ineq_y1_min; ineq_y1_max; ineq_u1_min; ineq_u1_max;...
%      ineq_y2_min; ineq_y2_max; ineq_u2_min; ineq_u2_max;...
%      ineq_y3_min; ineq_y3_max; ineq_u3_min; ineq_u3_max;...
%      Terminal_constraint]      ;
 

 

ceq = [];

end


function xk=yprediction(xo,uo,du,hstep,Nsteps,Np,cof,dx);
    %     xk=yprediction(xo,uo,n,hstep,Nsteps,Np,cof,dx, a, N);
%disp(size(dx))
xo=xo(1:4)  ;

du1 = zeros(Np,1); 
du2 = zeros(Np,1);  


for i = 1:Np
du1(i,:) = du(i);
end

for i = 1:Np
du2(i,:) = du(Np+i);
end



%%%%%%%%%%%%%%%%%%%%%%%%%
    u1(1)=du1(1)+uo(1)  ;   u2(1)=du2(1)+uo(2)  ;  u(:,1)=[u1(1),u2(1)] ;
for k=2:Np
    u1(k)=du1(k)+u1(k-1)  ;  u2(k)=du2(k)+u2(k-1)  ; u(:,k)=[u1(k),u2(k)] ;
end




%STEP 1
%k=0;
    x(:,1) = RK4(xo,uo,cof,dx,hstep,Nsteps)  ;
    
    
for k=1:Np-1
    
     x(:,k+1) = RK4(x(:,k),u(:,k),cof,dx,hstep,Nsteps)  ;
   % x(:,k+1)=x(:,k)+h*odefun(x(:,k),u(:,k));  %u(:,k)
end



xk=x';



end 





function [xpred]=RK4(xo,uo,cof,dx,hstep,Nsteps);


h=hstep ;
xk1=xo ;
    
for k=1:Nsteps
    
    %xk1=xk1 + delta*odefun(xk1,uo,cof);  %u(:,k)
    
    k1=h*odefun(xk1,uo,cof,dx);
    k2=h*odefun(xk1+k1/2,uo,cof,dx);
    k3=h*odefun(xk1+k2/2,uo,cof,dx);
    k4=h*odefun(xk1+k3,uo,cof,dx);
    
   xk1 = xk1 + 1/6*(k1+2*k2+2*k3+k4);
  
    
end
 xpred = xk1 ;

end
 


function xdot=odefun(x,u,cof,dx);

%values of parameters
A1=35; A3=35; A2=39  ;A4=39  ;
a1=0.97  ;a3=0.42;  a2=0.74  ;a4=0.426   ;g=981  ;

% alpha1=0.5588*0.9; alpha2=0.7504*1.1;  alpha3=0.6221*0.9; alpha4=0.6196*1.15; 
alpha1=0.8542*0.95; alpha3=0.7126*0.95;  alpha2=0.7722*0.95; alpha4=0.5149*0.95; 

gamma1=0.23   ;
gamma2=0.19   ;
h1=x(1)  ;  h2=x(2)  ;  h3=x(3)  ;  h4=x(4)  ;
 v1=u(1);
 v2=u(2);

dh1dt=-(a1*alpha1/A1)*(2*g*h1)^0.5+(a3*alpha3/A1)*(2*g*h3)^0.5+(gamma1/A1)*v1 + dx(1);
dh2dt=-(a2*alpha2/A2)*(2*g*h2)^0.5+(a4*alpha4/A2)*(2*g*h4)^0.5+(gamma2/A2)*v2 + dx(2)  ;
dh3dt=-(a3*alpha3/A3)*(2*g*h3)^0.5+((1-gamma2)/A3)*v2 + dx(3)  ;
dh4dt=-(a4*alpha4/A4)*(2*g*h4)^0.5+((1-gamma1)/A4)*v1 + dx(4)  ;


xdot=[dh1dt;dh2dt;dh3dt;dh4dt];


end