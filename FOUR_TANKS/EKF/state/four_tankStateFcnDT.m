function hk1 = four_tankStateFcnDT(hk,uk)


hstep = 1;
hk1 = hk(:);
Nsteps = 2;  % Number of integration time steps for Euler method


uk1 = [uk(:);0;0;0;0];  % The manipulated input + white noise of size 2
for i = 1:Nsteps
    hk1 = hk1 + hstep*four_tankStateFcnCT(hk1,uk1);
end




