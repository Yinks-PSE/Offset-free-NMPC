function yk1 = cstrStateFcnDT(yk,uk)

hstep = 1; %0.01; 0.05; 0.1; 0.5
yk1 = yk(:);
Nsteps = 2;

uk1 = [uk(:);0];
for i = 1:Nsteps
    yk1 = yk1 + hstep*cstrStateFcnCT(yk1,uk1);
end