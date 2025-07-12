function yk1 = conFerStateFcnDT(yk,uk)

hstep = 0.05; %0.01; 0.05; 0.1; 0.5
yk1 = yk(:);
Nsteps = 4;

uk1 = [uk(:);0;0;0];
for i = 1:Nsteps
    yk1 = yk1 + hstep*conFerStateFcnCT(yk1,uk1);
end