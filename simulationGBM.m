function [xx,u0,u1,u2] = simulationGBM(param,experiment)

%%%%% Model parameters
Dn = param.Dn; %6.6e-10; %cm2/s
csat = param.csat; %5e7; % cell/mL
chi = param.chi; %cm2/mmHg/s
Tg = param.Tg;
Td = param.Td; %1.7e5; %s
DO2 = param.DO2; % 1e-5 %cm2/s
alpha = param.alpha; %1e-9 %mmHg mL/ cell/s
O2_T  = param.O2_T; %2.5 %mmHg
O2_H  = param.O2_H; %7 %mmHg
O2_A = param.O2_A; % 1.6; %mmHg
dO2_A = param.dO2_A; % 0.1; % mmHg

%%%% Scenario considered
switch experiment
    case 'nc'
        % Data related with the scenario
        load('data/core','XX','YYa','N');
        T = N*24*3600; % total time (s)
        L = XX(end) - XX(1);
        
        % Initial conditions
        C0 = @(x) interp1(XX,YYa(:,1),x); % initial condition interpolant
        
        % Boundary conditions
        J = 1e6; %s/cm
        Jcell_r = J;
        Jcell_l = J;
        Joxygen_r = 0;
        Joxygen_l = 0;
        Koxygen_r = 1;
        Koxygen_l = 1;
        Or = 7;
        Ol = 7;
    case '1p'
        % Data related with the scenario
        load('data/palisade','XX','YY','N');
        %load('data/core','XX','YYa','YYd','N'); YY = YYa;
        T = N*24*3600; % total time (s)
        L = XX(end) - XX(1);
        
        % Initial conditions
        C0 = @(x) interp1(XX,YY(:,1),x); % initial condition interpolant
        
        % Boundary conditions
        J = 1e9; %s/cm 1e9
        Jcell_r = J;
        Jcell_l = J;
        Joxygen_r = 0; 
        Joxygen_l = 1; 
        Koxygen_r = 1; 
        Koxygen_l = 0; 
        Or = 2; % mmHg
        Ol = 0; % 
    case '2p'
        % Data related with the scenario
        load('data/double_palisade','XX','YY','N');
        T = N*24*3600; % total time (s)
        L = XX(end) - XX(1);
        
        % Initial conditions
        C0 = @(x) interp1(XX,YY(:,1),x); % initial condition interpolant
        
        % Boundary conditions
        J = 1.2e7; %s/cm
        Jcell_r = J;
        Jcell_l = J;
        Joxygen_r = 0;
        Joxygen_l = 0;
        Koxygen_r = 1;
        Koxygen_l = 1;
        Or = 7;
        Ol = 7;
end


%--------------------------------------------------------------------------

m = 0;
Nt = 100; % discretization in time
Nx = 664; % discretization in space
tspan = linspace(0,T,Nt);
xspan = linspace(0,L,Nx);

SOL = pdepe(m,@pdefun,@initial,@boundary,xspan,tspan);
xx = xspan;
u0 = SOL(:,:,3);
u1 = SOL(:,:,1);
u2 = SOL(:,:,2);

% --------------------------------------------------------------
function [c,f,s] = pdefun(x,t,u,DuDx)

%% Functions
s = u(3); th = O2_H;
PI_go =  1*(s<=0) + (1-s/th)*(s>0)*(s<=th)+0*(s>th);
PI_gr =  0*(s<=0) + (s/th)*(s>0)*(s<=th)+1*(s>th);

F =  1-u(1)/csat;
Snd =  1/2*(1-tanh((u(3)-O2_A)/dO2_A));
G =  1-(u(1)+u(2))/csat;
H = u(3)/(u(3) + O2_T);
    
% Function c
c = ones(3,1);

% Function f
fdif = zeros(3,1);
fdif(1) = Dn.*DuDx(1);
fdif(2) = 0;
fdif(3) = DO2*DuDx(3);

fchemo = zeros(3,1);
fchemo(1) = chi*PI_go*F*u(1)*DuDx(3);
fchemo(2) = 0;
fchemo(3) = 0;

f = fdif - fchemo;

% Function s
s = zeros(3,1);
s(1) = 1/Tg*PI_gr*G*u(1) -1/Td*Snd*u(1);
s(2) = 1/Td*Snd*u(1);
s(3) = -alpha*H*u(1);

end
% --------------------------------------------------------------
function u0 = initial(x)
u0 = zeros(3,1);
u0(1) = C0(x);
u0(2) = 0;
u0(3) = Or;
end
% --------------------------------------------------------------
function [pl,ql,pr,qr] = boundary(xl,ul,xr,ur,t)
pl = zeros(3,1);
ql = ones(3,1);
pr = zeros(3,1);
qr = ones(3,1);

% Alive cells
pr(1) = ur(1); qr(1) = Jcell_r; % right
pl(1) = ul(1); ql(1) = -Jcell_l; % left

% Dead cells
pr(2) = 0; qr(2) = 1; % right
pl(2) = 0; ql(2) = 1; % left

% Oxygen right
pr(3) = Koxygen_r*(ur(3) - Or); qr(3) = Joxygen_r; % right
pl(3) = Koxygen_l*(ul(3) - Ol); ql(3) = -Joxygen_l; % left

end
end