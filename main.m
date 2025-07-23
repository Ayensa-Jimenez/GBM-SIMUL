%%%%%%%%%%%%%%%%%%%%%%%%%% GBM SIMUL CHIP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       This script recreates GBM progression in microfluidic devices and
%       compares the results obtained with three different experiments. The
%       mathematical model and computational approach is based on
%       Ayensa-Jiménez et al. (2020). There are three different
%       experiments:
%           1. Necrotic core experiment, from Ayuso et al. (2017).
%           2. Pseudopalisade experiment, from Ayuso et al. (2016).
%           3. Double pseudopalisade experiment, from Ayensa-Jimenez et al.
%           (2020).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

%%%% Parameters for figures
color1 = [0.4667    0.6745    0.1882];
color2 = [0.8510    0.3255    0.0980];
lw = 1.5;
s = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Model parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dn =  6.6e-10; %cm2/s
csat  = 5e7; % cell/mL
chi = 7.5e-9; %cm2/mmHg/s
Tg = 7.2e5; %s
Td = 1.7e5; %s
DO2 = 1e-5; %cm2/s
alpha = 1e-9; %mmHg mL/ cell/s
O2_T  = 2.5; %mmHg
O2_H  = 7; %mmHg
O2_A = 1.6; %mmHg
dO2_A = 0.1; % mmHg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading parameters
param.Dn = Dn;
param.csat = csat;
param.chi = chi;
param.Tg = Tg;
param.Td = Td;
param.DO2 = DO2;
param.alpha = alpha;
param.O2_T = O2_T;
param.O2_H = O2_H;
param.O2_A = O2_A;
param.dO2_A = dO2_A;

%%%%%%%%%%%%%%%%%%%%%% Necrotic core experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation
exp = 'nc';
tic;
[xx,u0,u1,u2] = simulationGBM(param,exp);
toc;

%%%% Comparing profiles
[nt,nx] = size(u0);
load('data/core','XX','YYa','YYd')

%%
% Days: [0*,1*,2*,3,6];
fig = figure(1);
subplot(1,2,1); % Day 3
n1 = floor(nt/2);
xxe = XX;
yy1e = YYa(:,4); yy2e = YYd(:,4);
yy1 = u1(n1,:); yy2 = u2(n1,:);
plot(xxe,yy1e,'linewidth',lw,'linestyle','--','color',color1,'DisplayName','Exp: live cells'); hold on; grid on;
plot(xx,yy1,'linewidth',lw,'linestyle','-','color',color1,'DisplayName','Sim: live cells');
plot(xxe,yy2e,'linewidth',lw,'linestyle','--','color',color2,'DisplayName','Exp: dead cells');
plot(xx,yy2,'linewidth',lw,'linestyle','-','color',color2,'DisplayName','Sim: dead cells');
ylabel('$u$ (cell/mL)','fontsize',s,'interpreter','latex')
xlabel('$x$ (cm)','fontsize',s,'interpreter','latex')
legend('fontsize',s,'interpreter','latex','location','north')
ylim([0 1.2*csat]);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
title('Day 3','fontsize',s,'interpreter','latex')

subplot(1,2,2); % Day 6
n2 = nt;
yy1e = YYa(:,5); yy2e = YYd(:,5);
yy1 = u1(n2,:); yy2 = u2(n2,:);
plot(xxe,yy1e,'linewidth',lw,'linestyle','--','color',color1,'DisplayName','Exp: live cells'); hold on; grid on;
plot(xx,yy1,'linewidth',lw,'linestyle','-','color',color1,'DisplayName','Sim: live cells');
plot(xxe,yy2e,'linewidth',lw,'linestyle','--','color',color2,'DisplayName','Exp: dead cells');
plot(xx,yy2,'linewidth',lw,'linestyle','-','color',color2,'DisplayName','Sim: dead cells');
ylabel('$u$ (cell/mL)','fontsize',s,'interpreter','latex')
xlabel('$x$ (cm)','fontsize',s,'interpreter','latex')
legend('fontsize',s,'interpreter','latex','location','north')
ylim([0 1.2*csat]);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
title('Day 6','fontsize',s,'interpreter','latex')

%%%%%%%%%%%%%%%%%%%%%% Pseudopalisade experiment %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation
exp = '1p';
tic;
[xx,u0,u1,u2] = simulationGBM(param,exp);
toc;

%%%% Comparing profiles
[nt,nx] = size(u0);
load('data/palisade','XX','YY')

%%
% Days: [1*,3,6];
fig = figure(2);
subplot(1,2,1); % Day 6
n1 = floor(nt/2); 
xxe = XX;
yy1e = YY(:,2);
yy1 = u1(n1,:); yy2 = u2(n1,:);
plot(xxe,yy1e,'linewidth',lw,'linestyle','--','color',color1,'DisplayName','Exp: live cells'); hold on; grid on;
plot(xx,yy1,'linewidth',lw,'linestyle','-','color',color1,'DisplayName','Sim: live cells');
plot(xx,yy2,'linewidth',lw,'linestyle','-','color',color2,'DisplayName','Sim: dead cells');
ylabel('$u$ (cell/mL)','fontsize',s,'interpreter','latex')
xlabel('$x$ (cm)','fontsize',s,'interpreter','latex')
legend('fontsize',s,'interpreter','latex','location','north')
ylim([0 1.2*csat]);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
title('Day 6','fontsize',s,'interpreter','latex')

subplot(1,2,2); % Day 9
n2 = nt;
yy1e = YY(:,3);
yy1 = u1(n2,:); yy2 = u2(n2,:);
plot(xxe,yy1e,'linewidth',lw,'linestyle','--','color',color1,'DisplayName','Exp: live cells'); hold on; grid on;
plot(xx,yy1,'linewidth',lw,'linestyle','-','color',color1,'DisplayName','Sim: live cells');
plot(xx,yy2,'linewidth',lw,'linestyle','-','color',color2,'DisplayName','Sim: dead cells');
ylabel('$u$ (cell/mL)','fontsize',s,'interpreter','latex')
xlabel('$x$ (cm)','fontsize',s,'interpreter','latex')
legend('fontsize',s,'interpreter','latex','location','north')
ylim([0 1.2*csat]);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
title('Day 9','fontsize',s,'interpreter','latex')

%%%%%%%%%%%%%%%%% Double pseudopalisade experiment %%%%%%%%%%%%%%%%%%%%%%%%
%%%% Simulation
exp = '2p';
tic;
[xx,u0,u1,u2] = simulationGBM(param,exp);
toc;

%%%% Comparing profiles
[nt,nx] = size(u0);
load('data/double_palisade','XX','YY')

%% 
% Days: [1*,3*,5*,7,11,17,21];
fig = figure(3);
subplot(2,2,1); % Day 6
n1 = floor(6/20*nt); 
xxe = XX;
yy1e = YY(:,4);
yy1 = u1(n1,:); yy2 = u2(n1,:);
plot(xxe,yy1e,'linewidth',lw,'linestyle','--','color',color1,'DisplayName','Exp: live cells'); hold on; grid on;
plot(xx,yy1,'linewidth',lw,'linestyle','-','color',color1,'DisplayName','Sim: live cells');
plot(xx,yy2,'linewidth',lw,'linestyle','-','color',color2,'DisplayName','Sim: dead cells');
ylabel('$u$ (cell/mL)','fontsize',s,'interpreter','latex')
xlabel('$x$ (cm)','fontsize',s,'interpreter','latex')
legend('fontsize',s,'interpreter','latex','location','north')
ylim([0 1.2*csat]);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
title('Day 7','fontsize',s,'interpreter','latex')

subplot(2,2,2); % Day 9
n2 = floor(10/20*nt);
yy1e = YY(:,5);
yy1 = u1(n2,:); yy2 = u2(n2,:);
plot(xxe,yy1e,'linewidth',lw,'linestyle','--','color',color1,'DisplayName','Exp: live cells'); hold on; grid on;
plot(xx,yy1,'linewidth',lw,'linestyle','-','color',color1,'DisplayName','Sim: live cells');
plot(xx,yy2,'linewidth',lw,'linestyle','-','color',color2,'DisplayName','Sim: dead cells');
ylabel('$u$ (cell/mL)','fontsize',s,'interpreter','latex')
xlabel('$x$ (cm)','fontsize',s,'interpreter','latex')
legend('fontsize',s,'interpreter','latex','location','north')
ylim([0 1.2*csat]);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
title('Day 11','fontsize',s,'interpreter','latex')

subplot(2,2,3); % Day 6
n3 = floor(16/20*nt); 
xxe = XX;
yy1e = YY(:,6);
yy1 = u1(n3,:); yy2 = u2(n3,:);
plot(xxe,yy1e,'linewidth',lw,'linestyle','--','color',color1,'DisplayName','Exp: live cells'); hold on; grid on;
plot(xx,yy1,'linewidth',lw,'linestyle','-','color',color1,'DisplayName','Sim: live cells');;
plot(xx,yy2,'linewidth',lw,'linestyle','-','color',color2,'DisplayName','Sim: dead cells');
ylabel('$u$ (cell/mL)','fontsize',s,'interpreter','latex')
xlabel('$x$ (cm)','fontsize',s,'interpreter','latex')
legend('fontsize',s,'interpreter','latex','location','north')
ylim([0 1.2*csat]);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
title('Day 17','fontsize',s,'interpreter','latex')

subplot(2,2,4); % Day 6
n4 = floor(nt); 
xxe = XX;
yy1e = YY(:,7);
yy1 = u1(n4,:); yy2 = u2(n4,:);
plot(xxe,yy1e,'linewidth',lw,'linestyle','--','color',color1,'DisplayName','Exp: live cells'); hold on; grid on;
plot(xx,yy1,'linewidth',lw,'linestyle','-','color',color1,'DisplayName','Sim: live cells');;
plot(xx,yy2,'linewidth',lw,'linestyle','-','color',color2,'DisplayName','Sim: dead cells');
ylabel('$u$ (cell/mL)','fontsize',s,'interpreter','latex')
xlabel('$x$ (cm)','fontsize',s,'interpreter','latex')
legend('fontsize',s,'interpreter','latex','location','north')
ylim([0 1.2*csat]);
ax = gca;
ax.FontSize = floor(0.9*s);
ax.TickLabelInterpreter = 'latex';
title('Day 21','fontsize',s,'interpreter','latex')