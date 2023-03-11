%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all
close all
clc

%% Parameters

% plant parameters
m = .5;         % mass 
k = 25;         % spring stiffness (pre-charged to cancel gravity)
c = 1.5;        % damping coefficient
zMax = 3 ;      % impact height (vertical axis is inverted)
lambda = 0.95;  % restitution coefficient

A = [0 1; -k/m -c/m]; 
B = [0; 1/m] ; 
C = eye(2);
D = [0;0];

% controller parameters
Va = 3.5;       % reference high
Vb = 2.5;       % reference low
Ts = 0.25;      % time to maintain Vb after jump

K = lqr(A,B,diag([25 2]),1);        % controller gain matrix
Acl = A - B*K ;                     % closed-loop system
ssGain = (D + C/(eye(2) - Acl)*B);  % steady state gain
Fw = 1/(ssGain(1)) ;                % feed forward with unitary dc gain

% estimator parameters
s = 0.02;               % sample period during flows
gammac = 1/(s*19^2);    % \psi_M for this example is 18.6
gammad = 1;
%gammad = 0;
theta = K.';

parameters.K = K;
parameters.Fw = Fw;
parameters.zMax = zMax;
parameters.lambda = lambda;
parameters.Va = Va;
parameters.Vb = Vb;
parameters.Ts = Ts;
parameters.A = A;
parameters.B = B;
parameters.theta = theta;
parameters.gammac = gammac;
parameters.gammad = gammad;
parameters.s = s;

%% Create hybrid system
sys = hybridDiscretizedGD(parameters);

tspan = [0, 5];
jspan = [0, 500];
config = HybridSolverConfig('RelTol', 1e-4, 'AbsTol', 1e-4,  ...
                            'MaxStep', s/4, 'Refine', 4);

z0 = [0; 0];        % pressure mounter state
tau0 = Ts;          % timer for reference command
thetahat0 = [0;0];  % parameter estimate
tauS0 = s;          % timer for samples during flows
x0 = [z0;tau0;thetahat0;tauS0];
sol = sys.solve(x0, tspan, jspan, config);

% compute reference command
v = Fw*Vb*(sol.x(:,3) < Ts) + Fw*Va*(sol.x(:,3) >= Ts);
v = HybridArc(sol.t, sol.j, v);

%% Plots
HybridPlotBuilder.defaults.set('flow line width', 3, ...
                               'jump line width', 2,...
                               'label size', 34,...
                               'tick label size', 20,...
                               't_label', '$t$ [s]')

figure; clf;
tlt = tiledlayout(2, 1);

nexttile(1)
hold on; grid on;
HybridPlotBuilder()....
    .flowColor('blue')...
    .jumpMarker('none')...
    .labels('$z_1$')...
    .tLabel('')...
    .plotFlows(sol.select(1))
xlim(tspan)
ylim([0 zMax])
set(gca,'Box','on');

nexttile(2)
hold on; grid on;
HybridPlotBuilder()....
    .flowColor('blue')...
    .jumpColor('blue')...
    .jumpMarker('none')...
    .jumpLineStyle('-')...
    .labels('$v$')...
    .plotFlows(v)
xlim(tspan)
ylim([66 102])
set(gca,'Box','on');

tlt.Padding = "none";
tlt.TileSpacing = "compact";
pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), 1.1*pos(4)])
movegui(gcf,'north')

%%
inD1 = (sol.x(:,1) >= zMax) & (sol.x(:,2) >= 0); % get points in the jump set
isFlowing = ~inD1;                               % get points during flows
isJumping = inD1 + cat(1,0,inD1(1:end-1));       % get points either side of each jump

figure; clf;
hold on; grid on;
HybridPlotBuilder().... % plot only the flows
    .flowColor('blue')...
    .jumpColor('none')...
    .labels('$|\tilde\theta_s|$')...
    .plotFlows(sol.select(4:5),@(x) norm(theta - x))
HybridPlotBuilder().... % plot the jumps due to samples during flows
    .flowColor('none')...
    .jumpColor('blue')...
    .jumpLineStyle('-')...
    .jumpMarker('none')...
    .labels('$|\tilde\theta_s|$')...
    .filter(isFlowing)...
    .plotFlows(sol.select(4:5),@(x) norm(theta - x))
HybridPlotBuilder().... % plot the jumps due to the jump set
    .flowColor('none')...
    .jumpColor('red')...
    .jumpMarker('none')...
    .labels('$|\tilde\theta_s|$')...
    .filter(isJumping)...
    .plotFlows(sol.select(4:5),@(x) norm(theta - x))%
xlim(tspan)
set(gca,'Box','on');

tlt.Padding = "none";
tlt.TileSpacing = "compact";
pos = get(gcf, 'Position');
set(gcf, 'Position',  [pos(1), pos(2), 1.8*pos(3), pos(4)])
movegui(gcf,'north')
set(gca, 'LooseInset', get(gca,'TightInset'))
