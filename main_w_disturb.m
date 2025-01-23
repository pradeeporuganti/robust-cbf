% Author: Pradeep Sharma Oruganti
clc; close all; clear all;

%% Import all folders for TIRA
%% addpath(genpath(''))

% TIRA stuff (making things simple)
global system_choice
global u
global g_sensitivity_bounds_method
bool_discrete_time = 0;
g_sensitivity_bounds_method = nan;
system_choice = 1;

%% Import CasADi
%% ADD CASADI PATH HERE
%% addpath('')
import casadi.*

%% General
dt = 0.1;         % Time step [s]
Tfinal = 25;      % Final simulation time (s)
t = 0:dt:Tfinal;  % Simulation time steps (s)
T = 0.1;          % sampling step  
% rng(127)        % For disturbance reproducibility
% rng(156)
rng(456)

% Time interval
t_init = 0;         % Initial time
t_final = T;        % Final time

% final_step = 10;

%% Robot params
u1min = -1;         % min steering rate [rad/s]
u1max = 1;          % max steering rate [rad/s]
u2min = -2;         % min acceleration  [m/s]
u2max = 2;          % max acceleration  [m/s]
D = 5;              % safety distance requirements

dmin = -0.3; mmin = -0.5;
dmax = 0.3; mmax = 0.5;

ulims = [u1min, u1max; u2min, u2max];
dlims = [dmin; dmax];
mlims = [mmin; mmax];

%% Sensor params
variance = 10;        % sensor parameters currently chosen at random
alpha = -3.86;         % thresholds
beta = 3.89;
P_sns = eye(4);

%% Environment params
xo = 32.5;          % obstacle x-position [m]
yo = 25;            % obstacle y-position [m]

%% Nonlinear system model

% Model states
s1 = MX.sym('s1');          % x-position [m]
s2 = MX.sym('s2');          % y-position [m]
v = MX.sym('v');            % velocity [m/s]
theta = MX.sym('theta');    % heading [rad]
u1 = MX.sym('u1');          % steering [rad/s] 
u2 = MX.sym('u2');          % accel [m/s^2]    

x = [s1; s2; theta; v];     % state vector
uc = [u1; u2];

% Nonlinear dynamics
xdot = [v*cos(theta);...
    v*sin(theta);...
    u1;...
    u2];

C = eye(4);

f = Function('f', {x ,uc}, {xdot});

%% Integrator

% Fixed step Runge-Kutta 4 integrator
M = 4;                   % RK4 steps per interval
DT = dt/4;
X0 = MX.sym('X0', 4);
U = MX.sym('U', 2);
X = X0;

for j=1:M
    
     k1 = f(X, U);
     k2 = f(X + DT/2 * k1, U);
     k3 = f(X + DT/2 * k2, U);
     k4 = f(X + DT * k3, U);
     
     X = X + DT/6*(k1 +2*k2 +2*k3 +k4);
     
end

F = Function('F', {X0, U}, {X}, {'x0', 'u_in'}, {'x_next'});

%% MPC problem setup
Q = 1e4*eye(4);
R = eye(2);
N = 10;

% Start with an empty NLP
w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];
P = {};

% "Lift" initial conditions
Xk = MX.sym('X0', 4);
w = {w{:}, Xk};
lbw = [lbw; 0; 0; -pi; 0];
ubw = [ubw; 50;  50; pi; 2];
w0 = [w0; 0; 0; 0; 0];

% Add parameters for initial state
Pk = MX.sym('P0',4);
P = {P{:}, Pk};
g = {g{:}, Pk-Xk};
lbg = [lbg; 0; 0; 0; 0];
ubg = [ubg; 0; 0; 0; 0];

% Formulate the NLP
for k=0:N-1
    
    % New NLP variable for control
    Uk = MX.sym(['U_' num2str(k)], 2);
    w = {w{:}, Uk};
    lbw = [lbw; u1min; u2min];
    ubw = [ubw;  u1max; u2max];
    w0 = [w0;  0; 0];
    
    % New NLP variable for parameters
    J = J+(Xk-[45;21;0;0])'*Q*(Xk-[45;21;0;0])+Uk'*R*Uk;
   
    % Next point
    Fk =  F('x0', Xk, 'u_in', Uk);
    Xk_next = Fk.x_next;

    % New NLP variable for state at end of interval
    Xk = MX.sym(['X_' num2str(k+1)], 4);
    w = {w{:}, Xk};
    lbw = [lbw; 0; 0; -pi; 0];
    ubw = [ubw; 50; 50; pi; 2];
    w0 = [w0; 0; 0; 0; 0];
    
    % Add equality constraint
    g = {g{:}, Xk_next-Xk};
    lbg = [lbg; 0; 0; 0; 0];
    ubg = [ubg; 0; 0; 0; 0];
    
end

% Create an NLP solver
opts=struct;
opts.print_time = 0;
opts.ipopt.print_level=0;
prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}), 'p', vertcat(P{:}));
solver = nlpsol('solver', 'ipopt', prob, opts);

%% Safety filter 
% Start with an empty NLP
w_sf={};
w0_sf = [];
lbw_sf = [];
ubw_sf= [];
J_sf = 0;
g_sf={};
lbg_sf = [];
ubg_sf = [];
P_sf = {};

% optim variable
u_sf = MX.sym('u_sf', 2);      % safe control input
w_sf = {w_sf{:}, u_sf};
lbw_sf = [lbw_sf; u1min; u2min];
ubw_sf = [ubw_sf; u1max; u2max];
w0_sf = [w0_sf; 0; 0];

% working combo 
p1 = 1; p2 = 1; % HOCBF penalties
q1 = 1; q2 = 1; 

% create constraint function
% h
h = (s1-xo)^2 + (s2-yo)^2 - D^2; 

% \psi_1 = Lfh + p1 h^q1
% psi_1 = (2*v*cos(theta)*(s1-xo) ) + (2*v*sin(theta)*(s2-yo)) + (p1*h^q1);
% 
% % Lf\psi_1
% Lf_psi_1 = (v*cos(theta)*(2*v*cos(theta) + p1*q1*h^(q1-1)*(2*(s1-xo)))) + ...
%     (v*sin(theta)*(2*v*sin(theta) + p1*q1*h^(q1-1)*(2*(s2-yo))));
% 
% % Lg\psi_1
% Lg_psi_1 = [-2*v*sin(theta)*(s1-xo) + 2*v*cos(theta)*(s2-yo);...
%     2*cos(theta)*(s1-xo) + 2*sin(theta)*(s2-yo)];

% psi_2 = Lf\psi_1 + Lg\psi_1 + p2 \psi_1^q2 >= 0;
% con = Lf_psi_1 + Lg_psi_1'*u_sf + p2*psi_1^q2;
con = 2*v^2 + 4*v*cos(theta)*(s1-xo) + 4*v*sin(theta)*(s2-yo) + h + ...
    [2*v*cos(theta)*(s2-yo) - 2*v*sin(theta)*(s1-xo); ...
    2*cos(theta)*(s1-xo) + 2*sin(theta)*(s2-yo)]'*u_sf;

% create constraint function 
f_con = Function('f_con', {x, u_sf}, {con});

%% Main sim with receding horizon control

t0 = 0;             % initial time
% x_init = [5; 25; 0; 0];
% x0 = [5; 25; 0; 0];
x_init = [20; 25; 0; 5];
x0 = [20; 25; 0; 5];
u0 = zeros(2,1);    % initial control
tsim = [];          % sim time
tsim = [tsim t0];   % first time instant

xx(:,1) = x0;       % nominal state trajectory
xx_true(:,1) = x0;  % true state trajectory
xx_meas(:,1) = x0;  % true state trajectory
xxOpt = [];         % solution trajectory during each optimization window
u_cl=[];            % final control input
u_cl(:, 1) = u0;
x0_true = x0;
n(:, 1) = zeros(4,1);
mpciter = 0;

u_safe=[];          % final safe control input
u_safe(:, 1) = u0;
didnt = 0;

args.w0 = w0;       % initil guess for optimization problem

% margin stuff
margins = [];

% Reachable sets
succ_lows = []; succ_ups = [];

% special noises
spec_d = zeros(1,Tfinal/dt+1);
spec_m = zeros(1,Tfinal/dt+1);

% % disturbance set 
W_vertex = [1, 1; 1, -1; -1, -1; -1, 1].*0.3; 
sys.W = Polyhedron(W_vertex);
% optim_output = [];

% noise
noises = [];
meas_noises = [];
trajs_x = []; trajs_y = [];

% disturbance
pd = makedist('Normal','mu', 0,'sigma', 0.5);
t = truncate(pd,dmin,dmax);
m = truncate(pd,mmin,mmax);


main_loop = tic;
while mpciter < Tfinal/dt
% while mpciter < 10
    
    args.p = x0;
    
%     mpc_tic = tic;
    % Solve N-step optim problem to get u_des N-step policy
    sol = solver('x0', args.w0, 'lbx', lbw, 'ubx', ubw,...
        'lbg', lbg, 'ubg', ubg,'p', args.p);
%     mpc_toc = toc(mpc_tic);
    
    % N*size(w) optim solution used as init for next optim cycle
    xxOpt(:,mpciter+1) = full(sol.x);
    u_cl(:, mpciter+2) = xxOpt(5:6,mpciter+1);   % receding horizon control input

    %% TIRA test
    % Interval of initial states (defined by two column vectors)
    x_low = [xx_meas(1, mpciter+1) + mlims(1); xx_meas(2, mpciter+1) + mlims(1); xx_meas(3, mpciter+1)];
    x_up = [xx_meas(1, mpciter+1) + mlims(2); xx_meas(2, mpciter+1) + mlims(2); xx_meas(3, mpciter+1)];

       % matched
%     u = [xx_meas(4, mpciter+1) + (0.5*(u2max+dmax)*T); u1max];       
    % mis-matched
    u = [xx_meas(4, mpciter+1) + (0.5*u2max*T); u1max];

    p_up = [dmax;dmax;0]; p_low = [dmin;dmin;0];              % Interval of allowed input values (disturbances)
%     p_up = [0;0;dmax]; p_low = [0;0;dmin];

    % Call of the main over-approximation function
    [succ_low_umax, succ_up_umax] = TIRA([t_init, t_final],x_low,x_up,p_low,p_up);
    
    % matched
%     u = [xx_meas(4, mpciter+1) + (0.5*(u2min+dmin)*T); u1min];
    % mismatched
    u = [xx_meas(4, mpciter+1) + (0.5*u2min*T); u1min];

    [succ_low_umin,succ_up_umin] = TIRA([t_init, t_final],x_low,x_up,p_low,p_up);

    scl = [succ_low_umin  succ_low_umax];
    scu = [succ_up_umin succ_up_umax];

    succ_low = min(scl, [], 2);
    succ_up = max(scu, [], 2);
    
    succ_lows = [succ_lows succ_low]; succ_ups = [succ_ups succ_up]; % Reachable set overapproximation

%     norm_LY_psi_1 = norm_max(succ_low_umax, succ_up_umax, xx_true(:, mpciter+1));
    margin = drd_reach_sup_mis(succ_low, succ_up, ulims, dlims, xx_true(:, mpciter+1));   %   mismatch
%     margin = reach_sup(succ_low, succ_up, ulims, dlims, xx_true(:, mpciter+1));     % match
    margins = [margins margin];

    %% Safety filter
    % generate constraint based on current states
    f_con_val = full(f_con(x0, u_sf));
    
    % constraint for optimization
    g_sf = {f_con_val + margin};
    lbg_sf = 0;
    ubg_sf = inf; 
    
    % objective for safety filter optim problem
    J_sf = 0.5*(u_sf - u_cl(:, mpciter+2))'*eye(2)*(u_sf - u_cl(:, mpciter+2));

    % Create an NLP solver
    opts2=struct;
    opts2.print_time = 0;
    opts2.printLevel = 'none';
    prob_safe = struct('f', J_sf, 'x', vertcat(w_sf{:}), 'g', vertcat(g_sf{:}));
    safe_solver = qpsol('safe_solver', 'qpoases', prob_safe, opts2);
    
    safe_sol = safe_solver('x0', w0_sf, 'lbx', lbw_sf, 'ubx', ubw_sf,...
        'lbg', lbg_sf, 'ubg', ubg_sf);
    u_sf_Opt = full(safe_sol.x);
    u_safe(:, mpciter+2) = u_sf_Opt;
    w0_sf = u_safe(:, mpciter+2);
 
    %% Control application
    % simulate dynamics
%     xnext = F('x0', x0, 'u_in', u_safe(:, mpciter+2));
    xnext = F('x0', xx_true(:, mpciter+1), 'u_in', u_safe(:, mpciter+2));
    xx(:,mpciter+2) = full(xnext.x_next);

    % external disturbance
    noise = random(t,4,1);
    apps = [1 0 0 0;...
        0 1 0 0;...
        0 0 0 0;...
        0 0 0 0]*noise;       % mismatch
%         apps = [0 0 0 0;...
%         0 0 0 0;...
%         0 0 1 0;...
%         0 0 0 1]*noise;       % match
    noises = [noises apps];

    meas_noise = random(m,4,1);
    meas_apps = [1 0 0 0;...
        0 1 0 0;...
        0 0 0 0;...
        0 0 0 0]*meas_noise;
    meas_noises = [meas_noises meas_apps];
    
    xx_true(:, mpciter+2) = xx(:,mpciter+2) + apps;
    xx_meas(:, mpciter+2) = xx_true(:,mpciter+2) + meas_apps;

    % guess for next time step
    args.w0 = xxOpt(:,mpciter+1);
    % initial condition for next time step
    x0 = xx_meas(:, mpciter+2);

    % sim time
    tsim(mpciter+2) = tsim(mpciter+1)+dt;
    mpciter = mpciter + 1;
    
end
main_loop = toc(main_loop);

%% Randomize 

% figure(2)
% for i = 1:length(trajs_x(1,:))
%     hold on
%     plot(trajs_x(:,i), trajs_y(:,i));
% end 
% xlabel('X [m]');ylabel('Y [m]');
% xlim([0; 50]); ylim([0; 50]); 
% viscircles([xo yo], D, 'Color','k');

%% plots

barrier = [];

% x1 = 0:0.5:50;
% y1 = 0:0.5:50;
% [X1, Y1] = meshgrid(x1, y1);

% plot trajectory and obstacle
figure(1)
plot(xx_true(1,:), xx_true(2,:), 'b-', 'MarkerSize', 7,'LineWidth',1.5);
hold on
plot(xx_true(1,:), xx_true(2,:), 'k-', 'MarkerSize', 7,'LineWidth',1.5);
for i = 1:length(xx_true)
%     plot([xx_true(1,i); xx_true(2,i)] + sys.W, 'alpha', 0.05);
    barrier = [barrier ((xx_true(1,i) - xo)^2 + (xx_true(2,i)-yo)^2 - D^2)];
end
grid on; box on; hold on
xlabel('X [m]');ylabel('Y [m]');
xlim([0; 50]); ylim([0; 50]); 
viscircles([xo yo], D, 'Color','k');
for i = 2:length(succ_lows)
    
    succ_low = succ_lows(:, i-1);
    succ_up = succ_ups(:, i-1); 
    figure(1)   
    hold on
    plot([succ_low(1) succ_up(1) succ_up(1) succ_low(1) succ_low(1)], [succ_low(2) succ_low(2) succ_up(2) succ_up(2) succ_low(2)], 'k-');

end

% plot reachable sets
figure(3)
plot3(tsim, xx_true(1,:), xx_true(2,:),'lineWidth', 3)
box on;

for i = 2:length(succ_lows)
    
    succ_low = succ_lows(:, i-1);
    succ_up = succ_ups(:, i-1); 
    figure(3)   
    hold on
    plot3(tsim(i).*[1 1 1 1 1], [succ_low(1) succ_up(1) succ_up(1) succ_low(1) succ_low(1)], [succ_low(2) succ_low(2) succ_up(2) succ_up(2) succ_low(2)], 'k-');

end
xlabel('Simulation time [s]', 'FontSize',45); ylabel('$x$ [m]', 'interpreter', 'latex', 'FontSize',45); 
zlabel('$y$ [m]','interpreter', 'latex', 'FontSize',45);
legend('Robot trajectory', 'Estimated reachable set');

% plot barrier function
c = barrier > 0;

safe_barr = barrier;
unsafe_barr = barrier;

safe_barr(~c) = nan;
unsafe_barr(c) = nan;

figure(2)
hold on
plot(tsim, safe_barr, 'k', tsim, unsafe_barr, 'r', 'lineWidth', 1.5);
grid on; box on;
xlabel('Simulation time [s]'); ylabel('$h(x)$','Interpreter','latex')
%%


%% plots

% barrier = [];
% 
% % x1 = 0:0.5:50;
% % y1 = 0:0.5:50;
% % [X1, Y1] = meshgrid(x1, y1);
% 
% % plot trajectory and obstacle
% figure(1)
% plot(xx_true(1,:), xx_true(2,:), 'b-', 'MarkerSize', 7,'LineWidth',1.5);
% hold on
% plot(unsafe.xx_true(1,:), unsafe.xx_true(2,:), 'k-', 'MarkerSize', 7,'LineWidth',1.5);
% for i = 1:length(xx_true)
% %     plot([xx_true(1,i); xx_true(2,i)] + sys.W, 'alpha', 0.05);
%     barrier = [barrier ((xx_true(1,i) - xo)^2 + (xx_true(2,i)-yo)^2 - D^2)];
% end
% grid on; box on; hold on
% xlabel('X [m]');ylabel('Y [m]');
% xlim([0; 50]); ylim([0; 50]); 
% viscircles([xo yo], D, 'Color','k');
% 
% % plot reachable sets
% figure(3)
% plot3(tsim, xx_true(1,:), xx_true(2,:),'lineWidth', 3)
% box on;
% 
% for i = 2:length(succ_lows)
%     
%     succ_low = succ_lows(:, i-1);
%     succ_up = succ_ups(:, i-1); 
%     figure(3)   
%     hold on
%     plot3(tsim(i).*[1 1 1 1 1], [succ_low(1) succ_up(1) succ_up(1) succ_low(1) succ_low(1)], [succ_low(2) succ_low(2) succ_up(2) succ_up(2) succ_low(2)], 'k-');
% 
% end
% xlabel('Simulation time [s]', 'FontSize',31); ylabel('$x$ [m]', 'interpreter', 'latex', 'FontSize',31); 
% zlabel('$y$ [m]','interpreter', 'latex', 'FontSize',31);
% 
% % plot barrier function
% c = unsafe.barrier > 0;
% 
% safe_barr = unsafe.barrier;
% unsafe_barr = unsafe.barrier;
% 
% safe_barr(~c) = nan;
% unsafe_barr(c) = nan;
% 
% figure(2)
% hold on
% plot(tsim, safe_barr, 'k', tsim, unsafe_barr, 'r', 'lineWidth', 1.5);
% grid on; box on;
% xlabel('Simulation time [s]'); ylabel('$h(x)$','Interpreter','latex')
%%
