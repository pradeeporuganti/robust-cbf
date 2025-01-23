function margin = drd_reach_sup_mis(succ_low, succ_up, ulims, dlims, x_hat)

%% Import CasADi
%% addpath('')
import casadi.*

%% Environment params
xo = 32.5;          % obstacle x-position [m]
yo = 25;            % obstacle y-position [m]

%% Nonlinear system model

% Model states
s1 = MX.sym('s1');          % x-position [m]
s2 = MX.sym('s2');          % y-position [m]
v = MX.sym('v');            % velocity [m/s]
theta = MX.sym('theta');    % heading [rad]
ucon = MX.sym('ucon', 2);   % control inputs
d = MX.sym('d', 2);         % distrubance

x = [s1; s2; theta; v];     % state vector

u1min = ulims(1,1);         % min steering rate [rad/s]
u1max = ulims(1,2);         % max steering rate [rad/s]
u2min = ulims(2,1);         % min acceleration  [m/s]
u2max = ulims(2,2);         % max acceleration  [m/s]
D = 5;              % safety distance requirements
T = 0.1;


%% Tests
% % % working paramters combo 
% p1 = 1; p2 = 1; % HOCBF penalties
% q1 = 1; q2 = 1; 

% % h
% h = (s1-xo)^2 + (s2-yo)^2 - D^2; 
% 
% % Lfh
% Lfh = (2*v*cos(theta)*(s1-xo) ) + (2*v*sin(theta)*(s2-yo));
% 
% % \psi_1
% psi_1 = Lfh + h + [2*(s1-xo); 2*(s2-yo)]'*d;
% 
% % Lf2h
% Lf2h = 2*v^2;
% 
% %LgLfh
% LgLfh = [-2*v*sin(theta)*(s1-x0) + 2*v*cos(theta)*(s2-yo); ...
%     2*cos(theta)*(s1-xo) + 2*sin(theta)*(s2-yo)]';
% 
% %LfLYh
% LfLYh = [2*v*cos(theta);2*v*sin(theta)]';

% % \psi_1
% psi_1 = (2*v*cos(theta)*(s1-xo)) + (2*v*sin(theta)*(s2-yo)) + h;
% 
% %Lfpsi_1
% Lfpsi_1 = (v*cos(theta)*(2*v*cos(theta) + 2*(s1-xo))) + (v*sin(theta)*(2*v*sin(theta) + 2*(s2-yo)));
% 
% % Lgpsi_1
% Lgpsi_1 = [-2*v*sin(theta)*(s1-xo) + 2*v*cos(theta)*(s2-yo);...
%     2*cos(theta)*(s1-xo) + 2*sin(theta)*(s2-yo)];
% 
% % LYpsi_1
% LYpsi_1 = [2*v*cos(theta)+ 2*(s1-xo); 
%     2*v*sin(theta)+2*(s2-yo)];

% f_psi_1 = Function('f_psi_1', {x}, {psi_1});
% f_Lfpsi_1 = Function('f_Lfpsi_1', {x}, {Lfpsi_1});
% f_Lgpsi_1 = Function('f_Lgpsi_1', {x}, {Lgpsi_1});
% f_LYpsi_1 = Function('f_LYpsi_1', {x}, {LYpsi_1});
%% 

% psi_1 = 2*v*(s1 - xo)*cos(theta) + 2*v*(s2 - yo)*sin(theta) + 2*(s1 - xo)*d(1) + 2*(s2 - yo)*d(2);
% 
% con = (2*v*cos(theta) + 2*d(1))*v*cos(theta) + (2*v*sin(theta) + 2*d(2))*v*sin(theta) + ....
%     (2*v*(s2 - yo)*cos(theta) - 2*v*(s1 - xo)*sin(theta))*ucon(1) + ....
%     (2*(s2 - yo)*sin(theta) + 2*(s1 - xo)*cos(theta))*ucon(2) + ....
%     (2*v*cos(theta) + 2*d(1))*d(1) + (2*v*sin(theta) + 2*d(2))*d(2)  + psi_1;

h = (s1-xo)^2 + (s2-yo)^2 - D^2;

con = 2*v^2 + 4*v*cos(theta)*(s1-xo) + 4*v*sin(theta)*(s2-yo) + h + ...
    (2*v*cos(theta)*(s2-yo) - 2*v*sin(theta)*(s1-xo))*ucon(1) + ...
    (2*cos(theta)*(s1-xo) + 2*sin(theta)*(s2-yo))*ucon(2);

th = 2*v*d(2)*sin(theta) + 2*v*d(1)*cos(theta) + 2*d(1)^2 + 2*d(2)^2 + 4*d(1)*(s1-xo) + 4*d(2)*(s2-yo);

f_con = Function('f_con', {x, ucon, d}, {con});
f_th = Function('f_th', {x,d}, {th});

w = {};
w0 = [];
lbw = [];
ubw = [];

w = {w{:}, x, ucon, d};
lbw = [lbw; succ_low(1); succ_low(2); succ_low(3); x_hat(4) + (0.5*u2min*T); u1min; u2min; dlims(1); dlims(1)];
ubw = [ubw; succ_up(1); succ_up(2); succ_up(3); x_hat(4) + (0.5*u2max*T); u1max; u2max; dlims(2); dlims(2)];  
w0 = [w0; 0; 0; 0; 0; 0; 0; 0; 0];

% w2 = {};
% w20 = [];
% lbw2 = [];
% ubw2 = [];
% 
% w2 = {w2{:}, x, d};
% lbw2 = [lbw2; succ_low(1); succ_low(2); succ_low(3); x_hat(4) + (0.5*u2min*T); dlims(1); dlims(1)];
% ubw2 = [ubw2; succ_up(1); succ_up(2); succ_up(3); x_hat(4) + (0.5*u2max*T); dlims(2); dlims(2)];  
% w20 = [w20; 0; 0; 0; 0; 0; 0];

% J = -norm(LYpsi_1);
% J = -( (-full(f_Lfpsi_1(x_hat)) + Lfpsi_1) + ...
%     (-full(f_Lgpsi_1(x_hat)) + Lgpsi_1)'*ucon + ...
%     (-f_LYpsi_1(x_hat) + LYpsi_1)'*d + ...
%     (-full(f_psi_1(x_hat)) + psi_1) );

J = - (f_con(x, ucon, d) - f_con(x_hat, ucon, d) + f_th(x, d));

% J2 = - norm(f_th(x, d));

opts=struct;
opts.print_time = 0; 
opts.ipopt.print_level=0;

prob = struct('f', J, 'x', vertcat(w{:}));
solver = nlpsol('solver', 'ipopt', prob, opts);

% prob2 = struct('f', J2, 'x', vertcat(w2{:}));
% solver2 = nlpsol('solver', 'ipopt', prob2, opts);

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw);
% sol2 = solver2('x0', w20, 'lbx', lbw2, 'ubx', ubw2);

margin = full(sol.f);
% margin2 = full(sol2.f);

end

