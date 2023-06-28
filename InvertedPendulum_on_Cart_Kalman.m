clear; close all; clc
%% Initialization
m = 1;
M = 5;
Kf = 2;
g = -10;
d = 1;

s = -1;

% y = [x; dx; theta; dtheta];
A = [0 1 0 0;
    0 -d/M -m*g/M 0;
    0 0 0 1;
    0 -s*d/(M*Kf) -s*(m+M)*g/(M*Kf) 0];

B = [0; 1/M; 0; s*1/(M*Kf)];

C = [1 0 0 0];  

D = zeros(size(C,1),size(B,2));

Q = eye(4)*10;
R = 0.1;

%% Stability and Controllability
disp('Eigen Values:')
disp(eig(A))
disp('Rank of Controllability Matrix:')
disp(rank(ctrb(A,B)))

%% Feedback Gain
% K = lqr(A,B,Q,R);

% [X1,K,L1] = icare(A,B,Q,R);
K = [0 0 0 0];
% [X1,K,L1] = icare(A,B,Q,R,'anti');

% desired_poles = [-1; -2; -3; -4];
% K = place(A,B,desired_poles);

% K = [-10.0000  -24.4555  284.5793  123.0397];

%% Initial Conditions and Desired States (x, x_dot, theta, theta_dot)
Initial_States = [-2; 0; pi-5*pi/6; 0];
Desired_States = [2; 0; pi; 0];

%% Simulation
% sampling_time = 0.1;
% total_time = 0:sampling_time:10;
% 
% [t,state] = ode45(@(t,y)cart_pend_diff(y,m,M,Kf,g,d,-K*(y-Desired_States)),total_time,Initial_States);
% 
% plot(t,state(:,1))
% hold on
% plot(t,state(:,2))
% plot(t,state(:,3))
% plot(t,state(:,4))
% title('LQR Controller for an Inverted Pendulum')
% legend('Position (m)','Velocity (m/s)','Angle (deg)','Angular Velocity (deg/s)')
% hold off
% pause(1.5)
% 
% figure
% 
% time = numel(total_time);
% draw_cart_pendulum(state,time,sampling_time)

%%  Augment system with disturbances and noise
Vd = .1*eye(4);  % disturbance covariance
Vn = 1;       % noise covariance

BF = [B Vd 0*B];  

mysys = ss(A,BF,C,[0 0 0 0 0 Vn]); 

mysys_full = ss(A,BF,eye(4),zeros(4,size(BF,2)));  

%%  Build Kalman filter
[Kf,P,E] = lqe(A,Vd,C,Vd,Vn); 
% Kf = (lqr(A',C',Vd,Vn))';   

mysys_KF = ss(A-Kf*C,[B Kf],eye(4),0*[B Kf]);  % Kalman filter estimator

%%  Estimate linearized system in down position
dt = .01;
t = dt:dt:50;

uDIST = randn(4,size(t,2));
uNOISE = randn(size(t));
u = 0*t;
u(100:120) = 100;     
u(1500:1520) = -100;  

uAUG = [u; Vd*Vd*uDIST; uNOISE];

[y,t] = lsim(mysys,uAUG,t);
[xtrue,t] = lsim(mysys_full,uAUG,t);
[x,t] = lsim(mysys_KF,[u; y'],t);
figure
plot(t,xtrue,'-',t,x,'--','LineWidth',2)
legend('x', 'dx', 'theta', 'dtheta', 'x(KF)', 'dx(KF)', 'theta(KF)', 'dtheta(KF)')

figure
plot(t,y)
hold on
plot(t,xtrue(:,1),'r')
plot(t,x(:,1),'k--')
legend("Noisy Measurement","True Measurement","Kalman Filter Estimation")

%% Differential Equations
function dy = cart_pend_diff(y,m,M,L,g,d,u)

Sy = sin(y(3));
Cy = cos(y(3));
D = m*L*L*(M+m*(1-Cy^2));

dy(1,1) = y(2);
dy(2,1) = (1/D)*(-m^2*L^2*g*Cy*Sy + m*L^2*(m*L*y(4)^2*Sy - d*y(2))) + m*L*L*(1/D)*u;
dy(3,1) = y(4);
dy(4,1) = (1/D)*((m+M)*m*g*L*Sy - m*L*Cy*(m*L*y(4)^2*Sy - d*y(2))) - m*L*Cy*(1/D)*u +.01*randn;
end

