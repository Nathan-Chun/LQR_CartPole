clear all; clc;
nx = 4;
nu = 1;

% Target state
xf = [0; pi; 0; 0];
% xf = [0; 0; 0; 0];

uf = [0];

% IC
x0 = [0 0 0 0] + [1.5, deg2rad(-20), .3, 0];
x0 = x0(:);

% Simulation params
dt = 0.1;
tf = 5.0;
t_vec = 0: dt: tf;
N = length(t_vec);
X = cell(1,N);
for i = 1:N
    X{i} = zeros(nx);
end
X{1} = x0;

% Physical system params [M, m, l]
params = [1.2, 0.16, 0.55];

Q = diag([1, 1, .05, .1]);
R = 0.05*diag(ones(nu));

% Dynamics function
f = @(x,u) dynamics(params, x, u);

% Jacobians
n = length(xf);
A = zeros(n);
for i = 1:n
    dx = zeros(n,1);
    dx(i) = 1e-6;  % small perturbation
    f1 = rk4(f, xf + dx, uf, dt);
    f0 = rk4(f, xf - dx, uf, dt);
    A(:,i) = (f1 - f0) / (2e-6);  % central diff
end

m = length(uf);
B = zeros(n, m);
for i = 1:m
    du = zeros(m,1);
    du(i) = 1e-6;
    f1 = rk4(f, xf, uf + du, dt);
    f0 = rk4(f, xf, uf - du, dt);
    B(:,i) = (f1 - f0) / (2e-6);
end

[S, K] = ihlqr(A, B, Q, R);

% Dynamics
% M = cart mass; m = ball mass
function dx = dynamics(params, x, u)
    M = params(1);
    m = params(2);
    l = params(3);
    g = 9.81;
    q = x(1:2);
    qd = x(3:4);

    % Trig shorthand
    s = sin(q(2));
    c = cos(q(2));

    % Cart pole manipulator dynamics
    H = [M+m, m*l*c; m*l*c, m*l^2];
    C = [0, -m*qd(2)*l*s; 0, 0];
    G = [0; m*g*l*s];
    B = [1; 0];

    qdd = -H\(C*qd + G - B*u(1));
        size(qd)
    dx = [qd(1); qd(2); qdd(1); qdd(2)];

end


% Runga Kutta integration for sim
function xnew = rk4(f, x, u, dt)
    k1 = dt*f(x,u);
    k2 = dt * f(x + k1/2, u);
    k3 = dt * f(x + k2/2, u);
    k4 = dt * f(x + k3, u);
    xnew = x + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    % disp(xnew)
end

% Infinite horizon LQR
function [S, K] = ihlqr(A, B, Q, R, max_iter, tol)
    % Default iter and tol
    if nargin<6
        tol = 1e-6;
    end
    if nargin<5
        max_iter = 1000;
    end
    % x and u sizes based on B
    [nx, nu] = size(B);
    
    % Declare S 
    S = Q*1.0;
    
    for ricatti_iter = 1:max_iter
        K = (R+B'*S*B)\(B'*S*A);
        newS = Q + K'*R*K + (A - B*K)'*S*(A - B*K);
        if norm(S-newS, Inf) < tol
            fprintf("Converged after %d iterations\n", ricatti_iter);
            S = newS; 
            return;
        end
        S = newS;
    end
    error("Did not converge")
end



%State Space
% A=[0 1 0 0; 
%     0 0 -m2*g/m1 0; 
%     0 0 0 1; 
%     0 0 (m1+m2)*g/(l*m1) 0];
% B=[0; 1/m1; 0; -1/(l*m1)];
% % A=[0 1 0 0; 
% %     0 0 0 0; 
% %     0 0 0 1; 
% %     0 0 2*g/l 0];
% % B=[0; F/(m1+m2); 0; 0];
% C=[0 0 1 0];
% D=0;
% sys=ss(A,B,C,D)
% %initial(sys,[0 0 pi 0])
% %step(sys)
% 
% %Check Controllability
% rank(ctrb(A, B))
% 
% %eigenvals
% eig(A)
% 
% %lqr
% Q=[1 0 0 0;
%     0 1 0 0;
%     0 0 10 0;
%     0 0 0 100];
% R=.0001;
% 
% K= lqr(A, B, Q, R)
% Acl=A-B*K
% sys1 = ss(Acl, B, C, D)
% step(sys1)
% 
