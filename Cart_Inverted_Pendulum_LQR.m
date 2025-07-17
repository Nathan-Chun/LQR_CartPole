clear all; clc;
nx = 4;
nu = 1;

% Target state
xf = [0; 0; 0; 0]; % p, q, dp, dq

% IC
x0 = [0 0 0 0] + [0, 35 0, 0];

% state space matrices
M = 0.5;
m=0.5;
l=0.3;
l=0.3;
g=9.81;

A = [0 0 1 0;
    0 0 0 1;
    0 m*g/M 0 0;
    0 (m*g + M*g)/(M*l) 0 0];
B = [0;0;1/M;1/(M*l)];
C = eye(4);
D = 0;

% Check controllability (rank = # of states)
controllability = rank(ctrb(A, B))

% Stability (has one positive and one negative, so system is unstable)
eig(A)

% ss system
cartpole = ss(A,B,C,D);

Q = diag([100 100 10 10]);
R = 1; 

K = lqr(cartpole, Q, R)
% K = ihlqr(A,B,Q,R);
% K = [K(1), K(2), K(3), K(4)]

%% Infinite horizon LQR
% Demonstrates Ricatti Recursion, but doesn't actually work, because
% of linearization of dynamics and difference between internals of MATLAB
% lqr() function and this ihlqr function
function [S, K] = ihlqr(A, B, Q, R, max_iter, tol)
    % Default iter and tol
    if nargin<6
        tol = 1e-6;
    end
    if nargin<5
        max_iter = 1000;
    end

    % Declare S 
    S = Q*1.0;

    for ricatti_iter = 1:max_iter
        K = (R+B'*S*B)\(B'*S*A);
        newS = Q + K'*R*K + (A - B*K)'*S*(A - B*K);
        if norm(S-newS, Inf) < tol
            fprintf("Converged after %d iterations\n", ricatti_iter);
            S = newS; 
            K = reshape(K,1,[]);
            return;
        end
        S = newS;
    end

    error("Did not converge")
end
%% Disregard 

% % Dynamics
% % M = cart mass; m = ball mass
% function dx = dynamics(params, x, u)
%     M = params(1);
%     m = params(2);
%     l = params(3);
%     g = 9.81;
%     q = x(1:2);
%     qd = x(3:4);
% 
%     % Trig shorthand
%     s = sin(q(2));
%     c = cos(q(2));
% 
%     % Cart pole manipulator dynamics
%     H = [M+m, m*l*c; m*l*c, m*l^2];
%     C = [0, -m*qd(2)*l*s; 0, 0];
%     G = [0; m*g*l*s];
%     B = [1; 0];
% 
%     qdd = -H\(C*qd + G - B*u(1));
%         size(qd);
%     dx = [qd(1); qd(2); qdd(1); qdd(2)];
% 
% end

% % Runga Kutta integration for sim
% function xnew = rk4(f, x, u, dt)
%     k1 = dt*f(x,u);
%     k2 = dt * f(x + k1/2, u);
%     k3 = dt * f(x + k2/2, u);
%     k4 = dt * f(x + k3, u);
%     xnew = x + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
%     % disp(xnew)
% end
