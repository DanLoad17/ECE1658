clearvars
close all
clc

%2.1
l1 = 0.5;
l2 = 0.5;
l3 = 0.3;
m1 = 0.05;
m5 = 0.05;
m2 = 0.5;
m4 = 0.5;
m6 = 0.5;
m3 = 0.3;
g = 9.81;

%2.2
syms t q1 q2 q3 q4 q5 x1 x2 q1dot q2dot q3dot q4dot q5dot x1dot x2dot real
syms tau1 tau2 tau3 tau4 real

q = [q1; q2; q3; q4; q5];
qdot = [q1dot; q2dot; q3dot; q4dot; q5dot];
tau = [tau1; tau2; tau3; tau4];

qbar = [q; x1; x2];
qbardot = [qdot; x1dot; x2dot];

%2.3a - assume to find center
theta1 = q1;                     
theta2 = q1 + q2;                
theta3 = q1 + q2 + q3;           
theta4 = q1 + q2 + q3 + q4;      
theta5 = q1 + q2 + q3 + q4 + q5; 

% r1: stance foot (m1)
r1 = [x1; x2];
% r2: end of stance shank (m2), distance l1 from r1 along theta1
r2 = r1 + l1*[cos(theta1); sin(theta1)];
% r3: end of stance thigh / hip (m3), distance l2 from r2 along theta2
r3 = r2 + l2*[cos(theta2); sin(theta2)];
% r6: torso top (m6), distance l3 from r3 along theta3
r6 = r3 + l3*[cos(theta3); sin(theta3)];
% r4: end of swing thigh (m4), distance l2 from hip r3 along theta4
r4 = r3 + l2*[cos(theta4); sin(theta4)];
% r5: swing foot (m5), distance l1 from r4 along theta5
r5 = r4 + l1*[cos(theta5); sin(theta5)];

%2.3b
r1dot = jacobian(r1, qbar) * qbardot;
r2dot = jacobian(r2, qbar) * qbardot;
r3dot = jacobian(r3, qbar) * qbardot;
r4dot = jacobian(r4, qbar) * qbardot;
r5dot = jacobian(r5, qbar) * qbardot;
r6dot = jacobian(r6, qbar) * qbardot;

%2.3c Kinetic energy
K = (1/2) * ( ...
    m1 * (r1dot.' * r1dot) + ...
    m2 * (r2dot.' * r2dot) + ...
    m3 * (r3dot.' * r3dot) + ...
    m4 * (r4dot.' * r4dot) + ...
    m5 * (r5dot.' * r5dot) + ...
    m6 * (r6dot.' * r6dot) );

%2.3d
Dbar = hessian(K, qbardot);

%2.3e
D = Dbar(1:5, 1:5);

%2.4
P = m1*g*r1(2) + m2*g*r2(2) + m3*g*r3(2) + ...
    m4*g*r4(2) + m5*g*r5(2) + m6*g*r6(2);

G = jacobian(P,q)';

%2.5
B = [ 0 0 0 0;   % q1 unactuated
      1 0 0 0;   % q2 actuated
      0 1 0 0;   % q3 actuated
      0 0 1 0;   % q4 actuated
      0 0 0 1 ]; % q5 actuated

%2.6
C = sym(zeros(5,5));

for i = 1:5
    for j = 1:5
        for k = 1:5
            C(i,j) = C(i,j) + 1/2 * ( ...
                diff(D(i,j), q(k)) + ...
                diff(D(i,k), q(j)) - ...
                diff(D(j,k), q(i)) ) * qdot(k);
        end
    end
end

%2.7
px = r5 - [x1; x2];  

dpx_q = jacobian(px, qbar);   
E = [dpx_q; [zeros(2,5) eye(2)]];

%2.8
Dfun = matlabFunction(D, 'Vars', {q});
Dbarfun = matlabFunction(Dbar, 'Vars', {qbar});
pxfun = matlabFunction(px, 'Vars', {qbar});
Efun = matlabFunction(E, 'Vars', {qbar});
Cfun = matlabFunction(C, 'Vars', {q, qdot});
Gfun = matlabFunction(G, 'Vars', {q});

%2.9
data = struct();

data.D     = Dfun;      % Inertia matrix of pinned robot
data.Dbar  = Dbarfun;   % Inertia matrix of unpinned robot
data.E     = Efun;      % Impact map helper matrix
data.C     = Cfun;      % Coriolis matrix
data.G     = Gfun;      % Gravity vector

% Store the input matrix B
data.B     = B;

%3.1
