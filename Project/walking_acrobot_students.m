%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ECE1658
%% FINAL PROJECT
%% VHC for walking acrobot
%% PARTIAL CODE TO GET STUDENTS STARTED ON THE PROJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

number_steps=10;
symbolic_computations=1;
%% Physical parameters
l=1;
lc=0.5;
m=1;
Iz=1/12*m*l^2;
g=9.81;
psi=deg2rad(2);
%% Control parameters
Kp=20^2;
Kd=20*2;

%% Symbolic computations
if symbolic_computations
    % Define symbolic variables
    fprintf('\n Initializing symbolic variables...\n')
    % syms l lc Iz m g real
    syms t q1 q2 x1 x2 q1dot q2dot x1dot x2dot tau real
    
    q=[q1;q2];x=[x1;x2]; qbar=[q;x];
    qdot=[q1dot;q2dot];  xdot=[x1dot;x2dot]; qbardot=[qdot;xdot];
    
    fprintf('\n Symbolic model computation...\n')
    
    % Define centres of mass of two links
    rc1=x+lc*[cos(q1);sin(q1)];
    rc2=x+l*[cos(q1);sin(q1)]+lc*[cos(q1+q2);sin(q1+q2)];
    
    % Compute time derivatives of centres of mass
    rc1dot=jacobian(rc1,qbar)*qbardot;
    rc2dot=jacobian(rc2,qbar)*qbardot;
    
    % Define the total kinetic energy of the robot
    K=1/2*m*(rc1dot'*rc1dot+rc2dot'*rc2dot)+1/2*Iz*(q1dot^2+(q1dot+q2dot)^2);
    K=simplify(K);
    
    % Extract the square symmetric matrix of the kinetic energy
    Dbar=simplify(hessian(K,qbardot));
    
    % Extract the matrix of the kinetic energy of the pinned robot
    D = Dbar(1:2,1:2);
    
    % Define the potential energy of the pinnedrobot
    P = m*g*(lc*sin(q1)+l*sin(q1)+lc*sin(q1+q2));
    psi=deg2rad(incline_degrees);
    gvector=[cos(psi) -sin(psi);sin(psi) cos(psi)]*[0; -1];
    P=-(m*g*lc*[cos(q1);sin(q1)]+m*g*[l*cos(q1)+lc*cos(q1+q2);l*sin(q1)+lc*sin(q1+q2)])'*gvector;
    P=simplify(P);
    % Input matrix of pinned robot
    B=[0;1];
    
    % Computation of matrix C(q,qdot) of pinned robot
    C = sym(zeros(2,2));
    for i=1:2
        for j=1:2
            for k=1:2
                C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*qdot(i);
            end
        end
    end
    
    % Computation of gradient of the potential for the pinned robot
    
    G = jacobian(P,q)';
    
    % Computation of qddot (pinned robot)
    
    %     qddot = simplify(inv(D)*(-C*qdot-G+B*tau));
    Dfun=matlabFunction(D,'Vars',{q});
    Gfun=matlabFunction(G,'Vars',{q});
    Cfun=matlabFunction(C,'Vars',{[q;qdot]});
    fprintf('\n Impact map computation...\n')
    
    %% Impact map
    
    % Relabelling map
    
    T=-eye(2)*q + [pi;2*pi];
    Delta2=[T;-eye(2)*qdot];
    
    % First impact map
    
    px=l*[cos(q1)+cos(q1+q2);sin(q1)+sin(q1+q2)];
    xo=sym(zeros(2,1));
    
    E=[jacobian(px,q), eye(2)];
    Deltaqdot=[eye(2),zeros(2,4)]*inv([Dbar -E';E zeros(2,2)])*...
        [Dbar*[eye(2);jacobian(xo,q)];zeros(2,2)];
    Delta1=eval([q;Deltaqdot*qdot]);
    Delta1=simplify(Delta1);
    
    % Composition of two maps
    Delta=simplify(subs(Delta2,[q;qdot],Delta1));
    Deltafun=matlabFunction(Delta,'Vars',{[q;qdot]});
    save('walking_acrobot_model','D','Dfun','Cfun','Gfun','Deltafun','B');
else
    fprintf('\nLoading robot model...\n')
    load('walking_acrobot_model');
end
%% HERE WRITE YOUR CODE FOR THE VHC DESIGN
% The outcome of this part should be a parameter vector a, whose components
% a_1,ldots,a_k define the polynomial phi_a(theta)=a_1 + a_2 theta + ... +
% a_k theta^(k-1)

fprintf('\nDetermining vhc...\n')

%% ===================== Step 1 (Section 2.4) ============================
% Compute q+, q-, and q_bar from β

beta = 0.316637;  % initial leg aperture (can be modified later)

% Pre-impact configuration q-
q_minus = [(pi - beta)/2 ;
           pi + beta];

% Post-impact configuration q+
q_plus  = [(pi + beta)/2 ;
           pi - beta];

% Scuffing point midpoint configuration q_bar
q_bar   = [pi/2 ;
           pi];

% Difference in stance leg angles (used later for sigma and phi)
q_tilde1 = q_plus(1) - q_minus(1);

%% =============== Compute I and check v1, v2 feasibility ==============
% Requires: q_minus (2x1), q_tilde1 (scalar), and Deltafun available

% Compute I by evaluating Deltafun at basis velocities
% Deltafun expects a 4x1 vector [q; qdot] and returns a 4x1 [q_plus; post_qdot]
% We evaluate at q_minus and qdot = [1;0] and [0;1] to get columns of I
d1 = Deltafun([q_minus; [1;0]]);
d2 = Deltafun([q_minus; [0;1]]);
I = [d1(3:4), d2(3:4)];   % I is 2x2, columns = response to basis qdot

% Extract components for clarity
I11 = I(1,1); I12 = I(1,2);
I21 = I(2,1); I22 = I(2,2);

fprintf('\nMatrix I (dT*J) at q_minus:\n');
disp(I);

% Simple bounds from (11a) and (11b)
v1_upper = 2 * q_tilde1;   % v1 <= 2 q_tilde1
v2_lower = 2 * q_tilde1;   % v2 >  2 q_tilde1

fprintf('Simple bounds: v1 <= %.6f, v2 > %.6f\n', v1_upper, v2_lower);

% Example candidate values (from PDF) -- you can change these
v1_candidate = -0.894373;
v2_candidate = 1.9;

% Check basic inequalities
cond11a = (v1_candidate <= v1_upper);
cond11b = (v2_candidate > v2_lower);

fprintf('Candidate v1=%.6f : satisfies v1 <= 2*q_tilde1 ? %d\n', v1_candidate, cond11a);
fprintf('Candidate v2=%.6f : satisfies v2 >  2*q_tilde1 ? %d\n', v2_candidate, cond11b);

% Check the either/or conditions in (11c)
% Avoid divisions by zero: check I12 and (2*I12 + I22)
tol = 1e-12;
ok11c = false;
if abs(I12) > tol
    thresh1 = (I11 / I12) * q_tilde1;
    if abs(2*I12 + I22) > tol
        thresh2 = ((I21 + 2*I11) / (2*I12 + I22)) * q_tilde1;
    else
        thresh2 = NaN; % indicate can't compute
    end

    % First branch
    branch1 = (v1_candidate > thresh1) && (v1_candidate >= thresh2);
    % Second branch
    branch2 = (v1_candidate < thresh1) && (v1_candidate <= thresh2);

    ok11c = branch1 || branch2;
    fprintf('11c thresholds: thresh1=%.6f, thresh2=%.6f\n', thresh1, thresh2);
    fprintf('11c branch1? %d, branch2? %d --> 11c satisfied? %d\n', branch1, branch2, ok11c);
else
    fprintf('I12 is too small (or zero). Cannot evaluate thresholds for (11c) reliably.\n');
end

% Compute f(v1) per equation (7) if denominator valid
denom = -I11*q_tilde1 + I12*v1_candidate;
if abs(denom) < 1e-12
    error('Denominator for f(v1) is (near) zero; chosen v1 not allowed (pole).');
else
    f_of_v1 = - q_tilde1 * (-I21*q_tilde1 + I22*v1_candidate) / denom;
    fprintf('Computed f(v1) = %.6f\n', f_of_v1);
end

% Final feasibility flag: all conditions
feasible = cond11a && cond11b && ok11c && (abs(denom) >= 1e-12);
fprintf('Final feasibility for candidate (v1,v2) = (%.6f,%.6f): %d\n', v1_candidate, v2_candidate, feasible);

%% ================== Step 2: Solve 6x6 system for a ====================

% Construct matrix for φ and φ' boundary conditions
M = [ ...
    1      0       0^2       0^3       0^4       0^5 ;       % phi(0)   = q2+
    1      1       1^2       1^3       1^4       1^5 ;       % phi(1)   = q2-
    1    0.5   (0.5)^2   (0.5)^3   (0.5)^4   (0.5)^5 ;       % phi(.5) = pi
    0      1   2*0       3*0^2    4*0^3    5*0^4  ;          % phi'(0) = f(v1)
    0      1   2*1       3*1^2    4*1^3    5*1^4  ;          % phi'(1) = v1
    0      1   2*0.5 3*(0.5)^2 4*(0.5)^3 5*(0.5)^4 ];        % phi'(.5)= v2

% Right-hand side vector
b = [ ...
    q_plus(2) ;      % = q2+
    q_minus(2) ;     % = q2-
    pi ;             % midpoint
    f_of_v1 ;        % from equation (7)
    v1_candidate ;
    v2_candidate ];

% Solve for polynomial coefficients
a = M \ b;           % <= this is the result of Step 2

disp('Solved VHC polynomial parameters a =')
disp(a)

%% ================== Step 3: Define phi and sigma functions ==================

% Ensure a exists (computed in Step 2). Force column
a = a(:);

% MATLAB polyval expects coefficients from highest degree to lowest.
% Flip the a vector to match polyval's ordering.
a_poly = flip(a)';   % row vector for polyval / polyder

% Define phi and its derivatives (w.r.t. theta)
phi = @(theta) polyval(a_poly, theta);
phi_prime_coeffs = polyder(a_poly);
phi_doubleprime_coeffs = polyder(phi_prime_coeffs);
phiprime = @(theta) polyval(phi_prime_coeffs, theta);
phipprime = @(theta) polyval(phi_doubleprime_coeffs, theta);

% Define sigma(theta) and its derivatives (w.r.t. theta)
% First component: q_plus(1) - theta * q_tilde1
sigma = @(theta) [ q_plus(1) - theta * q_tilde1 ; phi(theta) ];
sigmaprime = @(theta) [ -q_tilde1 ; phiprime(theta) ];
sigmapprime = @(theta) [ 0 ; phipprime(theta) ];

% Register into data struct (used later by controller / simulations)
data.Kp=Kp;
data.Kd=Kd;
data.D=Dfun;
data.C=Cfun;
data.G=Gfun;
data.B=B;
data.phi = phi;
data.phiprime = phiprime;
data.phipprime = phipprime;
data.sigma = sigma;
data.sigmaprime = sigmaprime;
data.sigmapprime = sigmapprime;
data.q_plus = q_plus;
data.q_tilde1 = q_tilde1;


% (Optional) quick sanity prints
fprintf('VHC polynomial (coeffs a, lowest->highest):\n');
disp(a');
fprintf('VHC polynomial (for polyval, highest->lowest a_poly):\n');
disp(a_poly);

%% ================== Step 4 — VHC Verification (Section 2.4) ==================

fprintf("\n================ VHC Verification ================\n");

% Create grid over θ ∈ [0,1]
theta_grid = linspace(0,1,200);

% Evaluate σ(θ)
sig_eval = arrayfun(@(th) sigma(th), theta_grid, 'UniformOutput', false);
sig_eval = cell2mat(sig_eval); % 2xN matrix of q1,q2

q1_vals = sig_eval(1,:);
q2_vals = sig_eval(2,:);

%% 1. Check interpolation constraints (pass through required points)
fprintf("\nChecking σ(θ) boundary conditions:\n");

fprintf("σ(0)   should equal q+    →   "); disp(sigma(0)');
fprintf("σ(1)   should equal q-    →   "); disp(sigma(1)');
fprintf("σ(0.5) should equal q_bar →   "); disp(sigma(0.5)');

%% 2. Safe set check  (Required by Section 2.4)
% Safe walking region: 0 < q1 < π,   0 < q2 < 2π
in_bounds = (q1_vals > 0 & q1_vals < pi) & (q2_vals > 0 & q2_vals < 2*pi);

if all(in_bounds)
    fprintf("\n✓ VHC stays inside safe set W for all θ.\n");
else
    fprintf("\n✗ WARNING: VHC leaves safe set! Check polynomial or β.\n");
end

%% 3. Regularity check — must satisfy:  B⊥ Dσ'(θ) ≠ 0  ∀θ ∈ [0,1]
% For acrobot B = [0;1] so B⊥ = [1 0], making:
%   B⊥Dσ' = first row of D(q) * σ'(θ)
regularity = zeros(size(theta_grid));

for k = 1:length(theta_grid)
    th = theta_grid(k);
    q = sigma(th);              % q = (q1,q2)
    Dq = data.D(q);             % inertia matrix at q
    s_prime = sigmaprime(th);   % σ'(θ)
    regularity(k) = Dq(1,:) * s_prime;
end

if all(abs(regularity) > 1e-4)
    fprintf("✓ Regularity condition satisfied: B⊥Dσ'(θ) ≠ 0 everywhere.\n");
else
    fprintf("✗ Regularity violation detected — redesign VHC required.\n");
end

%% 4. Virtual mass & Potential (Section 2.4 end → leads into 2.5 later)

M_vals = zeros(size(theta_grid));
V_vals = zeros(size(theta_grid));

for k = 1:length(theta_grid)
    th = theta_grid(k);
    q = sigma(th);
    Dq = data.D(q);
    s_prime = sigmaprime(th);
    s_dprime = sigmapprime(th);

    % Virtual mass M(θ)
    M_vals(k) = s_prime' * (Dq * s_prime);

    % Virtual potential V(θ) approximated as G(q)' * sigma derivative (used for plotting)
    V_vals(k) = data.G(q)' * s_prime;
end

% Plot results
figure; subplot(2,2,1);
plot(theta_grid,q1_vals,'LineWidth',2); hold on
plot(theta_grid,q2_vals,'LineWidth',2);
title('σ(θ) trajectory'); legend('q1(θ)','q2(θ)'); grid on

subplot(2,2,2);
plot(theta_grid,regularity,'LineWidth',2);
title('Regularity B^\perp Dσ''(θ)'); grid on

subplot(2,2,3);
plot(theta_grid,M_vals,'LineWidth',2);
title('Virtual Mass M(θ)'); grid on

subplot(2,2,4);
plot(theta_grid,V_vals,'LineWidth',2);
title('Virtual Potential V(θ)'); grid on

fprintf("\n========== Verification complete ==========\n");

%% NOW YOU CAN SIMULATE THE ROBOT. PLACE YOUR CONTROLLER INSIDE THE FUNCTION acrobot AT THE END OF THIS SCRIPT
ops= odeset('reltol',1e-7,'abstol',1e-7,'Events',@ground_impact);
dt=1/60; % 60 fps; time increment in simulations and animations

fprintf('\n Simulating...\n')
%% DEFINE THE INITIAL CONDITION [q0;qdot0];

% NOTE: this is a placeholder initial condition set on the constraint manifold at theta=0.
% Later you should set qdot0 according to the computed limit-cycle theta_dot_a in the document.
    % ===================== VHC feedback-linearizing controller =====================
    % Output y = q2 - phi(theta) with theta = (q_plus(1) - q1)/q_tilde1
    % Note: q_plus and q_tilde1 are in the main workspace; ensure they are visible
    % through nested functions or store them in 'data' (we put them in data).
    
    % Get parameters and VHC functions from data
    phi = data.phi;
    phiprime = data.phiprime;
    phipprime = data.phipprime;
    sigma = data.sigma;
    sigmaprime = data.sigmaprime;
    sigmapprime = data.sigmapprime;
    Dfun = data.D;
    Cfun = data.C;
    Gfun = data.G;
    B = data.B;
    Kp = data.Kp;
    Kd = data.Kd;
    
    % Make sure q_plus(1) and q_tilde1 are available in data
    if isfield(data,'q_plus')
        q_plus_local = data.q_plus;
    else
        error('data.q_plus not found. Please set data.q_plus earlier.');
    end
    if isfield(data,'q_tilde1')
        q_tilde1_local = data.q_tilde1;
    else
        error('data.q_tilde1 not found. Please set data.q_tilde1 earlier.');
    end
    
    % Compute theta and its derivative
    theta = (q_plus_local(1) - q1) / q_tilde1_local;
    theta_dot = - q1dot / q_tilde1_local;  % d/dt[(q+1 - q1)/qtilde1] = -q1dot/qtilde1
    
    % Output and its derivative
    y = q2 - phi(theta);
    ydot = q2dot - phiprime(theta) * theta_dot;
    
    % Compute dynamics terms
    Dq = Dfun(q);
    Cq = Cfun([q;qdot]);
    Gq = Gfun(q);
    Dinv = inv(Dq);
    
    % Build S = [ phiprime/q_tilde1 , 1 ] to combine qddot into yddot
    S = [ phiprime(theta)/q_tilde1_local , 1 ]; % 1x2
    
    % extra term from phi''(theta)*(theta_dot)^2 (see derivation)
    extra = - phipprime(theta) * (theta_dot)^2;
    
    % Compute L = S * Dinv * (-C*qdot - G) + extra
    L = S * ( Dinv * ( - Cq * qdot - Gq ) ) + extra;
    
    % Compute B2 = S * Dinv * B  (scalar multiplier of tau in yddot)
    B2 = S * ( Dinv * B );
    
    % Safety: if B2 is too small, avoid division by (near) zero
    if abs(B2) < 1e-8
        % fallback: saturate/regularize or set tau to simple PD on q2
        % We choose a safe small-control fallback to avoid NaNs
        warning('B2 is very small in controller (possible singularity). Using PD fallback on q2.');
        % simple PD directly on q2 (not ideal, but prevents blowup)
        tau = - (Kp * sin(y) + Kd * ydot);
    else
        % Desired linearizing input (v) with PD on y, using sin(y) for nonlinearity
        v = - Kp * sin(y) - Kd * ydot;
        % Solve for tau :  yddot = L + B2 * tau  =>  tau = (v - L) / B2
        tau = (v - L) / B2;
    end
    % =============== end of VHC controller ============================
    
    % Compute accelerations after control
    qddot = Dinv * ( - Cq * qdot - Gq + B * tau );
    
    xdot = [qdot; qddot];


T=[];
X=[];
Te=[];
Ie=[];
Xe=[];
% Simulate number_steps steps
disp("------------------------------------------------------------");
disp("Running Closed-Loop Walking Simulation...");
disp("------------------------------------------------------------");

post_impact_state = [q0 ; qdot0];

for step=1:number_steps
    fprintf('\n...step %d...\n',step);
    [t,x,te,xe,ie]=ode45(@(t,x) acrobot(t,x,data),0:dt:10,post_impact_state,ops);
    % Application of the impact map
    impact_state=xe(end,:)';
    post_impact_state=Deltafun(impact_state);
    T{step}=t;
    X{step}=x;
    Ie{step}=ie;
    Te{step}=te;
    Xe{step}=xe;
end

fprintf('\n Setting up animation...\n')
figure(1);

%% Animation of the simulation results
ref=0;time_passed=0;step=1;
Axis=[-1 4 0 2];
Time=text(-1+2,1.8,['time= ','0',' secs,',' step= ',num2str(step)]);
axis(Axis);
stance_leg=line([ref l*cos(q0(1))],[0 l*sin(q0(1))],'color','red','linewidth',2);
swing_leg=line([ref+l*cos(q0(1)) ref+l*cos(q0(1))+l*cos(q0(1)+q0(2))],...
    [l*sin(q0(1)) l*sin(q0(1))+l*sin(q0(1)+q0(2))],'linewidth',2);
fprintf('\n Animation is ready...\n')

animation_slowdown_factor=1; % >1 means slow down
for step=1:length(Ie)
    t=T{step};
    x=X{step};
    xe=Xe{step};
    xe=xe(end,:);
    for k=2:length(t)
        t0=clock;
        drawnow;
        q=x(k,1:2)';
        xdata1=[ref ref+l*cos(q(1))];
        xdata2=[ref+l*cos(q(1)) ref+l*cos(q(1))+l*cos(q(1)+q(2))];
        ydata1=[0 l*sin(q(1))];
        ydata2=[l*sin(q(1)) l*sin(q(1))+l*sin(q(1)+q(2))];
        set(stance_leg,'xdata',xdata1,'ydata',ydata1);
        set(swing_leg,'xdata',xdata2,'ydata',ydata2);
        set(Time,'String',['time= ',num2str(round(time_passed+t(k),1)),' secs,',' step= ',num2str(step)]);
        current_axis=gca;
        if ref>.95*current_axis.XLim(end)
            current_axis.XLim=[.95*ref .95*ref+5];
            Time.Position=[.95*ref+2 1.8 0];
            Axis=axis;
        else
            axis(Axis)
        end
        while etime(clock,t0)<animation_slowdown_factor*(t(k)-t(k-1))
        end
    end
    time_passed=time_passed+t(end);
    ref=ref+l*(cos(xe(1))+cos(xe(1)+xe(2)));
end


%% FUNCTIONS
function xdot=acrobot(t,x,data)
q1=x(1);
q2=x(2);
q1dot=x(3);
q2dot=x(4);
q=x(1:2);
qdot=x(3:4);
Kp=data.Kp;
Kd=data.Kd;
D=data.D;
C=data.C;
G=data.G;
B=data.B;
phi=data.phi;
phiprime=data.phiprime;
phipprime=data.phipprime;
% DEFINE YOUR CONTROLLER HERE
% tau=....
% For now set tau = 0 as placeholder. Replace with VHC controller later.
    % ===================== VHC feedback-linearizing controller =====================
    % Output y = q2 - phi(theta) with theta = (q_plus(1) - q1)/q_tilde1
    % Note: q_plus and q_tilde1 are in the main workspace; ensure they are visible
    % through nested functions or store them in 'data' (we put them in data).
    
    % Get parameters and VHC functions from data
    phi = data.phi;
    phiprime = data.phiprime;
    phipprime = data.phipprime;
    sigma = data.sigma;
    sigmaprime = data.sigmaprime;
    sigmapprime = data.sigmapprime;
    Dfun = data.D;
    Cfun = data.C;
    Gfun = data.G;
    B = data.B;
    Kp = data.Kp;
    Kd = data.Kd;
    
    % Make sure q_plus(1) and q_tilde1 are available in data
    if isfield(data,'q_plus')
        q_plus_local = data.q_plus;
    else
        error('data.q_plus not found. Please set data.q_plus earlier.');
    end
    if isfield(data,'q_tilde1')
        q_tilde1_local = data.q_tilde1;
    else
        error('data.q_tilde1 not found. Please set data.q_tilde1 earlier.');
    end
    
    % Compute theta and its derivative
    theta = (q_plus_local(1) - q1) / q_tilde1_local;
    theta_dot = - q1dot / q_tilde1_local;  % d/dt[(q+1 - q1)/qtilde1] = -q1dot/qtilde1
    
    % Output and its derivative
    y = q2 - phi(theta);
    ydot = q2dot - phiprime(theta) * theta_dot;
    
    % Compute dynamics terms
    Dq = Dfun(q);
    Cq = Cfun([q;qdot]);
    Gq = Gfun(q);
    Dinv = inv(Dq);
    
    % Build S = [ phiprime/q_tilde1 , 1 ] to combine qddot into yddot
    S = [ phiprime(theta)/q_tilde1_local , 1 ]; % 1x2
    
    % extra term from phi''(theta)*(theta_dot)^2 (see derivation)
    extra = - phipprime(theta) * (theta_dot)^2;
    
    % Compute L = S * Dinv * (-C*qdot - G) + extra
    L = S * ( Dinv * ( - Cq * qdot - Gq ) ) + extra;
    
    % Compute B2 = S * Dinv * B  (scalar multiplier of tau in yddot)
    B2 = S * ( Dinv * B );
    
    % Safety: if B2 is too small, avoid division by (near) zero
    if abs(B2) < 1e-8
        % fallback: saturate/regularize or set tau to simple PD on q2
        % We choose a safe small-control fallback to avoid NaNs
        warning('B2 is very small in controller (possible singularity). Using PD fallback on q2.');
        % simple PD directly on q2 (not ideal, but prevents blowup)
        tau = - (Kp * sin(y) + Kd * ydot);
    else
        % Desired linearizing input (v) with PD on y, using sin(y) for nonlinearity
        v = - Kp * sin(y) - Kd * ydot;
        % Solve for tau :  yddot = L + B2 * tau  =>  tau = (v - L) / B2
        tau = (v - L) / B2;
    end
    % =============== end of VHC controller ============================
    
    % Compute accelerations after control
    qddot = Dinv * ( - Cq * qdot - Gq + B * tau );
    
    xdot = [qdot; qddot];

end

function [value,isterminal,direction]=ground_impact(t,x)
q1=x(1);
q2=x(2);
% impact occurs when q2 = -2*q1+2*pi
value=q2+2*q1-2*pi;

% We exclude the scuffing point from the impact conditions
if abs(q1-pi/2)<0.01
    isterminal=0;
else
    isterminal=1;
end

% We distinguish between impact on S^+ or S^- by changing the way in which
% ode45 monitors the zero crossing
if abs(q1)<pi/2
    direction=-1;
else
    direction=1;
end
end
