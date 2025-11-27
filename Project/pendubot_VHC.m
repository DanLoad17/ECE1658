%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ECE1658
%% Lecture X
%% Pendubot with a VHC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc
symbolic_computations=1;
vhc_search=1;
simulation=1;
plot_errors=1;
find_constrained_dyn=0;
animation=1;

% Flag consistency
animation=animation*simulation;
plot_errors=plot_errors*simulation;
% Initial conditions
q10=0;
q1dot0=10;
percent_error_from_manifold=1.01;
vhc_degree=-1;
%% Set 1
% q0=[0;pi-pi/12];
% qdot0=[0;0];
%% Set 2
% q0=[0;pi-pi/6];
% qdot0=[0;0];
%% Set 3
% q0=[pi;-pi/60];
% qdot0=[0;.1];
%% Set 4
% q0=[pi;-pi/60];
% qdot0=[0;0];
ops= odeset('reltol',1e-8,'abstol',1e-8);
%%

if symbolic_computations
    fprintf('\nSymbolic computation of the model...\n')
    % Define symbolic variables
    syms l lc Iz m g real
    syms t q1 q2 q1dot q2dot tau real
    q=[q1;q2];
    qdot=[q1dot;q2dot];
    % Define centres of mass of two links
    rc1=lc*[sin(q1);-cos(q1)];
    rc2=l*[sin(q1);-cos(q1)]+lc*[sin(q1+q2);-cos(q1+q2)];
    
    % Compute time derivatives of centres of mass
    rc1dot=jacobian(rc1,q)*qdot;
    rc2dot=jacobian(rc2,q)*qdot;
    
    % Define the total kinetic energy of the robot
    K=1/2*m*(rc1dot'*rc1dot+rc2dot'*rc2dot)+1/2*Iz*(q1dot^2+(q1dot+q2dot)^2);
    K=simplify(K);
    
    % Extract the square symmetric matrix of the kinetic energy
    D=simplify(hessian(K,qdot));
    
    % Define the potential energy of the robot
    
    P = m*g*(-lc*cos(q1)-l*cos(q1)-lc*cos(q1+q2));
    
    % Input matrix
    
    B=[1;0]; Bperp=[0 1];
    
    % Computation of matrix C(q,qdot)
    C = sym(zeros(2,2));
    for i=1:2
        for j=1:2
            for k=1:2
                C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+diff(D(k,i),q(j))-diff(D(i,j),q(k)))*qdot(i);
            end
        end
    end
    
    % Computation of gradient of the potential
    
    G = jacobian(P,q)';
    
    % Computation of qddot
    
    % qddot = simplify(inv(D)*(-C*qdot-G+B*tau));
    
    % Physical parameters
    l=1;
    lc=0.5;
    m=1;
    Iz=1/12*m*l^2;
    g=9.81;
    save('lec22_pendubot_VHC.mat','D','C','G','l','lc','m','Iz','g','B',...
        'Bperp','q1','q2','q1dot','q2dot','q','qdot')
else
    load('lec22_pendubot_VHC.mat')
end
syms theta delta phi phiprime phipprime real
sigma=[theta;phi];

%% Search for virtual constraints
if vhc_search
    fprintf('\nSearching for a virtual constraint of degree %d...\n',vhc_degree)
    LHS=subs(Bperp*subs(D,q,sigma)*[1;phiprime]-delta);
    Phiprime = simplify(solve(LHS,phiprime));
    Phiprime=subs(subs(Phiprime));
    [delta,fval,flag] =fsolve(@constraint_search,0.1,optimset('Display','off'),vhc_degree);
    close
    if flag==1
        fprintf('\n...found it! delta = %f\n',delta)
    else
        error('problem finding the VHC')
    end
    figure
    hold on;
    vhc_gen=matlabFunction(eval(Phiprime),'Vars',{theta,phi});
    [Theta,Phi]=ode45(vhc_gen,linspace(0,2*pi,100),pi,ops);
    plot(Theta,Phi)
    xlabel('\theta')
    ylabel('\phi(\theta)')
    title('VHC')
    drawnow
    
    phi_spline=spline(Theta,Phi);
    phiprime_spline=fnder(phi_spline);
    phipprime_spline=fnder(phiprime_spline);    
    save('lec22_splines.mat','phi_spline','phiprime_spline','phipprime_spline');
else
    load('lec22_splines.mat');
end
figure;hold on;
axis([-2 2 -2.5 2.5]);axis equal;
s1=linspace(0,2*pi,15);
s2=ppval(phi_spline,s1);
for j=1:length(s1)
    Q=[s1(j);s2(j)];
    link1=line([0 l*sin(Q(1))],[0 -l*cos(Q(1))],'color','red','linewidth',2);
    link2=line([l*sin(Q(1)) l*sin(Q(1))+l*sin(Q(1)+Q(2))],...
        [-l*cos(Q(1)) -l*cos(Q(1))-l*cos(Q(1)+Q(2))],'linewidth',1);
end
s1=linspace(0,2*pi,1000);
s2=ppval(phi_spline,s1);
plot(l*sin(s1)+l*sin(s1+s2),-l*cos(s1)-l*cos(s1+s2),'g','linewidth',2)
title('acrobot configurations on the constraint set')
drawnow;
%% Numerical simulation
if simulation
    fprintf('\nSimulating the closed-loop system\n')
    D=subs(D);
    Dinv=inv(D);
    C=subs(C);
    G=subs(G);
    Dfun = matlabFunction(D,'Vars',{q});
    Dinvfun = matlabFunction(Dinv,'Vars',{[q1;q2]});
    Cfun= matlabFunction(C,'Vars',{[q1;q2;q1dot;q2dot]});
    Gfun = matlabFunction(G,'Vars',{[q1;q2]});    
    K=[1 2];
    data.D=Dfun;
    data.Dinv=Dinvfun;
    data.C=Cfun;
    data.G=Gfun;
    data.phi=phi_spline;
    data.phiprime=phiprime_spline;
    data.phipprime=phipprime_spline;
    data.Kp=K(1);
    data.Kd=K(2);
    dt=1/60; % 60 fps; time increment in simulations and animations
    
    % Initial conditions
    q0=[q10;ppval(phi_spline,q10)*percent_error_from_manifold];
    qdot0=[q1dot0;ppval(phiprime_spline,q10)*q1dot0];
    %
    
    [t,x]=ode45(@pendubot,0:dt:15,[q0;qdot0],ops,data);
    q1=x(:,1);
    q2=x(:,2);
    q1dot=x(:,3);
    q2dot=x(:,4);
end

%% Plotting manifold errors
if plot_errors
    fprintf('\nPlotting constraint manifold errors...\n')
    e=q2-ppval(phi_spline,wrapTo2Pi(q1));
    edot=q2dot-ppval(phiprime_spline,wrapTo2Pi(q1)).*q1dot;
    
    figure
    hold on
    plot(t,wrapToPi(e));
    plot(t,edot);
    xlabel('t')
    title('Manifold errors, e and edot')
    drawnow;
end
%% Constrained dynamics
if find_constrained_dyn==1
    fprintf('\nDerivation of constrained dynamics...\n')
    denominator=simplify(Bperp*D*[1;phiprime]);
    Psi1=-Bperp*G/denominator;
    Psi2=-Bperp*(D*[0;phipprime]+subs(C,qdot,[1;phiprime])*[1;phiprime])/denominator;    
    Psi1=simplify(subs(Psi1,q,sigma));
    Psi2=simplify(subs(Psi2,q,sigma));
    Psi1=subs(Psi1);
    Psi2=subs(Psi2);
    Psi1fun=@(theta,phi,phiprime) eval(Psi1);
    Psi2fun=@(phi,phiprime,phipprime) eval(Psi2);
    data2.phi=phi_spline;
    data2.phiprime=phiprime_spline;
    data2.phipprime=phipprime_spline;
    data2.Psi1=Psi1fun;
    data2.Psi2=Psi2fun;
    fprintf('\n...finding virtual mass and virtual potential...\n')
    ops2=odeset('RelTol',1e-4,'AbsTol',1e-4);
    [Theta,X]=ode45(@mass_potential,linspace(0,2*pi,1000),[1;0],ops2,data2);
    M=spline(Theta,X(:,1));
    V=spline(Theta,X(:,2));
    figure
    plot(Theta,X(:,1));
    hold on
    plot(Theta,X(:,2));
    legend('M','V')
    drawnow;
    save('lec22_Lagrangian_structure','M','V');
else
    load('lec22_Lagrangian_structure');
end
figure;hold on;
fprintf('\n...plotting phase portrait of constrained dynamics...\n')
[Theta,Thetadot]=meshgrid(linspace(-pi,pi,100),linspace(-20,20,100));
E=1/2*ppval(M,wrapTo2Pi(Theta)).*Thetadot.^2+ppval(V,wrapTo2Pi(Theta));
VV=ppval(V,linspace(0,2*pi,1000));
Emax=max(VV);
Emin=min(VV);
contour(Theta,Thetadot,E,linspace(Emin,Emax,10));
contour(Theta,Thetadot,E,linspace(sqrt(Emax),sqrt(Emax)+10,10).^2);
contour(Theta,Thetadot,E,[Emax Emax],'color','k','linewidth',2);
xlabel('\theta')
ylabel('thetadot')
axis tight
title('Orbits of the constrained dynamics')
if simulation
    hold on
    plot(wrapToPi(q1),q1dot,'r')
    plot(wrapToPi(q1(1)),q1dot(1),'ro')
end
drawnow;
%% Animation of the simulation results
if animation
    fprintf('\nSetting up animation...\n')
    figure
    Axis=[-2 2 -2 2];
    axis equal;
    Time=text(1,1.8,['time= ',num2str(t(1)),' secs']);
    axis(Axis);
    q0(1)=q0(1);
    link1=line([0 l*sin(q0(1))],[0 -l*cos(q0(1))],'color','red','linewidth',2);
    link2=line([l*sin(q0(1)) l*sin(q0(1))+l*sin(q0(1)+q0(2))],...
        [-l*cos(q0(1)) -l*cos(q0(1))-l*cos(q0(1)+q0(2))],'linewidth',2);
    fprintf('\nAnimation is ready...\n')
    pause
    tic
    animation_slowdown_factor=1; % >1 means slow down
    for k=2:length(t)
        t0=clock;
        drawnow;
        q=x(k,1:2)';
        q(1)=q(1);
        xdata1=[0 l*sin(q(1))];
        xdata2=[l*sin(q(1)) l*sin(q(1))+l*sin(q(1)+q(2))];
        ydata1=[0 -l*cos(q(1))];
        ydata2=[-l*cos(q(1)) -l*cos(q(1))-l*cos(q(1)+q(2))];
        set(link1,'xdata',xdata1,'ydata',ydata1);
        set(link2,'xdata',xdata2,'ydata',ydata2);
        set(Time,'String',['time= ',num2str(round(t(k),1)),' secs']);
        axis(Axis);
        while etime(clock,t0)<animation_slowdown_factor*(t(k)-t(k-1))
        end
    end
    toc
end


function endpoint=constraint_search(delta,k)
syms theta phi real
vhc_generator=matlabFunction(3*delta-3/2*cos(phi)-1,'Vars',{theta,phi});
ops= odeset('reltol',1e-8,'abstol',1e-8);
[Theta,Phi]=ode45(vhc_generator,[0,2*pi],pi,ops);
% plot(Theta,Phi,'g')
% title(['delta= ',num2str(delta)])
% drawnow
% pause(.3)

% degree of phi
% k=0;
% k=-1;
% k=1;
% k=-2;
% k=2;

endpoint=Phi(end)-pi-2*k*pi;
end

function xdot=pendubot(t,x,data);

q=x(1:2);
q1=wrapTo2Pi(x(1));q1dot=x(3);
qdot=x(3:4);
phi=ppval(data.phi,q1);
phiprime=ppval(data.phiprime,q1);
phipprime=ppval(data.phipprime,q1);
D=data.D(q);
Dinv=data.Dinv(q);
C=data.C(x);
G=data.G(q);
Kp=data.Kp;
Kd=data.Kd;
B=[1;0];
e=(q(2)-phi);
edot=qdot(2)-phiprime*qdot(1);

% The implicit form of the constraint is
% h(q) = q2 - phi(q1)
tau=1/([-phiprime 1]*Dinv*B)*([-phiprime 1]*Dinv*(C*qdot+G)+phipprime*q1dot^2-Kp*sin(e)-Kd*edot);
qddot=Dinv*(-C*qdot-G+B*tau);
xdot=[qdot;qddot];
end

function xdot=mass_potential(theta,x,data2)
M=x(1);
V=x(2);
phi=ppval(data2.phi,theta);
phiprime=ppval(data2.phiprime,theta);
phipprime=ppval(data2.phipprime,theta);
Mdot = -2*M*data2.Psi2(phi,phiprime,phipprime);
Vdot = -data2.Psi1(theta,phi,phiprime)*M;
xdot=[Mdot;Vdot];
end
