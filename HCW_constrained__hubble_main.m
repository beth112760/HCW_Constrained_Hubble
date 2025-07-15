format long
clear all; close all; clc;
%% Code to solve the constrained HCW Equations (NFZ)
% Written by Bethany Hintz - 15 July 2025
% This code solves the constrained HCW for rendezvous and docking
% with the Hubble Space Telescope

% Initial Values
epsilon = 10;
N = 50;% Number of nodes
M = 3; % Number of elements
n = 1.1313658 * 10^(-3); %orbital frequency
tf = pi/n; % Final time

% Initial Parameteres
Neq = 12;
mu    = 398600.4;
a     = 6778.14+515;%;    % radius of the orbit
n_var = sqrt(mu/(a^3)); 

% Initial conditions for x_bar
% x1 x2 y1 y2 z1 z2
x0 = [50,0.025,50,0.025,50,0.025];
%Final conditions for x_bar
% x1 x2 y1 y2 z1 z2
xf = [0,0,0,0,0,0];


BCtype = 'fixed';%'fixed','free','P0-Pf','P0-Vf','V0-Pf'
BC = [x0,xf]';

%No fly zone parameters
xcenter = 33.8463703737795;
ycenter = -12.1462049693612;
zcenter = 34.0072223233482;
radii = 4.6716;

initial_guess = [ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1);ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1);ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1); ones(1*N*M,1)]; 

%Time guess
t1=1187.22410564302;
t2 = 1376.67230936008;
%% Approximation matrix
D = zeros(N,N,M);
phi = zeros(N,N,M);
phid = zeros(N,N,M);

for k = 1
    t0e            = (k-1)*tf/M; % Initial time of each segment (element)
    tfe            = t1;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N); 
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end

for k = 2
    t0e            = t1; % Initial time of each segment (element)
    tfe            = t2;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N); 
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end

for k = M
    t0e            = t2; % Initial time of each segment (element)
    tfe            = k*tf/M;     % Final time of each segment (element)
    te(:,:,k)      = linspace(t0e,tfe,N); 
    %     tinterp(:,:,k) = linspace(t0e,tfe,5*n);
    for i = 1 : N
        phi(i,:,k)        = rbf0(epsilon,te(:,i,k),te(:,:,k));
        phid(i,:,k)       = rbf1(epsilon,te(:,i,k),te(:,:,k));
        %         phi_interp(i,:,k) = rbf0(c,te(:,i,k),tinterp(:,:,k));
    end
        D(:,:,k) = phid(:,:,k)/(phi(:,:,k));
end
 
%% CRBF Solution using fsolve
ip=0;
Jlocal_pattern = ones(Neq*N);
J_pattern = kron(eye(M),Jlocal_pattern);

J_local = zeros(Neq*N,Neq*N,M);
J    = zeros(Neq*N*M, Neq*N*M); % Local matrices are to be concatenated sequentially

ns = 0;
Elb   = zeros(1,M);
Erb   = Elb;
SLb   = zeros(Neq,1,M);
%% 

SRb   = SLb;
R = zeros(N,M);
sizeJ_local = size(J_local);
LJ_local    = sizeJ_local(1);

f = @(X) HCWNAE_local_constrained_WORKS(X,BC,n_var,D,N,Neq,ns,M,ip,J,J_local,Elb,Erb,SLb,SRb,R,LJ_local,BCtype,xcenter,ycenter,zcenter);
options = optimoptions(@fsolve,'Display','iter','Algorithm','levenberg-marquardt','SpecifyObjectiveGradient',false,'JacobPattern',J_pattern,'StepTolerance',1e-20,'FunctionTolerance',1e-20,'UseParallel',true,'FiniteDifferenceType','central');
[xx,fval1,exitflag1,output1] = fsolve(f,initial_guess,options);

x = reshape(xx,[Neq*N,M]);
X1 = x(1:N,:);
X2 = x(N+1:2*N,:);
Y1 = x(2*N+1:3*N,:);
Y2 = x(3*N+1:4*N,:);
Z1 = x(4*N+1:5*N,:);
Z2 = x(5*N+1:6*N,:);
L1 = x(6*N+1:7*N,:);
L2 = x(7*N+1:8*N,:);
L3 = x(8*N+1:9*N,:);
L4 = x(9*N+1:10*N,:);
L5 = x(10*N+1:11*N,:);
L6 = x(11*N+1:12*N,:);

%% Clean Data %%
T = te(:);
[~, ind] = unique(T);
duplicate_ind = flip(setdiff(1:size(T), ind));
duplicate_value = T(duplicate_ind);
for i = 1 : length(duplicate_ind)
    T(duplicate_ind(i))  = [];
    X1(duplicate_ind(i)) = [];
    X2(duplicate_ind(i)) = [];
    Y1(duplicate_ind(i)) = [];
    Y2(duplicate_ind(i)) = [];
    Z1(duplicate_ind(i)) = [];
    Z2(duplicate_ind(i)) = [];
    L1(duplicate_ind(i)) = [];
    L2(duplicate_ind(i)) = [];

end
x1 = X1(:);
x2 = X2(:);
y1 = Y1(:);
y2 = Y2(:);
z1 = Z1(:);
z2 = Z2(:);

%% Plot Results %%
figure(1)
plot(T,x1)
hold on
plot(T,y1)
plot(T,z1)
xlabel('Time')
ylabel('Position')
title('Unconstrained HCW Equations')
legend('x','y','z')

figure(2)
plot(T,x2)
hold on
plot(T,y2)
plot(T,z2)
xlabel('Time')
ylabel('Velocity')
title('Unconstrained HCW Equations')
legend('x','y','z')

figure(3)
plot3(x1,y1,z1,'*b','LineWidth',5)
xlabel('X Position')
ylabel('Y Position')
zlabel('Z position')
title('HCW Equations')

hold on
[X,Y,Z] = sphere(14);
xplot = X*radii;
yplot = Y*radii;
zplot = Z*radii;
surf(xplot+xcenter,yplot+ycenter,zplot+zcenter)
legend('Trajectory','NFZ')
plot3(x1(1,1),y1(1,1),z1(1,1),'.k','MarkerSize',40)
plot3(x1(end,1),y1(end,1),z1(end,1),'.g','MarkerSize',40)

% Note: load the unconstrained data to plot against constrained data.
% Comment out the k = 2 section of the function