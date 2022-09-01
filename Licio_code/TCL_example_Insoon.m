clear all 
close all
clc


%% Model parameters 

R = 2; % Thermal resistence in Celsius/kW
C = 2; % Thermal capacity ini kW/Celsius
P = 14; % Range of energy transfer to or from the thermal mass in kW
eta = 0.7; % Control efficiency
h = 5/60; % Discretization time in hours

alpha = exp(-h/(C*R));

T = 90/60; % simulation time in hours

X_scale = 0:h:T;
X_scale = 60*X_scale;

N_steps = T/h;
N_trajectories = 1000;

Al = 19; % Lower bound on the safe set
Ah = 22; % Upper bound on the safe set

mu = 0; sigma = 0.25^2; % Estimate on the mean and variance
W = 50*[-0.5*sqrt(sigma/12),0.5*sqrt(sigma/12)];  % Support of the distribution

%% Dynamical simulation

% Simulation parameters
x0 = 19.5; % Initial temperature
theta = 32; % Environment temperature
p = 0.05; % success probability

param.alpha = alpha;
param.theta = theta;
param.eta = eta;
param.R = R;
param.P = P;
param.p = p;
param.SafeSet = [Al,Ah];

pd = makedist('Normal','mu',mu,'sigma',sigma);
pd_truncated = truncate(pd,W(1),W(2));

param.noise = pd_truncated;

x = zeros(N_trajectories,N_steps + 1); 
input = zeros(N_trajectories,N_steps);
noise = zeros(N_trajectories,N_steps);

for ell = 1:N_trajectories
     temp = generateTrajec(@f,x0,N_steps,param);
     x(ell,:) = temp.x;
     input(ell,:) = temp.u;
     noise(ell,:) = temp.noise;
end


plot(X_scale,x) % Plot the generated trajectories

%% Value function computation

part_size = 0.2;
Partition = generatePartition(Ah+2,Al-2,part_size);
SafeSet_index = intersect(find(Al <= Partition),find(Partition <= Ah));
partLength = length(Partition); % Size of the partition;
L = 90/5; % Number of value iteration steps;

valueFunc = zeros(partLength,L);
valueFunc(SafeSet_index,L) = 1;

% Parameters defining the ambiguity set
param.Partition = Partition;
param.meanCenter = 0;
param.Sigma = sigma;
param.meanError = .2;
param.gamma = 1;

m = 50; % Number of samples to update the value function

final = L-2;

for i =L-1:-1:final
    S = ((W(2) - W(1))*rand(m,1) - W(2));
    S = S';
    temp_ValueFunc = valueFunc(:,i+1);
    parfor j = 1:partLength
        currentState = Partition(j);
        vf0 = iterateValueFunction(temp_ValueFunc,@f,currentState,0,S,param);
        vf1 = iterateValueFunction(temp_ValueFunc,@f,currentState,1,S,param);
        valueFunc(j,i) = max(vf0,vf1);
        fprintf('Iteration missing: %d \n',((i-final)*partLength+ (partLength -j)));
    end
    fprintf('Values function after %d iterations: \n',L-i);
    fprintf('%.2g \n',valueFunc(:,i));
end

plot(Partition,valueFunc)

%save 01_06_2022_v2

%% Auxiliary functions

% Function used to generate the trajectory of the TLC
function out = generateTrajec(f,x0,L,param)
    out.x = zeros(L + 1,1);
    out.u = zeros(L,1);
    out.noise = zeros(L,1);
    
    out.x(1) = x0;
    
    for i =2:L+1
        u = binornd(1,param.p);
        omega = random(param.noise);
        previous = out.x(i-1);
        next = f(previous,u,omega,param);
        out.u(i-1) = u;
        out.noise(i-1) = omega;
        out.x(i) = next;
    end
end

%%

% Vector field of the paper

function out = f(x,u,omega,param)
    out = param.alpha*x + (1-param.alpha)*(param.theta - param.eta*param.R*param.P*u) + omega; % Vector field
end

%%

% Defining the partition

function out = generatePartition(arg1,arg2,h)

    upper = max(arg1,arg2);
    low = min(arg1,arg2);
    
    out = low:h:upper;    
    out = out';
end

% Computing the element of the partition

function out = computeElementPartition(func,currentState,input,noise,Partition,param)

out.nextState = [];
out.elementPartition = [];

if ~isempty(func)
    nextState = func(currentState,input,noise,param);
    out.elementPartition = find(nextState > Partition,1,'last');
    out.nextState = nextState;
else
    out.elementPartition = find(currentState > Partition,1,'last');
end

end

function out = iterateValueFunction(ValueFunc,func,currentState,input,Samples,param)

if currentState < param.SafeSet(1) || currentState > param.SafeSet(2)
    out = 0;
else
    N = size(Samples,2);
    Partition = param.Partition;
    c = zeros(N,1);
    for i = 1:N
        transition_next = computeElementPartition(func,currentState,input,Samples(i),Partition,param);
        c(i) = ValueFunc(transition_next.elementPartition);
    end
    mu = param.meanCenter;
    Sigma = param.Sigma;
    
    ambiguityOpt = MomentBasedAmbiguity(c,mu,Sigma,Samples);
    ambiguityOpt.b = param.meanError;
    ambiguityOpt.gamma = param.gamma;
    ambiguityOpt = ambiguityOpt.MomentBasedLinProg;
    
    out = ambiguityOpt.OptRes.opt_value;
    
end

end





