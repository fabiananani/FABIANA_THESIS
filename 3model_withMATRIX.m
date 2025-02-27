clear all
close all
clc
tic

%% Set up Gurobi
addpath('C:\gurobi1103\win64\examples\matlab')
gurobi_setup

%% Input data
path_main = 'C:\Users\nanf\Desktop\H2-districts';
path_input = fullfile(path_main, 'Input');

% Read input files
Input1 = readtable(fullfile(path_input, 'Inputs_BoundaryLoads_2022.xlsx'));
Input2 = readtable(fullfile(path_input, 'nestel_demand.csv'));
Input3 = readtable(fullfile(path_input, 'PriceImpExp_ele.xlsx'));

%% Pre-processing
nHours = 8760;
irradiance = Input1.G;            
P_load = Input2.P_el_nest_h;      
c_gridimp = Input3.Priceimpgreen; 
c_gridexp = Input3.PriceexpPV;    

% Constants
bigM = 600;    
MaxSimTime = 1000;  
deltat = 3600;  
HHV = 39.4 * 3.6 * 10^3;
LHV = 33.3 * 3.6 * 10^3;
eff_PV = 0.17;
eff_FC = 0.51 * HHV / LHV;
eff_ch = 0.95;  
eff_disch = 0.95;

% Fuel Cell degradation parameters
nH_lifetime = 40000;     
Delta_eff_FC = 0.1 * eff_FC;
Const1 = Delta_eff_FC / nH_lifetime;
Const2 = Const1 / HHV / eff_FC^2;

% Hydrogen Cost
c_h2 = 4;  % CHF/kg

% Decision variable limits
S_FC_min = 0;
S_FC_max = 60;
Area_PV_min = 0;
Area_PV_max = 200;
C_b_max = 96 * 3.6 * 10^3;
C_b_min = C_b_max / 3;
P_b_max = C_b_max / 3600;

%% Define decision variables for Gurobi (MIQP)
numVars = 3 + 8 * nHours;  

lb = zeros(numVars, 1);  
ub = inf(numVars, 1);    

vtype = repmat('C', numVars, 1);
vtype(nHours*3+1 : nHours*3+nHours) = 'B';  % Binary for FC operation

% Objective function coefficients
f = zeros(numVars, 1);

% Cost function (Installation)
UP_PV = 1300; 
UP_FC = 950; 
UP_b = 1000;
ann_PV = 0.04 / (1 - (1 + 0.04)^(-25));
ann_FC = 0.04 / (1 - (1 + 0.04)^(-7.5));
ann_b = 0.04 / (1 - (1 + 0.04)^(-12));

f(1) = UP_PV * ann_PV / 1000;
f(2) = UP_FC * ann_FC / 1000;
f(3) = UP_b * ann_b / 1000;

% Power import/export costs
for t = 1:nHours
    f(3 + t) = c_gridimp(t) / 1000;
    f(3 + nHours + t) = -c_gridexp(t) / 1000;
    f(3 + 7 * nHours + t) = c_h2 / 1000; % Hydrogen cost
end

%% Constraints (Sparse Matrices)
Aeq = sparse(nHours, numVars);
beq = P_load;

for t = 1:nHours
    Aeq(t, 3 + t) = 1;    
    Aeq(t, 3 + nHours + t) = -1;  
    Aeq(t, 1) = -irradiance(t) * eff_PV / 1000;  
    Aeq(t, 2) = -1;   
end

% Inequality Constraints (Use Sparse Matrices)
A = sparse([], [], [], 8 * nHours, numVars);
b = zeros(8 * nHours, 1);

for t = 1:nHours
    % Battery charge/discharge limits
    A(t, 3 + 2 * nHours + t) = 1;
    A(t, 3 + 3 * nHours + t) = -1;
    b(t) = P_b_max;

    % FC power constraints
    A(2 * nHours + t, 2) = 1;
    A(2 * nHours + t, 3 + nHours + t) = -1;
    A(2 * nHours + t, 3 + 5 * nHours + t) = -bigM;
    b(2 * nHours + t) = 0;
    
    % w_hours constraints
    if t == 1
        A(4 * nHours + t, 3 + 6 * nHours + t) = 1;  
        A(4 * nHours + t, 3 + 5 * nHours + t) = -1; 
    else
        A(4 * nHours + t, 3 + 6 * nHours + t) = 1;  
        A(4 * nHours + t, 3 + 6 * nHours + t - 1) = -1;  
        A(4 * nHours + t, 3 + 5 * nHours + t) = -1; 
    end
    b(4 * nHours + t) = 0;
end

%% Quadratic Constraints for m_flow_H2
Q = sparse(numVars, numVars);
for t = 1:nHours
    Q(2, 3 + 6 * nHours + t) = Const2 * 3600;  
end

%% **Fixing the `model.sense` Issue**
numIneqConstraints = size(A,1);    % Number of inequality constraints
numEqConstraints = size(Aeq,1);    % Number of equality constraints

model.sense = [repmat('<', numIneqConstraints, 1); repmat('=', numEqConstraints, 1)]; 
% **Ensures `model.sense` matches total constraint count**

%% Solve MIQP problem using Gurobi

model.A = sparse([A; Aeq]);  % Combine inequality and equality constraints
model.Q = Q;  % Quadratic term (bilinear constraint)
model.rhs = [b; beq];  % Ensure `rhs` matches the number of constraints
model.vtype = vtype;
model.modelsense = 'min';
model.obj = f;
model.lb = lb;
model.ub = ub;

params.TimeLimit = MaxSimTime;
params.NonConvex = 2;  % Allow non-convex quadratic constraints

result = gurobi(model, params);

