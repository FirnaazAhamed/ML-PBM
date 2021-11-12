%% ML-PBM v11.01

% Developed and written by Firnaaz Ahamed (2019)
% Contact: firnaaz.ahamed@monash.edu; firnaaz.ahamed@gmail.com

%% Simulation Code

close all; clear all; clc;

% Setting input parameters

% Cellulose and enzymes

mss = 10;           % Initial concentration of cellulose, g/L

L_CBH = 11.2;       % CBH concentration, mg enzyme/g cellulose
L_EG = 0;           % EG concentration, mg enzyme/g cellulose
L_BG = 0;           % BG concentration, mg enzyme/g cellulose

MW_CBHi = 70000;    % CBH molecular weight, g/mol 
MW_EGi = 50000;     % EG molecular weight, g/mol
MW_BGi = 135000;    % BG molecular weight, g/mol

% Mesh setting

p = 20;             % Number of discrete pivots 

% Time span 

t_end = 1.8e5;      % Time span of simulation, seconds

% Initial distribution - Base fitted cellulose distribution.

load('Raw_Avicel_calibrated_initialDist_Engel_2012_final','MnL','MwL',...
    'MnH','MwH','N','r_massL');

% List of model paramaters

% CBH-related parameters

% 1  kp_h_CBH - Rate contant of CBH hydrolysis, 1/s
% 2  kp_f_CBH - Rate constant of CBH complexation, DP.L/mol.s
% 3  kp_e_CBH - Rate constant of CBH decomplexation, 1/s
% 4  tau_CBH - CBH enzyme footprint/area coverage, mol/m2
% 5  k_ads_CBH - Rate constant of CBH adsorption, L/mol.s
% 6  k_des_CBH - Rate constant of CBH desoprtion, 1/s
% 7  k_If_CBH(1) - Forward glucose inhibition rate constant, L/mol.s
% 8  k_If_CBH(2) - Forward cellobiose inhibition rate constant, L/mol.s
% 9  k_Ir_CBH(1) - Reverse glucose inhibition rate constant, 1/s
% 10 k_Ir_CBH(2) - Reverse cellobiose inhibition rate constant, 1/s

% Substrate-related parameters

% 11 pnt_zone - number of layers in penetration zone
% 12 r_massL - mass ratio of cellulose in penetration zone

% BG-related parameters

% 13 k_h_BG - Rate constant of BG hydrolysis, L/mol.s
% 14 k_If_BG - Forward glucose inhibition rate constant, L/mol.s
% 15 k_Ir_BG - Reverse glucose inhibition rate constant, 1/s

% EG-related parameters

% 16 kp_h_EG - Rate constant of EG hydrolysis (insoluble cellulose), DP/s 
% 17 kp_hs_EG - Rate constant of EG hydrolysis (soluble cellulose), L/mol.s
% 18 kp_f_EG - Rate constant of EG complexation, DP.L/mol.s
% 19 kp_e_EG - Rate constant of EG decomplexation, 1/s
% 20 tau_EG - EG enzyme footprint/area coverage, mol/m2
% 21 k_ads_EG - Rate constant of EG adsorption, L/mol.s
% 22 k_des_EG - Rate constant of EG desorption, 1/s
% 23 k_If_EG(1) - Forward glucose inhibition rate constant, L/mol.s
% 24 k_If_EG(2) - Forward cellobiose inhibition rate constant, L/mol.s
% 25 k_Ir_EG(1) - Reverse glucose inhibition rate constant, 1/s
% 26 k_Ir_EG(2) - Reverse cellobiose inhibition rate constant, 1/s

% Input the values for model parameters below:

prm(1) = 4.491;
prm(2) = 4.32e4;
prm(3) = 10;
prm(4) = 2.22e-8;
prm(5) = 7.3722e5;
prm(6) = 0.01;
prm(7) = 0;
prm(8) = 0;
prm(9) = 0;
prm(10) = 0;

prm(11) = 1;
prm(12) = 0.028;

prm(13) = 3550;
prm(14) = 0;
prm(15) = 0;

prm(16) = 21;
prm(17) = 3;
prm(18) = 4e5;
prm(19) = 10;
prm(20) = 8.04e-9;
prm(21) = 7.11e5;
prm(22) = 0.01;
prm(23) = 0;
prm(24) = 0;
prm(25) = 0;
prm(26) = 0;

% Run simulation

[x_piv,R_in,layers,C_T,C_S,C_INT,MnT,MwT,pdiT,C,CXS,CNS,CNXS,CFS,CS,CI,...
    CT,CC,E_S_EG,E_S_CBH,E_F_EG,E_F_CBH,E_F_BG,R,conv_rem,t] = MLPBM...
    (N,p,prm,mss,MnL,MwL,MnH,MwH,L_CBH,L_EG,L_BG,MW_CBHi,MW_EGi,MW_BGi,...
    t_end);

% List of Outputs

% x_piv     % DP mesh
% R_in      % Initial particle radius
% layers    % Total no. of layers in cellulose particles
% C_T       % Total initial cellulose distribution
% C_S       % Surface initial cellulose distribution, mol/L
% C_INT     % Internal initial cellulose distribution, mol/L
% MnT       % Overall initial number-average DP
% MwT       % Overall initial weight-average DP
% pdiT      % Overall initial polydispersity index

% t         % Time, s

% Time,DP-dependent variables

% C         % Soluble products, mol/L
% CXS       % CBH-bound complex, mol/L
% CNS       % EG-bound complex, mol/L
% CNXS      % CBH-EG-bound complex, mol/L
% CFS       % Free un-bound surface polymers, mol/L
% CS        % Total surface polymers, mol/L
% CI        % Internal polymers, mol/L 
% CT        % Total polymers, mol/L
% CC        % Insoluble polymers, mol/L

% Time-dependent variables

% E_S_EG    % Surface-adsorbed EG, mol/L
% E_S_CBH   % Surface-adsorbed CBH, mol/L
% E_F_EG    % Free EG, mol/L
% E_F_CBH   % Free CBH, mol/L
% E_F_BG    % Free BG, mol/L

% R         % Transient of particle radius, m

% conv_rem  % Overall conversion