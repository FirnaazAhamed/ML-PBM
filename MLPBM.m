function [x_piv,R_in,layers,C_T,C_S,C_INT,MnT,MwT,pdiT,C,CXS,CNS,CNXS,...
    CFS,CS,CI,CT,CC,E_S_EG,E_S_CBH,E_F_EG,E_F_CBH,E_F_BG,R,conv_rem,t]...
    = MLPBM(N,p,prm,mss,MnL,MwL,MnH,MwH,L_CBH,L_EG,L_BG,MW_CBHi,MW_EGi,...
    MW_BGi,t_end)

[x_piv,R_in,layers,C_T,C_S,C_INT,MnT,MwT,pdiT,q,vm,R0,ro,L,n,exp_ratio,...
    E_T_CBH,E_T_EG,E_T_BG] = MLPBM_setting(p,N,prm,mss,MnL,MwL,MnH,MwH,...
    L_CBH,L_EG,L_BG,MW_CBHi,MW_EGi,MW_BGi);

[t,y] = ODE_setting(p,q,C_S,C_INT,prm,x_piv,vm,t_end,E_T_EG,E_T_CBH,...
    E_T_BG,R0,ro,L,n,exp_ratio);

[C,CXS,CNS,CNXS,CFS,CS,CI,CT,CC,E_S_EG,E_S_CBH,E_F_EG,E_F_CBH,E_F_BG,R,...
    conv_rem] = post_processing(t,y,p,q,x_piv,ro,L,n,mss,E_T_EG,E_T_CBH,...
    E_T_BG);

end

%%-------------------------------------------------------------------------

function [x_piv,R_in,layers,C_T,C_S,C_INT,MnT,MwT,pdiT,q,vm,R0,ro,L,n,...
    exp_ratio,E_T_CBH,E_T_EG,E_T_BG] = MLPBM_setting(p,N,prm,mss,MnL,MwL,...
    MnH,MwH,L_CBH,L_EG,L_BG,MW_CBHi,MW_EGi,MW_BGi)

% General settings

vm = 1;   % Size of monomer

qmax = 1+(log(N/(p+1))/log(1+(vm/(p+1))));  % Maximum pivots in cont region
q = floor(qmax);                            % No. of pivots in cont region
ratio = (N/(p+1))^(1/(q-1));                % Ratio of geometric progr

% Pivots for discrete-continuous mesh

x_piv = zeros(p+q,1);
x_piv(1:p) = 1:p;

for i = p+1:p+q
    x_piv(i) = (p+1)*ratio^(i-(p+1));       % Geometric mesh
end

% Boundary points for discrete-continuous mesh

x_bound = zeros(p+q+1,1);
x_bound(1) = 0.5;

for i=2:p+q
    x_bound(i) = (x_piv(i)+x_piv(i-1))/2;
end

x_bound(p+q+1) = x_piv(end);

% Initial distribution

c_in = @(x,alpha,beta,pin) pin*gampdf(x,alpha,beta); % Gamma distribution

% Physical properties

layers = round(2/(((mss*prm(12))/prm(11))/mss));     % Total no. of layers
ro = 1500;                   % Cellulose density, g/L
L = (N/2)*10.38e-10;         % Length of microfibril, m
R0 = 1e-9;                   % Width of single layer, m
R_in = layers*R0;            % Total radius of particle, m
n = mss/(ro*pi*(R_in^2)*L);  % No. of particles per unit vol, 1/m3
p_zone = prm(11);            % No. of layers in penetration zone
i_zone = (R_in/R0)-p_zone;   % No. of layers in internal zone

% Radiuses of discrete layers

R_set = sort(0:R0:R_in,'descend')';

% Mass of polymers in each layer

mass_l = zeros(length(R_set)-1,1); 

for i = 2:length(R_set)
    mass_l(i-1) = (n*ro*pi*(R_set(i-1)^2)*L)-(n*ro*pi*(R_set(i)^2)*L);
end

% Initial cellulose distribution in each layer

pin_l = zeros(length(mass_l),1);
c_in_tl = zeros(length(mass_l),p+q); % Molar concentration density
C_in_tl = zeros(length(mass_l),p+q); % Molar concentration

optionsfmincon = optimoptions('fsolve','display','none');

for i = 1:p_zone
    pin_l(i) = fsolve(@(x) sum(c_in(x_piv,MnL/(MwL-MnL),MwL-MnL,x)...
        .*(x_bound(2:end)-x_bound(1:p+q)).*(162*x_piv+18))-mass_l(i),0,...
        optionsfmincon);
    c_in_tl(i,:) = c_in(x_piv,MnL/(MwL-MnL),MwL-MnL,pin_l(i));
end

for i = p_zone+1:p_zone+i_zone
    pin_l(i) = fsolve(@(x) sum(c_in(x_piv,MnH/(MwH-MnH),MwH-MnH,x)...
        .*(x_bound(2:end)-x_bound(1:p+q)).*(162*x_piv+18))-mass_l(i),0,...
        optionsfmincon);
    c_in_tl(i,:) = c_in(x_piv,MnH/(MwH-MnH),MwH-MnH,pin_l(i));
end

for i = 1:p+q
    C_in_tl(:,i) = c_in_tl(:,i).*(x_bound(i+1)-x_bound(i));
end

% Overall initial cellulose distribution

C_T = zeros(p+q,1);

for i = 1:p+q
    C_T(i) = sum(C_in_tl(:,i));
end

% Surface initial cellulose distribution

C_S = zeros(p+q,1);

for i = 1:p+q
    C_S(i) = sum(C_in_tl(1,i));
end

% Internal initial cellulose distribution

C_INT = zeros(p+q,1);

for i = 1:p+q
    C_INT(i) = sum(C_in_tl(2:end,i));
end

% Ratio of molar concentration to total mass in each layer 

C_ratio = zeros(length(R_set)-1,p+q);

for i = 1:length(R_set)-1
    for j = 1:p+q
        C_ratio(i,j) = C_in_tl(i,j)/sum(C_in_tl(i,:)'.*(162*x_piv+18));
    end
end

% Piece-wise function for ratio of molar concentration to total mass as a
% function of particle radius

exp_ratio = cell(p+q,1);

for i = 1:p+q
    exp_ratio{i} = mkpp(flipud(R_set),flipud(C_ratio(:,i)));
end

% General properties of initial cellulose distribution

MnT = sum(C_T.*x_piv)/sum(C_T);             % Overall number-averaged DP
MwT = sum(C_T.*(x_piv.^2))/sum(C_T.*x_piv); % Overall mass-averaged DP
pdiT = MwT/MnT;                             % Overall poly-dispersity index

% Enzyme loadings, mol/L

E_T_CBH = ((L_CBH/1000)*sum(C_T.*(162*x_piv+18)))/MW_CBHi;
E_T_EG = ((L_EG/1000)*sum(C_T.*(162*x_piv+18)))/MW_EGi;     
E_T_BG = ((L_BG/1000)*sum(C_T.*(162*x_piv+18)))/MW_BGi;     

end

%%-------------------------------------------------------------------------

function [t,y] = ODE_setting(p,q,C_S,C_INT,prm,x_piv,vm,t_end,E_T_EG,...
    E_T_CBH,E_T_BG,R0,ro,L,n,exp_ratio)

% Initial conditions

init = zeros(6*p+6*q+9,1);
init(1:p+q) = 0;                    % 1 to p+q --> free liberated polymers
init(p+q+1:2*p+2*q) = 0;            % p+q+1 to 2p+2q --> CBH-bound polymers
init(2*p+2*q+1:3*p+3*q) = 0;        % 2p+2q+1 to 3p+3q --> EG-bound polymers
init(3*p+3*q+1:4*p+4*q) = 0;        % 3p+3q+1 to 4p+4q --> EG-CBH-bound polymers
init(4*p+4*q+1:5*p+5*q) = C_S;      % 4p+4q+1 to 5p+5q --> surface-accessible polymer
init(5*p+5*q+1:6*p+6*q) = C_INT;    % 5p+5q+1 to 6p+6q --> internal inaccessible polymer
init(6*p+6*q+1) = 0;                % 6p+6q+1 --> surface adsorbed EG
init(6*p+6*q+2) = 0;                % 6p+6q+2 --> surface adsorbed CBH
init(6*p+6*q+3:6*p+6*q+4) = 0;      % 6p+6q+3 to 6p+6q+4 --> inhibited EG
init(6*p+6*q+5:6*p+6*q+6) = 0;      % 6p+6q+5 to 6p+6q+6 --> inhibited free CBH
init(6*p+6*q+7:6*p+6*q+8) = 0;      % 6p+6q+7 to 6p+6q+8 --> inhibited adsorbed CBH    
init(6*p+6*q+9) = 0;                % 6p+6q+9 --> inhibited BG

% Rate kernels

% CBH

kp_h_CBH = prm(1);                      % Rate constant of hydrolysis, 1/s     
m_h_CBH = 0;
k_h_CBH = kp_h_CBH*(x_piv.^m_h_CBH);    % Rate coefficients of hydrolysis, 1/s

kp_f_CBH = prm(2);                      % Rate constant of complexation, DP.L/mol.s
m_f_CBH = -1;
k_f_CBH = kp_f_CBH*(x_piv.^m_f_CBH);    % Rate coefficients of complexation (L/mol.s)

kp_e_CBH = prm(3);                      % Rate constant of decomplexation, 1/s
m_e_CBH = 0;
k_e_CBH = kp_e_CBH*(x_piv.^m_e_CBH);    % Rate coefficients of decomplexation, 1/s

tau_CBH = prm(4);                       % Enzyme footprint, mol/m2  
k_ads_CBH = prm(5);                     % Adsorption rate constant, L/mol.s
k_des_CBH = prm(6);                     % Desorption rate constant, 1/s

k_If_CBH = [prm(7) prm(8)];             % Forward glucose/cellobiose inhibition, L/mol.s 
k_Ir_CBH = [prm(9) prm(10)];            % Reverse glucose/cellobiose inhibition, 1/s 

% EG

kp_h_EG = prm(16);                      % Rate constant of hydrolysis (insoluble cellulose), DP/s 
m_h_EG = -1;
k_h_EG = kp_h_EG*(x_piv.^m_h_EG);       % Rate coefficients of hydrolysis (insoluble cellulose),1/s

kp_hs_EG = prm(17);                     % Rate constant of hydrolysis (soluble), L/mol.s 
m_hs_EG = 0;
k_hs_EG = kp_hs_EG*(x_piv.^m_hs_EG);    % Rate coefficients of hydrolysis (soluble), L/mol.s
k_hs_EG(2) = 0;                         % EG unable to hydrolyze cellobiose

kp_f_EG = prm(18);                      % Rate constant of complexation, DP.L/mol.s
m_f_EG = -1;
k_f_EG = kp_f_EG*(x_piv.^m_f_EG);       % Rate coefficients of complexation, L/mol.s

kp_e_EG = prm(19);                      % Rate constant of decomplexation, 1/s
m_e_EG = 0;
k_e_EG = kp_e_EG*(x_piv.^m_e_EG);       % Rate coefficients of decomplexation, 1/s

tau_EG = prm(20);                       % Enzyme footprint, mol/m2 
k_ads_EG = prm(21);                     % Adsorption rate constant, L/mol.s  
k_des_EG = prm(22);                     % Desorption rate constant, 1/s

k_If_EG = [prm(23) prm(24)];            % Forward glucose/cellobiose inhibition, L/mol.s 
k_Ir_EG = [prm(25) prm(26)];            % Reverse glucose/cellobiose inhibition, 1/s 

% BG

k_h_BG = prm(13);                       % Rate constant of hydrolysis, L/mol.s

k_If_BG = prm(14);                      % Forward glucose inhibition, L/mol.s 
k_Ir_BG = prm(15);                      % Reverse glucose inhibition, 1/s

% Fixed pivot technique

[n_exo,n_endo] = fp(p,q,vm,x_piv);      % Particle allocation functions

% ODE solver

tstart = 0;                                     
tfinal = t_end;                   

tic;

options = odeset('AbsTol',1e-8,'RelTol',1e-6,...
    'NonNegative',1:(6*p+6*q+9),'OutputFcn',@odeprog,'Events',@odeabort);

[t,y] = ode15s(@(t,y) MLPBM_ode(p,q,n_exo,n_endo,E_T_EG,...
    E_T_CBH,E_T_BG,k_h_CBH,k_h_EG,k_hs_EG,k_h_BG,k_f_CBH,k_f_EG,k_e_CBH,...
    k_e_EG,k_ads_EG,k_des_EG,k_ads_CBH,k_des_CBH,tau_EG,tau_CBH,k_If_CBH,...
    k_Ir_CBH,k_If_EG,k_Ir_EG,k_If_BG,k_Ir_BG,x_piv,R0,ro,L,n,...
    exp_ratio,y,t),[tstart tfinal],init,options);

time_fp_sim =...
    '\nTotal time elapsed for ML-PBM simulation is %.1f seconds.\n';

fprintf(time_fp_sim,toc)

end

%%-------------------------------------------------------------------------

function [n_exo,n_endo] = fp(p,q,vm,x_piv)

% Particle allocation function - fixed pivot technique

% Chain-end dimer scission

n_exo = zeros(p+q,p+q);
n_exo(1,3) = 1;
n_exo(2,3) = 1; n_exo(2,4) = 2; n_exo(2,5:end) = 1;

for i = 3:p+q
    for j = i:p+q
        if j==i && x_piv(j)-2*vm>x_piv(i-1)
            n_exo(i,j) = ((x_piv(j)-2*vm)-x_piv(i-1))/(x_piv(i)-x_piv(i-1));
        elseif j==i && x_piv(j)-2*vm==x_piv(i-1)
            n_exo(i,j) = 0;
        elseif j~=i && x_piv(j)-2*vm>x_piv(i-1) && x_piv(j)-2*vm<x_piv(i)
            n_exo(i,j) = ((x_piv(j)-2*vm)-x_piv(i-1))/(x_piv(i)-x_piv(i-1));
        elseif j~=i && x_piv(j)-2*vm>x_piv(i) && x_piv(j)-2*vm<x_piv(i+1)
            n_exo(i,j) = (x_piv(i+1)-(x_piv(j)-2*vm))/(x_piv(i+1)-x_piv(i));
        elseif j~=i && x_piv(j)-2*vm==x_piv(i)
            n_exo(i,j) = 1;
        elseif j~=i && x_piv(j)-2*vm==x_piv(i-1)
            n_exo(i,j) = 0;
        else
            n_exo(i,j) = 0;
        end
    end
    n_exo(i,i) = n_exo(i,i)-1;
end

% Random scission

n_endo = zeros(p+q,p+q);

for j = 2:p+q
    n_endo(1,j) = 2/(x_piv(j-1));
end

for i = 2:p
    for j = i+1:p+q
        n_endo(i,i) = -1;
        n_endo(i,j) = 2/(x_piv(j-1));
    end
end

for i = p+1:p+q
    for j = i+1:p+q
        n_endo(i,j) = (x_piv(i+1)-x_piv(i-1))/(x_piv(j-1));
    end
    n_endo(i,i) = -1;
end

end

%%-------------------------------------------------------------------------

function ode = MLPBM_ode(p,q,n_exo,n_endo,E_T_EG,E_T_CBH,...
    E_T_BG,k_h_CBH,k_h_EG,k_hs_EG,k_h_BG,k_f_CBH,k_f_EG,k_e_CBH,k_e_EG,...
    k_ads_EG,k_des_EG,k_ads_CBH,k_des_CBH,tau_EG,tau_CBH,k_If_CBH,...
    k_Ir_CBH,k_If_EG,k_Ir_EG,k_If_BG,k_Ir_BG,x_piv,R0,ro,L,n,...
    exp_ratio,y,t)

dydt = zeros(6*p+6*q+9,1);

% Enzymes

E_F_EG = E_T_EG-y(6*p+6*q+1)-sum(y(2*p+2*q+1:4*p+4*q))...
    -sum(y(6*p+6*q+3:6*p+6*q+4));
E_F_CBH = E_T_CBH-y(6*p+6*q+2)-sum(y(p+q+1:2*p+2*q))...
    -sum(y(3*p+3*q+1:4*p+4*q))-sum(y(6*p+6*q+5:6*p+6*q+8));
E_F_BG = E_T_BG - y(6*p+6*q+9);

% Particle radius

x3 = [x_piv;x_piv;x_piv];
x5 = [x_piv;x_piv;x_piv;x_piv;x_piv];
R = sqrt(sum(y(p+q+1:6*p+6*q).*(162*x5+18))/(ro*pi*L*n));

if R-R0<=0;
    R_ratio = 0;
else
    R_ratio = (R-R0)/R;
end

% Total surface area

As = (n*2*pi*R*L)/1000;

% Soluble polymers

dydt(1) = sum(n_exo(1,:)'.*k_h_CBH.*y(p+q+1:2*p+2*q))...
    +sum(n_endo(1,7:end)'.*k_h_EG(7:end).*y(2*p+2*q+7:3*p+3*q))...
    +sum(n_endo(1,4:end)'.*k_h_EG(4:end).*y(3*p+3*q+4:4*p+4*q))...
    +sum(E_F_EG*n_endo(1,2:6)'.*k_hs_EG(2:6).*y(2:6))...
    +2*k_h_BG*E_F_BG*y(2)+k_Ir_CBH(1)*y(6*p+6*q+5)+k_Ir_CBH(1)*y(6*p+6*q+7)...
    +k_Ir_EG(1)*y(6*p+6*q+3)+k_Ir_BG*y(6*p+6*q+9)...
    -k_If_CBH(1)*E_F_CBH*y(1)-k_If_CBH(1)*y(6*p+6*q+2)*y(1)...
    -k_If_EG(1)*E_F_EG*y(1)-k_If_BG*E_F_BG*y(1);

dydt(2) = sum(n_exo(2,:)'.*k_h_CBH.*y(p+q+1:2*p+2*q))...
    +sum(n_endo(2,7:end)'.*k_h_EG(7:end).*y(2*p+2*q+7:3*p+3*q))...
    +sum(n_endo(2,4:end)'.*k_h_EG(4:end).*y(3*p+3*q+4:4*p+4*q))...
    +sum(E_F_EG*n_endo(2,3:6)'.*k_hs_EG(3:6).*y(3:6))...
    -k_hs_EG(2)*E_F_EG*y(2)-k_h_BG*E_F_BG*y(2)...
    +k_Ir_CBH(2)*y(6*p+6*q+6)+k_Ir_CBH(2)*y(6*p+6*q+8)+k_Ir_EG(2)*y(6*p+6*q+4)...
    -k_If_CBH(2)*E_F_CBH*y(2)-k_If_CBH(2)*y(6*p+6*q+2)*y(2)...
    -k_If_EG(2)*E_F_EG*y(2);

for i = 3:5;
    dydt(i) = sum(n_endo(i,7:end)'.*k_h_EG(7:end).*y(2*p+2*q+7:3*p+3*q))...
        +sum(n_endo(i,i+1:end)'.*k_h_EG(i+1:end).*y(3*p+3*q+i+1:4*p+4*q))...
        +sum(E_F_EG*n_endo(i,i+1:6)'.*k_hs_EG(i+1:6).*y(i+1:6))...
        -k_hs_EG(i)*E_F_EG*y(i);
end

for i = 6;
    dydt(i) = sum(n_endo(i,i+1:end)'.*k_h_EG(i+1:end).*y(2*p+2*q+i+1:3*p+3*q))...
        +sum(n_endo(i,i+1:end)'.*k_h_EG(i+1:end).*y(3*p+3*q+i+1:4*p+4*q))...
        -k_hs_EG(i)*E_F_EG*y(i);
end

% CBH-bound surface polymers

dydt(p+q+3) = sum(n_exo(3,3:end)'.*k_h_CBH(3:end).*y(p+q+3:2*p+2*q));


for i = p+q+4:p+q+6;
    dydt(i) = sum(n_exo(i-p-q,:)'.*k_h_CBH.*y(p+q+1:2*p+2*q))...
        -(k_f_EG(i-p-q)*y(6*p+6*q+1)*y(i))...
        +(k_e_EG(i-p-q)*y(i-p-q+3*p+3*q));

end

for i = p+q+7:2*p+2*q-1;
    dydt(i) = sum(n_exo(i-p-q,:)'.*k_h_CBH.*y(p+q+1:2*p+2*q))...
        +(y(6*p+6*q+2)*k_f_CBH(i-p-q)*y(i-p-q+4*p+4*q))...
        -(k_e_CBH(i-p-q)*y(i))...
        -(k_f_EG(i-p-q)*y(6*p+6*q+1)*y(i))...
        +(k_e_EG(i-p-q)*y(i-p-q+3*p+3*q));

end

for i = 2*p+2*q;
    dydt(i) = sum(n_exo(i-p-q,:)'.*k_h_CBH.*y(p+q+1:2*p+2*q))...
        +(y(6*p+6*q+2)*k_f_CBH(i-p-q)*y(i-p-q+4*p+4*q))...
        -(k_e_CBH(i-p-q)*y(i))...
        -(k_f_EG(i-p-q)*y(6*p+6*q+1)*y(i))...
        +(k_e_EG(i-p-q)*y(i-p-q+3*p+3*q));
end

% EG-bound surface polymers

for i = 2*p+2*q+7:3*p+3*q;
    dydt(i) = (k_f_EG(i-2*p-2*q)*y(6*p+6*q+1)*y(i-2*p-2*q+4*p+4*q))...
        -(k_e_EG(i-2*p-2*q)*y(i))...
        -(k_h_EG(i-2*p-2*q)*y(i));
end

% CBH-EG-bound surface polymers

for i = 3*p+3*q+4:4*p+4*q;
    dydt(i) = (k_f_EG(i-3*p-3*q)*y(6*p+6*q+1)*y(i-3*p-3*q+p+q))...
        -(k_e_EG(i-3*p-3*q)*y(i))...
        -(k_h_EG(i-3*p-3*q)*y(i));
end

% Surface polymers

m = zeros(p+q,1);

for i = 4*p+4*q+7:5*p+5*q-1;
    m(i-4*p-4*q) = ...
        (sum(n_endo(i-4*p-4*q,i-4*p-4*q+1:end)'.*k_h_EG(i-4*p-4*q+1:end)...
        .*y(i-4*p-4*q+2*p+2*q+1:3*p+3*q))...
        +sum(n_endo(i-4*p-4*q,i-4*p-4*q+1:end)'.*k_h_EG(i-4*p-4*q+1:end)...
        .*y(i-4*p-4*q+1+3*p+3*q:4*p+4*q))...
        -(k_f_CBH(i-4*p-4*q)*y(6*p+6*q+2)*y(i))...
        +(k_e_CBH(i-4*p-4*q)*y(i-4*p-4*q+p+q))...
        -(k_f_EG(i-4*p-4*q)*y(6*p+6*q+1)*y(i))...
        +(k_e_EG(i-4*p-4*q)*y(i-4*p-4*q+2*p+2*q)));
end

m(p+q) = (-(k_f_CBH(p+q)*y(6*p+6*q+2)*y(5*p+5*q))...
    +(k_e_CBH(p+q)*y(2*p+2*q))...
    -(k_f_EG(p+q)*y(6*p+6*q+1)*y(5*p+5*q))...
    +(k_e_EG(p+q)*y(3*p+3*q)));

for i = 4*p+4*q+7:5*p+5*q;
    dydt(i) = m(i-4*p-4*q)...
        -(sum(m.*(162*x_piv+18))+sum(dydt(p+q+1:4*p+4*q).*(162*x3+18)))...
        *R_ratio*ppval(exp_ratio{i-4*p-4*q},R-R0);
end

% Internal inaccessible polymer

for i = 5*p+5*q+7:6*p+6*q;
    dydt(i) = (sum(m.*(162*x_piv+18))+sum(dydt(p+q+1:4*p+4*q)...
        .*(162*x3+18)))*R_ratio*ppval(exp_ratio{i-5*p-5*q},R-R0);
end

% Surface-adsorbed EG

dydt(6*p+6*q+1) = (k_ads_EG*E_F_EG*(As*tau_EG-y(6*p+6*q+1)))...
    -(k_des_EG*y(6*p+6*q+1))...
    -sum(y(6*p+6*q+1)*k_f_EG(7:end).*y(4*p+4*q+7:5*p+5*q))...
    +sum(k_e_EG(7:end).*y(2*p+2*q+7:3*p+3*q))...
    -sum(y(6*p+6*q+1)*k_f_EG(4:end).*y(p+q+4:2*p+2*q))...
    +sum(k_e_EG(4:end).*y(3*p+3*q+4:4*p+4*q));

% Surface-adsorbed CBH

dydt(6*p+6*q+2) = (k_ads_CBH*E_F_CBH*(As*tau_CBH-y(6*p+6*q+2)))...
    -(k_des_CBH*y(6*p+6*q+2))...
    -sum(y(6*p+6*q+2)*k_f_CBH(7:end).*y(4*p+4*q+7:5*p+5*q))...
    +sum(k_e_CBH(7:end).*y(p+q+7:2*p+2*q))...
    +k_Ir_CBH(1)*y(6*p+6*q+7)+k_Ir_CBH(2)*y(6*p+6*q+8)...
    -y(6*p+6*q+2)*sum(k_If_CBH'.*y(1:2));

% Inhibited EG

for i = 6*p+6*q+3:6*p+6*q+4;
    dydt(i) = k_If_EG(i-6*p-6*q-3+1)*E_F_EG*y(i-6*p-6*q-3+1)...
        -k_Ir_EG(i-6*p-6*q-3+1)*y(i);
end

% Inhibited CBH

for i = 6*p+6*q+5:6*p+6*q+6;
    dydt(i) = k_If_CBH(i-6*p-6*q-5+1)*E_F_CBH*y(i-6*p-6*q-5+1)...
        -k_Ir_CBH(i-6*p-6*q-5+1)*y(i);
end

for i = 6*p+6*q+7:6*p+6*q+8;
    dydt(i) = k_If_CBH(i-6*p-6*q-7+1)*y(6*p+6*q+2)*y(i-6*p-6*q-7+1)...
        -k_Ir_CBH(i-6*p-6*q-7+1)*y(i);
end

% Inhibited BG

for i = 6*p+6*q+9;
    dydt(i) = k_If_BG*E_F_BG*y(i-6*p-6*q-9+1)...
        -k_Ir_BG*y(i);
end

% Output

ode = dydt;

end

%%-------------------------------------------------------------------------

function [C,CXS,CNS,CNXS,CFS,CS,CI,CT,CC,E_S_EG,E_S_CBH,E_F_EG,E_F_CBH,...
    E_F_BG,R,conv_rem] = post_processing(t,y,p,q,x_piv,ro,L,n,mss,E_T_EG,...
    E_T_CBH,E_T_BG)

% Molar concentrations, mol/L

C = y(:,1:p+q);                     % Soluble products
CXS = y(:,p+q+1:2*p+2*q);           % CBH-bound complex
CNS = y(:,2*p+2*q+1:3*p+3*q);       % EG-bound complex
CNXS = y(:,3*p+3*q+1:4*p+4*q);      % CBH-EG-bound complex
CFS = y(:,4*p+4*q+1:5*p+5*q);       % Free un-bound surface polymers
CS = CXS+CNS+CNXS+CFS;              % Total surface polymers
CI = y(:,5*p+5*q+1:6*p+6*q);        % Internal polymers 
CT = C+CS+CI;                       % Total polymers
CC = CS+CI;                         % Insoluble polymers

E_S_EG = y(:,6*p+6*q+1);            % Surface-adsorbed EG
E_S_CBH = y(:,6*p+6*q+2);           % Surface-adsorbed CBH

CNI(:,1) = y(:,6*p+6*q+3);          % Inhibited EG
CNI(:,2) = y(:,6*p+6*q+4);

CXI(:,1) = y(:,6*p+6*q+5)+y(:,6*p+6*q+7);   % Inhibited CBH
CXI(:,2) = y(:,6*p+6*q+6)+y(:,6*p+6*q+8);

CBI = y(:,6*p+6*q+9);               % Inhibited BG
   
E_F_EG = zeros(length(t),1);        % Free EG
E_F_CBH = zeros(length(t),1);       % Free CBH
E_F_BG = zeros(length(t),1);        % Free BG

for i = 1:length(t)
    E_F_EG(i) = E_T_EG-E_S_EG(i)-sum(CNS(i,:))-sum(CNXS(i,:))-sum(CNI(i,:));
    E_F_CBH(i) = E_T_CBH-E_S_CBH(i)-sum(CXS(i,:))-sum(CXI(i,:))-sum(CNXS(i,:));
    E_F_BG(i) = E_T_BG-sum(CBI(i));
end

% Transient particle radius, m

R = zeros(length(t),1);

for i = 1:length(t)
    R(i) = sqrt(sum(CC(i,:)'.*(162*x_piv+18))/(ro*pi*L*n));
end

% Mass of insoluble polymers, g/L

mss_CC = zeros(length(t),1);

for i = 1:length(t)
    mss_CC(i) = sum(CC(i,:)'.*(162*x_piv+18));
end

% Conversion

conv_rem = zeros(length(t),1);

for i = 1:length(t)
    conv_rem(i) = (mss-mss_CC(i))/mss;
end

end