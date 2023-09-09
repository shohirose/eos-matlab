import eos.SoaveRedlichKwongEos

%% Test 1: Methane
% Methane
Pc = 4.6e6;     % Critical pressure [Pa]
Tc = 190.6;     % Critical temperature [K]
omega = 0.008;  % Acentric factor
Mw = 16.0425;   % Molecular weight [g/mol]
P = 3.0e6;      % Pressure [Pa]
T = 180;        % Temperature [K]
srk = SoaveRedlichKwongEos(Pc,Tc,omega,Mw);
[z,result] = srk.zFactors(P,T);
phi = srk.fugacityCoeff(z,result);
assert(max(abs((z - [0.1393, 0.2132, 0.6475]')./z)) < 1e-3);
assert(max(abs((phi - [0.7808, 0.7863, 0.7441]')./phi)) < 1e-3);

%% Test 2: Mixture of Methane and Propane
% CH4, C3H8
Pc = [4.6e6, 4.246e6]';     % Critical pressure [Pa]
Tc = [190.6, 369.8]';       % Critical temperature [K]
omega = [0.008, 0.152]';    % Acentric factor
Mw = [16.0425, 44.1]';      % Molecular weight [g/mol]
K = [0, 0.009;
     0.009, 0];             % Binary interaction parameters
P = 3.0e6;                  % Pressure [Pa]
T = 180;                    % Temperature [K]
x = [0.85, 0.15]';          % Overall composition
srk = SoaveRedlichKwongEos(Pc,Tc,omega,Mw,K);
[z,result] = srk.zFactors(P,T,x);
phi = srk.fugacityCoeff(z,result);
assert(max(abs((z - 0.109951)./z)) < 1e-3);
assert(max(abs((phi - [0.835510, 0.003603]')./phi)) < 1e-3);