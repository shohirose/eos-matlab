import eos.VanDerWaalsEos

%% Test 1: Methane
% Methane
Pc = 4.6e6;     % Critical pressure [Pa]
Tc = 190.6;     % Critical temperature [K]
omega = 0.008;  % Acentric factor
Mw = 16.0425;   % Molecular weight [g/mol]
vdw = VanDerWaalsEos(Pc,Tc,Mw);

P = 3.0e6;      % Pressure [Pa]
T = 180;        % Temperature [K]
[z,result] = vdw.zFactors(P,T);
phi = vdw.fugacityCoeff(z,result);
assert(abs((z - 0.6999)/z) < 1e-3);
assert(abs((phi - 0.7769)/phi) < 1e-3);

%% Test 2: Mixture of Methane and Propane
% CH4, C3H8
Pc = [4.6e6, 4.246e6]';     % Critical pressure [Pa]
Tc = [190.6, 369.8]';       % Critical temperature [K]
omega = [0.008, 0.152]';    % Acentric factor
Mw = [16.0425, 44.1]';      % Molecular weight [g/mol]
K = [0, 0.009;
     0.009, 0];             % Binary interaction parameters
x = [0.85, 0.15]';          % Overall composition

vdw = VanDerWaalsEos(Pc,Tc,Mw,K);

P = 3.0e6;                  % Pressure [Pa]
T = 180;                    % Temperature [K]
[z,result] = vdw.zFactors(P,T,x);
phi = vdw.fugacityCoeff(z,result);
assert(max(abs((z - 0.158531)./z)) < 1e-3);
assert(max(abs((phi - [0.870881, 0.048737]')./phi)) < 1e-3);
