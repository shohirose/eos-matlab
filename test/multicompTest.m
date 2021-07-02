import eos.VanDerWaalsEos
import eos.SoaveRedlichKwongEos
import eos.PengRobinsonEos

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

%% Test 1: VanDerWaalsEos
eos = VanDerWaalsEos(Pc,Tc,Mw,K);
[z,s] = eos.zFactors(P,T,x);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((z - 0.158531)./z)) < 1e-3);
assert(max(abs((phi - [0.870881, 0.048737]')./phi)) < 1e-3);

%% Test 2: SoaveRedlichKwongEos
eos = SoaveRedlichKwongEos(Pc,Tc,omega,Mw,K);
[z,s] = eos.zFactors(P,T,x);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((z - 0.109951)./z)) < 1e-3);
assert(max(abs((phi - [0.835510, 0.003603]')./phi)) < 1e-3);

%% Test 3: PengRobinsonEos
eos = PengRobinsonEos(Pc,Tc,omega,Mw,K);
[z,s] = eos.zFactors(P,T,x);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((z - 0.097246)./z)) < 1e-3);
assert(max(abs((phi - [0.809473, 0.003713]')./phi)) < 1e-3);