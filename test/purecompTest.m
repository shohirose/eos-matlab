import eos.VanDerWaalsEos
import eos.SoaveRedlichKwongEos
import eos.PengRobinsonEos

% Methane
Pc = 4.6e6;     % Critical pressure [Pa]
Tc = 190.6;     % Critical temperature [K]
omega = 0.008;  % Acentric factor
Mw = 16.0425;   % Molecular weight [g/mol]

P = 3.0e6;      % Pressure [Pa]
T = 180;        % Temperature [K]

%% VanDerWaalsEos test:
eos = VanDerWaalsEos(Pc,Tc,Mw);
[z,s] = eos.zFactors(P,T);
phi = eos.fugacityCoeff(z,s);
assert(abs((z - 0.6999)/z) < 1e-3);
assert(abs((phi - 0.7769)/phi) < 1e-3);

%% SoaveRedlichKwongEos test:
eos = SoaveRedlichKwongEos(Pc,Tc,omega,Mw);
[z,s] = eos.zFactors(P,T);
z = sort(z);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((z - [0.1393, 0.2132, 0.6475]')./z)) < 1e-3);
assert(max(abs((phi - [0.7808, 0.7863, 0.7441]')./phi)) < 1e-3);

%% PengRobinsonEos test:
eos = PengRobinsonEos(Pc,Tc,omega,Mw);
[z,s] = eos.zFactors(P,T);
z = sort(z);
phi = eos.fugacityCoeff(z,s);
assert(max(abs((z - [0.1235, 0.1988, 0.6240]')./z)) < 1e-3);
assert(max(abs((phi - [0.7591, 0.7654, 0.7256]')./phi)) < 1e-3);
