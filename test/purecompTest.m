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
result = eos.zFactors(P,T);
phi = eos.fugacityCoeff(result);
z = result.z;
assert(abs((z - 0.6999)/z) < 1e-3);
assert(abs((phi - 0.7769)/phi) < 1e-3);

%% SoaveRedlichKwongEos test:
eos = SoaveRedlichKwongEos(Pc,Tc,omega,Mw);
result = eos.zFactors(P,T);
z = result.z;
phi = eos.fugacityCoeff(result);
assert(max(abs((z - [0.1393, 0.2132, 0.6475]')./z)) < 1e-3);
assert(max(abs((phi - [0.7808, 0.7863, 0.7441]')./phi)) < 1e-3);

%% PengRobinsonEos test:
eos = PengRobinsonEos(Pc,Tc,omega,Mw);
result = eos.zFactors(P,T);
z = result.z;
phi = eos.fugacityCoeff(result);
assert(max(abs((z - [0.1237, 0.1984, 0.6242]')./z)) < 1e-3);
assert(max(abs((phi - [0.7595, 0.7657, 0.7257]')./phi)) < 1e-3);
