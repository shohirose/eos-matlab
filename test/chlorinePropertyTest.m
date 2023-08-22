import eos.VanDerWaalsEos
import eos.SoaveRedlichKwongEos
import eos.PengRobinsonEos

% Cl2
Pc = 76.1*101325;   % Pa
Tc = 417;           % K
omega = 0.090;
Mw = 70.91;         % g/mol

P = 15.5e6; % Pa
T = 521;    % K

%% Test 1: VanDerWaalsEos
eos = VanDerWaalsEos(Pc,Tc,Mw);
Pr = P/Pc;
Tr = T/Tc;
A = eos.reducedAttractionParam(Pr,Tr,1);
B = eos.reducedRepulsionParam(Pr,Tr);
assert(abs(A - 0.54326)/A < 1e-3);
assert(abs(B - 0.2011)/B < 1e-3);

[z,s] = eos.zFactors(P,T);
assert(length(z) == 1);
assert(abs(z - 0.5983)/z < 1e-3);

phi = eos.fugacityCoeff(z,s);
assert(abs(phi - 0.6795)/phi < 1e-3);

%% Test 2: SoaveRedlichKwongEos
eos = SoaveRedlichKwongEos(Pc,Tc,omega,Mw);
[z,s] = eos.zFactors(P,T);
assert(length(z) == 1);
assert(abs(z - 0.6808)/z < 1e-3);

phi = eos.fugacityCoeff(z,s);
assert(abs(phi - 0.7134)/phi < 1e-3);

%% Test 3: PengRobinsonEos
eos = PengRobinsonEos(Pc,Tc,omega,Mw);
[z,s] = eos.zFactors(P,T);
assert(length(z) == 1);
assert(abs(z - 0.6434)/z < 1e-3);

phi = eos.fugacityCoeff(z,s);
assert(abs(phi - 0.6782)/phi < 1e-3);