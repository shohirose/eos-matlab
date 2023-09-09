import eos.VanDerWaalsEos
import eos.SoaveRedlichKwongEos
import eos.PengRobinsonEos
import eos.BubblePointPressureCalculator

Pc = 4.6e6; % [Pa]
Tc = 190.6; % [K]
omega = 0.008;
Mw = 16.0425; % [g-mole]
vdw = VanDerWaalsEos(Pc,Tc,Mw);
srk = SoaveRedlichKwongEos(Pc,Tc,omega,Mw);
pr = PengRobinsonEos(Pc,Tc,omega,Mw);

%% Test 1: Vapor pressure of methane with vdW EoS
T = 150;
calculator = BubblePointPressureCalculator(vdw,1e-6,50);
[Pbub,report] = calculator.compute(T,0);
assert(abs(Pbub - 1.634e6)/Pbub < 1e-3);

%% Test 2: Vapor pressure of methane with SRK EoS
T = 150;
calculator = BubblePointPressureCalculator(srk,1e-6,50);
[Pbub,report] = calculator.compute(T,0);
assert(abs(Pbub - 1.0550e6)/Pbub < 1e-3);


%% Test 3: Vapor pressure of methane with PR EoS
T = 150;
calculator = BubblePointPressureCalculator(pr,1e-6,50);
[Pbub,report] = calculator.compute(T,0);
assert(abs(Pbub - 1.0511e6)/Pbub < 1e-3);

%% Test 4: Bubble point pressure of CH4-N2 system with PR EoS
Pc = [4.604e6, 3.394e6];
Tc = [190.6, 126.1];
omega = [0.011, 0.04];
Mw = [16.04, 28.01];
K = [0.0, 0.03; 0.03, 0.0];
T = 120;
x = [0.5, 0.5];
pr = PengRobinsonEos(Pc,Tc,omega,Mw,K);
calculator = BubblePointPressureCalculator(pr,1e-6,50);
[Pbub,report] = calculator.compute(T,0,0,x);
assert(abs(Pbub - 1.4041e6)/Pbub < 1e-3);