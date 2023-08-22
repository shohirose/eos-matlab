import eos.VanDerWaalsEos
import eos.SoaveRedlichKwongEos
import eos.PengRobinsonEos

% Chlorine
Pc = 76.1*101325;   % Pa
Tc = 417;           % K
omega = 0.090;
Mw = 70.91;         % g/mol

vdw = VanDerWaalsEos(Pc,Tc,Mw);
srk = SoaveRedlichKwongEos(Pc,Tc,omega,Mw);
pr = PengRobinsonEos(Pc,Tc,omega,Mw);

volrange = @(eos) logspace(log10(1.05*eos.RepulsionParam),-1,200); % [m3]
V1 = volrange(vdw);
V2 = volrange(srk);
V3 = volrange(pr);

%% T = 521 K
T = 521;    % K
P1 = vdw.pressure(T,V1);
P2 = srk.pressure(T,V2);
P3 = pr.pressure(T,V3);
figure;
loglog(V1,P1,V2,P2,V3,P3);
axis([1e-5,0.1,1e5,1e9]);
legend('VDW','SRK','PR');
xlabel('Volume [m3/mol]');
ylabel('Pressure [Pa]');
s = sprintf('T = %.1f K',T);
title(s);

%% T = 417 K
T = 417;    % K
P1 = vdw.pressure(T,V1);
P2 = srk.pressure(T,V2);
P3 = pr.pressure(T,V3);
figure;
loglog(V1,P1,V2,P2,V3,P3);
axis([1e-5,0.1,1e5,1e9]);
legend('VDW','SRK','PR');
xlabel('Volume [m3/mol]');
ylabel('Pressure [Pa]');
s = sprintf('T = %.1f K',T);
title(s);

%% T = 293 K
T = 293;    % K
P1 = vdw.pressure(T,V1);
P2 = srk.pressure(T,V2);
P3 = pr.pressure(T,V3);
figure;
semilogx(V1,P1,V2,P2,V3,P3);
axis([1e-5,0.1,-5e7,5e7]);
legend('VDW','SRK','PR');
xlabel('Volume [m3/mol]');
ylabel('Pressure [Pa]');
s = sprintf('T = %.1f K',T);
title(s);