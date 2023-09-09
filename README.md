# eos-matlab

This is a package for cubic equations of state (EoS) for MATLAB.

# Overview

- `eos` : EoS package
  - `VanDerWaalsEos` : Van der Waals EoS.
  - `SoaveRedlichKwongEos` : Soave-Redlich-Kwong EoS.
  - `PengRobinsonEos` : Peng-Robinson EoS.
  - `BubblePointPressureCalculator` : Bubble point pressure calculator.

# How to Use

Please copy the `+eos` folder wherever you want in your computer, and add the path to the folder where `+eos` exists using `addpath` function in the MATLAB terminal.

```matlab
> % Add the path to the package folder
> addpath('path-to-the-folder-containing-eos-package')
```

An example of usage is:

```matlab
import eos.PengRobinsonEos

Pc = 4.6e6; % Critical pressure [Pa]
Tc = 190.6; % Critical temperature [K]
Mw = 16.0425; % Molecular weight [g-mole]
omega = 0.008; % Acentric factor

% Create an object of PR EoS
preos = PengRobinson(Pc,Tc,omega,Mw);
P = 3e6; % pressure [Pa]
T = 180; % temperature [K]

% Compute Z-factor
[z,params] = preos.zFactors(P,T);
% Compute fugacity coefficients
phi = preos.fugacityCoeff(z,params);
```

If you need a isothermal line, you can calculate:

```matlab
% Compute isothermal line
V = linspace(1e-3,1e-2);
Pnew = preos.pressure(T,V);
```

# License

MIT License. Please refer to the `LICENSE` file under the project root directory.