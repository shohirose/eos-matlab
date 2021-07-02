# eos-matlab

This is a package for cubic equations of state (EoS) for MATLAB.

# Overview

- `eos` : EoS package
  - `VanDerWaalsEos` : Van der Waals EoS.
  - `SoaveRedlichKwongEos` : Soave-Redlich-Kwong EoS.
  - `PengRobinsonEos` : Peng-Robinson EoS.

# How to Use

Please copy the `+eos` folder wherever you want in your computer, and add the path to the folder where `+eos` exists using `addpath` function in the MATLAB terminal.

```matlab
> % Add the path to the package folder
> addpath('path-to-the-folder-containing-eos-package')
> % Then you can import anything in the package.
> import eos.VanDerWaalsEos
> Pc = 4.6e6;
> Tc = 190.6;
> Mw = 16.0425;
> eos = VanDerWaalsEos(Pc,Tc,Mw);
```

# License

MIT License. Please refer to the `LICENSE` file under the project root directory.