classdef VaporLiquidFlash
    % VAPORLIQUIDFLASH Flash calculation for vapor-liquid equilibria.
    
    properties
        Eos
        RootSearchMethod
        Tolerance
        MaxIter
    end
    
    methods
        function obj = VaporLiquidFlash(eos, method, tol, maxiter)
            % Constructor
            %
            % Parameters
            % ----------
            % eos : Cubic EoS
            % method : str
            %   Root search method, 'SS' or 'Newton'
            % tol : double
            %   Tolerance for convergence
            % maxiter : int
            %   Maximum iteration number
            arguments
                eos {mustBeA(eos, 'eos.CubicEosBase')}
                method string
                tol doubleint{mustBeInteger}
            end
            obj.Eos = eos;
            obj.RootSearchMethod = method;
            obj.Tolerance = tol;
            obj.MaxIter = maxiter;
        end
        
        function [Pvap, report] = vaporPressure(obj, T, Pinit)
            % Calculate vapor pressure
            %
            % Parameters
            % ----------
            % T : double
            %   Temperature [K]
            % Pinit : double
            %   Initial pressure [Pa]
            %
            % Returns
            % -------
            % pvap : double
            %   Vapor pressure
            % report : struct
            %   Iteraction report
            arguments
                obj {mustBeA(obj,'eos.VaporLiquidFlash')}
                T double
                Pinit double = 1
            end
            
            Pvap = Pinit;
            
            report.P = zeros(obj.MaxIter);
            report.Iter = zeros(obj.MaxIter);
            report.Eps = zeros(obj.MaxIter);
            
            for i = 1:obj.MaxIter
                [z,s] = obj.Eos.zFactor(Pvap,T);
                zV = max(z);
                zL = min(z);
                phiV = obj.Eos.fugacityCoeff(zV,s);
                phiL = obj.Eos.fugacityCoeff(zL,s);
                Pvap = Pvap*phiL/phiV;
                eps = abs(1 - phiL/phiV);
                
                report.P(i) = Pvap;
                report.Iter(i) = i;
                report.Eps(i) = eps;
                if eps < obj.Tolerance
                    report.P((i+1):obj.MaxIter) = [];
                    report.Iter((i+1):obj.MaxIter) = [];
                    report.Eps((i+1):obj.MaxIter) = [];
                    break;
                end
            end
        end
    end
end