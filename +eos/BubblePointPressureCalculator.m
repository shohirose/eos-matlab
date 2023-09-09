classdef BubblePointPressureCalculator
    %BUBBLEPOINTPRESSURECALCULATOR Calculate bubble point pressure    
    properties
        Eos
        Tolerance
        MaxIter
    end

    methods (Access = private)
        function Pbub = estimateBubblePointPressure(obj,T,xall)
            %ESTIMATEBUBBLEPOINTPRESSURE Estimate bubble point pressure
            % using the Wilson equation.
            %
            % Parameters
            % ----------
            % T    : Temperature
            % xall : Overall composition
            %
            % Returns
            % -------
            % Pbub : Bubble point pressure
            arguments
                obj
                T (1,1) {mustBeNumeric}
                xall (:,1) {mustBeNumeric} = 1
            end
            Pc = obj.Eos.CriticalPressure;
            Tc = obj.Eos.CriticalTemperature;
            omega = obj.Eos.AcentricFactor;
            Pbub = Pc.*10.^(7.0/3.0*(1 + omega).*(1 - Tc/T));
            if nargin == 3
                Pbub = dot(xall,Pbub);
            end
        end

        function K = estimateEquilibriumRatio(obj,P,T)
            % ESTIMATEEQUILIBRIUMRATIO Estimate equilibrium ratio using the
            % Wilson equation.
            %
            % Parameters
            % ----------
            % P  : Pressure
            % T  : Temperature
            %
            % Returns
            % -------
            % K : Equilibrium ratio
            Pc = obj.Eos.CriticalPressure;
            Tc = obj.Eos.CriticalTemperature;
            omega = obj.Eos.AcentricFactor;
            K = Pc/P.*exp(5.373*(1 + omega).*(1 - Tc/T));
        end

        function [Pb,report] = computeForPureComponent(obj,T,Pinit)
            %COMPUTEFORPURECOMPONENT Compute bubble point pressure for pure
            % component systems.
            %
            % Parameters
            % ----------
            % T     : Temperature [K]
            % Pinit : Initial pressure [Pa]
            %
            % Returns
            % -------
            % Pb     : Bubble point pressure
            % report : Iteration report
            arguments
                obj
                T (1,1) {mustBeNumeric}
                Pinit (1,1) {mustBeNumeric}
            end
            Pb = Pinit;
            report.P = zeros(obj.MaxIter,1);
            report.Eps = zeros(obj.MaxIter,1);

            for i = 1:obj.MaxIter
                [zL,resultL] = obj.Eos.zFactors(Pb,T);
                [zV,resultV] = obj.Eos.zFactors(Pb,T);
                zV = max(zV);
                zL = min(zL);
                if abs(zL - zV)/zL < 1e-6
                    warning('The same Z-factors are obtained.');
                    break;
                end
                phiV = obj.Eos.fugacityCoeff(zV,resultV);
                phiL = obj.Eos.fugacityCoeff(zL,resultL);
                f = phiL/phiV;
                Pb = Pb*f;
                eps = abs(1 - f);
                
                report.P(i) = Pb;
                report.Eps(i) = eps;
                report.NumIter = i;
                if eps < obj.Tolerance
                    report.P((i+1):obj.MaxIter) = [];
                    report.Eps((i+1):obj.MaxIter) = [];
                    break;
                elseif i == obj.MaxIter
                    warning('Max iteration reached.');
                end
            end
        end

        function [Pb,report] = computeForMultiComponents(...
                obj,T,Pinit,Kinit,xall)
            %COMPUTEFORMULTICOMPONENTS Compute bubble point pressure for
            % multi-component systems.
            %
            % Parameters
            % ----------
            % T     : Temperature [K]
            % Pinit : Initial pressure [Pa]
            % Kinit : Initial equilibrium ratio
            % xall  : Overall composition
            %
            % Returns
            % -------
            % Pb     : Bubble point pressure
            % report : Iteration report
            arguments
                obj
                T (1,1) {mustBeNumeric}
                Pinit (1,1) {mustBeNumeric}
                Kinit (:,1) {mustBeNumeric}
                xall (:,1) {mustBeNumeric}
            end

            % Multi-components
            Pb = Pinit;
            K = Kinit;
            xL = xall;
            xV = K.*xL;
            report.P = zeros(obj.MaxIter,1);
            report.Eps = zeros(obj.MaxIter,1);
            
            for i = 1:obj.MaxIter
                [zL,resultL] = obj.Eos.zFactors(Pb,T,xL);
                [zV,resultV] = obj.Eos.zFactors(Pb,T,xV);
                phiV = obj.Eos.fugacityCoeff(max(zV),resultV);
                phiL = obj.Eos.fugacityCoeff(min(zL),resultL);
                f = sum(xL.*phiL./phiV);
                Pb = Pb*f;
                eps = abs(1 - f);
                
                report.P(i) = Pb;
                report.Eps(i) = eps;
                report.NumIter = i;
                if eps < obj.Tolerance
                    report.P((i+1):obj.MaxIter) = [];
                    report.Eps((i+1):obj.MaxIter) = [];
                    break;
                elseif i == obj.MaxIter
                    warning('Max iteration reached.');
                end
            end
        end
    end
    
    methods
        function obj = BubblePointPressureCalculator(eos,tol,maxiter)
            %BUBBLEPOINTPRESSURECALCULATOR Create an object of this class
            %
            % Parameters
            % ----------
            % eos    : Cubic Eos
            % tol    : Tolerance for convergence
            % maxiter: Maximum iterations
            obj.Eos = eos;
            obj.Tolerance = tol;
            obj.MaxIter = maxiter;
        end
        
        function [Pb, report] = compute(obj,T,Pinit,Kinit,xall)
            %COMPUTE Calculate bubble point pressure using successive
            % substitution method.
            %
            % Parameters
            % ----------
            % T     : Temperature [K]
            % Pinit : Initial pressure [Pa]. If zero is passed, it will be
            %         automatically calculated.
            % Kinit : Initial equilibrium ratio. If zero is passed, it will
            %         be automatically calculated.
            % xall  : Overall composition
            %
            % Returns
            % -------
            % Pb     : Bubble point pressure
            % report : Iteraction report
            arguments
                obj
                T {mustBeNumeric}
                Pinit {mustBeNumeric}
                Kinit (:,1) {mustBeNumeric} = 1
                xall (:,1) {mustBeNumeric} = 1
            end

            if nargin == 3
                if Pinit == 0
                    Pinit = obj.estimateBubblePointPressure(T);
                end
                [Pb,report] = obj.computeForPureComponent(T,Pinit);
            else
                if Pinit == 0
                    Pinit = obj.estimateBubblePointPressure(T,xall);
                end
                if Kinit == 0
                    Kinit = obj.estimateEquilibriumRatio(Pinit,T);
                end
                [Pb,report] = obj.computeForMultiComponents(T,Pinit,...
                    Kinit,xall);
            end
        end
    end
end

