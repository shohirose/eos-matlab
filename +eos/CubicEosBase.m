classdef CubicEosBase
    % Base class for two-parameter cubic equations of state.
    
    properties (SetAccess = private)
        CriticalPressure    % Critical pressure [Pa]
        CriticalTemperature % Critical temperature [K]
        MolecularWeight     % Molecular weight [g/mol]
        AcentricFactor      % Acentric factor
        CorrectionFactor    % Correction factor for the attraction
                            % parameter
        OmegaA              % Coefficient for attraction parameter
        OmegaB              % Coefficient for repulsion parameter
        AttractionParam     % Attraction parameter
        RepulsionParam      % Repulsion parameter
        MixingRule          % Mixing rule
    end

    methods
        function obj = CubicEosBase(OmegaA,OmegaB,Pc,Tc,Mw,omega,alpha,K)
            % Construct cubic EOS
            %
            % Parameters
            % ----------
            % OmegaA : Coefficient for attraction parameter
            % OmegaB : Coefficient for repulsion parameter
            % Pc     : Critical pressure [Pa]
            % Tc     : Critical temperature [K]
            % Mw     : Molecular weight [g/mol]
            % omega  : Acentric factor
            % alpha  : Correction factor
            % K      : Binary interaction parameter (optional)
            arguments
                OmegaA (1,1) {mustBeNumeric}
                OmegaB (1,1) {mustBeNumeric}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                alpha {eos.mustBeCorrectionFactor}
                K (:,:) {mustBeNumeric} = 1
            end
            obj.OmegaA = OmegaA;
            obj.OmegaB = OmegaB;
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            obj.AcentricFactor = omega;
            obj.CorrectionFactor = alpha;
            if nargin == 8
                obj.MixingRule = eos.MixingRule(K);
            end
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = OmegaA*(R*Tc).^2./Pc;
            obj.RepulsionParam = OmegaB*R*Tc./Pc;
        end
        function obj = setParams(obj,Pc,Tc,Mw,omega,alpha,K)
            % Set parameters
            %
            % Parameters
            % ----------
            % Pc    : Critical pressure [Pa]
            % Tc    : Critical temperature [K]
            % Mw    : Molecular weight [g/mol]
            % omega : Acentric factor
            % alpha : Correction factor
            % K     : Binary interaction parameters (Only for
            %         multi-component systems)
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                alpha {eos.mustBeCorrectionFactor}
                K (:,:) {mustBeNumeric} = 1
            end
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            obj.AcentricFactor = omega;
            obj.CorrectionFactor = alpha;
            if nargin == 7
                obj.MixingRule = eos.MixingRule(K);
            end
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = obj.OmegaA*(R*Tc).^2./Pc;
            obj.RepulsionParam = obj.OmegaB*R*Tc./Pc;
        end
        function Pr = reducedPressure(obj,P)
            % Compute reduced pressure
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            %
            % Returns
            % -------
            % Pr : Reduced pressure
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                P {mustBeNumeric} 
            end
            Pr = P./obj.CriticalPressure;
        end
        function Tr = reducedTemperature(obj,T)
            % Compute reduced temperature
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            %
            % Returns
            % -------
            % Tr : Reduced temperature
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                T {mustBeNumeric} 
            end
            Tr = T./obj.CriticalTemperature;
        end
        function A = reducedAttractionParam(obj,Pr,Tr)
            % Compute reduced attraction parameter
            %
            % Parameters
            % ----------
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % A : Reduced attraction parameter
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                Pr (:,1) {mustBeNumeric} 
                Tr (:,1) {mustBeNumeric} 
            end
            alpha = obj.CorrectionFactor.calcAlpha(Tr);
            A = obj.OmegaA*alpha.*Pr./Tr.^2;
        end
        function B = reducedRepulsionParam(obj,Pr,Tr)
            % Compute reduced repulsion parameter
            %
            % Parameters
            % ----------
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % B : Reduced repulsion parameter
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                Pr (:,1) {mustBeNumeric} 
                Tr (:,1) {mustBeNumeric} 
            end
            B = obj.OmegaB*Pr./Tr;
        end
        function P = pressure(obj,T,V,x)
            % Compute pressure
            %
            % P = obj.PRESSURE(T,V)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % V : Molar volume [m3/mole]
            % x : Composition (Only for multi-component systems)
            %
            % Returns
            % -------
            % P : Pressure [Pa]
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                T (1,1) {mustBeNumeric}
                V (:,1) {mustBeNumeric}
                x (:,1) {mustBeNumeric} = 1
            end
            Tr = obj.reducedTemperature(T);
            alpha = obj.CorrectionFactor.calcAlpha(Tr);
            if nargin == 3
                % Pure component
                a = alpha*obj.AttractionParam;
                b = obj.RepulsionParam;
            else
                % Multi-components
                ai = alpha.*obj.AttractionParam;
                bi = obj.RepulsionParam;
                [a,b] = obj.MixingRule.apply(x,ai,bi);
            end
            P = obj.pressureImpl(T,V,a,b);
        end
        function lnPhi = lnFugacityCoeff(obj,z,params)
            % Compute fugacity coefficients
            %
            % Parameters
            % ----------
            % z      : Z-factor
            % params : struct returned by zFactors function
            %
            % Returns
            % -------
            % lnPhi : Fugacity coefficients
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                z (:,1) {mustBeNumeric}
                params struct
            end
            if isfield(params,'x')
                % Multi-component
                A = params.A;
                B = params.B;
                Ai = params.Aij*params.x;
                Bi = params.Bi;
                lnPhi = obj.lnFugacityCoeffImpl(z,A,B,Ai,Bi);
            else
                % Pure component
                A = params.A;
                B = params.B;
                lnPhi = obj.lnFugacityCoeffImpl(z,A,B);
            end
        end
        function phi = fugacityCoeff(obj,z,params)
            % Compute fugacity coefficients
            %
            % Parameters
            % ----------
            % z      : Z-factors
            % params : struct returned by zFactors function
            %
            % Returns
            % -------
            % phi : Fugacity coefficients
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                z (:,1) {mustBeNumeric}
                params struct
            end
            phi = exp(obj.lnFugacityCoeff(z,params));
        end
        function [z,params] = zFactors(obj,P,T,x)
            % Compute Z-factors
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            % x : Composition (Only for multi-component systems)
            %
            % Returns
            % -------
            % z      : Z-factors
            % params : struct
            %     P   : Pressure
            %     T   : Temperature
            %     x   : Composition
            %     A   : Mixed attraction parameter
            %     B   : Mixed repulsion parameter
            %     Ai  : i-th component attraction parameter
            %     Bi  : i-th component repulsion parameter
            %     Aij : Combined attraction parameter between i and j
            %     components
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                P (1,1) {mustBeNumeric}
                T (1,1) {mustBeNumeric}
                x (:,1) {mustBeNumeric} = 1
            end
            Pr = obj.reducedPressure(P);
            Tr = obj.reducedTemperature(T);
            if nargin == 3
                % Pure component
                A = obj.reducedAttractionParam(Pr,Tr);
                B = obj.reducedRepulsionParam(Pr,Tr);
            else
                % Multi-components
                Ai = obj.reducedAttractionParam(Pr,Tr);
                Bi = obj.reducedRepulsionParam(Pr,Tr);
                [A,B,Aij] = obj.MixingRule.apply(x,Ai,Bi);
            end
            y = roots(obj.zFactorCubicEq(A,B));
            z = sort(y(imag(y) == 0));
            params.P = P;
            params.T = T;
            params.A = A;
            params.B = B;
            if nargin > 3
                params.x = x;
                params.Ai = Ai;
                params.Bi = Bi;
                params.Aij = Aij;
            end
        end
    end
end