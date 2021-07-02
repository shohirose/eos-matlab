classdef CubicEosBase
    % Base class for two-parameter cubic equations of state.
    
    properties (SetAccess = private)
        CriticalPressure    % Critical pressure [Pa]
        CriticalTemperature % Critical temperature [K]
        MolecularWeight     % Molecular weight [g/mol]
        OmegaA              % Coefficient for attraction parameter
        OmegaB              % Coefficient for repulsion parameter
        AttractionParam     % Attraction parameter
        RepulsionParam      % Repulsion parameter
        BinaryInteractionParams % Binary interaction parameters
    end
    properties (Dependent)
        NumComponents % Number of components
    end
    methods
        function obj = CubicEosBase(OmegaA,OmegaB,Pc,Tc,Mw,K)
            % Construct cubic EOS
            %
            % Parameters
            % ----------
            % OmegaA : Coefficient for attraction parameter
            % OmegaB : Coefficient for repulsion parameter
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            % K  : Binary interaction parameter
            arguments
                OmegaA (1,1) {mustBeNumeric}
                OmegaB (1,1) {mustBeNumeric}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric}
            end
            obj.OmegaA = OmegaA;
            obj.OmegaB = OmegaB;
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            obj.BinaryInteractionParams = K;
            R = eos.ThermodynamicConstants.Gas;
            obj.AttractionParam = OmegaA*(R*Tc).^2./Pc;
            obj.RepulsionParam = OmegaB*R*Tc./Pc;
        end
        function ncomp = get.NumComponents(obj)
            ncomp = len(obj.CriticalPressure);
        end
        function obj = setParams(obj,Pc,Tc,Mw,K)
            % Set parameters
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            % K  : Binary interaction parameters
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric} = 1
            end
            obj.CriticalPressure = Pc;
            obj.CriticalTemperature = Tc;
            obj.MolecularWeight = Mw;
            obj.BinaryInteractionParams = K;
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
        function A = reducedAttractionParam(obj,Pr,Tr,alpha)
            % Compute reduced attraction parameter
            %
            % Parameters
            % ----------
            % Pr : Reduced pressure
            % Tr : Reduced temperature
            % alpha : Temperature correction factor
            %
            % Returns
            % -------
            % A : Reduced attraction parameter
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                Pr {mustBeNumeric} 
                Tr {mustBeNumeric} 
                alpha {mustBeNumeric} 
            end
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
                Pr {mustBeNumeric} 
                Tr {mustBeNumeric} 
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
            % V : Volume [m3]
            % x : Composition
            %
            % Returns
            % -------
            % P : Pressure [Pa]
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                T (1,1) {mustBeNumeric}
                V (1,1) {mustBeNumeric}
                x (:,1) {mustBeNumeric} = 1
            end
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            if nargin ~= 3
                a = alpha.*obj.AttractionParam;
                b = obj.RepulsionParam;
            else
                ai = alpha.*obj.AttractionParam;
                bi = obj.RepulsionParam;
                [a,b] = obj.applyMixingRule(x,ai,bi);
            end
            P = obj.pressureImpl(T,V,a,b);
        end
        function lnPhi = lnFugacityCoeff(obj,z,s)
            % Compute fugacity coefficients
            %
            % Parameters
            % ----------
            % z : Z-factors
            % s : struct containing parameters
            %
            % Returns
            % -------
            % lnPhi : Fugacity coefficients
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                z (:,1) {mustBeNumeric}
                s struct
            end
            if isfield(s,'x')
                lnPhi = obj.lnFugacityCoeffImpl(z,s.A,s.B,s.x,s.Aij,s.Bi);
            else
                lnPhi = obj.lnFugacityCoeffImpl(z,s.A,s.B);
            end
        end
        function phi = fugacityCoeff(obj,z,s)
            % Compute fugacity coefficients
            %
            % Parameters
            % ----------
            % z : Z-factors
            % s : struct containing parameters
            %
            % Returns
            % -------
            % phi : Fugacity coefficients
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                z (:,1) {mustBeNumeric}
                s struct
            end
            phi = exp(obj.lnFugacityCoeff(z,s));
        end
        function [z,s] = zFactors(obj,P,T,x)
            % Compute Z-factors
            %
            % Parameters
            % ----------
            % P : Pressure [Pa]
            % T : Temperature [K]
            % x : Composition
            %
            % Returns
            % -------
            % z : Z-factors
            % s : struct containing parameters
            arguments
                obj {mustBeA(obj,'eos.CubicEosBase')}
                P (1,1) {mustBeNumeric}
                T (1,1) {mustBeNumeric}
                x (:,1) {mustBeNumeric} = 1
            end
            Pr = obj.reducedPressure(P);
            Tr = obj.reducedTemperature(T);
            alpha = obj.temperatureCorrectionFactor(Tr);
            if nargin == 3
                A = obj.reducedAttractionParam(Pr,Tr,alpha);
                B = obj.reducedRepulsionParam(Pr,Tr);
            else
                Ai = obj.reducedAttractionParam(Pr,Tr,alpha);
                Bi = obj.reducedRepulsionParam(Pr,Tr);
                [A,B,Aij] = obj.applyMixingRule(x,Ai,Bi);
            end
            y = roots(obj.zFactorCubicEq(A,B));
            z = y(imag(y) == 0);
            if nargout > 1
                s.P = P;
                s.T = T;
                s.A = A;
                s.B = B;
                if nargin > 3
                    s.x = x;
                    s.Ai = Ai;
                    s.Bi = Bi;
                    s.Aij = Aij;
                end
            end
        end
        function [a,b,aij] = applyMixingRule(obj,x,ai,bi)
            % Apply mixing rule to attraction and repulsion parameters.
            %
            % [a,b,aij] = obj.APPLYMIXINGRULE(x,ai,bi)
            %
            % Parameters
            % ----------
            % x  : Composition
            % ai : Attraction parameter of each component
            % bi : Repulsion parameter of each component
            %
            % Returns
            % -------
            % a   : Attraction parameter of the mixture
            % b   : Repulsion parameter of the mixture
            % aij : Combined attraction parameter between i and j
            % components
            K = obj.BinaryInteractionParams;
            % Combining rule with correction parameters
            aij = (1 - K).*sqrt(kron(ai,ai'));
            % Quadratic mixing
            a = x'*aij*x;
            % Linear mixing
            b = x'*bi;
        end
    end
end