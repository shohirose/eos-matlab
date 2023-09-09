classdef SoaveRedlichKwongEos < eos.CubicEosBase
    % SoaveRedlichKwongEos Soave-Redlich-Kwong equation of state
    %
    %  This class provides methods to calculate thermodynamic properties
    %  based on Soave-Redlich-Kwong equation of state.
    
    properties (SetAccess = private)
        AcentricFactor % Acentric factor
    end
    methods (Static)
        function coeffs = zFactorCubicEq(A,B)
            % Compute coefficients of Z-factor cubic equation
            %
            % coeffs = ZFACTORCUBICEQ(A,B)
            %
            % Parameters
            % ----------
            % A : Reduced attraction parameter
            % B : Reduced repulsion parameter
            %
            % Returns
            % -------
            % coeffs : Coefficients of the cubic equation of Z-factor
            arguments
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
            end
            coeffs = [1, -1, A - B - B^2, -A*B];
        end
        function lnPhi = lnFugacityCoeffImpl(z,A,B,Ai,Bi)
            % Compute natural log of fugacity coefficients
            %
            % lnPhi = LNFUGACITYCOEFFIMPL(z,A,B,x,Aij,Bi)
            %
            % Parameters
            % ----------
            % z : Z-factor
            % A : Attraction parameter
            % B : Repulsion parameter
            % x : Composition
            % Ai : = Aij*xj
            % Bi : Repulstion parameter of i component
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
            arguments
                z (:,1) {mustBeNumeric}
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
                Ai (:,:) {mustBeNumeric} = 1
                Bi (:,1) {mustBeNumeric} = 1
            end
            Q = A/B*log(B./z + 1);
            if nargin > 3
                Ak = repmat(Ai/A,1,length(z));
                Bk = repmat(Bi/B,1,length(z));
                Qk = repmat(Q',length(Ai),1);
                zk = repmat(z',length(Ai),1);
                lnPhi = Bk.*(zk - 1) - log(zk - B) - Qk.*(2*Ak - Bk);
            else
                lnPhi = z - 1 - log(z - B) - Q;
            end
        end
        function P = pressureImpl(T,V,a,b)
            % Compute pressure.
            %
            % P = PRESSUREIMPL(T,V,a,b)
            %
            % Parameters
            % ----------
            % T : Temperature [K]
            % V : Molar volume [m3/mole]
            % a : Attraction parameter
            % b : Repulsion parameter
            %
            % Returns
            % -------
            % P : Pressure [Pa]
            arguments
                T {mustBeNumeric}
                V {mustBeNumeric}
                a {mustBeNumeric}
                b {mustBeNumeric}
            end
            R = eos.ThermodynamicConstants.Gas;
            P = R*T./(V - b) - a./(V.*(V + b));
        end
    end
    methods
        function obj = SoaveRedlichKwongEos(Pc,Tc,omega,Mw,K)
            % Construct SRK EoS
            %
            % obj = SOAVEREDLICHKWONGEOS(Pc,Tc,omega,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            % K : Binary interaction parameters
            %
            % Returns
            % -------
            % obj : SOAVEREDLICHKWONGEOS
            arguments
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric} = 1
            end
            if nargin > 4
                args = {0.42747,0.08664,Pc,Tc,Mw,K};
            else
                args = {0.42747,0.08664,Pc,Tc,Mw};
            end
            obj@eos.CubicEosBase(args{:});
            obj.AcentricFactor = omega;
        end
        function obj = setParams(obj,Pc,Tc,omega,Mw,K)
            % Set parameters
            %
            % obj = obj.SETPARAMS(Pc,Tc,omega,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % omega : Acentric factor
            % Mw : Molecular weight [g/mol]
            % K : Binary interaction parameters
            %
            % Returns
            % -------
            % obj : SOAVEREDLICHKWONGEOS
            arguments
                obj {mustBeA(obj,'eos.SoaveRedlichKwongEos')}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric} = 1
            end
            if nargin > 5
                obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw,K);
            else
                obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw);
            end
            obj.AcentricFactor = omega;
        end
        function m = calcm(obj)
            % Calc m using the correlation of Soave(1972).
            %
            % m = obj.CALCM();
            %
            % Returns
            % -------
            % m : Parameter m
            omega = obj.AcentricFactor;
            m = 0.48 + 1.574*omega - 0.176*omega.^2;
        end
        function alpha = temperatureCorrectionFactor(obj,Tr)
            % Compute temperature correction factor.
            %
            % alpha = obj.TEMPERATURECORRECTIONFACTOR(Tr)
            %
            % Parameters
            % ----------
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % alpha : Temperature correction factor
            arguments
                obj {mustBeA(obj, 'eos.SoaveRedlichKwongEos')}
                Tr (:,1) {mustBeNumeric}
            end
            m = obj.calcm();
            alpha = (1 + m.*(1 - sqrt(Tr))).^2;
        end
    end
end