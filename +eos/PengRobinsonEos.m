classdef PengRobinsonEos < eos.CubicEosBase
    % PengRobinsonEos Peng-Robinson equation of state
    %
    %  This class provides methods to calculate thermodynamic properties
    %  based on Peng-Robinson equation of state.
    
    properties (Constant, Access = private)
        Sqrt2 = sqrt(2)
        Delta1 = 1 + sqrt(2)
        Delta2 = 1 - sqrt(2)
    end
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
            coeffs = [1, B - 1, A - 2*B - 3*B^2, -A*B + B^2 + B^3];
        end
        function lnPhi = lnFugacityCoeffImpl(z,A,B,x,Aij,Bi)
            % Compute the natural log of fugacity coefficients.
            %
            % lnPhi = LNFUGACITYCOEFFIMPL(z,s)
            %
            % Parameters
            % ----------
            % z : Z-factor
            % A : Attraction parameter
            % B : Repulsion parameter
            % x : Composition
            % Aij : Attraction parameter between i and j components
            % Bi : Repulstion parameter of i component
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
            arguments
                z (:,1) {mustBeNumeric}
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
                x (:,1) {mustBeNumeric} = 1
                Aij (:,:) {mustBeNumeric} = 1
                Bi (:,1) {mustBeNumeric} = 1
            end
            Sqrt2 = eos.PengRobinsonEos.Sqrt2;
            Delta1 = eos.PengRobinsonEos.Delta1;
            Delta2 = eos.PengRobinsonEos.Delta2;
            if nargin > 3
                Ai = Aij*x;
                lnPhi = zeros(length(x),length(z));
                for i = 1:length(z)
                    lnPhi(:,i) = Bi/B*(z(i) - 1) - log(z(i) - B) ...
                        - A/(2*Sqrt2*B)*log((z(i) + Delta1*B)/(z(i) + Delta2*B)) ...
                        *(2*Ai/A - Bi/B);
                end
            else
                lnPhi = z - 1 - log(z - B) ...
                    - A/(2*Sqrt2*B)*log((z + Delta1*B)./(z + Delta2*B));
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
            % V : Volume [m3]
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
            P = R*T./(V - b) - a./((V - b).*(V + b) + 2*b*V);
        end
    end
    methods
        function obj = PengRobinsonEos(Pc,Tc,omega,Mw,K)
            % Construct PR EoS
            %
            % obj = PENGROBINSONEOS(Pc,Tc,omega,Mw,K)
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
            % obj : PengRobinsonEos
            arguments
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric} = 1
            end
            if nargin == 4
                K = zeros(length(Pc));
            end
            obj@eos.CubicEosBase(0.45724,0.07780,Pc,Tc,Mw,K);
            obj.AcentricFactor = omega;
        end
        function obj = setParams(obj,Pc,Tc,omega,Mw,K)
            % Set parameters.
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
            % obj : PengRobinsonEos
            arguments
                obj {mustBeA(obj,'eos.PengRobinsonEos')}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                omega (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric} = 1
            end
            if nargin == 5
                K = zeros(length(Pc));
            end
            obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw,K);
            obj.AcentricFactor = omega;
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
                obj {mustBeA(obj,'eos.PengRobinsonEos')}
                Tr (:,1) {mustBeNumeric}
            end
            omega = obj.AcentricFactor;
            m = 0.3796 + 1.485*omega - 0.1644*omega.^2 + 0.01667*omega.^3;
            alpha = (1 + m.*(1 - sqrt(Tr))).^2;
        end
    end
end