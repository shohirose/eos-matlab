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
        function lnPhi = lnFugacityCoeffImpl(z,A,B,Ai,Bi)
            % Compute the natural log of fugacity coefficients.
            %
            % lnPhi = LNFUGACITYCOEFFIMPL(z,s)
            %
            % Parameters
            % ----------
            % z : Z-factor
            % A : Attraction parameter
            % B : Repulsion parameter
            % Ai : Ai = Aij*xj
            % Bi : Repulstion parameter of i component
            %
            % Returns
            % -------
            % lnPhi : Natural log of fugacity coefficients
            arguments
                z (:,1) {mustBeNumeric}
                A (1,1) {mustBeNumeric}
                B (1,1) {mustBeNumeric}
                Ai (:,1) {mustBeNumeric} = 1
                Bi (:,1) {mustBeNumeric} = 1
            end
            Sqrt2 = eos.PengRobinsonEos.Sqrt2;
            Delta1 = eos.PengRobinsonEos.Delta1;
            Delta2 = eos.PengRobinsonEos.Delta2;
            Q = A/(2*Sqrt2*B)*log((z + Delta1*B)./(z + Delta2*B));
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
            P = R*T./(V - b) - a./((V - b).*(V + b) + 2*b*V);
        end
        
        function m = calcM(omega)
            % Compute parameter m.
            %
            % m = obj.CALCM()
            %
            % Returns
            % ----------
            % m : Parameter m
            arguments
                omega (:,1) {mustBeNumeric}
            end
            m = 0.37464 + 1.54226*omega - 0.26992*omega.^2;
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
            % K : Binary interaction parameters (optional)
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
            m = eos.PengRobinsonEos.calcM(omega);
            alpha = eos.SoaveCorrectionFactor(m);
            if nargin > 4
                args = {0.45724,0.07780,Pc,Tc,Mw,omega,alpha,K};
            else
                args = {0.45724,0.07780,Pc,Tc,Mw,omega,alpha};
            end
            obj@eos.CubicEosBase(args{:});
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
            m = eos.PengRobinsonEos.calcM(omega);
            alpha = eos.SoaveCorrectionFactor(m);
            if nargin > 5
                obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw,omega,alpha,K);
            else
                obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw,omega,alpha);
            end
        end
    end
end