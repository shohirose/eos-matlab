classdef VanDerWaalsEos < eos.CubicEosBase
    % Van der Waals equation of state.
    
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
            coeffs = [1, -B - 1, A, -A*B];
        end
        function lnPhi = lnFugacityCoeffImpl(z,A,B,Ai,Bi)
            % Compute the natural log of fugacity coeffcients
            %
            % lnPhi = LNGUGACITYCOEFF(z,s)
            %
            % Parameters
            % ----------
            % z : Z-factors
            % A : Attraction parameter
            % B : Repulsion parameter
            % Ai : = Aij*xj
            % Bi : Repulsion parameter for i component
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
            if nargin > 3
                % Multi-component systems
                Ak = repmat(Ai/A,1,length(z));
                Bk = repmat(Bi/B,1,length(z));
                zk = repmat(z',length(Ai),1);
                lnPhi = Bk.*(zk - 1) - log(zk - B) - A./zk.*(2*Ak - Bk);
            else
                % Pure component systems
                lnPhi = z - 1 - log(z - B) - A./z;
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
            P = R*T./(V - b) - a./V.^2;
        end
    end
    methods
        function obj = VanDerWaalsEos(Pc,Tc,Mw,K)
            % Construct VDW EoS
            %
            % obj = VANDERWAALSEOS(Pc,Tc,Mw)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            % K  : Binary interaction parameters (optional)
            %
            % Returns
            % -------
            % obj : VANDERWAALSEOS
            arguments
               Pc (:,1) {mustBeNumeric}
               Tc (:,1) {mustBeNumeric}
               Mw (:,1) {mustBeNumeric}
               K (:,:) {mustBeNumeric} = 1
            end
            omega = 0;
            alpha = eos.UnitCorrectionFactor();
            if nargin > 3
                args = {0.421875,0.125,Pc,Tc,Mw,omega,alpha,K};
            else
                args = {0.421875,0.125,Pc,Tc,Mw,omega,alpha};
            end
            obj@eos.CubicEosBase(args{:});
        end
        function obj = setParams(obj,Pc,Tc,Mw,K)
            % Set parameters
            %
            % obj = obj.SETPARAMS(Pc,Tc,Mw,K)
            %
            % Parameters
            % ----------
            % Pc : Critical pressure [Pa]
            % Tc : Critical temperature [K]
            % Mw : Molecular weight [g/mol]
            % K  : Binary interaction parameters
            %
            % Returns
            % -------
            % obj : VANDERWAALSEOS
            arguments
                obj {mustBeA(obj,'eos.VanDerWaalsEos')}
                Pc (:,1) {mustBeNumeric}
                Tc (:,1) {mustBeNumeric}
                Mw (:,1) {mustBeNumeric}
                K (:,:) {mustBeNumeric} = 1
            end
            omega = 0;
            alpha = eos.UnitCorrectionFactor();
            if nargin > 4
                obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw,omega,alpha,K);
            else
                obj = setParams@eos.CubicEosBase(obj,Pc,Tc,Mw,omega,alpha);
            end
        end
    end
end