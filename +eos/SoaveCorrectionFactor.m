classdef SoaveCorrectionFactor
    %SOAVECORRECTIONFACTOR Correction factor proposed by Soave (1972).
    
    properties
        M
    end
    
    methods
        function obj = SoaveCorrectionFactor(m)
            %SOAVECORRECTIONFACTOR Create an object of this class.
            %
            % Parameters
            % ----------
            % m : Coefficient calculated from acentric factor
            obj.M = m;
        end

        function m = get.M(obj)
            m = obj.M;
        end

        function obj = set.M(obj,m)
            obj.M = m;
        end
        
        function alpha = calcAlpha(obj,Tr)
            %CALCALPHA Compute the correction factor
            %
            % Parameters
            % ----------
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % alpha : Temperature correction factor for the attraction
            % parameter
            m = obj.M;
            alpha = (1 + m.*(1 - sqrt(Tr))).^2;
        end

        function beta = calcBeta(obj,Tr)
            %CALCBETA Compute the logarithmic derivative of the correction
            % factor
            %
            % Parameters
            % ----------
            % Tr : Reduced temperature
            %
            % Returns
            % -------
            % beta : The logarithmic derivative of temperature correction
            % factor for the attraction parameter
            m = obj.M;
            alpha = (1 + m.*(1 - sqrt(Tr))).^2;
            beta = -m*sqrt(Tr/alpha);
        end
    end
end

