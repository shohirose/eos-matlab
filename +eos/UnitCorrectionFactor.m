classdef UnitCorrectionFactor
    %UNITCORRECTIONFACTOR Correction factor euql to one.
    
    methods
        function obj = UnitCorrectionFactor()
            %UNITCORRECTIONFACTOR Create an object of this class
        end
        
        function alpha = calcAlpha(~,~)
            %CALCALPHA Compute the correction factor
            alpha = 1;
        end

        function beta = calcBeta(~,~)
            %CALCBETA Compute the logarithmic derivative of correction
            %factor
            beta = 0;
        end
    end
end

