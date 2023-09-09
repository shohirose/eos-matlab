classdef MixingRule
    %MIXINGRULE Mixing rule for cubic EoS
    %   
    
    properties
        BinaryInteractionParams % Binary interaction parameter
    end
    
    methods
        function obj = MixingRule(K)
            %MIXINGRULE Construct mixing rule
            %
            % Parameters
            % ----------
            % K : Binary interaction parameters
            %
            arguments
                K (:,:) {mustBeNumeric}
            end
            obj.BinaryInteractionParams = K;
        end

        function K = get.BinaryInteractionParams(obj)
            K = obj.BinaryInteractionParams;
        end

        function obj = set.BinaryInteractionParams(obj,K)
            arguments
                obj
                K (:,:) {mustBeNumeric}
            end
            obj.BinaryInteractionParams = K;
        end
        
        function [A,B,Aij] = apply(obj,x,Ai,Bi)
            %APPLY Apply mixing rule
            %
            % Paramters
            % ---------
            % x  : Phase composition
            % Ai : Reduced attraction parameters of each component
            % Bi : Reduced repulsion parameters of each component
            %
            % Returns
            % -------
            % A : Reduced attraction parameters of mixture
            % B : Reduced repulsion parameters of mixture
            % Aij : Average of reduced attraction parameters of i and j
            % components
            arguments
                obj {mustBeA(obj, 'eos.MixingRule')}
                x  (:,1) {mustBeNumeric}
                Ai (:,1) {mustBeNumeric}
                Bi (:,1) {mustBeNumeric}
            end
            K = obj.BinaryInteractionParams;
            Aij = (1 - K).*sqrt(kron(Ai,Ai'));
            A = x'*Aij*x;
            B = x'*Bi;
        end
    end
end

