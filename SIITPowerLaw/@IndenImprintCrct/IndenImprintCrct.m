classdef IndenImprintCrct < LhCurve
    %INDEN Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        elasticModulus
        yieldStress
        
        % the euqation applied to describe stress strain relation in most 
        % studies on uniaxial propeties identification with spherical 
        % indentation is the power law equation based on stress and full
        % true strain instead of plastic strain, as a result, the hardening
        % index here is the index in power law stress-strain (not plastic strain)
        % relation.
        hardeningIndex
    end
    
    methods
        function obj = IndenImprintCrct(LhPoints, indenterR,unloadingNum, creepOrNot, elasticModulus, ...
                                                                      yieldStress, hardeningIndex)
            
            obj = obj@LhCurve(LhPoints, indenterR, unloadingNum,creepOrNot);
            obj.elasticModulus = elasticModulus;
            obj.yieldStress = yieldStress;
            obj.hardeningIndex = hardeningIndex;
            
        end
        
        [truePene, trueContactRadius, trueCntctArea] = Kim2006(obj)
        
        [truePene, trueContactRadius, trueCntctArea] = Hill1989(obj)

        [truePene, trueContactRadius, trueCntctArea] =  Matthews1980(obj)    
        
    end
end

