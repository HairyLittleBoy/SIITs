classdef sphrInden < abaqusSimu
    %sphrInden generate Lh curve with interpolation, ANN and FEM methods
    %   
    
    properties 
        LhGenrMethod
    end
    
    methods
        function obj = sphrInden(mlPs,niu,indenterR,penetration,miu,cpuNum,meshSizeCoef,method)
            obj = obj@abaqusSimu(mlPs,niu,indenterR,penetration,miu,cpuNum,meshSizeCoef);
            obj.LhGenrMethod = method;
        end

% ----- LhPrdctWthAnn       
        [ LhAnn ] = LhCalcAnn(obj)

% ----- LhCalcNumr
        LhNumr = LhCalcNumr( obj )
       
% ----- LhGenr
        Lhdata = LhGenr(obj)
    end
end

