classdef methodsInPapers < LhCurve
    %METHODSINPAPERS contains the identification methds proposed in
    %published papers
    
    properties
        indenterElasModu
    end
    
    methods
        function obj = methodsInPapers(LhPoints, indenterR,unloadingNum,creepOrNot,indenterElasModu)            
            obj = obj@LhCurve(LhPoints, indenterR, unloadingNum,creepOrNot);
            obj.indenterElasModu = indenterElasModu;
        end
  
% below are all kinds of identification methods proposed for spherical indentation
%    in the published papers, standards, or some other materials

        matPara = GBT1(obj, contactR)   

        matPara = Cao2004(obj, hg2R)

        outputArg = SKKang2013(obj)           
        
        matPara = Haggag1993(obj)
        
        matPata = ChenXi2006(obj,relativePene)
        
        matPara = ZhangTH2009(obj)
    
        elasModu = OliverPharr(obj)   % hg2R is hg/R, where hg is the representative penetration used in Cao 2004
                                                                              % an array [hg2R1,hg2R2] should be provided here
    end
end

