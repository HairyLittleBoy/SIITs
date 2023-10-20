classdef LhCurve
    %LHCURVE 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties 
        LhPoints
        indenterR
        unloadingNum         % how many unloading between start and 
                                             % final unloading ( final unloading not included)
        creepOrNot           % a string indicating the Lh curve contains creep part or not
                             % 'creep' of 'nocreep'
    end

    methods
        function obj = LhCurve(LhP, indenterR, unloadingNum, creepOrNot)
            %--------------------------------------------------------------------------
            %   LhPoints is n x 2 matrix, the 1st colume is penetration in
            %   mm and the 2nd is load in N
            %--------------------------------------------------------------------------
            obj.LhPoints = LhP;
            obj.indenterR = indenterR;
            obj.unloadingNum = unloadingNum;
            obj.creepOrNot = creepOrNot;
        end
        
        [Pks, Lhloading, Lhunloading, creep] = Lhsplit ( obj )
         LhloadingWhole = multiCombine ( obj )
    end
end

