classdef NonLinearStructFEM
    %NonLinearStructFEM object contains all information of an element 
    %   
    
    properties
        NaturalCoords %   the natural coordinates of a material point, [2 x 1] vector [r; s]
        nodeCoords0   %   coordinates of the nodes of  the element at TIME 0, 
%                         [8 x 1] vector : [X11;X12;X21;X22;X31;X32;X41;X42]
        dispT         %   displacements of nodes at TIME t
%                         [4 x 2] matrix : [u11    u12
%                                           u21    u22
%                                           u31    u32    (disp of node 3 in first and second coords directions)
%                                           u41    u42]
    end
    
    methods

        function obj = NonLinearStructFEM(NaturalCoords,nodeCoords0,dispT)
            %--------------------------------------------------------------------------
            %   LhPoints is n x 2 matrix, the 1st colume is penetration in
            %   mm and the 2nd is load in N
            %--------------------------------------------------------------------------
            obj.NaturalCoords = NaturalCoords;
            obj.nodeCoords0 = nodeCoords0;
            obj.dispT = dispT;
        end
        
        J0_inv = C2D4JacobiTL( obj )
        
        [BL0, BL1, BNL] = strainDispMatrix(obj)
    end
    
end

