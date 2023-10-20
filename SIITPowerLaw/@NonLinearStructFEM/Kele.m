function [ Ke ] = Kele( obj )
%KELE this function calculate the element stiffness matrix 
%   this function is for axisymetric problem

%    input : NonLinearStructFEM object, 
%                 obj.NaturalCoords : the natural coordinates of a material point, [2 x 1] vector [r; s]

%                 obj.nodeCoords0 :  coordinates of the nodes of  the element at TIME 0, 
%                                                      [8 x 1] vector : [X11;X12;X21;X22;X31;X32;X41;X42]

%                       obj.dispT : displacements of nodes at TIME t
%                                   [4 x 2] matrix : [u11    u12
%                                                     u21    u22
%                                                     u31    u32         (disp of node 3 in first and second coords directions)
%                                                     u41    u42]

%   output : Ke, element stiffness matrix
%   the (r,s) coordinates of gauss integral points
%   (1) : sqrt(3)/3,sqrt(3)/3
%   (2) : -sqrt(3)/3,sqrt(3)/3
%   (3) : -sqrt(3)/3,-sqrt(3)/3
%   (4) : sqrt(3)/3,-sqrt(3)/3
%                       ^s
%                       |
%            -----------|-----------
%           |           |           |
%           |      *(2) |    *(1)   |
%           |           |           |
%   --------|-----------|-----------|----------->r
%           |           |           |
%           |      *(3) |    *(4)   |
%           |           |           |
%            -----------|-----------
%                       |
%                       |

[BL0, BL1, BNL, detJ0] = strainDispMatrix(obj);

ke = zeros(8,8);      
for i = 1:4
    yita = gaussIntPoints(i,1);
    kesi = gaussIntPoints(i,2);
    B = 0.5.*[-1+yita      0      1-yita      0       1+yita      0      -1-yita      0;
             0      -1+kesi      0     -1-kesi        0      1+kesi      0      1-kesi;
          -1+kesi   -1+yita   -1-kesi   1-yita    1+kesi     1+yita   1-kesi    -1-yita];
     ke = ke + B'*D*B;
end
ke = ke.*det(J);

end

