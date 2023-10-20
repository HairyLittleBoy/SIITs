function DefGrad_0tX = DefGrad( obj )
%DefGrad function generates the deformation gradient matrix
%    input : NonLinearStructFEM object, 
%                 obj.NaturalCoords : the natural coordinates of a material point, [2 x 1] vector [r; s]

%                 obj.nodeCoords0 :  coordinates of the nodes of  the element at TIME 0, 
%                                                      [8 x 1] vector : [X11;X12;X21;X22;X31;X32;X41;X42]

%                       obj.dispT : displacements of nodes at TIME t
%                                   [4 x 2] matrix : [u11    u12
%                                                     u21    u22
%                                                     u31    u32         (disp of node 3 in first and second coords directions)
%                                                     u41    u42]
%
% phpx is a [2 x4] matrix :  [ph1pr, ph2pr,ph3pr, ph4pr;
%                             ph1ps, ph2ps,ph3ps, ph4ps;]

% php0X is a [2 x4] matrix : [ph1p0X1, ph2p0X1,ph3p0X1, ph4p0X1;
%                             ph1p0X2, ph2p0X2,ph3p0X2, ph4p0X2;]

[J0_inv,~] = C2D4JacobiTL( obj );
phpx = 0.25.*[(1+obj.NaturalCoords(1,1))    -(1+obj.NaturalCoords(1,1))    -(1-obj.NaturalCoords(1,1))   (1-obj.NaturalCoords(1,1));
              (1+obj.NaturalCoords(2,1))     (1-obj.NaturalCoords(2,1))    -(1-obj.NaturalCoords(2,1))  -(1+obj.NaturalCoords(2,1))];
           
php0X = J0_inv * phpx;
DefGrad_0tX = php0X * obj.dispT + eye(2);
                               
end

