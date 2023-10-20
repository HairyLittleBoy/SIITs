function GL_0tepsilon = GreenLagrangeStrain( obj )
%GreenLagrangeStrain function generates the green-lagrange strain through
%deformation gradient matrix
%   Green-Lagrange strain (GL_0tepsilon) is calculated through Cauchy-Green
%   deformation tensor (CG_0tC)

%          obj.NaturalCoords : the natural coordinates of a material point, [2 x 1] vector [r; s]

%            obj.nodeCoords0 :  coordinates of the nodes of  the element at TIME 0, 
%                                                [8 x 1] vector : [X11;X12;X21;X22;X31;X32;X41;X42]

%                  obj.dispT : displacements of nodes at TIME t
%                                   [4 x 2] matrix : [u11    u12
%                                                     u21    u22
%                                                     u31    u32   (disp of node 3 in first and second coords directions)
%                                                     u41    u42]

DefGrad_0tX = DefGrad( obj );

CG_0tC = DefGrad_0tX' * DefGrad_0tX;
GL_0tepsilon_temp = 0.5.*(CG_0tC - eye(2));

h1 = 0.25*(1+obj.NaturalCoords(1,1))*(1+obj.NaturalCoords(2,1));
h2 = 0.25*(1-obj.NaturalCoords(1,1))*(1+obj.NaturalCoords(2,1));
h3 = 0.25*(1-obj.NaturalCoords(1,1))*(1-obj.NaturalCoords(2,1));
h4 = 0.25*(1+obj.NaturalCoords(1,1))*(1-obj.NaturalCoords(2,1));
X1 = [h1  h2  h3  h4]*[obj.nodeCoords0(1,1);  obj.nodeCoords0(3,1);  obj.nodeCoords0(5,1);  obj.nodeCoords0(7,1)];
tu1 = [h1  h2  h3  h4] * obj.dispT(:,1);
GL_0tepsilon33 = tu1/X1 + 0.5*(tu1/X1)^2;

GL_0tepsilon = [GL_0tepsilon_temp(1,1),GL_0tepsilon_temp(2,2),GL_0tepsilon_temp(1,2),GL_0tepsilon33];

end

