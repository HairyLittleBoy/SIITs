function [J0_inv,detJ0] = C2D4JacobiTL( obj )
%C3D4JACOBITL generate the element Jacobi matrix 
%   this function generates the element Jacobi matrix for 2D 4 node 
%   isoparametric element under total lagrange formulation

%    input : NonLinearStructFEM object, 
%                  where obj.nodeCoords0 is the coordinates of the nodes of
%                  the element at TIME 0, [8 x 1] vector : [X11;X12;X21;X22;X31;X32;X41;X42]
%                  for instance, the X31 present the X1 coordinate of node
%                  3 and X32 presents the X2 coordinate of node 3
%   output : J0_inv, inverse of the Jacobian matrix for continium 2 dimensional 4 node
%                  axisymmetric element

% the interpolation function of the element : 
%  h1(r,s) = 1/4(1+r)(1+s)        h2(r,s) = 1/4(1-r)(1+s)
%  h1(r,s) = 1/4(1-r)(1-s)          h4(r,s) = 1/4(1+r)(1-s)

% element Jacobian Matrix = [p0X1pr   p0X2pr
%                                                     p0X1ps   p0X2ps]
% where tX1 = h1*tX11 + h2*tX21 + h3*tX31 + h4*tX41
% where tX2 = h1*tX12 + h2*tX22 + h3*tX32 + h4*tX42

% interpFuncMatrix = 0.25.*[1 0 -1 0 -1 0 1 0; 
%                           1 0  1 0 -1 0 -1 0; 
%                           0 1 0 -1 0 -1 0 1; 
%                           0 1 0 1 0 -1 0 -1];
interpFuncMatrix = [0.25    0    -0.25    0    -0.25    0    0.25      0; 
                    0.25    0     0.25    0    -0.25    0    -0.25     0; 
                      0    0.25     0   -0.25     0   -0.25     0    0.25; 
                      0    0.25     0    0.25     0   -0.25     0    -0.25];                      
J0_temp = interpFuncMatrix * obj.nodeCoords0;
% J0(1,1) = J0_temp(1,1);
% J0(2,1) = J0_temp(2,1);
% J0(1,2) = J0_temp(3,1);
% J0(2,2) = J0_temp(4,1);
J0=reshape(J0_temp,[2,2]);
J0_inv = inv(J0);
detJ0 = det(J0);
end

