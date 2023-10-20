function [BL0, BL1, BNL, detJ0] = strainDispMatrix(obj)
%strainDispMatrix function generates the strain-displacement matrix
%   this function generates the  strain displacement matrix for 2D 4 node 
%   isoparametric element under total lagrange formulation

%    input : NonLinearStructFEM object, 
%                 obj.NaturalCoords : the natural coordinates of a material point, [2 x 1] vector [r; s]

%                 obj.nodeCoords0 :  coordinates of the nodes of  the element at TIME 0, 
%                                                      [8 x 1] vector : [X11;X12;X21;X22;X31;X32;X41;X42]

%                       obj.dispT : displacements of nodes at TIME t
%                                   [4 x 2] matrix : [u11    u12
%                                                     u21    u22
%                                                     u31    u32         (disp of node 3 in first and second coords directions)
%                                                     u41    u42]

%   output : BL0, linear strain-displacement matrix withiout initial
%                            displacement effect
%            BL1, linear strain-displacement matrix with initial
%                             displacement effect
%            BNL, nonlinear strian-displacment matrix

% the interpolation function of the element : 
%  h1(r,s) = 1/4(1+r)(1+s)        h2(r,s) = 1/4(1-r)(1+s)
%  h3(r,s) = 1/4(1-r)(1-s)        h4(r,s) = 1/4(1+r)(1-s)



% phpx is a [2 x4] matrix :  [ph1pr, ph2pr,ph3pr, ph4pr;
%                             ph1ps, ph2ps,ph3ps, ph4ps;]

% php0X is a [2 x4] matrix : [ph1p0X1, ph2p0X1,ph3p0X1, ph4p0X1;
%                             ph1p0X2, ph2p0X2,ph3p0X2, ph4p0X2;]


h1 = 0.25*(1+obj.NaturalCoords(1,1))*(1+obj.NaturalCoords(2,1));
h2 = 0.25*(1-obj.NaturalCoords(1,1))*(1+obj.NaturalCoords(2,1));
h3 = 0.25*(1-obj.NaturalCoords(1,1))*(1-obj.NaturalCoords(2,1));
h4 = 0.25*(1+obj.NaturalCoords(1,1))*(1-obj.NaturalCoords(2,1));

[J0_inv,detJ0] = C2D4JacobiTL( obj );
phpx = 0.25.*[(1+obj.NaturalCoords(1,1))    -(1+obj.NaturalCoords(1,1))    -(1-obj.NaturalCoords(1,1))   (1-obj.NaturalCoords(1,1));
              (1+obj.NaturalCoords(2,1))     (1-obj.NaturalCoords(2,1))    -(1-obj.NaturalCoords(2,1))  -(1+obj.NaturalCoords(2,1))];
           
 php0X = J0_inv * phpx;
 X1 = [h1  h2  h3  h4]*[obj.nodeCoords0(1,1);  obj.nodeCoords0(3,1);  obj.nodeCoords0(5,1);  obj.nodeCoords0(7,1)];
 
 BL0 = [php0X(1,1)            0          php0X(1,2)          0           php0X(1,3)         0          php0X(1,4)           0;
            0             php0X(2,1)          0          php0X(2,2)          0          php0X(2,3)         0             php0X(2,4);
        php0X(2,1)        php0X(1,1)     php0X(2,2)      php0X(1,2)      php0X(2,3)     php0X(1,3)     php0X(2,4)        php0X(1,4);
          h1/X1               0            h2/X1             0              h3/X1           0             h4/X1              0];
             
phuMatrix =   php0X *  obj.dispT;

BL1_1 = [phuMatrix(1,1)        0             phuMatrix(1,2)         0;
               0         phuMatrix(2,1)             0         phuMatrix(2,2)];
                  
BL1_2 = [php0X(1,1)           0            php0X(2,1)           0            php0X(1,3)          0           php0X(1,4)          0;
         php0X(1,2)           0            php0X(2,2)           0            php0X(2,3)          0           php0X(2,4)          0;
             0            php0X(1,1)           0           php0X(2,1)             0          php0X(1,3)          0            php0X(1,4);
             0            php0X(1,2)           0           php0X(2,2)             0          php0X(2,3)          0            php0X(2,4)];     
                      
 BL1_row12 =  BL1_1 * BL1_2;
 BL1_3 = [phuMatrix(1,1)      phuMatrix(1,2)      phuMatrix(2,1)      phuMatrix(2,2)];
 
 BL1_4 = [    BL1_2(2,:);
              BL1_2(4,:);
              BL1_2(1,:);
              BL1_2(3,:)    ];
                   
BL1_row3 = BL1_3 * BL1_4;

hu = [h1  h2  h3  h4] * obj.dispT(:,1);

BL1 = [                     BL1_row12;
                             BL1_row3;
       hu*h1/X1^2  0  hu*h2/X1^2   0  hu*h3/X1^2   0   hu*h4/X1^2   0];
BNL = [  BL1_2;
        BL0(4,:)];
end

