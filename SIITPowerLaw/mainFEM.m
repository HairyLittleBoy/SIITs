clear all
NaturalCoords = [0;0];
nodeCoords0 = [0.3;0.3;0.1;0.3;0.1;0.1;0.3;0.1];
dispT = [0.1,0.1;0.1,0.1;0,0;0,0];
objFEM = NonLinearStructFEM(NaturalCoords,nodeCoords0,dispT);
J0_inv = C2D4JacobiTL( objFEM );
[BL0, BL1, BNL] = strainDispMatrix(objFEM);

%%
clear all
nodeCoords0(:,1) = [1.5;1;0;1;0;0;1.5;0];
nodeCoords0(:,2) =  [3;1;1.5;1;1.5;0;3;0];
nodeCoords0(:,3) =  [3;2;1.5;2;1.5;1;3;1];
nodeCoords0(:,4) =  [1.5;2;0;2;0;1;1.5;1];
NaturalCoords = [0 0];
dispT(:,:,1) = [0  0
                0  0
                0  0
                0  0];
            
               
for i = 1 :4
    obj = NonLinearStructFEM(NaturalCoords,nodeCoords0(:,i),dispT(:,:,1));
    J0_inv(:,:,i) = C2D4JacobiTL( obj );
end