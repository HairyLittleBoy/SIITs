function matPara = ChenXi2006(obj)

% -----------------------------------------------------------------------------
% this method is from the paper "Zhao M., Ogasawara N., Chiba N., et al. 
% A new approach to measure the elastic¨Cplastic properties of bulk 
% materials using spherical indentation [J]. Acta Mater, 2006, 54(1): 23-32."
% hg2R : [h1/R, h2/R]
% the sigY and hdnIndx are for the power law : sig = sigY*(1+E/sigY*epsP)^n
% -----------------------------------------------------------------------------

h2R = [0.13, 0.3]';

[Pks, ~, Lhunloading, ~] = Lhsplit ( obj );
if Pks.High1(end,1)/obj.indenterR < h2R(2)
    warning('for the method ChenXi2006, the maximum h/R should be no smaller than 0.3')
    warning('the unloading slope and C will be extrapolated and the accuracy is not assured')
end
ft1 = fittype( 'poly1' );
unloadingSlope = zeros(length(Lhunloading),1);
for i = 1 : length(Lhunloading)
    [flt1, ~] = fit( Lhunloading{i}(:,1), Lhunloading{i}(:,3), ft1 );
    unloadingSlope(i) = flt1.p1;
end
if Pks.High1(end,1)/obj.indenterR >= h2R(2)
    S = interp1(Pks.High1(:,1)./obj.indenterR, unloadingSlope, h2R(2));
    C = interp1(Pks.High1(:,1)./obj.indenterR, Pks.High1(:,2)./Pks.High1(:,1).^2, h2R);
else
    S = interp1(Pks.High1(:,1)./obj.indenterR, unloadingSlope, h2R(2),'linear','extrap');
    ft2 = fittype( 'exp1' );
    [flt2,~] = fit( Pks.High1(:,1)./obj.indenterR, Pks.High1(:,2)./Pks.High1(:,1).^2, ft2 );
    C = flt2.a.*exp(flt2.b.*h2R);
end
CS.C = C;
CS.S = S;

sig0 = 500;

sp = [1000, 2, 0.2];
lb = [2, 0.01, 0];
ub = [3000, 10, 0.6];

options = optimset('Algorithm','interior-point','Disp','iter');
problem = createOptimProblem('fmincon','objective',@(x)ChenXi2006error(x, obj, CS),'x0',sp,...
                              'lb',lb,'ub',ub,'options',options);
ms = MultiStart;
ms.UseParallel = 'always';
[x,funcValue] = run(ms,problem,500);

E2sigY = x(1);
sigY2sig0 = x(2);
n = x(3);

matPara.sigY = sigY2sig0 * sig0;
matPara.elasModu = E2sigY * matPara.sigY;
matPara.hdnIndx = n;

end

