function matPara = ZhangTH2009(obj)

% -------------------------------------------------------------------------
% this method is from 姜鹏,张泰华,杨荣,梁乃刚.基于球形压入法提取材料的塑性
% 力学参数[J].力学学报,2009,41(05):730-738.
% Pneg J., Taihua Z., Yang R., et al. A New Spherical Indentation-based 
% Method to Extract Plastic Material Parameters [J]. Chinese Journal of 
% Theoretical and Applied Mechanics, 2009, 41(05): 730-738.
% -------------------------------------------------------------------------

% the stress-strain relation applied in this method:
%                   sig = E*eps                              eps < epsY
%                   sig = K*eps^n = E*epsy^(1-n) * eps^n     eps > epsY 

h2R = [0.05;0.1];
[Pks, Lhloading, Lhunloading, ~] = Lhsplit ( obj );
Wu = zeros(length(Lhloading), 1);
Wt = zeros(length(Lhloading), 1);
for i = 1 : length(Lhloading)
    Wu(i) = -trapz(Lhunloading{i}(:,1),Lhunloading{i}(:,3));
    if i == 1
        Wt(i) = trapz(Lhloading{i}(:,1),Lhloading{i}(:,3));
    else
        Wt(i) = trapz(Lhloading{i}(:,1),Lhloading{i}(:,3)) + Wt(i-1);
    end
end
ft = fittype( 'power2' );
[flt1, ~] = fit( Pks.High1(:,1)./obj.indenterR, Wt, ft );
[flt2, ~] = fit( Pks.High1(:,1)./obj.indenterR, Wu, ft );
Wth2R = flt1.a * h2R.^flt1.b + flt1.c;
Wuh2R = flt2.a * h2R.^flt2.b + flt2.c;
W = Wuh2R ./ Wth2R;

sp = [0.005, 0.2];
lb = [0, 0.005];
ub = [0.01, 0.4];

options = optimset('Algorithm','interior-point','Disp','iter');
problem = createOptimProblem('fmincon','objective',@(x)ZhangTH2009error(x,W,obj),'x0',sp,...
                              'lb',lb,'ub',ub,'options',options);
ms = MultiStart;
ms.UseParallel = 'always';
[x,funcValue] = run(ms,problem,200);
matPara.epsY = x(1);
matPara.hdnIndx = x(2);

% evaluate the elastic mudulus according to the Eq(24) Eq(25) in the paper
LhloadingWhole = multiCombine ( obj );
pene01 = h2R(2)*obj.indenterR;
if LhloadingWhole(end,1) > pene01
    Load01 = interp1(LhloadingWhole(:,1),LhloadingWhole(:,3),pene01);     % the load when h=0.1R
else
    ft = fittype( 'poly3' );
    [flt3, ~] = fit( LhloadingWhole(:,1), LhloadingWhole(:,3), ft );
    Load01 = flt3.p1*pene01^3 + flt3.p2*pene01^2 + flt3.p3*pene01 + flt3.p4;    % the load when h=0.1R, extrp may lead to large discrepency
                                                                                % a good method should be applied here to evaluate the Lh data
end
epsrc = 0.0418;       % the epsrc at h/R = 0.1;
temp1 = matPara.epsY^(1-matPara.hdnIndx)*(matPara.epsY + epsrc)^matPara.hdnIndx;    % Eq(25)
fs1 = @(E)-1.356*(log(E/((1-0.3^2)*(E*temp1))))^3+11.5*(log(E/((1-0.3^2)*(E*temp1))))^2+ ...
           25.055*log(E/((1-0.3^2)*(E*temp1)))-64.419 - Load01/(E*temp1*(h2R(2)*obj.indenterR)^2);
[matPara.elasModu, funcValue] = fsolve(fs1,300000);
matPara.sigY = matPara.elasModu * matPara.epsY;

end




