function matPara = Cao2004(obj, hg2R)
% -------------------------------------------------------------------------
% Cao Y P, Lu J. A new method to extract the plastic properties of metal
% materials from an instrumented spherical indentation loading curve[J].
% Acta Materialia. 2004, 52(13): 4023-4032.
% prompt = ' input the elastic modulus of the tested material : ';
% matPara.elasModu = input(prompt);
% if isempty(matPara.elasModu)
%     warning( 'an Elastic modulus will be evaluated based on the Oliver-Pharr method');
% the sigY and hdnIndx are for the power law : sig = sigY*(1+E/sigY*epsP)^n

% -------------------------------------------------------------------------
matPara.elasModu = OliverPharr(obj);
% end
hg2R = linspace(hg2R(1),hg2R(2),5);
hg2R_fix = [0.01;0.02;0.03;0.04;0.05;0.06;0.07;0.08;0.09;0.1];
C1_fix = [-21.835;-5.263;-2.021;-2.717;-2.740;-3.236;-2.999;-2.007;-1.609;-1.356];
C2_fix = [331.125;54.511;16.173;25.896;26.723;36.951;32.912;19.869;14.623;11.5];
C3_fix = [-1097.3520;126.801;166.247;71.876;35.29;-48.125;-39.816;4.926;19.0950;25.055];
C4_fix = [1155.349;-615.371;-468.866;-260.281;-165.404;7.134;1.682;-49.791;-62.449;-64.419];
if sum(hg2R>0.1)>0||sum(hg2R<0.01)>0
    error('the value of hg2R should be in the range of [0.01,0.1]')
end
hg = hg2R.*obj.indenterR;
if obj.unloadingNum == 0    
    LhloadingCao(:,1) = Lhloading.part1(:,1);
    LhloadingCao(:,2) = Lhloading.part1(:,3);
else
    % if the Lh curve contains multiple loading-unloading
    % cycles, then the loading part should be combined together
    % for Cao2004
    LhloadingCao = multiCombine ( obj );

end
Pg = interp1(LhloadingCao(:,1),LhloadingCao(:,3),hg);
epsR = 0.00939 + 0.435.*(hg2R)-1.106.*(hg2R).^2;
C1_interp = interp1(hg2R_fix,C1_fix,hg2R);
C2_interp = interp1(hg2R_fix,C2_fix,hg2R);
C3_interp = interp1(hg2R_fix,C3_fix,hg2R);
C4_interp = interp1(hg2R_fix,C4_fix,hg2R);
C = [C1_interp;C2_interp;C3_interp;C4_interp];     % C为4行,每行是C1,C2,C3,C4,每列对应一个hg2R值
% the following part is within function dimlsCao2004, to
% solve the sigR
sigR = zeros(1,length(hg));
for i = 1:length(hg)
    fs = @(x)Pg(i)-x*hg(i)^2*(C(1,i)*(log(matPara.elasModu/x))^3+ ...
        C(2,i)*(log(matPara.elasModu/x))^2+C(3,i)*log(matPara.elasModu/x)+C(4,i));
    sigR(i) = fsolve(fs,300);
end
% the following part is within function YNsolve, to
% solve yield stress and hardening index
[epsRx, sigRy] = prepareCurveData( epsR, sigR );
ft = fittype( ['Y*(1+',num2str(matPara.elasModu),'/Y*x)^N'], 'independent', 'x', 'dependent', 'y' );
% Fit model to data.
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0];
opts.MaxFunEvals = 5000;
opts.MaxIter = 500;
opts.StartPoint = [0.25 500];
opts.Upper = [0.5 2000];
[fitresult, ~] = fit( epsRx, sigRy, ft, opts );
matPara.sigY= fitresult.Y;
matPara.hdnIndx = fitresult.N;
matPara.KInPower = matPara.sigY*(matPara.sigY / matPara.elasModu)^(-matPara.hdnIndx);

end