function matPara = GBT1(obj, contactR)

%--------------------------------------------------------------------------------------------------
% This GBT1 function is the methods applied in the standard GB/T 39635-2020,
% which is drafted leadingly by Zhengyue Yu in Jiaotong University, Shanghai.
% The methods in this standard are actually from  work of Kang, S.-K., et al.
% (2013). IJP 49: 1-15, and the works of Tabor
%
% According to the GB/T 39635-2020, this method requires L-h
% curve with multiple loading-unloading cycles, so LhDivideMtpl
% should be applied.

% the sigY and hdnIndx are for the power law sig = K * eps ^hdnIndx
%--------------------------------------------------------------------------------------------------

% Step-1: the evaluation of hardening index n ---------

% With the Eq.(45) and Eq.(55) in Kang, S.-K., et al.
% (2013). IJP 49: 1-15, hardening index can be estimated

warning ([' GBT1@methodsInPapers: As far as can be seen, the GBT1 method applys better for R = 0.25mm, ', ...
    'and the hardening index will be evaluated at a1 = 0.15mm and a2 = 0.2mm,', ...
    ' also the a1 and a2 should not positioned within the range of first two loading', ...
    ' cycles '])

contactR1 = contactR(1); contactR2 = contactR(2);

if obj.unloadingNum > 0
    [ Pks, ~, ~ ] = Lhsplit( obj );
else
    error(['GBT1@methodsInPapers: The GBT1 method  only support Lh curve with mutiple loading-unloading cycles,', ...
        'however, one with single cycle is provided, error@GBT1 < methodsInPapers'])
end
LhloadingWhole = multiCombine ( obj );
ft = fittype( 'power1' );
[flt, ~] = fit( LhloadingWhole(2:end,2), LhloadingWhole(2:end,3), ft );
loadingRate1 = flt.a*flt.b*contactR1^(flt.b-1);
loadingRate2 = flt.a*flt.b*contactR2^(flt.b-1);

s =  loadingRate2 / loadingRate1;

% To apply to various indenter radius, the coefficients in Eq.(55) in Kang, S.-K., et al.
% (2013). IJP 49: 1-15 should be evaluated beforehand
hdnIndx = linspace(0.005,0.45,5);
epsY = linspace(0.001,0.01,10);
kk = 0.5098 + 0.0048.*exp(hdnIndx./0.0598);      % Eq.(40) in Kang, S.-K., et al.(2013). IJP 49: 1-15

funcPI3 = (2/3).*(1./hdnIndx+1.5.*kk);

loadingSlope = zeros(length(epsY),length(hdnIndx));    % the loading slope p in Eq.(45)
fitX = zeros(length(epsY)*length(hdnIndx),1);

for i = 1:length(epsY)
funcPI2 = (2/3).*(1-1./hdnIndx).*(1/4.*(1/epsY(i))*(1/obj.indenterR)).^(-hdnIndx);
loadingSlope(i,:) = (2.*funcPI2.*contactR(2)+funcPI3.*(hdnIndx+2).*contactR(2).^(hdnIndx+1))./ ...
                    (2.*funcPI2.*contactR(1)+funcPI3.*(hdnIndx+2).*contactR(1).^(hdnIndx+1));
        
fitX((i-1)*length(hdnIndx)+1:i*length(hdnIndx),1) = loadingSlope(i,:)';
% plot(loadingSlope(i,:),hdnIndx,'d')
% hold on
end
fitY = repmat(hdnIndx',length(epsY),1);
matPara.dispersionPN = mean([(loadingSlope(end,1) - loadingSlope(1,1)),(loadingSlope(end,end) - loadingSlope(1,end))]);
if matPara.dispersionPN > 0.01
    warning(' GBT1@methodsInPapers: the dispersion of p-n data points is large, please change the contactRs and try again')
end
ft = fittype( 'poly3' );
[fitresult, ~] = fit( fitX, fitY, ft);
coefsInEq55.p1 = fitresult.p1;
coefsInEq55.p2 = fitresult.p2;
coefsInEq55.p3 = fitresult.p3;
coefsInEq55.p4 = fitresult.p4;

% calculate the hardenning index
matPara.hdnIndx = coefsInEq55.p1 * s^3 + coefsInEq55.p2 * s^2 + coefsInEq55.p3 * s + coefsInEq55.p4;


if (matPara.hdnIndx < 0 || matPara.hdnIndx > 0.5)
    error (' GBT1@methodsInPapers: the hardeining index evaluated is out of normal range ')
elseif (matPara.hdnIndx > 0.3)
    warning(' GBT1@methodsInPapers: the hardeining index evaluated is maybe too large ')
end


% Step-2: estimate the true contact radius ---------

% The method given in GB/T 39635-2020 to estimate the true
% contact radius or contact penetration(hc) can be found in the
% paper "Kim, S. H., et al. (2006). MSEA 415(1-2): 59-65".

% To use the the function in IndenInprintCrct to estimate the
% true penetration and contact area, an IndenInprintCrct object
% should be built first.

% For GBT1 method, the Kim2006 method should be applied, which
% only take hardening index into considerration, so the elastic
% modulus and yield stress can be set to zero when define the
% IndenInprintCrct object.

pileUpSnkInObj = IndenImprintCrct(obj.LhPoints, obj.indenterR, obj.unloadingNum, ...
    obj.creepOrNot,0,0, matPara.hdnIndx);

% although Kim2006 is applied in GB/T 39635-2020, it is found
% that the method is hard to use, the unloading slope from Lh
% curve obtaind by FEM simulation is increasing along the loading-
% unloading cycles, resulting in a large ture contact area,
% which leads to the early decrease of the sigIT points.

[~, trueContactR, trueCntctArea] = Hill1989(pileUpSnkInObj);    % the true contact area at each loading peak

% Step-3: true stress sigIT and true strain epsIT ---------
constraintFactorPhi = 3;
ksi = 0.14;
PmTrue =  Pks.High1(:,2) ./ trueCntctArea;
sigIT = PmTrue ./ constraintFactorPhi;           % Eq. (4) in GB/T 39635-2020
epsIT =  ksi ./ sqrt(1 - (trueContactR ./ obj.indenterR) .^ 2) .*  trueContactR ./ obj.indenterR; % Eq. (5) in GB/T 39635-2020

% Step-4: solve the yield stress
% K in sig=K*eps^n should be first fitted

[xData, yData] = prepareCurveData( epsIT, sigIT );
ft = fittype( 'power1' );
[ft, ~] = fit( xData, yData, ft );
matPara.KInPower = ft.a;

matPara.elasModu = OliverPharr(obj);
% end

% with K in sig=K*eps^n and elastic modulus, we can evaluate
% the yield stress with the yield point, according to Eq. (8)
% in GB/T 39635-2020, yield stress can be calculated with from
% following nonlinear equation.
fs = @(epsY)matPara.KInPower*epsY^matPara.hdnIndx - matPara.elasModu*(epsY - 0.002);
epsY = fsolve(fs, 0.002);
matPara.sigY = matPara.elasModu * (epsY - 0.002);

matPara.sigU =  matPara.KInPower *  matPara.hdnIndx ^ matPara.hdnIndx;

end