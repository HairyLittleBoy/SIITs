function [ bestYN ,localYN, lnrCoef, errorOpt] = LOFMethod(sphrIndenObj,LhExpObj,numOfLocal)
%LOFMethod identify the Y and N in power law with LOF method
%   sphrIndenObj: obj of sphrInden class
%   LhExpObj: obj of LhCurve class
%   numOfLocal: num Of Local domains
%   maxIter: maxIter in optimization
%   MatParaIni: initial guess of Y and N
MatParaIni = sphrIndenObj.modelParas(2:3);
Nlocal = linspace(0.5,0,numOfLocal);

Ylocal = linspace(0.2*MatParaIni(1),5*MatParaIni(1),numOfLocal);
errorOpt = zeros(numOfLocal,1);

for i = 1:length(Nlocal)
    sp = [Ylocal(i), Nlocal(i)];
    lb = [min(Ylocal),min(Nlocal)];
    ub = [max(Ylocal),max(Nlocal)];
    
    options = optimoptions('patternsearch');
    options = optimoptions(options,'AccelerateMesh', true);
    options = optimoptions(options,'UseCompletePoll', true);
    options = optimoptions(options,'UseCompleteSearch', true);
    options = optimoptions(options,'Display', 'Iter');
    options = optimoptions(options,'UseVectorized', false);
    options = optimoptions(options,'UseParallel', true);
    options = optimoptions(options,'MaxFunctionEvaluations',100);


    [x,funcValue,~,~] = patternsearch(@(x)errorLh( x, sphrIndenObj, LhExpObj),sp, ...
        [],[],[],[],lb,ub,[],options);

    localYN.YOpt(i) = x(1);
    localYN.nOpt(i) = x(2); 
    errorOpt(i) = funcValue;
    disp(['the ',num2str(i),'th local domain is done'])
end
[localYN.YOpt,Index] = sort(localYN.YOpt);
localYN.nOpt = localYN.nOpt(Index);
errorOpt = errorOpt(Index);
[~,minErrorPosition] = min(errorOpt);
bestYN.nBest = localYN.nOpt(minErrorPosition);
bestYN.YBest = localYN.YOpt(minErrorPosition);

% ----------------------------------------------------------------
ft = fittype( 'poly1' );
[xData, yData] = prepareCurveData( localYN.nOpt', localYN.YOpt' );
yData(find(xData>0.5)) = [];
xData(find(xData>0.5)) = [];
[rlt, ~] = fit( xData, yData, ft );
lnrCoef.p1 = rlt.p1;
lnrCoef.p2 = rlt.p2;


end

