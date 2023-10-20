clear all

LhdataExp= -LhImport('Lh_fit');

obj = LhCurve(LhdataExp);
sphrIndenObj=sphrInden([120000,200,0.1],0.3,0.5,0.17,0,1,8,'interp');
numOfLocal = 20;
maxIterOpt = 100;
MatParaIni = [358,0.1];
for i=1:2
[ bestYN(i) ,localYN(i), lnrCoef(i), errorOpt(:,i)] = LOFMethod(sphrIndenObj,obj,numOfLocal,maxIterOpt,MatParaIni)
MatParaIni = [bestYN.YBest,bestYN.nBest];
end
[ bestYN ,localYN, lnrCoef, errorOpt] = LOFMethod(sphrIndenObj,obj,numOfLocal,maxIterOpt,MatParaIni);
% a =        4596;
% b =      0.9872;
% x=linspace(0,0.172,1000);
% y =  a.*x.^b;

Lhbest = LhImport('Lh_R0p5E120Y323p1130501N0p1516852546fc0');
Lhbest=LhDivide(Lhbest);
plot(LhdataExp(:,1),LhdataExp(:,2),'r')
hold on
plot(Lhbest(:,1),Lhbest(:,2),'gd')
z = smooth(Lhbest(:,2),'rloess');
plot(Lhbest(:,1),z,'k')

x = LhdataSimu(:,1);
y = LhdataSimu(:,2);
y = localYN.nOpt;

%%
clear all
indenterR = 0.25;
unloadingNum = 8;
% E = 210000, sigY = 453, n=0.137
LhPoints= -LhImport('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\FEM\LhR025');
obj = LhCurve(LhPoints, indenterR, unloadingNum);
[Pks, Lhloading, Lhunloading ] = Lhsplit ( obj );
plot(obj.LhPoints(:,1),obj.LhPoints(:,2),'d')
hold on
plot(pksLocLow,pksLow,'rd')

plot(pksLocTemp1,pksTemp1,'gd','MarkerSize',20)
plot(pksLocHigh2,pksHigh2,'kd')

figure
for i  = 1:8
% eval(['plot(Lhloading.part',num2str(i),'(:,1),Lhloading.part',num2str(i),'(:,3))'])
% hold on
eval(['plot(Lhunloading.part',num2str(i),'(:,1),Lhunloading.part',num2str(i),'(:,3))'])
hold on
end

figure
for i  = 1:9
% eval(['plot(Lhloading.part',num2str(i),'(:,2),Lhloading.part',num2str(i),'(:,3))'])
hold on
eval(['plot(Lhunloading.part',num2str(i),'(:,2),Lhunloading.part',num2str(i),'(:,3))'])
end

%% let us test the functions @methodsInPapers
clear all
indenterR = 0.25;
unloadingNum = 8;
indenterElasModu = 1000000;
% E = 210000, sigY = 453, n=0.137
LhPoints=-LhImport('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\FEM\LhR025');
% LhPoints= -LhImport('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\v3\Lh_fit');

obj = methodsInPapers(LhPoints, indenterR,unloadingNum,indenterElasModu);
matPara1 = GBT1(obj,[0.15,0.2]);
hg2R = [0.0138,0.0815];
matPara2 = Cao2004(obj,hg2R);

contact1 = linspace(0.149,0.162,30);
for i = 1:length(contact1)
    matPara1 = GBT1(obj,[contact1(i),0.2]);
    x(i) = matPara1.sigY;
    y(i) = matPara1.hdnIndx;
    dPN(i) = matPara1.dispersionPN;
    plot(matPara1.sigY,matPara1.hdnIndx,'d')
    hold on
end


%%
for i = 8
    eval(['tempPartx = Lhunloading.part',num2str(i),'(:,1);']);
    eval(['tempParty = Lhunloading.part',num2str(i),'(:,3);']);
    tempPartx = tempPartx - min(tempPartx);
    
    tempParty(find(tempPartx == 0)) = [];
    tempPartx(find(tempPartx == 0)) = [];
    
    [xS, yS] = prepareCurveData( tempPartx, tempParty );
    plot(xS,yS,'d')
    hold on
end

%%
clear all
indenterR = 0.25;
unloadingNum = 8;
indenterElasModu = 1000000;
LhdataExp=-LhImport('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\FEM\LhR025');

obj = LhCurve(LhdataExp, indenterR, unloadingNum);
sphrIndenObj=sphrInden([210000,300,0.1],0.3,indenterR,0.17,0,1,8,'interp');
numOfLocal = 20;
maxIterOpt = 100;
MatParaIni = [358,0.1];
[ bestYN(i) ,localYN(i), lnrCoef(i), errorOpt(:,i)] = LOFMethod(sphrIndenObj,obj,numOfLocal,maxIterOpt,MatParaIni)

%% 
N=0.15;
sigY=400;
E=210000;
K = 1.023317427622841e+03;
n = 0.151002237771716;

[powerParaOut,powerTypeOut] = powerParaTrans([E,sigY,N], 'plStrn')
[powerParaOut,powerTypeOut] = powerParaTrans([E,K,n], 'strn')



%%
clear all
indenterR = 1.588/2;
unloadingNum = 7;
indenterElasModu = 600000;

LhPoints=-LhImport('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\v3\NaiBo20230223');
obj = methodsInPapers(LhPoints, indenterR,unloadingNum,indenterElasModu);

LhCurveMatrixU = unique(obj.LhPoints,'rows','stable') ;
LhCurveMatrixU(:,2) = smooth(LhCurveMatrixU(:,2),20);
LhCurveMatrixU(:,1) = smooth(LhCurveMatrixU(:,1),20);

[maxLoad,maxLoadNum] = max(LhCurveMatrixU(:,2)) ;
maxLoadLoc = LhCurveMatrixU(maxLoadNum,1);
peneNum=[1:length(LhCurveMatrixU)];

% pkNumLow is the number of the unloading bottom points in the whole LhCurveMatrixU data
[pksLow, pkNumLow] = findpeaks(-LhCurveMatrixU(:,2), peneNum,  ...
    'MinPeakHeight',-10, 'SortStr','ascend', 'NPeaks',obj.unloadingNum) ;
pksLow = -pksLow;

[pkNumLow, index1]= sort(pkNumLow);   % have to be sorted
pksLow = pksLow(index1);
pksLocLow = LhCurveMatrixU(pkNumLow, 1);

plot(LhCurveMatrixU(:,1),LhCurveMatrixU(:,2),'d')
hold on
plot(pksLocLow,pksLow,'rd')
%----------------------------
[pksTemp1, pkNumTemp1] = findpeaks(LhCurveMatrixU(:,2), peneNum, 'MinPeakHeight', 10) ;
pksLocTemp1 = LhCurveMatrixU(pkNumTemp1, 1);

% the maximum creep displacement is  set to 5 micrometer
pksLocTemp1Change = find((pksLocTemp1(2:end) - pksLocTemp1(1:length(pksLocTemp1)-1)) > 5e-3);

creepNum = cell(1,obj.unloadingNum);
creepPart = cell(1,obj.unloadingNum);
for i = 1 : obj.unloadingNum+1
    if i == 1
    creepNum{i} = find(abs(LhCurveMatrixU(:,2) - mean(pksTemp1(1:pksLocTemp1Change(i)))) < 0.5);
    elseif i == obj.unloadingNum+1
    creepNum{i} = find(abs(LhCurveMatrixU(:,2) - mean(pksTemp1(pksLocTemp1Change(i-1)+1:end))) < 0.5);
    else
    creepNum{i} = find(abs(LhCurveMatrixU(:,2) - mean(pksTemp1(pksLocTemp1Change(i-1)+1:pksLocTemp1Change(i)))) < 0.5);    
    end
    for j = 1:length(creepNum{i})-1
        temp1 = creepNum{i}(j+1) - creepNum{i}(j);
        if temp1 ~= 1
           creepNum{i}(j:end) = [];
           break
        end
    end
    creepPart{i} = LhCurveMatrixU(creepNum{i},:);
    plot(creepPart{i}(:,1),creepPart{i}(:,2),'yd')
    hold on
end

%----------------------------------------
plot(pksLocTemp1,pksTemp1,'kd')

plot(pksLocHigh2,pksHigh2,'gd')

kkk = diff(LhCurveMatrixU(:,2))./diff(LhCurveMatrixU(:,1));
plot(LhCurveMatrixU(2:end,1),kkk,'d')

%%
clear all
indenterR = 1.588/2;
unloadingNum = 7;
indenterElasModu = 600000;

LhPoints=-LhImport('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\v3\NaiBo20230223');
obj = methodsInPapers(LhPoints, indenterR,unloadingNum,indenterElasModu);
[Pks, Lhloading, Lhunloading, creep] = Lhsplit ( obj ,'creep');

plot(LhCurveMatrixU(:,1),LhCurveMatrixU(:,2),'rd')
hold on
plot(pksLocTemp1,pksTemp1,'yd')

for i = 1:8
plot(mean(Lhloading{i}(:,1)),mean(Lhloading{i}(:,5)),'bd')
hold on
end

for i = 1:8
plot(Lhloading{i}(:,1),Lhloading{i}(:,3),'r')
hold on
end


for i = 1:8
plot(Lhunloading{i}(:,1),Lhunloading{i}(:,3),'b')
hold on
end
for i = 1:8
plot(Lhunloading{i}(:,1),Lhunloading{i}(:,3))
hold on
end
for i = 1:8
plot(creep.creepPart{i}(:,1),creep.creepPart{i}(:,2),'g')
hold on
end

plot(pksLocHigh2,pksHigh2,'kd')
plot(pksLocHigh1,pksHigh1,'kd')

plot(LhloadingWhole(2:end,1),smooth(diff(LhloadingWhole(:,3))./diff(LhloadingWhole(:,1)),100))

plot(LhloadingWhole(:,1),LhloadingWhole(:,5),'d')


TSSimu = diff(LSimu)./diff(pene);
plot(pene(2:end),TSSimu,'d')


%%

plot(LhloadingWhole(:,1),LhloadingWhole(:,3),'gd');
hold on


for i = 1:length(errorOpt)
    x=[localYN.YOpt(i),localYN.nOpt(i)];
    sphrIndenObj.modelParas(2) = x(1);
    sphrIndenObj.modelParas(3) = x(2);
    LhdataSimu = LhGenr(sphrIndenObj);
    tt = find(LhdataSimu(:,1)>=0);
    pene = LhdataSimu(tt,1);
    LSimu = LhdataSimu(tt,2);
    if i==1
    plot(pene,LSimu,'r')
    else
    plot(pene,LSimu,'k')
    end
    hold on
end

    
    
    
