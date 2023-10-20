%% example for Lh curve with multiple loading-unloading cycles
% E = 210000MPa, yield stress = 453MPa, hardening index = 0.137 in 
% sig = K * eps^n type of power law
clear all
indenterR = 0.25;
unloadingNum = 8;
indenterElasModu = 1000000;
creepOrNot = 'nocreep';
LhPoints=-LhImport('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\FEM\LhR025');

obj = methodsInPapers(LhPoints, indenterR,unloadingNum,creepOrNot,indenterElasModu);
matParaGBT1 = GBT1(obj,[0.166,0.2]);

hg2R = [0.042,0.0815];
matParaCao = Cao2004(obj,hg2R);

% define the sphrInden object, mlPs is the modelParas in strsPlstrnCrv
% class, the power law is sig = sigY * (1 + E/sigY * epsP)^N. It has to be
% noticed that, most non-numerical identification methods are based on the
% sig = K * eps^n type of power law, it means that we have to TRANSFORM the
% results given by GBT1 or Cao2004 mothod to the parameters in 
% sig = sigY * (1 + E/sigY * epsP)^N
[powerParaOutGBT1,powerTypeOut] = powerParaTrans([matParaGBT1.elasModu, matParaGBT1.KInPower, matParaGBT1.hdnIndx], 'strn');
[powerParaOutCao,powerTypeOut] = powerParaTrans([matParaCao.elasModu, matParaCao.KInPower, matParaCao.hdnIndx], 'strn');

% use only the parameters provided by GBT1 method in LOF 
sphrIndenObj = sphrInden(powerParaOutCao, 0.3, indenterR, 0.6*indenterR, 0, 1, 5, 'interp');
LhExpObj = LhCurve(LhPoints, indenterR, 8);
[ bestYN ,localYN, lnrCoef, errorOpt] = LOFMethod(sphrIndenObj,LhExpObj,8,10);

% the parameters generated are for sig = K * eps^n type of power law, let's
% transform them to sig = K * eps^n
[powerParaGBT1Final,powerTypeOut] = powerParaTrans([matParaCao.elasModu, bestYN.YBest, bestYN.nBest], 'plStrn');


%% Now, data from spherical indentation test experiment will be tried,
% however, the Lh curve data is not good enough, means we have to do some
% doing to fix the data

clear all
indenterR = 1.5875/2;
unloadingNum = 7;
indenterElasModu = 600000;
creepOrNot = 'creep';

LhPoints=-LhImport('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\v3\NaiBo20230223.txt');
obj = methodsInPapers(LhPoints, indenterR,unloadingNum,creepOrNot,indenterElasModu);

% matParaGBT1 = GBT1(obj,[0.2,0.25]);     % matParaGBT1 gives ridiculus result
hg2R = [0.042,0.0815];
matParaCao = Cao2004(obj,hg2R);
matParaHag = Haggag1993(obj);
matParaGBT1 =  GBT1(obj, [0.15,0.2]);
matParaChenXi = ChenXi2006(obj);
matParaZhangTH = ZhangTH2009(obj);


penetration = 0.8*max(LhPoints(:,1));
sphrIndenObj = sphrInden([matParaHag.elasModu,matParaHag.sigY,matParaHag.hdnIndx], ...
                        0.3, indenterR, 0.95*max(LhPoints(:,1)), 0, 1, 5, 'interp');
[ bestYN ,localYN, lnrCoef, errorOpt] = LOFMethod(sphrIndenObj,obj,100);

finalSigY = mean([matParaHag.sigY,matParaCao.sigY]);
finalN = (finalSigY - lnrCoef.p2) / lnrCoef.p1;
[powerParaOutfinal,powerTypeOut] = powerParaTrans([matParaCao.elasModu, finalSigY, finalN],'plStrn');








