function [ LhAnn ] = LhCalcAnn(obj)
% predict Lh with ANN based on power law model
% indenterR: indenter radius in mm
netRand1 = randperm(100,2);
netRand2 = randperm(400,2);

netPath = 'E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\MatLabMatFiles\PrLawNets\';

for i=1:length(netRand1)
    numNet1 = netRand1(i);
    
    load([netPath,'nets_fp1_hdn10_',num2str(numNet1)])
    load([netPath,'nets_fp2_hdn10_',num2str(numNet1)])
    load([netPath,'nets_fp3_hdn10_',num2str(numNet1)])
    load([netPath,'nets_fp4_hdn10_',num2str(numNet1)])
    load([netPath,'nets_fp5_hdn10_',num2str(numNet1)])
    load([netPath,'nets_fp1_hdn40_',num2str(numNet1)])
    load([netPath,'nets_fp2_hdn40_',num2str(numNet1)])
    load([netPath,'nets_fp3_hdn40_',num2str(numNet1)])
    load([netPath,'nets_fp4_hdn40_',num2str(numNet1)])
    load([netPath,'nets_fp5_hdn40_',num2str(numNet1)])
end
for i=1:length(netRand2)
    numNet2 = netRand2(i);
    load([netPath,'nets_fp1_hdn30_',num2str(numNet2)])
    load([netPath,'nets_fp2_hdn30_',num2str(numNet2)])
    load([netPath,'nets_fp3_hdn30_',num2str(numNet2)])
    load([netPath,'nets_fp4_hdn30_',num2str(numNet2)])
    load([netPath,'nets_fp5_hdn30_',num2str(numNet2)])
    load([netPath,'nets_fp1_hdn20_',num2str(numNet2)])
    load([netPath,'nets_fp2_hdn20_',num2str(numNet2)])
    load([netPath,'nets_fp3_hdn20_',num2str(numNet2)])
    load([netPath,'nets_fp4_hdn20_',num2str(numNet2)])
    load([netPath,'nets_fp5_hdn20_',num2str(numNet2)])
    
end


fp1Total = 0;
fp2Total = 0;
fp3Total = 0;
fp4Total = 0;
fp5Total = 0;
x3 = [obj.modelParas(1)/obj.modelParas(2);obj.modelParas(3)];
numNN = 1000;               % 用4次多项式拟合h/R-pm曲线，系数作为神经网络的输出。
% 为了提高神经网络的泛化能力，预测每个系数用1000个小网络

for i = 1:length(netRand1)
    numNet1 = netRand1(i);
    eval(['net1i = nets_fp1_hdn10_',num2str(numNet1),';'])
    eval(['net2i = nets_fp2_hdn10_',num2str(numNet1),';'])
    eval(['net3i = nets_fp3_hdn10_',num2str(numNet1),';'])
    eval(['net4i = nets_fp4_hdn10_',num2str(numNet1),';'])
    eval(['net5i = nets_fp5_hdn10_',num2str(numNet1),';'])
    
    fp1 = net1i(x3);
    fp2 = net2i(x3);
    fp3 = net3i(x3);
    fp4 = net4i(x3);
    fp5 = net5i(x3);
    
    fp1Total = fp1Total + fp1;
    fp2Total = fp2Total + fp2;
    fp3Total = fp3Total + fp3;
    fp4Total = fp4Total + fp4;
    fp5Total = fp5Total + fp5;
end

for i = 1:length(netRand2)
    numNet2 = netRand2(i);
    eval(['net1i = nets_fp1_hdn30_',num2str(numNet2),';'])
    eval(['net2i = nets_fp2_hdn30_',num2str(numNet2),';'])
    eval(['net3i = nets_fp3_hdn30_',num2str(numNet2),';'])
    eval(['net4i = nets_fp4_hdn30_',num2str(numNet2),';'])
    eval(['net5i = nets_fp5_hdn30_',num2str(numNet2),';'])
    
    fp1 = net1i(x3);
    fp2 = net2i(x3);
    fp3 = net3i(x3);
    fp4 = net4i(x3);
    fp5 = net5i(x3);
    
    fp1Total = fp1Total + fp1;
    fp2Total = fp2Total + fp2;
    fp3Total = fp3Total + fp3;
    fp4Total = fp4Total + fp4;
    fp5Total = fp5Total + fp5;
    
end

for i = 1:length(netRand2)
    numNet2 = netRand2(i);
    eval(['net1i = nets_fp1_hdn20_',num2str(numNet2),';'])
    eval(['net2i = nets_fp2_hdn20_',num2str(numNet2),';'])
    eval(['net3i = nets_fp3_hdn20_',num2str(numNet2),';'])
    eval(['net4i = nets_fp4_hdn20_',num2str(numNet2),';'])
    eval(['net5i = nets_fp5_hdn20_',num2str(numNet2),';'])
    
    
    fp1 = net1i(x3);
    fp2 = net2i(x3);
    fp3 = net3i(x3);
    fp4 = net4i(x3);
    fp5 = net5i(x3);
    
    fp1Total = fp1Total + fp1;
    fp2Total = fp2Total + fp2;
    fp3Total = fp3Total + fp3;
    fp4Total = fp4Total + fp4;
    fp5Total = fp5Total + fp5;
end
for i = 1:length(netRand1)
    numNet1 = netRand1(i);
    
    eval(['net1i = nets_fp1_hdn40_',num2str(numNet1),';'])
    eval(['net2i = nets_fp2_hdn40_',num2str(numNet1),';'])
    eval(['net3i = nets_fp3_hdn40_',num2str(numNet1),';'])
    eval(['net4i = nets_fp4_hdn40_',num2str(numNet1),';'])
    eval(['net5i = nets_fp5_hdn40_',num2str(numNet1),';'])
    
    fp1 = net1i(x3);
    fp2 = net2i(x3);
    fp3 = net3i(x3);
    fp4 = net4i(x3);
    fp5 = net5i(x3);
    
    fp1Total = fp1Total + fp1;
    fp2Total = fp2Total + fp2;
    fp3Total = fp3Total + fp3;
    fp4Total = fp4Total + fp4;
    fp5Total = fp5Total + fp5;
    
end
fpAverageOutput(1) = fp1Total / (2*length(netRand1)+2*length(netRand2));
fpAverageOutput(2) = fp2Total / (2*length(netRand1)+2*length(netRand2));
fpAverageOutput(3) = fp3Total / (2*length(netRand1)+2*length(netRand2));
fpAverageOutput(4) = fp4Total / (2*length(netRand1)+2*length(netRand2));
fpAverageOutput(5) = fp5Total / (2*length(netRand1)+2*length(netRand2));

pm(:,1) = linspace(0,0.3,1000);
pm(:,2) = fpAverageOutput(1).*pm(:,1).^4 + fpAverageOutput(2).*pm(:,1).^3 + ...
    fpAverageOutput(3).*pm(:,1).^2 + fpAverageOutput(4).*pm(:,1) + ...
    fpAverageOutput(5);
LhAnn(:,1) = pm(:,1).*obj.indenterR;
LhAnn(:,2) = (obj.modelParas(2)/200).*pm(:,2).*pi.*(obj.indenterR^2 - (obj.indenterR - LhAnn(:,1)).^2);  %在神经网络的训练中，输入是
%[E/Y,N],而Y被设定为200MPa
%所以这里要比上200
end