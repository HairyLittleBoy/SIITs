function error = errorLh( x, sphrIndenObj, LhExpObj)
% 此函数为优化目标函数，对比exact目标ph曲线与预测的ph曲线之间的差异
% 预测的ph曲线可以由FEM生成也可以由解析方（ANA）方法生成
% 最后优化的目标量为C的差异和ph曲线点的差异相加组成
if length(sphrIndenObj.modelParas) == 3
    sphrIndenObj.modelParas(2) = x(1);
    sphrIndenObj.modelParas(3) = x(2);
elseif length(sphrIndenObj.modelParas) == 4
    sphrIndenObj.modelParas(2) = x(1);
    sphrIndenObj.modelParas(3) = x(2);
    sphrIndenObj.modelParas(4) = x(3);
end
    
LhdataSimu = LhGenr(sphrIndenObj);

% ft = fittype( 'poly3' );
% [fltExp, ~] = fit( LhExpObj.LhPoints(:,1), LhExpObj.LhPoints(:,2), ft );
% LExp = fltExp.p1.*LhExpObj.LhPoints(:,1).^3 + fltExp.p2.*LhExpObj.LhPoints(:,1).^2 + ...
%     fltExp.p3.*LhExpObj.LhPoints(:,1) + fltExp.p4;
% [fltSimu, ~] = fit( LhdataSimu(:,1), LhdataSimu(:,2), ft );
% LSimu = fltSimu.p1.*LhExpObj.LhPoints(:,1).^3 + fltSimu.p2.*LhExpObj.LhPoints(:,1).^2 + ...
%     fltSimu.p3.*LhExpObj.LhPoints(:,1) + fltSimu.p4;

tt = find(LhdataSimu(:,1)>0);
pene = LhdataSimu(tt,1);
LSimu = LhdataSimu(tt,2);

[ ~, Lhloading, ~ ] = Lhsplit( LhExpObj );
if LhExpObj.unloadingNum == 0    
    LhloadingWhole(:,1) = Lhloading.part1(:,1);
    LhloadingWhole(:,2) = Lhloading.part1(:,3);
    LhloadingWhole(:,3) = Lhloading.part1(:,5);
else
    % if the Lh curve contains multiple loading-unloading
    % cycles, then the loading part should be combined together

    LhloadingWhole = multiCombine ( LhExpObj );

end
pene(pene>=max(LhloadingWhole(:,1))) = [];
LExp = interp1(LhloadingWhole(:,1),LhloadingWhole(:,3),pene);

C_Simu = LSimu(1:length(pene))./pene.^2;
C_Exp = LExp./pene.^2;
error = mean(abs(C_Simu-C_Exp)./C_Exp) + mean(abs(LSimu(1:length(pene)) - LExp)./LExp); 

% mtd = sphrIndenObj.LhGenrMethod;
% fid = fopen('optimLog.txt','a');
% newline = {['x:',num2str(x(1)),'  ',num2str(x(2)),'  error:',num2str(error),'  ',datestr(clock),'  ',mtd]};
% fprintf(fid,'%s\t\n',char(newline));
% fclose(fid);
end

