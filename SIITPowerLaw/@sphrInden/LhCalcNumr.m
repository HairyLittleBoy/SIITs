function LhNumr = LhCalcNumr( obj )
% this function use E, Y, N, R to predict Lh curve，
% with pm=p1_dim*(a/R)^3+p2_dim*(a/R)^2+p3_dim*(a/R)+p4_dim
% formula to calculate pm(pm=L/(pi*a^2)),coefficients
% p1_dim-p4_dim are obtained with FEM and interpolation method
% it is usually more accurate than ANN
pene = linspace(0,obj.penetration,1000)';
EY = obj.modelParas(1)/obj.modelParas(2);
aR = sqrt(pene.*(2*obj.indenterR-pene))./obj.indenterR;

if length(obj.modelParas) == 3
    load('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\MatLabMatFiles\p1TOp4PrwLaw.mat', ...
        'p1_dim','p2_dim','p3_dim','p4_dim','XX','YY');
    p1_intrp = interp2(XX,YY,p1_dim',obj.modelParas(3),EY,'spline');
    p2_intrp = interp2(XX,YY,p2_dim',obj.modelParas(3),EY,'spline');
    p3_intrp = interp2(XX,YY,p3_dim',obj.modelParas(3),EY,'spline');
    p4_intrp = interp2(XX,YY,p4_dim',obj.modelParas(3),EY,'spline');
    pm_intrp = p1_intrp.*aR.^3+p2_intrp.*aR.^2+p3_intrp.*aR+p4_intrp;
    pm_intrp = pm_intrp.*(obj.modelParas(2)/200);  % 有限元模拟中全部用200MPa为屈服应力，只改变
    % E从而改变E/Y的值，所以这里要用Y/200对pm进行
    % 成倍数缩放
elseif length(obj.modelParas) == 4
    AY = obj.modelParas(3)/obj.modelParas(2);
    method = 'interp';
    switch method
        case 'interp'
            load('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\MatLabMatFiles\p1TOp3Voce.mat', ...
                'p1_dim','p2_dim','p3_dim','XX','YY','ZZ');
            p1_intrp = interp3(XX,YY,ZZ,p1_dim,EY,obj.modelParas(4),AY,'spline');
            p2_intrp = interp3(XX,YY,ZZ,p2_dim,EY,obj.modelParas(4),AY,'spline');
            p3_intrp = interp3(XX,YY,ZZ,p3_dim,EY,obj.modelParas(4),AY,'spline');
            pm_intrp = p1_intrp.*aR.^2+p2_intrp.*aR.^1+p3_intrp;
            pm_intrp = pm_intrp.*((obj.modelParas(2))/200);
        case 'Ann'
            load('E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\MatLabMatFiles\p1TOp3VoceWithAnn.mat', ...
                'netp1','netp2','netp3')
            p1_intrp = netp1([obj.modelParas(4);EY;AY]);
            p2_intrp = netp2([obj.modelParas(4);EY;AY]);
            p3_intrp = netp3([obj.modelParas(4);EY;AY]);
            pm_intrp = p1_intrp.*aR.^2+p2_intrp.*aR.^1+p3_intrp;
            pm_intrp = pm_intrp.*((obj.modelParas(2))/200);
    end
end
LhNumr(:,1) = pene';
LhNumr(:,2) = pm_intrp.*(pi.*aR.^2)./(1/obj.indenterR^2);   % 有限元模拟中全部用R=1mm，所
% 以这里用1/R^2成倍数缩放

end