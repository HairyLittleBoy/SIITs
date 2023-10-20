function [powerParaOut,powerTypeOut] = powerParaTrans(powerParaIn, powerTypeIn)
%POWERPARATRANS transforms the parameters between two types of power law
%stress - (plastic) strain model
%   powerType : plStrn -- sig = sigY * (1 + E/sigY * epsP)^N
%                 Strn -- sig = K * eps^n (sig > yield stress)
%   Input : powerParaIn -- the parameters in strss-(pl.) strain model needs
%                          to be transformed;
%           powerTypeIn -- the power model type needs to be transformed

%   powerPara for sig = sigY * (1 + E/sigY * epsP)^N : 
%                 powerPara(1) -- E
%                 powerPara(2) -- sigY
%                 powerPara(3) -- N

%   powerPara for sig = K * eps^n (sig > yield stress) : 
%                 powerPara(1) -- E
%                 powerPara(2) -- K
%                 powerPara(3) -- n
if ~strcmp(powerTypeIn,'plStrn') && ~strcmp(powerTypeIn,'strn')
    error('error @powerParaTrans : the 2nd input should be ''plStrn'' or ''strn''')
end
switch powerTypeIn
    case 'plStrn'
        % the paras in sig = sigY * (1 + E/sigY * epsP)^N will be
        % transformed to paras in sig = K * eps^n
        E = powerParaIn(1);
        sigY = powerParaIn(2);
        N = powerParaIn(3);
        epsP=linspace(0,1,1000);
        sig=sigY .* (1 + E/sigY .* epsP).^N;
        epsE = sig./E;
        eps=epsP+epsE;
        ft = fittype( 'power1' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [sigY 0.2];
        [fitresult, ~] = fit( eps', sig', ft ,opts);
        powerParaOut(1) = E;
        powerParaOut(2) = fitresult.a;
        powerParaOut(3) = fitresult.b;
        powerTypeOut = 'sig = K * eps^n ';
    case 'strn'
        % the paras in sig = K * eps^n will be transformed to paras in 
        % sig = sigY * (1 + E/sigY * epsP)^N
        E = powerParaIn(1);
        K = powerParaIn(2);
        n = powerParaIn(3);
        powerParaOut(1) = E;
        eps=linspace(0,1,1000);
        sig = K .* eps.^n;

        powerParaOut(1) = E;
        sigYtemp = (K/(E^n))^(1/(-n+1));
        epsP = eps(sig > sigYtemp) - sig(sig > sigYtemp) ./ E;
        epsP(1) = 0;
        xData = epsP';
        yData = sig(sig > sigYtemp)';
        yData(1) = sigYtemp;
        ft = fittype( ['sigY * (1 + ',num2str(E),'/sigY * x)^N'], 'independent', 'x', 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
        opts.Display = 'Off';
        opts.StartPoint = [0.2 sigYtemp];

        [fitresult, ~] = fit( xData, yData, ft, opts);
        powerParaOut(2) = fitresult.sigY;
        powerParaOut(3) = fitresult.N;
        powerTypeOut = 'sig = sigY * (1 + E/sigY * epsP)^N';
end

end






