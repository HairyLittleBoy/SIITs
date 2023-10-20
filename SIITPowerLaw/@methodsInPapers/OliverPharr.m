function elasModu = OliverPharr(obj)
    %  estimate the elastic modulue the Oliver-Pharr method
    [ Pks, ~, Lhunloading ] = Lhsplit( obj );
    omega = 0.75;
    unloadingSlope = zeros(length(Lhunloading), 1);
    ft = fittype( 'power1' );
    for i = 1 : length(Lhunloading)
        tempPartx = Lhunloading{i}(:,1);
        tempParty = Lhunloading{i}(:,3);
        tempPartx = tempPartx - min(tempPartx);

        tempParty(find(tempPartx == 0)) = [];
        tempPartx(find(tempPartx == 0)) = [];

        [xS, yS] = prepareCurveData( tempPartx, tempParty );
        [fsS, ~] = fit( xS, yS, ft );
        unloadingSlope(i) = fsS.a*fsS.b*max(tempPartx)^(fsS.b - 1);
    end
    hs = omega .* Pks.High1(:,2) ./ unloadingSlope;
    hc =  Pks.High1(:,1) - hs;
    trueCntctArea = sqrt(hc.*(2*obj.indenterR - hc));
    Er = sqrt(pi)/2.*unloadingSlope./sqrt(trueCntctArea);
    elasModuAll = (1 - 0.3^2)./(1./Er - (1 - 0.22^2)/obj.indenterElasModu);
    % the beta is a parameter in GB¨MT 37782-2019(Cai) appendix A
    hu2D = Pks.High2(:,1) ./ (2*obj.indenterR);
    beta = 1383.*hu2D.^2 - 39.69.*hu2D + 1.387;
    elasModu = mean(elasModuAll(1:obj.unloadingNum) ./ beta);
end