function [truePene, trueContactRadius, trueCntctArea] = Kim2006(obj)
    % (1)
    %        The method given in "Kim, S. H., et al. (2006). MSEA 415(1-2): 59-65"
    %        where the influence of elastic deflection on penetration (hd) is also
    %        considerred, which is not in the paper of Taljat, B. and G. M. Pharr (2004).
    %        IJSS 41(14): 3891-3904".
    %        With Eq.(1), Eq.(4) and Eq.(10) in the paepr, the hc can be estimated
    omega = 0.75;
    [ Pks, ~, Lhunloading ] = Lhsplit( obj );
    S = zeros(length(fieldnames(Lhunloading)), 1);
    ft = fittype( 'poly1' );

    for i = 1 : length(fieldnames(Lhunloading))
        tempPart = ['Lhunloading.part',num2str(i)];
        eval(['upperPortionUnloadingNum = find(',tempPart,'(:,3) > 0.7 * max(',tempPart,'(:,3)));'])
        eval(['upperPortionUnloading = ',tempPart,'(upperPortionUnloadingNum,:);'])

        [xS, yS] = prepareCurveData( upperPortionUnloading(:,1), upperPortionUnloading(:,3) );
        [fsS, ~] = fit( xS, yS, ft );
        S(i) = fsS.p1;
    end

    hd = omega .* Pks.High1(:,2) ./ S;
    hcstar =  Pks.High1(:,1) - hd;
    coef_hPileStar2hcstar = 0.131 .* (1 - 3.423 *  obj.hardeningIndex + 0.079 *  obj.hardeningIndex ^ 2) .* ...
        (1 + 6.258 .* Pks.High1(:,1) ./ obj.indenterR - 8.072 .* (Pks.High1(:,1) ./ obj.indenterR) .^2 );
    hPileStar = coef_hPileStar2hcstar .* hcstar;
    hc = Pks.High1(:,1) - hd + hPileStar;
    trueContactRadius = sqrt(hc.*(2*obj.indenterR));
    trueContactArea = pi .* trueContactRadius .^ 2 ;

    truePene = hc;
    trueCntctArea = trueContactArea;

end