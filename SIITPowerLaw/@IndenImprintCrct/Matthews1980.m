function [truePene, trueContactRadius, trueCntctArea] =  Matthews1980(obj)
    % (3) 
    %        The methods mentioned in Taljat and Pharr (2004) IJSS
    %        41 (2004) 3891�C3904. These two methods proposed by Hill
    %        and Matthews estimate the true contact penetration h
    %        based on only hardening indenx n, while the method
    %        given by Taljat and Pharr (2004) took the E/Y (yield
    %        strain) into consideration 
    %       

    c2Hill = 0.5 * ((2 + obj.hardeningIndex) / 2) ^ (2*(1 - obj.hardeningIndex) / obj.hardeningIndex);
    [ Pks, ~, ~ ] = Lhsplit( obj );

    hc = c2Hill .* Pks.High1(:,1);
    ac = sqrt( hc .* ( 2*obj.indenterR - hc ));

    trueContactRadius = ac;
    trueContactArea = pi .* ac .^ 2 ;           
    truePene = hc;
    trueCntctArea = trueContactArea;

end   