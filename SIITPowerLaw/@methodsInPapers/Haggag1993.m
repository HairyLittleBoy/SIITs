function matPara = Haggag1993(obj)

% -------------------------------------------------------------------------
% the method in this function is from the paper "In-situ measurements of 
% mechanical properties using novel automated ball indentation system, 
% ASTM Special Technical Publication"

%  the relation describing stress-strain is 
%  sig = K * epsP ^ n
%  Haggag1993 needs several loading-unloading cycles
% -------------------------------------------------------------------------

matPara.elasModu = OliverPharr(obj);
if obj.unloadingNum == 0    
    error('Haggag1993 needs several loading-unloading cycles')
else
    % if the Lh curve contains multiple loading-unloading
    % cycles, then the loading part should be combined together
    % for Haggag1993
    [Pks, ~, Lhunloading, ~] = Lhsplit ( obj );
end
hp = zeros(length(Lhunloading),1);           % penetration after unloading (plastic penetration)
dp = zeros(length(Lhunloading),1);           % imprint diameter after unloading (plastic penetration)
dp0 = zeros(length(Lhunloading),1);          % initial point for solving dp (plastic penetration)
C = zeros(length(Lhunloading),1);            % the C in the paper (Eq.(5))
matPara.epsP = zeros(length(Lhunloading),1); % plastic strain 
temp1 = zeros(length(Lhunloading),1);        % temp1 * cf = phi in the paper
cf = zeros(length(Lhunloading),1);           % the constraint factor
ht = zeros(length(Lhunloading),1);           % total penetration before unloading (max penetration in each cycle)
dt = zeros(length(Lhunloading),1);           % total imprint diameter before unloading (max imprint diameter in each cycle)
matPara.epsP = zeros(length(Lhunloading)+1,1);
matPara.sig = zeros(length(Lhunloading)+1,1);
D = 2 * obj.indenterR;
alfam = 1;                                   % alfam in the paper(Eq.(8))
cfMax = 2.87 * alfam;                        % maximum constraint factor(Eq.(8) in the paper)
tau = (cfMax -1.12)/log(27);                 % Eq.(9) in the paper
for i = 1 : length(Lhunloading)
    hp(i) = Lhunloading{i}(end,1);
    ht(i) = Pks.High1(i,1);
    dt(i) = 2*sqrt(ht(i)*D-ht(i)^2);
    C(i) = 5.47 * Pks.High1(i,2)*(1/obj.indenterElasModu - 1/matPara.elasModu);
    
    fs1 = @(dp)nthroot((0.5*C(i)*D*(hp(i)^2+(dp/2)^2)/(hp(i)^2+(dp/2)^2-hp(i)*D)),3) - dp;    % Eq(4) in the paper
    dp0(i) = 2*sqrt(obj.indenterR^2 - (obj.indenterR - hp(i))^2)*0.99;
    dp(i) = fzero(fs1,dp0(i));    
    
    matPara.epsP(i+1) = 0.2*dp(i)/D;
    
    temp1(i) = matPara.epsP(i+1)*matPara.elasModu*pi*dp(i)^2/1.72/Pks.High1(i,2);
    fs2 = @(cf)cf - ((temp1(i)*cf<=1)*1.12 +(temp1(i)*cf>1 && temp1(i)*cf<=27)*(1.12+tau*log(temp1(i)*cf))+(temp1(i)*cf>27)*cfMax);  % Eq(6) in the paper
    cf(i) = fzero(fs2,cfMax);    
    matPara.sig(i+1) = 4*Pks.High1(i,2)/(pi*dp(i)^2*cf(i));
end
phi = cf.*temp1;

betam = 0.2285;                         % Eq(12) in the paper (betam is Bm in the paper)
Pdt2 = Pks.High1(:,2)./dt.^2;
dtD = dt./D;
ft = fittype( 'power1' );

[flt1, ~] = fit( dtD, Pdt2, ft);         % fit the Eq(11) in the paper for the A
matPara.sigY = betam*flt1.a;             % yield stress is betam * A

matPara.epsP(1) = 1e-5;
matPara.sig(1) = matPara.sigY;

[flt2, ~] = fit( matPara.epsP, matPara.sig, ft); 
matPara.Kp = flt2.a;
matPara.hdnIndx = flt2.b;
end














