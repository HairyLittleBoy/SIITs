function error = ChenXi2006error(x, obj, CS)

% -------------------------------------------------------------------------
% ChenXi2006error is the error = error1 + error2 + error3 in the 
% paper of ChenXi 2006, where error1 = C1/sigR1 - f1(m1,n)
%                             error2 = C1/sigR2 - f2(m2,n)
%                             error3 = S1/h2/E - g(m2,n)
% -------------------------------------------------------------------------

epsPR = [0.0374, 0.0674]';
h2R = [0.13, 0.3]';
sig0 = 500;                   % the reference yield stress(paper of ChenXi 2006, pp29 left,above Fig.5)

E2sigY = x(1);
sigY2sig0 = x(2);
n = x(3);

sigY = sigY2sig0 * sig0;     % yield stress
E = E2sigY * sigY;           % elastic modulus
Ebar = E / (1 - 0.3^2);      % the poission's ratio is set to be 0.3

sigR = zeros(2,1);
for i = 1:2
    fs1 = @(sigR)sigR-sigY*(E/sigY*(sigR/E + epsPR(i)))^n;
    sigR = fzero(fs1,sigY);
    sigR(i) = sigR;
end

m1 = E/sigR(1);m2 =  E/sigR(1); 
h1 = 32.77 - 52.59*log(m1) + 33.46*(log(m1))^2 - 4.8*(log(m1))^3 + 0.2147*(log(m1))^4;      % Eq(16)
h2 = 3.817 - 12.73*log(m2) + 11.99*(log(m2))^2 - 2.032*(log(m2))^3 + 0.1049*(log(m2))^4;    % Eq(17)
k1 = 1.001 + 0.261*n - 0.5217*n^2 + 0.1547*n^3;   % Eq(18)
k2 = 1.002 + 0.7637*n - 1.92*n^2 + 1.255*n^3;     % Eq(19)
A1 = 3.66556 + 0.0244179*n;    % Eq(21)
A2 = 6.06122 - 2.15891*n;      % Eq(22)
q = 29.0856 - 24.3547*n;       % Eq(23)
p = 1.31861 - 0.154675*n;      % Eq(24)

f1 = h1 * k1;
f2 = h2 * k2;
g = A2 + (A1 - A2)/(1 + (m2/q)^p);

error1 = CS.C(1)/sigR(1) - f1;
error2 = CS.C(2)/sigR(2) - f2;
error3 = CS.S/(obj.indenterR*h2R(2)*Ebar) - g;
error = abs(error1) + abs(error2) + abs(error3);

end









