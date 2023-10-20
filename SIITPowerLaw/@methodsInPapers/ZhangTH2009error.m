function error = ZhangTH2009error(x,W,obj)
% ------------------------------------------------------------------------
% this function is the Eq(17) and Eq(21) in the paper, from which epsy and 
% n are solved

 % x(1): epsy,
 % x(2):n
% ------------------------------------------------------------------------
error1 = (0.87*0.5*(x(1)/(0.7781*x(1)+0.0155))^(1-x(2)))/((1/2-(1/(x(2)+1)))*(x(1)/(0.7781*x(1)+0.0155))^(1+x(2))+1/(x(2)+1)) - W(1);
error2 = (-48.755*x(2)^2 + 50.885*x(2) - 12.615)*x(1) + 0.845*x(2)^2 - 1.085*x(2) +1.44 - W(1)/W(2);
error = abs(error1) + abs(error2);

end

