function strsPlstrnPoints = strsPlstrnVoce(obj)
    if length(obj.modelParas) ~= 4
        error('the length of modelParas is not matched ')
    else
        strsPlstrnPoints(:,1) = linspace(0,1,400);
        disp('voce model:[E,Y,A,m]')
        strsPlstrnPoints(:,2) = obj.modelParas(2)+ ...
            obj.modelParas(3).*(1-exp(-obj.modelParas(4).*strsPlstrnPoints(:,1)));
    end
end