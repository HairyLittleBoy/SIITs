classdef classifier < LhCurve
    %CLASSIFIER  contains the methods which classfies the given L-h curve, provide
    %            the possible range of parameter in the equation describing the stress
    %            plastic strain curve
    %            此类所包含的成员函数对L-h曲线进行分类，给出对应描述应力应变曲线
    %            方程参数的可能取值范围
    % based on the jornal paper:
    % Zhankun, Sun, et al. "Study on concavity-convexity 
    % transition of loading curve for spherical indentation." 
    % Mechanics of Materials 114(2017): 107-118.
    
    properties
        classfyMethod
    end
    
    methods
        function obj = classifier(LhP,clsMd)
            obj = obj@LhCurve(LhP);
            obj.classfyMethod = clsMd;
        end
        
        function outputArg = ccTranstn(obj,inputArg)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            outputArg = obj.Property1 + inputArg;
        end
    end
end

