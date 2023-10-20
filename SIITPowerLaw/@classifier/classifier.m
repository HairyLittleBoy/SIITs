classdef classifier < LhCurve
    %CLASSIFIER  contains the methods which classfies the given L-h curve, provide
    %            the possible range of parameter in the equation describing the stress
    %            plastic strain curve
    %            �����������ĳ�Ա������L-h���߽��з��࣬������Ӧ����Ӧ��Ӧ������
    %            ���̲����Ŀ���ȡֵ��Χ
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
            %METHOD1 �˴���ʾ�йش˷�����ժҪ
            %   �˴���ʾ��ϸ˵��
            outputArg = obj.Property1 + inputArg;
        end
    end
end

