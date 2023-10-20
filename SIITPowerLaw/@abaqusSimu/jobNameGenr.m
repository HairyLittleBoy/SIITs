function jobName = jobNameGenr(obj)
    % generate the jobName of FEM simulation
    Ename = num2str(obj.modelParas(1)/1000);
    Rname = num2str(obj.indenterR);
    Yname = num2str(obj.modelParas(2),10);
    miuname = num2str(obj.miu);

    Ename(Ename=='.')='p';
    Rname(Rname=='.')='p';
    Yname(Yname=='.')='p';
    miuname(miuname=='.')='p';


    if length(obj.modelParas) == 3
        Nname = num2str(obj.modelParas(3),10);
        Nname(Nname=='.')='p';
        jobName = ['R',Rname,'E',Ename,'Y',Yname,'N',Nname,'fc',miuname];
    elseif length(obj.modelParas) == 4
        Aname = num2str(obj.modelParas(3),10);
        mname = num2str(obj.modelParas(4),10);
        Aname(Aname=='.')='p';
        mname(mname=='.')='p';
        jobName = ['R',Rname,'E',Ename,'Y',Yname,'A',Aname,'m',mname,'fc',miuname];
    end
end