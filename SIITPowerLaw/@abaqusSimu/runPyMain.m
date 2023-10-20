function runPyMain(obj)
    % run pymain_jobName in abaqus

    jobName = jobNameGenr(obj);
    rewrtPyMain(obj);

    mkdir(jobName);
    copyfile( 'E:\MyPapers\MyMatlabFuncs\SIITPowerLaw\v3\SIIT.py', jobName);
    copyfile( ['pymain_',jobName,'.py'], jobName,'f');
    cd(jobName);

    system(['abaqus cae noGUI=','pymain_',jobName,'.py']);
    cd ..
    movefile( [jobName,'\Lh_',jobName,'.txt'])
%    movefile( [jobName,'\UnldSrfNdsCrd_',jobName,'.txt'])

    rmdir(jobName,'s');
end