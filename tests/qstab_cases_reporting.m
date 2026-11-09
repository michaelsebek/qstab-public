function cases=qstab_cases_reporting(writeReport)
%QSTAB_CASES_REPORTING Isolated I/O regressions; never alter a real archive.
% The private writer handle is supplied by test_qstab. Each case uses an
% artificial report and a separate temporary package root, then removes it.
    cases={ ...
      'Default reports and fixed ZIP ignore the current MATLAB folder',@()localDefault(writeReport); ...
      'Relative and absolute OutputDir never change the fixed ZIP path',@()localCustom(writeReport); ...
      'Repeated runs replace ZIP entries without touching historical copies',@()localReplace(writeReport); ...
      'Completed failed tests retain FAILED and the same run ID in the ZIP',@()localFailed(writeReport); ...
      'Invalid report paths are rejected instead of selecting pwd',@()localInvalid(writeReport)};
end

function localDefault(writer)
    [base,root,oldPwd]=localSetup();
    tidy=onCleanup(@()localTeardown(oldPwd,base)); %#ok<NASGU>
    cd(fullfile(base,'caller'));
    p=pwd;
    R=writer(localFixture(false),root,'qstab_test_results');
    assert(strcmp(pwd,p),'Report persistence must not change pwd.');
    assert(strcmp(R.archiveFile,fullfile(root,'qstab_test_results.zip')));
    assert(strcmp(R.outputDirectory,fullfile(root,'qstab_test_results')));
    assert(exist(fullfile(p,'qstab_test_results.zip'),'file')~=2);
    assert(exist(fullfile(R.outputDirectory,'qstab_test_results.zip'),'file')~=2);
    localCheckArchive(R);
end

function localCustom(writer)
    [base,root,oldPwd]=localSetup();
    tidy=onCleanup(@()localTeardown(oldPwd,base)); %#ok<NASGU>
    cd(fullfile(base,'caller'));
    R=writer(localFixture(false),root,fullfile('custom raw','results'));
    assert(strcmp(R.outputDirectory,fullfile(root,'custom raw','results')));
    fixed=fullfile(root,'qstab_test_results.zip');
    assert(strcmp(R.archiveFile,fixed));
    localCheckArchive(R);
    raw=fullfile(base,'external raw');
    R=writer(localFixture(false),root,raw);
    assert(strcmp(R.outputDirectory,raw)&&strcmp(R.archiveFile,fixed));
    assert(exist([raw '.zip'],'file')~=2,'OutputDir must not create a second archive name.');
    localCheckArchive(R);
end

function localReplace(writer)
    [base,root,oldPwd]=localSetup();
    tidy=onCleanup(@()localTeardown(oldPwd,base)); %#ok<NASGU>
    cd(fullfile(base,'caller'));
    old=fullfile(root,'qstab_report_0_1_2_fixture.zip');
    nestedDir=fullfile(root,'qstab_test_results');mkdir(nestedDir);
    nested=fullfile(nestedDir,'qstab_test_results.zip');
    localPut(old,'historical-manual-copy');
    localPut(nested,'historical-nested-copy');
    localPut(fullfile(base,'obsolete.txt'),'must not survive the replacement');
    fixed=fullfile(root,'qstab_test_results.zip');
    zip(fixed,{'obsolete.txt'},base);
    R=writer(localFixture(false),root,nestedDir);
    localCheckArchive(R);
    firstRun=R.runId;
    S=localFixture(false);S.runId='fixture-second-run';
    S.toolbox.version='fixture-second-version';S.elapsedSeconds=123;
    R=writer(S,root,nestedDir);
    assert(~strcmp(R.runId,firstRun));
    localCheckArchive(R);
    assert(strcmp(fileread(old),'historical-manual-copy'));
    assert(strcmp(fileread(nested),'historical-nested-copy'));
end

function localFailed(writer)
    [base,root,oldPwd]=localSetup();
    tidy=onCleanup(@()localTeardown(oldPwd,base)); %#ok<NASGU>
    cd(fullfile(base,'caller'));
    R=writer(localFixture(true),root,'qstab_test_results');
    assert(strcmp(R.status,'FAILED')&&R.failed==1&&~R.ok);
    localCheckArchive(R);
    T=fileread(fullfile(R.outputDirectory,'report.txt'));
    assert(~isempty(strfind(T,'FAILED'))); %#ok<STREMP>
    assert(~isempty(strfind(T,R.runId))); %#ok<STREMP>
end

function localInvalid(writer)
    [base,root,oldPwd]=localSetup();
    tidy=onCleanup(@()localTeardown(oldPwd,base)); %#ok<NASGU>
    cd(fullfile(base,'caller'));
    bad={'',[],{'folder'},char('aa','bb'),'   '};
    for k=1:numel(bad)
        rejected=false;
        try
            writer(localFixture(false),root,bad{k});
        catch ME
            rejected=strcmp(ME.identifier,'qstab:TestReportPath');
        end
        assert(rejected,'Invalid OutputDir must fail with the path error.');
    end
    assert(exist(fullfile(root,'qstab_test_results.zip'),'file')~=2);
end

function R=localFixture(failed)
    item=struct('group','fixture','name','Artificial persistence fixture', ...
        'status','PASS','seconds',0,'message','','identifier','');
    status='OVERALL_OK';
    if failed
        status='FAILED';item.status='FAIL';
        item.message='Deliberate artificial failure, not a numerical result.';
    end
    R=struct('toolbox',struct('version','reporting-fixture'), ...
        'environment',struct('version','artificial fixture','computer','fixture', ...
            'quaternionPath','','qstabPath',''), ...
        'timestamp','20000101T000000','runId','fixture-first-run', ...
        'nativeMATLABRun',false,'nativeQuaternionOK',false,'ok',~failed, ...
        'status',status,'passed',double(~failed),'failed',double(failed), ...
        'skipped',0,'elapsedSeconds',0,'cases',item);
end

function localCheckArchive(R)
    target=tempname(R.packageRoot);mkdir(target);
    tidy=onCleanup(@()rmdir(target,'s')); %#ok<NASGU>
    unzip(R.archiveFile,target);
    d=dir(target);d=d(~ismember({d.name},{'.','..'}));
    assert(numel(d)==2&&~any([d.isdir]));
    assert(isequal(sort({d.name}),sort({'report.mat','report.txt'})));
    S=load(fullfile(target,'report.mat'),'report');
    assert(isequaln(S.report,R),'Saved MAT report must equal the returned report.');
    for k=1:numel(R.archiveEntries)
        n=R.archiveEntries{k};
        assert(isequal(localBytes(fullfile(target,n)), ...
                       localBytes(fullfile(R.outputDirectory,n))));
    end
end

function [base,root,oldPwd]=localSetup()
    oldPwd=pwd;base=tempname;mkdir(base);
    [ok,info]=fileattrib(base);assert(ok);base=info.Name;
    root=fullfile(base,'isolated package root');mkdir(root);
    mkdir(fullfile(base,'caller'));
end
function localTeardown(oldPwd,base)
    cd(oldPwd);
    if exist(base,'dir')==7,rmdir(base,'s');end
end
function localPut(name,text)
    fid=fopen(name,'w');assert(fid>=0);
    tidy=onCleanup(@()fclose(fid)); %#ok<NASGU>
    fprintf(fid,'%s',text);
end
function b=localBytes(name)
    fid=fopen(name,'rb');assert(fid>=0);
    tidy=onCleanup(@()fclose(fid)); %#ok<NASGU>
    b=fread(fid,Inf,'*uint8');
end
