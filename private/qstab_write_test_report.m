function report = qstab_write_test_report(report,packageRoot,outputDir)
%QSTAB_WRITE_TEST_REPORT Internal report persistence with one fixed ZIP path.
%   Called by test_qstab; not a numerical solver or a new public data type.
%   The archive ALWAYS resides beside test_qstab.m and is ALWAYS named
%   qstab_test_results.zip. OutputDir affects only report.mat/report.txt.
%   Tests call this helper with isolated temporary roots via a function
%   handle supplied by test_qstab; no test touches the installed archive.

    packageRoot=localTextPath(packageRoot,'packageRoot');
    [rootOK,rootInfo]=fileattrib(packageRoot);
    if ~rootOK||~rootInfo.directory
        error('qstab:TestReportPath','Package root does not exist: %s',packageRoot);
    end
    packageRoot=rootInfo.Name;
    outputDir=localTextPath(outputDir,'OutputDir');
    if ~localAbsolutePath(outputDir)
        % Drive-relative paths (C:folder) and Windows current-drive-root
        % paths (\folder) are deliberately not resolved through pwd.
        if ~isempty(regexp(outputDir,'^[A-Za-z]:','once')) || ...
                (ispc&&outputDir(1)=='\')
            error('qstab:TestReportPath','Use a fully qualified drive/UNC path or a relative OutputDir.');
        end
        outputDir=fullfile(packageRoot,outputDir);
    end
    if exist(outputDir,'dir')~=7
        [made,msg]=mkdir(outputDir);
        if ~made,error('qstab:TestReport','Cannot create %s: %s',outputDir,msg);end
    end
    [outOK,outInfo]=fileattrib(outputDir);
    if ~outOK||~outInfo.directory
        error('qstab:TestReport','Cannot resolve report directory: %s',outputDir);
    end
    report.packageRoot=packageRoot;
    report.outputDirectory=outInfo.Name;
    report.archiveFile=fullfile(packageRoot,'qstab_test_results.zip');
    report.archivePolicy='fixedPackageRoot';
    report.reportSchemaVersion=1;
    report.archiveEntries={'report.mat','report.txt'};
    if ~isfield(report,'runId')||isempty(report.runId)
        [~,token]=fileparts(tempname);
        report.runId=[datestr(now,30) '_' token];
    end
    save(fullfile(report.outputDirectory,'report.mat'),'report','-v7');
    localWriteText(report); % fclose occurs before compression

    temporaryArchive=[tempname(packageRoot) '.zip'];
    verifyDir=tempname(packageRoot);
    tidy=onCleanup(@()localCleanup(temporaryArchive,verifyDir)); %#ok<NASGU>
    zip(temporaryArchive,report.archiveEntries,report.outputDirectory);
    [made,msg]=mkdir(verifyDir);
    if ~made,error('qstab:TestArchive','Cannot create archive verification directory: %s',msg);end
    unzip(temporaryArchive,verifyDir);
    contents=dir(verifyDir);
    contents=contents(~ismember({contents.name},{'.','..'}));
    if numel(contents)~=2||any([contents.isdir])|| ...
            ~isequal(sort({contents.name}),sort(report.archiveEntries))
        error('qstab:TestArchive','ZIP must contain exactly report.mat and report.txt at its root.');
    end
    for k=1:numel(report.archiveEntries)
        name=report.archiveEntries{k};
        if ~isequal(localBytes(fullfile(verifyDir,name)), ...
                    localBytes(fullfile(report.outputDirectory,name)))
            error('qstab:TestArchive','ZIP read-back differs from the saved %s.',name);
        end
    end
    % Replace rather than append: obsolete entries cannot survive.
    packedBytes=localBytes(temporaryArchive);
    [moved,msg]=movefile(temporaryArchive,report.archiveFile,'f');
    if ~moved
        error('qstab:TestArchive','Cannot replace %s: %s',report.archiveFile,msg);
    end
    if ~isequal(packedBytes,localBytes(report.archiveFile))
        error('qstab:TestArchive','Final archive differs from the verified temporary ZIP.');
    end
end

function p=localTextPath(x,label)
    valid=(ischar(x)&&isrow(x))||(isstring(x)&&isscalar(x));
    if ~valid,error('qstab:TestReportPath','%s must be a nonempty text path.',label);end
    p=char(x);
    if isempty(strtrim(p)),error('qstab:TestReportPath','%s must be a nonempty text path.',label);end
end

function tf=localAbsolutePath(p)
    if ispc
        tf=~isempty(regexp(p,'^[A-Za-z]:[\\/]','once'))|| ...
            (numel(p)>=2&&((p(1)=='\'&&p(2)=='\')||(p(1)=='/'&&p(2)=='/')));
    else
        tf=p(1)=='/';
    end
end

function bytes=localBytes(filename)
    fid=fopen(filename,'rb');
    if fid<0,error('qstab:TestArchive','Cannot read %s.',filename);end
    closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
    bytes=fread(fid,Inf,'*uint8');
end

function localWriteText(report)
    fid=fopen(fullfile(report.outputDirectory,'report.txt'),'w');
    if fid<0,error('qstab:TestReport','Cannot write report.txt.');end
    closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
    env=report.environment;
    fprintf(fid,'qstab %s\n%s\nRuntime: %s\nPlatform: %s\n', ...
        report.toolbox.version,report.status,env.version,env.computer);
    fprintf(fid,'Run ID: %s\nTimestamp: %s\n',report.runId,report.timestamp);
    fprintf(fid,'Native MATLAB run: %d\nNative quaternion tested: %d\n', ...
        report.nativeMATLABRun,report.nativeQuaternionOK);
    fprintf(fid,'quaternion: %s\nqstab.isstable: %s\n',env.quaternionPath,env.qstabPath);
    fprintf(fid,'Reports: %s\nArchive: %s\n',report.outputDirectory,report.archiveFile);
    fprintf(fid,'Archive policy: fixedPackageRoot (send this ZIP only)\n');
    fprintf(fid,'Passed: %d  Failed: %d  Skipped: %d  Time: %.3f s\n\n', ...
        report.passed,report.failed,report.skipped,report.elapsedSeconds);
    for k=1:numel(report.cases)
        e=report.cases(k);
        fprintf(fid,'%s [%s] %s (%.4f s)\n',e.status,e.group,e.name,e.seconds);
        if ~isempty(e.message),fprintf(fid,'  %s\n',e.message);end
    end
end

function localCleanup(archiveName,verifyDir)
    % Only temporary files made by this invocation are ever removed.
    if exist(archiveName,'file')==2,delete(archiveName);end
    if exist(verifyDir,'dir')==7,rmdir(verifyDir,'s');end
end
