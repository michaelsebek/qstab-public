function report = test_qstab(varargin)
%TEST_QSTAB Run the complete qstab regression and native integration suite.
%   REPORT=test_qstab runs core, MathWorks quaternion, and hidden-figure
%   tests. Results are saved to qstab_test_results/report.mat and report.txt.
%   Always send <package root>/qstab_test_results.zip. Its absolute path
%   is report.archiveFile. The ZIP is recreated and read back on each run.
%   OutputDir changes ONLY the raw report location, never the ZIP name/path.
%   A relative OutputDir is resolved against this package root, not pwd.
%   No personal quaternion or polynomial toolbox is required.
%
%   Options: RequireQuaternion=true, Plots=true, Verbose=true, OutputDir.
%   EigenComparison=false stores a reference-backed right-spectrum study.
%   EigenOptions=struct() can explicitly supply LegacyFunction=@qeig; no
%   private environment is initialized. Differences stay visible and cause
%   EIGEN_COMPARISON_ISSUES, without relabelling regression groups.
%   PaperExtensions=false stores 16 additional finite source checks.
%   PaperFigures=false computes eight figure datasets and exports PDF/PNG/FIG/CSV
%   to <package root>/qstab_paper_figures. Both payloads stay in report.mat.
%   PaperExamples=false, PaperOptions=struct() optionally store selected
%   source-linked qB01/qB02/qB03 cases as report.paperExamples.
%   BoxExamples=false optionally stores component-box radius and gauge examples.
%   DiscreteExamples=false optionally stores z/zinv one-parameter margin examples.
%   DiagnosticExamples=false optionally saves read-only diagnostic example data.
%   Benchmarks=false (default), BenchmarkOptions=struct() optionally run the
%   matched method study AFTER passing tests and BEFORE the same ZIP write.
%   Full data are report.benchmarks inside report.mat, not a third ZIP entry.
%   For auxiliary numeric execution without MATLAB quaternion:
%     test_qstab('RequireQuaternion',false,'Plots',false)
%   Such a run is explicitly labelled CORE_OK_NATIVE_PENDING, not a native
%   MATLAB validation. No clear classes/functions or startup changes occur.

    % Optional BallExamples retains full coefficient-ball data in report.mat.
    root=fileparts(mfilename('fullpath'));
    ip=inputParser;
    addParameter(ip,'RequireQuaternion',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'Plots',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'EigenComparison',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'EigenOptions',struct(),@(x)isstruct(x)&&isscalar(x));
    addParameter(ip,'PaperExtensions',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PaperFigures',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PaperExamples',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'PaperOptions',struct(),@(x)isstruct(x)&&isscalar(x));
    addParameter(ip,'BoxExamples',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'BallExamples',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'DiscreteExamples',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'DiagnosticExamples',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'Benchmarks',false,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'BenchmarkOptions',struct(),@(x)isstruct(x)&&isscalar(x));
    addParameter(ip,'OutputDir',fullfile(root,'qstab_test_results'),...
        @(x)ischar(x)||(isstring(x)&&isscalar(x)));
    parse(ip,varargin{:}); o=ip.Results;
    if o.PaperFigures && ~o.Plots
        error('qstab:TestOptions','PaperFigures=true requires Plots=true.');
    end
    oldPath=path; restorePath=onCleanup(@()path(oldPath)); %#ok<NASGU>
    addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
    addpath(fullfile(root,'examples'));
    % No private q-aliases are installed, enabled or tested by the public suite.
    wrapperConflicts=struct('name',{},'locations',{});wrappersReady=true;
    env=struct('version',version,'computer',computer,'isOctave',exist('OCTAVE_VERSION','builtin')~=0,...
        'quaternionPath',which('quaternion'),'qstabPath',which('qstab.isstable'));
    env.wrapperConflicts=wrapperConflicts;
    nativeReady=false; nativeReason='';
    if ~env.isOctave
        try
            q=quaternion(0,0,0,0);
            [qa,qb,qc,qd]=parts(q);
            nativeReady=strcmp(class(q),'quaternion')&&isequal([qa qb qc qd],[0 0 0 0]);
            if ~nativeReady, nativeReason='The quaternion constructor/parts contract is not the expected native API.'; end
        catch ME
            nativeReason=ME.message;
        end
    else
        nativeReason='Auxiliary Octave run; MathWorks quaternion integration is not exercised.';
    end
    [~,runToken]=fileparts(tempname);
    runId=[datestr(now,30) '_' runToken];
    reportOptions=o;
    if isfield(reportOptions.EigenOptions,'LegacyFunction') && ...
            isa(reportOptions.EigenOptions.LegacyFunction,'function_handle')
        % Never persist the supplied handle or its captured workspace.
        reportOptions.EigenOptions.LegacyFunction=func2str(reportOptions.EigenOptions.LegacyFunction);
    end
    report=struct('toolbox',qstab.about(),'environment',env,'timestamp',datestr(now,30),...
        'runId',runId,'nativeMATLABRun',~env.isOctave,'nativeQuaternionAvailable',nativeReady,...
        'nativeQuaternionReason',nativeReason,'options',reportOptions,'ok',false,'status','RUNNING',...
        'passed',0,'failed',0,'skipped',0,'elapsedSeconds',0,...
        'cases',struct('group',{},'name',{},'status',{},'seconds',{},'message',{},'identifier',{}));
    suites={ 'resultfirst',qstab_cases_resultfirst();...
             'families',test_qstab_families_public();...
             'inputerrors',qstab_cases_friendly();...
             'familyplots',test_qstab_exploration_public('numeric');...
             'guide',{'All public guide keys run headlessly and preserve caller state',@guideCases};...
             'nominal',qstab_cases_nominal();...
             'routhsingular',qstab_cases_routhsingular();...
             'robust',qstab_cases_robust();...
             'rays',qstab_cases_rays();...
             'radius',qstab_cases_radius();...
             'values',qstab_cases_values();...
             'zeroexclude',qstab_cases_zeroexclude();...
             'fourgenerator',qstab_cases_fourgenerator();...
             'leading',qstab_cases_leading();...
             'affinezero',qstab_cases_affinezero();...
             'intervalreductions',qstab_cases_intervalreductions();...
             'reporting',qstab_cases_reporting(@qstab_write_test_report);...
             'benchmarking',qstab_cases_benchmarks(@qstab_write_test_report);...
             'diagnostics',qstab_cases_diagnostics(@qstab_write_test_report);...
             'balls',qstab_cases_balls(@qstab_write_test_report);...
             'ballfix',qstab_cases_ballfix();...
             'discrete',qstab_cases_discrete(@qstab_write_test_report);...
             'boxes',qstab_cases_boxes(@qstab_write_test_report);...
             'boxfix',qstab_cases_boxfix();...
             'papers',qstab_cases_papers(@qstab_write_test_report);...
             'paperfix',qstab_cases_paperfix(@qstab_write_test_report);...
             'paperextensions',qstab_cases_paperextensions(@qstab_write_test_report);...
             'eigencomparison',qstab_cases_eigencomparison(@qstab_write_test_report);...
             'eigenfix',qstab_cases_eigenfix(@qstab_write_test_report);...
             'rightpairs',qstab_cases_rightpairs(@qstab_write_test_report);...
             'publicacceptance',qstab_cases_public(@qstab_public_manifest,@qstab_public_path_plan,@qstab_public_session,@qstab_write_test_report,root);...
             'manifestfix',qstab_cases_manifestfix(@qstab_public_manifest);...
             'publicartifacts',qstab_cases_public_artifacts(@qstab_public_figures,@qstab_write_test_report);...
             'provenance',qstab_cases_provenance(@qstab_source_manifest,root);...
             'namespace',qstab_cases_namespace(root);...
             'distribution',{'Public distribution has no qfiles aliases or private object dependencies',@publicSeparation}};
    started=tic;
    for s=1:size(suites,1)
        if strcmp(suites{s,1},'wrappers') && ~wrappersReady
            runCases(suites{s,1},suites{s,2},false, ...
                ['Foreign wrapper names on path: ',strjoin({wrapperConflicts.name},', '), ...
                 '. Not shadowed; run this suite in the standalone environment.']);
        else
            runCases(suites{s,1},suites{s,2},true,'');
        end
    end
    coreFailed=report.failed;
    runCases('native',[qstab_cases_native();qstab_cases_diagnostics_native();qstab_cases_balls_native();qstab_cases_ballfix_native();qstab_cases_discrete_native();qstab_cases_boxes_native();qstab_cases_boxfix_native();qstab_cases_papers_native();qstab_cases_paperfix_native();qstab_cases_paperextensions_native();qstab_cases_eigencomparison_native();qstab_cases_eigenfix_native();qstab_cases_rightpairs_native();test_qstab_exploration_public('native');qstab_cases_friendly('native')],nativeReady,nativeReason);
    if o.Plots
        runCases('graphics',[qstab_cases_graphics();qstab_cases_diagnostic_graphics();qstab_cases_paperfigure_graphics();test_qstab_exploration_public('graphics')],true,'');
    else
        runCases('graphics',[qstab_cases_graphics();qstab_cases_diagnostic_graphics();qstab_cases_paperfigure_graphics();test_qstab_exploration_public('graphics')],false,'Graphics tests explicitly disabled.');
    end
    report.elapsedSeconds=toc(started);
    report.coreOK=coreFailed==0;
    report.nativeQuaternionOK=nativeReady&&~any(strcmp({report.cases.group},'native') & strcmp({report.cases.status},'FAIL'));
    report.wrapperRequirementsSatisfied=wrappersReady;
    report.requirementsSatisfied=(~o.RequireQuaternion||nativeReady)&&wrappersReady;
    report.ok=report.failed==0&&report.requirementsSatisfied;
    if report.failed>0
        report.status='FAILED';
    elseif ~wrappersReady
        report.status='INCOMPLETE_WRAPPER_CONFLICT';
    elseif ~report.requirementsSatisfied
        report.status='INCOMPLETE_NATIVE_QUATERNION';
    elseif env.isOctave||~nativeReady
        report.status='CORE_OK_NATIVE_PENDING';
    elseif ~o.Plots
        report.status='NATIVE_CORE_OK_GRAPHICS_PENDING';
    else
        report.status='OVERALL_OK';
    end
    report.benchmarks=struct('status','NOT_REQUESTED','ok',true, ...
        'sourceRunId',runId,'message','Optional study not requested.');
    if o.Benchmarks
        if report.failed==0 && wrappersReady
            try
                names=fieldnames(o.BenchmarkOptions);args=cell(1,2*numel(names));
                for k=1:numel(names)
                    args{2*k-1}=names{k};args{2*k}=o.BenchmarkOptions.(names{k});
                end
                report.benchmarks=run_qstab_benchmarks(args{:});
                report.benchmarks.sourceRunId=runId;
                if ~report.benchmarks.ok
                    report.ok=false;report.status='BENCHMARK_ISSUES';
                end
            catch ME
                report.benchmarks=struct('status','ERROR','ok',false, ...
                    'sourceRunId',runId,'identifier',ME.identifier,'message',ME.message);
                report.ok=false;report.status='BENCHMARK_ERROR';
            end
        else
            report.benchmarks.status='SKIPPED_REGRESSION_FAILURE';
            report.benchmarks.ok=false;
            report.benchmarks.message='Resolve regression/path failures before a method comparison.';
        end
    end
    report.diagnosticExamples=struct('status','NOT_REQUESTED','ok',true,'sourceRunId',runId);
    if o.DiagnosticExamples
        if report.failed==0 && wrappersReady
            try
                report.diagnosticExamples=run_qstab_diagnostics('Plots',false,'Verbose',o.Verbose);
                report.diagnosticExamples.sourceRunId=runId;
                report.diagnosticExamples.status='COMPLETED';
                if ~report.diagnosticExamples.ok
                    report.ok=false;report.status='DIAGNOSTIC_ISSUES';
                end
            catch ME
                report.diagnosticExamples=struct('status','ERROR','ok',false, ...
                    'sourceRunId',runId,'identifier',ME.identifier,'message',ME.message);
                report.ok=false;report.status='DIAGNOSTIC_ERROR';
            end
        else
            report.diagnosticExamples.status='SKIPPED_REGRESSION_FAILURE';
            report.diagnosticExamples.ok=false;
        end
    end
    report.coefficientBallExamples=struct('status','NOT_REQUESTED','ok',true,'sourceRunId',runId);
    if o.BallExamples
        if report.failed==0 && wrappersReady
            try
                report.coefficientBallExamples=run_qb02_balls('Verbose',o.Verbose);
                report.coefficientBallExamples.sourceRunId=runId;
                report.coefficientBallExamples.status='COMPLETED';
                if ~report.coefficientBallExamples.ok
                    report.ok=false;report.status='COEFFICIENT_BALL_ISSUES';
                end
            catch ME
                report.coefficientBallExamples=struct('status','ERROR','ok',false, ...
                    'sourceRunId',runId,'identifier',ME.identifier,'message',ME.message);
                report.ok=false;report.status='COEFFICIENT_BALL_ERROR';
            end
        else
            report.coefficientBallExamples.status='SKIPPED_REGRESSION_FAILURE';
            report.coefficientBallExamples.ok=false;
        end
    end
    report.discreteMarginExamples=struct('status','NOT_REQUESTED','ok',true,'sourceRunId',runId);
    if o.DiscreteExamples
        if report.failed==0 && wrappersReady
            try
                report.discreteMarginExamples=run_qb01_discrete('Verbose',o.Verbose);
                report.discreteMarginExamples.sourceRunId=runId;
                report.discreteMarginExamples.status='COMPLETED';
                if ~report.discreteMarginExamples.ok
                    report.ok=false;report.status='DISCRETE_MARGIN_ISSUES';
                end
            catch ME
                report.discreteMarginExamples=struct('status','ERROR','ok',false, ...
                    'sourceRunId',runId,'identifier',ME.identifier,'message',ME.message);
                report.ok=false;report.status='DISCRETE_MARGIN_ERROR';
            end
        else
            report.discreteMarginExamples.status='SKIPPED_REGRESSION_FAILURE';
            report.discreteMarginExamples.ok=false;
        end
    end
    report.componentBoxExamples=struct('status','NOT_REQUESTED','ok',true,'sourceRunId',runId);
    if o.BoxExamples
        if report.failed==0 && wrappersReady
            try
                report.componentBoxExamples=run_qb02_boxes('Verbose',o.Verbose);
                report.componentBoxExamples.sourceRunId=runId;
                report.componentBoxExamples.status='COMPLETED';
                if ~report.componentBoxExamples.ok
                    report.ok=false;report.status='COMPONENT_BOX_ISSUES';
                end
            catch ME
                report.componentBoxExamples=struct('status','ERROR','ok',false, ...
                    'sourceRunId',runId,'identifier',ME.identifier,'message',ME.message);
                report.ok=false;report.status='COMPONENT_BOX_ERROR';
            end
        else
            report.componentBoxExamples.status='SKIPPED_REGRESSION_FAILURE';
            report.componentBoxExamples.ok=false;
        end
    end
    report.paperExamples=struct('status','NOT_REQUESTED','ok',true,'sourceRunId',runId);
    if o.PaperExamples
        if report.failed==0 && wrappersReady
            try
                paperOptions=o.PaperOptions;
                if ~isfield(paperOptions,'UseQuaternion'),paperOptions.UseQuaternion=nativeReady;end
                if ~isfield(paperOptions,'Verbose'),paperOptions.Verbose=o.Verbose;end
                paperOptionNames=fieldnames(paperOptions);paperArgs=cell(1,2*numel(paperOptionNames));
                for paperIndex=1:numel(paperOptionNames)
                    paperArgs{2*paperIndex-1}=paperOptionNames{paperIndex};
                    paperArgs{2*paperIndex}=paperOptions.(paperOptionNames{paperIndex});
                end
                report.paperExamples=run_qstab_papers(paperArgs{:});
                report.paperExamples.sourceRunId=runId;
                if ~report.paperExamples.ok
                    report.ok=false;report.status='PAPER_REFERENCE_ISSUES';
                end
            catch ME
                report.paperExamples=struct('status','ERROR','ok',false, ...
                    'sourceRunId',runId,'identifier',ME.identifier,'message',ME.message);
                report.ok=false;report.status='PAPER_REFERENCE_ERROR';
            end
        else
            report.paperExamples.status='SKIPPED_REGRESSION_FAILURE';
            report.paperExamples.ok=false;
        end
    end
    report.paperExtensions=struct('status','NOT_REQUESTED','ok',true,'sourceRunId',runId);
    if o.PaperExtensions
        if report.failed==0 && wrappersReady
            try
                report.paperExtensions=run_qstab_paper_extensions('UseQuaternion',nativeReady,'Verbose',o.Verbose);
                report.paperExtensions.sourceRunId=runId;
                if ~report.paperExtensions.ok
                    report.ok=false;report.status='PAPER_EXTENSION_ISSUES';
                end
            catch ME
                report.paperExtensions=struct('status','ERROR','ok',false,'sourceRunId',runId, ...
                    'identifier',ME.identifier,'message',ME.message);
                report.ok=false;report.status='PAPER_EXTENSION_ERROR';
            end
        else
            report.paperExtensions.status='SKIPPED_REGRESSION_FAILURE';report.paperExtensions.ok=false;
        end
    end
    report.paperFigures=struct('status','NOT_REQUESTED','ok',true,'sourceRunId',runId);
    if o.PaperFigures
        if report.failed==0 && wrappersReady
            try
                report.paperFigures=qstab_paper_figures('UseQuaternion',nativeReady);
                report.paperFigures.sourceRunId=runId;
                if report.paperFigures.ok
                    report.paperFigures.export=export_qstab_paper_figures(report.paperFigures);
                    report.paperFigures.ok=report.paperFigures.export.ok;
                end
                if ~report.paperFigures.ok
                    report.ok=false;report.status='PAPER_FIGURE_ISSUES';
                end
            catch ME
                % Keep computed ordinary data, if any, when an export fails.
                report.paperFigures.status='ERROR';report.paperFigures.ok=false;
                report.paperFigures.sourceRunId=runId;
                report.paperFigures.identifier=ME.identifier;report.paperFigures.message=ME.message;
                report.ok=false;report.status='PAPER_FIGURE_ERROR';
            end
        else
            report.paperFigures.status='SKIPPED_REGRESSION_FAILURE';report.paperFigures.ok=false;
        end
    end
    report.eigenComparison=struct('status','NOT_REQUESTED','ok',true,'sourceRunId',runId);
    if o.EigenComparison
        if report.failed==0 && wrappersReady
            try
                settings=o.EigenOptions;
                if ~isfield(settings,'UseQuaternion'),settings.UseQuaternion=nativeReady;end
                if ~isfield(settings,'Verbose'),settings.Verbose=o.Verbose;end
                names=fieldnames(settings);args=cell(1,2*numel(names));
                for k=1:numel(names)
                    args{2*k-1}=names{k};args{2*k}=settings.(names{k});
                end
                report.eigenComparison=run_qstab_eigen_comparison(args{:});
                report.eigenComparison.sourceRunId=runId;
                if ~report.eigenComparison.ok
                    report.ok=false;report.status='EIGEN_COMPARISON_ISSUES';
                end
            catch ME
                report.eigenComparison=struct('status','ERROR','ok',false,'sourceRunId',runId, ...
                    'identifier',ME.identifier,'message',ME.message);
                report.ok=false;report.status='EIGEN_COMPARISON_ERROR';
            end
        else
            report.eigenComparison.status='SKIPPED_REGRESSION_FAILURE';report.eigenComparison.ok=false;
        end
    end
    % elapsedSeconds retains its original meaning: regression-suite time.
    % totalExecutionSeconds additionally includes optional study/assessment.
    report.totalExecutionSeconds=toc(started);
    % Embed fingerprints inside the existing MAT report; ZIP name/location
    % and the two-file archive contract remain unchanged.
    report.sourceManifest=qstab_source_manifest(root);
    report=qstab_write_test_report(report,root,o.OutputDir);
    fprintf('\nqstab: %s | %d passed, %d failed, %d skipped | %.2f s\n',...
        report.status,report.passed,report.failed,report.skipped,report.elapsedSeconds);
    if o.Benchmarks
        fprintf('Benchmark: %s (full data in report.mat); total execution %.2f s\n', ...
            report.benchmarks.status,report.totalExecutionSeconds);
    end
    if o.DiagnosticExamples
        fprintf('Diagnostics: %s (ordinary data in report.mat)\n',report.diagnosticExamples.status);
    end
    if o.BallExamples
        fprintf('Coefficient balls: %s (ordinary data in report.mat)\n',report.coefficientBallExamples.status);
    end
    if o.DiscreteExamples
        fprintf('Discrete margins: %s (ordinary data in report.mat)\n',report.discreteMarginExamples.status);
    end
    if o.BoxExamples
        fprintf('Component boxes: %s (ordinary data in report.mat)\n',report.componentBoxExamples.status);
    end
    if o.PaperExamples
        fprintf('Selected paper examples: %s (ordinary data in report.mat)\n',report.paperExamples.status);
    end
    if o.PaperExtensions
        fprintf('Additional paper checks: %s (ordinary data in report.mat)\n',report.paperExtensions.status);
    end
    if o.PaperFigures
        fprintf('Selected paper figures: %s (data and export metadata in report.mat)\n',report.paperFigures.status);
        if isfield(report.paperFigures,'export')
            fprintf('Figure files: %s\n',report.paperFigures.export.outputDirectory);
        end
    end
    if o.EigenComparison
        fprintf('Right-spectrum comparison: %s (ordinary data in report.mat)\n',report.eigenComparison.status);
    end
    fprintf('Report: %s\n',fullfile(report.outputDirectory,'report.txt'));
    fprintf('Archive: %s\n',report.archiveFile);
    fprintf('SEND THIS ZIP: qstab %s | %s | run %s\n%s\n', ...
        report.toolbox.version,report.status,report.runId,report.archiveFile);
    % Old manual copies are never used, updated or silently deleted.
    oldCopies=dir(fullfile(root,'qstab_report_*.zip'));
    oldNested=fullfile(root,'qstab_test_results','qstab_test_results.zip');
    if ~isempty(oldCopies)||exist(oldNested,'file')==2
        fprintf(['NOTE: historical manual ZIP(s) were left unchanged. ', ...
            'Only the archive printed above belongs to this run.\n']);
    end

    function runCases(group,cases,enabled,reason)
        for c=1:size(cases,1)
            e=struct('group',group,'name',cases{c,1},'status','SKIP',...
                'seconds',0,'message',reason,'identifier','');
            if enabled
                tc=tic;
                try
                    cases{c,2}(); e.status='PASS'; e.message='';
                    report.passed=report.passed+1;
                catch ME
                    e.status='FAIL'; e.message=ME.message; e.identifier=ME.identifier;
                    if ~isempty(ME.stack)
                        e.message=sprintf('%s [%s:%d]',e.message,ME.stack(1).name,ME.stack(1).line);
                    end
                    report.failed=report.failed+1;
                end
                e.seconds=toc(tc);
            else
                report.skipped=report.skipped+1;
            end
            report.cases(end+1)=e;
            if o.Verbose
                fprintf('%-4s [%s] %s\n',e.status,e.group,e.name);
                if strcmp(e.status,'FAIL'), fprintf('     %s\n',e.message); end
            end
        end
    end
end

function guideCases()
before=rng; oldPath=path; oldDir=pwd; oldFigures=findall(0,'Type','figure');
guide=qstab_guide_example('all');
assert(iscell(guide)&&numel(guide)==11&&all(cellfun(@(entry) isstruct(entry)&&entry.ok,guide)));
assert(isequaln(before,rng)&&strcmp(oldPath,path)&&strcmp(oldDir,pwd));
assert(isequal(oldFigures,findall(0,'Type','figure')));
end

function publicSeparation()
    root=fileparts(mfilename('fullpath'));
    assert(exist(fullfile(root,'qwrappers'),'dir')~=7);
    names={'qeig','qroots','qisstable','qrobuststable','qraymargin','qstabilityradius'};
    for k=1:numel(names)
        assert(exist(fullfile(root,[names{k} '.m']),'file')~=2);
    end
    assert(exist(fullfile(root,'src','+qstab','isstable.m'),'file')==2);
    [~,A] = qstab.isstable([2 3 1]);assert(strcmp(A.status,'stable'));
    [~,B] = qstab.raymargin([-.5 1],1,1,'Variable','z');assert(max(abs(B.interval-[-.5 1.5]))<1e-8);
end
