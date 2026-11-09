function E=run_qstab_eigen_comparison(varargin)
%RUN_QSTAB_EIGEN_COMPARISON Reference-backed, read-only right-spectrum study.
% E=run_qstab_eigen_comparison checks qstab.eig on fourteen deterministic
% matrices with analytic reference multisets, preserving algebraic counts.
% E=run_qstab_eigen_comparison('LegacyFunction',@qeig,'Timing',true)
% additionally invokes the EXISTING function explicitly supplied by the user.
% Initialize that function's own environment before calling this runner.
% No qeig alias, path addition, global initialization or source edit occurs.
%
% Options: Cases (IDs or 'all'), UseQuaternion (true), Tolerance (1e-7),
% LegacyFunction ([]), CheckLegacyVectors (true), Timing (false),
% CheckPrimaryVectors (false) explicitly checks qstab.eigpairs as well,
% outside one-output timing. Its assessment is retained in primary.eigenpairs.
% Repetitions (3), WarmupRuns (1), MinSampleTime (.05), MaxBatchCalls (256),
% Verbose (true). With no legacy function, timing comparison is not applicable.
% Existing [V,D]=qeig(A) is checked separately from qstab's [lambda,info].
% A defective fixture does not require a full eigenvector basis.
%
% E.ok requires the selected checks, including requested legacy checks.
% COMPLETED_WITH_DIFFERENCES is a completed experiment, NOT an endorsement of
% a differing/failed legacy solver. All source data and failures are retained.
% Output is ordinary serializable data. This function itself writes no files.
ip=inputParser;
addParameter(ip,'LegacyFunction',[],@(x)isequal(x,[])||isa(x,'function_handle'));
addParameter(ip,'Cases','all',@(x)ischar(x)||isstring(x)||iscell(x));
addParameter(ip,'UseQuaternion',true,@boolean);
addParameter(ip,'Tolerance',1e-7,@positive);
addParameter(ip,'CheckLegacyVectors',true,@boolean);
addParameter(ip,'CheckPrimaryVectors',false,@boolean);
addParameter(ip,'Timing',false,@boolean);
addParameter(ip,'Repetitions',3,@posint);
addParameter(ip,'WarmupRuns',1,@nonnegint);
addParameter(ip,'MinSampleTime',.05,@nonnegative);
addParameter(ip,'MaxBatchCalls',256,@posint);
addParameter(ip,'Verbose',true,@boolean);
parse(ip,varargin{:});o=ip.Results;
C=select(qstab_eigen_cases(),o.Cases); % Reject selections before any solver call.
oldPath=path;oldDir=pwd;oldRng=rng;
cleanup=onCleanup(@()restore(oldPath,oldDir,oldRng)); %#ok<NASGU>
if o.UseQuaternion
    if exist('OCTAVE_VERSION','builtin')~=0
        error('qstab:eigencomparison:Native','UseQuaternion=true requires MathWorks quaternion.');
    end
    q=quaternion(0,0,0,0);[a,b,c,d]=parts(q);
    if ~strcmp(class(q),'quaternion')||~isequal([a b c d],[0 0 0 0])
        error('qstab:eigencomparison:Native','The native quaternion contract is not available.');
    end
end
legacyRequested=isa(o.LegacyFunction,'function_handle');
meta=source(o.LegacyFunction);savedOptions=o;savedOptions.LegacyFunction=meta.label;
% A fieldless struct([]) cannot receive a whole case by indexed assignment
% in MATLAB. Use ONE schema for both the allocated array and every row.
% Nested assessment payloads may still differ on success/error branches.
caseTemplate=struct('id','','primary',[],'legacy',[], ...
    'mutualMatch',[],'classMultisetsAgree',false,'timing', ...
    struct('status','NOT_REQUESTED','ok',true,'eligible',false,'secondOverFirst',NaN));
E=struct('schema','qstab.eigen-comparison.1','toolbox',qstab.about(), ...
    'runtime',version,'computer',computer,'isOctave',exist('OCTAVE_VERSION','builtin')~=0, ...
    'nativeQuaternionInputs',o.UseQuaternion,'options',savedOptions,'catalogue',C, ...
    'primaryEigenpairsRequested',o.CheckPrimaryVectors,'primaryEigenpairsOK',true, ...
    'legacyRequested',legacyRequested,'legacySourceBefore',meta,'legacySourceAfter',[], ...
    'legacySourceUnchanged',[],'cases',repmat(caseTemplate,1,numel(C)), ...
    'completed',false,'ok',false, ...
    'qstabOK',true,'legacyOK',~legacyRequested,'allClassMultisetsAgree',false, ...
    'status','RUNNING','matchedPrimary',0,'matchedLegacy',0,'legacyIssues',0, ...
    'timingEligible',0,'elapsedSeconds',0,'environmentChangedByCalls',[], ...
    'note',['Analytic references are independent of eig; finite fixtures are not an exhaustive proof. ', ...
        'One-output timing does not compare qstab diagnostics with a legacy eigenvector interface.']);
started=tic;
for k=1:numel(C)
    F=C(k);A=inputMatrix(F.matrix,o.UseQuaternion);
    row=caseTemplate;
    row.id=F.id;
    try
        [lambda,info]=qstab.eig(A);
        row.primary=qstab_eigen_assess(A,F.reference,lambda,'Info',info, ...
            'Tolerance',o.Tolerance,'Diagonalizable',F.diagonalizable);
    catch ME
        row.primary=failure(ME);
    end
    if o.CheckPrimaryVectors
        try
            [V,D,pairInfo]=qstab.eigpairs(A);
            if isfield(row.primary,'canonical') && ~isempty(row.primary.canonical)
                e=row.primary.canonical;
            else
                [da,db,dc,dd]=qstab.components(D);
                e=complex(diag(da),hypot(hypot(diag(db),diag(dc)),diag(dd)));
            end
            row.primary.eigenpairs=qstab_eigen_assess(A,F.reference,e, ...
                'Tolerance',o.Tolerance,'CheckVectors',true,'Vectors',V,'Diagonal',D, ...
                'Info',pairInfo,'Diagonalizable',F.diagonalizable);
        catch ME
            row.primary.eigenpairs=failure(ME);
        end
        E.primaryEigenpairsOK=E.primaryEigenpairsOK && row.primary.eigenpairs.ok;
        row.primary.ok=row.primary.ok && row.primary.eigenpairs.ok;
        if ~row.primary.eigenpairs.ok
            row.primary.status='RIGHT_EIGENVECTOR_ISSUE';
            row.primary.message=row.primary.eigenpairs.message;
        end
    end
    if row.primary.ok,E.matchedPrimary=E.matchedPrimary+1;end
    E.qstabOK=E.qstabOK&&row.primary.ok;
    row.legacy=struct('status','NOT_REQUESTED','ok',true,'spectrumAssessment',[], ...
        'vectorAssessment',[],'message','');
    if legacyRequested
        row.legacy.ok=false;
        try
            lambda=o.LegacyFunction(A);
            row.legacy.spectrumAssessment=qstab_eigen_assess(A,F.reference,lambda, ...
                'Tolerance',o.Tolerance,'Diagonalizable',F.diagonalizable);
            row.legacy.ok=row.legacy.spectrumAssessment.ok;
            row.legacy.status=row.legacy.spectrumAssessment.status;
            row.legacy.message=row.legacy.spectrumAssessment.message;
            if o.CheckLegacyVectors
                try
                    [V,D]=o.LegacyFunction(A);
                    row.legacy.vectorAssessment=qstab_eigen_assess(A,F.reference,lambda, ...
                        'Tolerance',o.Tolerance,'CheckVectors',true,'Vectors',V,'Diagonal',D, ...
                        'Diagonalizable',F.diagonalizable);
                catch ME
                    row.legacy.vectorAssessment=failure(ME);
                end
                row.legacy.ok=row.legacy.ok&&row.legacy.vectorAssessment.ok;
                if ~row.legacy.vectorAssessment.ok
                    row.legacy.status='EIGENVECTOR_ISSUE';
                    row.legacy.message=row.legacy.vectorAssessment.message;
                end
            end
            if row.primary.ok&&row.legacy.spectrumAssessment.ok
                row.mutualMatch=qstab_eigen_match(row.primary.canonical, ...
                    row.legacy.spectrumAssessment.canonical,o.Tolerance);
                row.classMultisetsAgree=row.mutualMatch.ok;
            end
        catch ME
            row.legacy.status='ERROR';row.legacy.message=ME.message;
            row.legacy.identifier=ME.identifier;row.legacy.ok=false;
        end
        if row.legacy.ok,E.matchedLegacy=E.matchedLegacy+1;
        else,E.legacyIssues=E.legacyIssues+1;end
    end
    if o.Timing
        if ~legacyRequested
            row.timing.status='NOT_APPLICABLE_NO_LEGACY';
        elseif row.primary.ok&&~isempty(row.legacy.spectrumAssessment)&& ...
                row.legacy.spectrumAssessment.ok&&row.classMultisetsAgree
            row.timing=qstab_eigen_timing(@qstab.eig,o.LegacyFunction,A,F.reference, ...
                'Tolerance',o.Tolerance,'Repetitions',o.Repetitions, ...
                'WarmupRuns',o.WarmupRuns,'MinSampleTime',o.MinSampleTime, ...
                'MaxBatchCalls',o.MaxBatchCalls);
            if row.timing.eligible,E.timingEligible=E.timingEligible+1;end
        else
            row.timing.status='NOT_RUN_UNVERIFIED_SPECTRUM';
        end
    end
    E.cases(k)=row;
    if o.Verbose
        fprintf('  eig %-21s qstab=%-12s legacy=%-19s timing=%s\n', ...
            F.id,row.primary.status,row.legacy.status,row.timing.status);
        if legacyRequested&&~row.legacy.ok
            fprintf('      %s\n',row.legacy.message);
            if ~isempty(row.legacy.spectrumAssessment)
                fprintf('      one-output spectrum: %s\n',row.legacy.spectrumAssessment.status);
            end
            if ~isempty(row.legacy.vectorAssessment)
                fprintf('      eigenvectors/diagonal: %s\n',row.legacy.vectorAssessment.status);
            end
        end
    end
end
E.legacyOK=~legacyRequested||E.matchedLegacy==numel(C);
E.allClassMultisetsAgree=legacyRequested&&all([E.cases.classMultisetsAgree]);
E.legacySourceAfter=source(o.LegacyFunction);
if ~isempty(meta.sha256)&&~isempty(E.legacySourceAfter.sha256)
    E.legacySourceUnchanged=strcmp(meta.sha256,E.legacySourceAfter.sha256);
end
E.environmentChangedByCalls=struct('path',~strcmp(path,oldPath), ...
    'directory',~strcmp(pwd,oldDir),'randomState',~isequaln(rng,oldRng));
E.completed=true;E.ok=E.qstabOK&&E.legacyOK;
for k=1:numel(E.cases),E.ok=E.ok&&E.cases(k).timing.ok;end
if ~isempty(E.legacySourceUnchanged),E.ok=E.ok&&E.legacySourceUnchanged;end
E.ok=E.ok&&~any(structfun(@(x)x,E.environmentChangedByCalls));
if E.ok,E.status='COMPLETED';
elseif ~E.qstabOK,E.status='QSTAB_ISSUES';
else,E.status='COMPLETED_WITH_DIFFERENCES';end
E.elapsedSeconds=toc(started);
if o.Verbose
    fprintf('qstab eigen comparison: %s | qstab %d/%d',E.status,E.matchedPrimary,numel(C));
    if legacyRequested,fprintf(' | legacy %d/%d | %d timing pairs',E.matchedLegacy,numel(C),E.timingEligible);
    else,fprintf(' | existing qeig NOT REQUESTED');end
    fprintf('\n');
end
end
function C=select(C,value)
if ischar(value)
    if ~isrow(value),error('qstab:eigencomparison:Cases','Cases must be text IDs.');end
    ids={value};
elseif isstring(value),ids=cellstr(value(:));
elseif iscell(value),ids=value;
else,error('qstab:eigencomparison:Cases','Cases must be text IDs.');end
if isempty(ids)||~all(cellfun(@(x)ischar(x)&&isrow(x)&&~isempty(x),ids))
    error('qstab:eigencomparison:Cases','Cases must be a nonempty list of IDs.');
end
if isequal(ids,{'all'}),return;end
if numel(unique(ids))~=numel(ids)||~all(ismember(ids,{C.id}))
    error('qstab:eigencomparison:Cases','Unknown or duplicate case ID.');
end
C=C(ismember({C.id},ids));
end
function A=inputMatrix(C,native)
if native,A=quaternion(C.a,C.b,C.c,C.d);else,A=C;end
end
function E=failure(ME)
E=struct('status','ERROR','ok',false,'message',ME.message,'identifier',ME.identifier);
end
function S=source(f)
S=struct('label','','file','','sha256','','hashingAvailable',false,'message','');
if ~isa(f,'function_handle'),S.message='External solver not requested.';return;end
S.label=func2str(f);
try
    % functions() is used for diagnostic metadata ONLY. Its workspace is
    % deliberately NOT copied, serialized or used to dispatch algorithms.
    d=functions(f);if isfield(d,'file')&&ischar(d.file),S.file=d.file;end
catch ME
    S.message=ME.message;
end
if ~isempty(S.file)&&exist(S.file,'file')==2
    fid=fopen(S.file,'rb');
    if fid>=0
        closeFile=onCleanup(@()fclose(fid)); %#ok<NASGU>
        raw=fread(fid,Inf,'*uint8');
        try
            md=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
            md.update(typecast(raw(:),'int8'));digest=typecast(md.digest(),'uint8');
            S.sha256=lower(reshape(dec2hex(digest,2).',1,[]));S.hashingAvailable=true;
        catch ME
            S.message=['Source hash unavailable: ' ME.message];
        end
    end
end
end
function restore(p,d,r)
if ~strcmp(path,p),path(p);end
if ~strcmp(pwd,d),cd(d);end
if ~isequaln(rng,r),rng(r);end
end
function tf=boolean(x),tf=islogical(x)&&isscalar(x);end
function tf=positive(x),tf=nonnegative(x)&&x>0;end
function tf=posint(x),tf=positive(x)&&x==floor(x);end
function tf=nonnegint(x),tf=nonnegative(x)&&x==floor(x);end
function tf=nonnegative(x),tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0;end
