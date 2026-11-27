function report=run_qstab_oneparameter(varargin)
%RUN_QSTAB_ONEPARAMETER Reproduce the scalar one-parameter article examples.
%   R=run_qstab_oneparameter runs E1--E7 and the directional diagnostics.
%   R=run_qstab_oneparameter('E1') runs the calibration example only.
%   R=run_qstab_oneparameter('Cases',{'E1','E5','E5V'}) selects examples.
%   C=run_qstab_oneparameter('list') lists all example keys.
%
%   Options (name/value pairs):
%     Cases             'all' or a cell of catalogue keys.
%     UseQuaternion     true: use MATLAB quaternion; false: components.
%     ExactCertificates 'auto' (default), true or false. Auto executes the
%                       optional exact certificates when Symbolic Math
%                       Toolbox is usable. Otherwise its absence is recorded.
%                       No Python process is started by this MATLAB runner.
%     Plots             true: save the E1 curves and E3 gain geometry.
%     Visible           'off' (default) or 'on'. Existing figures are kept.
%     MaxTime           5 seconds per numerical global-radius calculation.
%     MaxCells          1024 per global-radius calculation.
%     AbsTol, RelTol    1e-4 for numerical global-radius brackets.
%     Save              true: MAT/TXT/CSV and an automatic new result ZIP.
%     Verbose           true: print example summaries and result locations.
%     OutputDir         []: <package>/oneparameter_results/<unique run>.
%
%   The fixed ZIP is <package>/qstab_oneparameter_results.zip; its absolute
%   path is R.archiveFile. Errors within individual examples are retained
%   in the report and do not prevent the final archive from being written.
%   An invalid option or unavailable requested quaternion class is a
%   preflight error and creates no report.
%
%   R.ok reports that the executed checks passed. R.articleComplete is a
%   stricter field: all examples and all required exact proofs must have run.
%   A numerical bound, an unresolved ray contact, or an archived certificate
%   is never silently promoted to a newly verified exact result. Timing and
%   SVD witness errors are diagnostics, not portable reference constants.
%   E7 uses an explicit amplitude/frequency normalization, recorded in R.
%
%   Example:
%     R=run_qstab_oneparameter('all','Plots',true);
%     R.status
%     R.archiveFile
%     assert(R.ok,'Inspect the report: %s',R.archiveFile);
%
%   See also qstab_oneparameter_cases, qstab.freqmargin, qstab.raymargin,
%   qstab.stabilityradius, qstab.radiuscertificate, qstab.raycertificate.
root=fileparts(fileparts(mfilename('fullpath')));oldpath=path;
restore=onCleanup(@()path(oldpath)); %#ok<NASGU>
addpath(root,fullfile(root,'src'),fullfile(root,'examples'));
catalogue=qstab_oneparameter_cases();selection='all';
if ~isempty(varargin)
 first=varargin{1};if isstring(first)&&isscalar(first),first=char(first);end
 if ischar(first)&&strcmpi(first,'list')
  if numel(varargin)>1,error('qstab:oneparameter:Options','The list command takes no options.');end
  report=catalogue;for k=1:numel(report),fprintf('%-4s %s\n',report(k).id,report(k).title);end;return
 elseif iscell(first)||(ischar(first)&&(strcmpi(first,'all')||ismember(first,{catalogue.id})))
  selection=first;varargin(1)=[];
 end
end
ip=inputParser;ip.FunctionName='run_qstab_oneparameter';ip.PartialMatching=false;
addParameter(ip,'Cases',selection,@(x)ischar(x)||isstring(x)||iscell(x));
addParameter(ip,'UseQuaternion',true,@flag);addParameter(ip,'ExactCertificates','auto',@exactOption);
addParameter(ip,'Plots',true,@flag);addParameter(ip,'Visible','off',@(x)ischar(x)&&any(strcmp(x,{'on','off'})));
addParameter(ip,'MaxTime',5,@positive);addParameter(ip,'MaxCells',1024,@(x)positive(x)&&isfinite(x)&&x==fix(x));
addParameter(ip,'AbsTol',1e-4,@positiveFinite);addParameter(ip,'RelTol',1e-4,@positiveFinite);
addParameter(ip,'Save',true,@flag);addParameter(ip,'Verbose',true,@flag);
addParameter(ip,'OutputDir','',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});o=ip.Results;selection=o.Cases;
if isstring(selection)&&isscalar(selection),selection=char(selection);
elseif isstring(selection),selection=cellstr(selection);end
if ischar(selection),if strcmpi(selection,'all'),selection={catalogue.id};else,selection={selection};end;end
if isempty(selection)||~iscell(selection)||~all(cellfun(@ischar,selection))|| ...
 any(~ismember(selection,{catalogue.id}))||numel(unique(selection))~=numel(selection)
 error('qstab:oneparameter:Selection','Use unique keys from run_qstab_oneparameter(''list'').');
end
C=catalogue(ismember({catalogue.id},selection));
if o.Plots&&~o.Save,error('qstab:oneparameter:Options','Plots requires Save=true.');end
if o.UseQuaternion
 try,q=quaternion(0,0,0,0);[a,b,c,d]=parts(q); %#ok<ASGLU>
 catch err,error('qstab:oneparameter:Quaternion','Native quaternion unavailable. Use UseQuaternion=false for a component-only run. %s',err.message);end
end
dataDir=fullfile(root,'reproducibility','qb01');
source=jsondecode(fileread(fullfile(dataDir,'source_provenance.json')));
[~,token]=fileparts(tempname);runId=[datestr(now,30) '_' token];out=char(o.OutputDir);
if isempty(out),out=fullfile(root,'oneparameter_results');elseif ~absolute(out),out=fullfile(root,out);end
out=fullfile(out,runId);[runExact,exactState]=exactAvailable(o.ExactCertificates);
try
 installation=qstab_check_installation('PackageRoot',root,'CheckResolution',false);
 sourceManifest=installation.inventory;installationIntegrity=installation.integrity;
catch err
 sourceManifest=struct('status','UNAVAILABLE','message',err.message);installationIntegrity=struct('ok',false,'status','UNAVAILABLE');
end
report=struct('schema','qstab.oneparameter.study.1','toolbox',qstab.about(), ...
 'runId',runId,'runtime',version,'computer',computer,'isOctave',exist('OCTAVE_VERSION','builtin')~=0, ...
 'nativeQuaternionInputs',o.UseQuaternion,'options',o,'source',source,'sourceCatalogue',C, ...
 'sourceHashes',sourceHashes(dataDir),'sourceManifest',sourceManifest,'installationIntegrity',installationIntegrity, ...
 'resolvedFunctions',resolvedFunctions(),'ok',false,'articleComplete',false,'numericalAllResolved',true, ...
 'status','RUNNING','passed',0,'failed',0,'checksPassed',0,'checksFailed',0,'elapsedSeconds',0, ...
 'outputDirectory',out,'archiveFile',fullfile(root,'qstab_oneparameter_results.zip'), ...
 'cases',struct('id',{},'title',{},'status',{},'output',{},'checks',{},'seconds',{},'error',{}), ...
 'exact',struct('requested',o.ExactCertificates,'executed',runExact,'status',exactState, ...
 'ok',false,'radiusCertificates',0,'rayCertificates',0,'rayLiftingCertificates',0,'identityChecks',0,'identities',struct()), ...
 'graphics',struct('executed',false,'files',{{}},'error',''), ...
 'assurance',['Numerical comparisons, fresh exact proofs and archived reference certificates are separate. ' ...
 'No native MATLAB validation is inferred from an auxiliary execution.']);
started=tic;
for k=1:numel(C)
 t=tic;row=struct('id',C(k).id,'title',C(k).title,'status','ERROR','output',struct(), ...
  'checks',struct([]),'seconds',0,'error',struct());
 try
  [row.output,row.checks]=evaluateCase(C(k),o,runExact);
  if all([row.checks.pass])&&~isempty(row.checks),row.status='MATCHED';report.passed=report.passed+1;
  else,row.status='MISMATCH';report.failed=report.failed+1;end
  report.checksPassed=report.checksPassed+sum([row.checks.pass]);report.checksFailed=report.checksFailed+sum(~[row.checks.pass]);
  if isfield(row.output,'unresolvedRays')&&row.output.unresolvedRays>0
   report.numericalAllResolved=false;if strcmp(row.status,'MATCHED'),row.status='MATCHED_WITH_UNRESOLVED_RAY';end
  end
  if runExact
   report.exact.radiusCertificates=report.exact.radiusCertificates+row.output.exactRadiusCount;
   report.exact.rayCertificates=report.exact.rayCertificates+row.output.exactRayCount;
   report.exact.rayLiftingCertificates=report.exact.rayLiftingCertificates+row.output.exactLiftingCount;
  end
 catch err,row.error=struct('identifier',err.identifier,'message',err.message,'stack',err.stack);report.failed=report.failed+1;end
 row.seconds=toc(t);report.cases(end+1)=row;
 if o.Verbose
  fprintf('%-4s %-30s %s (%.3f s)\n',row.id,row.status,row.title,row.seconds);
  if strcmp(row.status,'ERROR'),fprintf('  %s\n',row.error.message);
  else,for j=find(~[row.checks.pass]),fprintf('  MISMATCH: %s\n',row.checks(j).name);end;end
 end
end
if runExact
 try
  if any(ismember({C.id},{'E1','E5V'}))
   report.exact.identities=verify_qstab_oneparameter_identities('Verbose',o.Verbose);
   report.exact.identityChecks=report.exact.identities.passed;
  else
   report.exact.identities=struct('ok',true,'executed',false,'status','NOT_APPLICABLE_TO_SELECTION','passed',0);
  end
  report.exact.ok=report.exact.identities.ok&&report.failed==0;
  if report.exact.ok,report.exact.status='EXECUTED_PASSED';else,report.exact.status='EXECUTED_REVIEW_REQUIRED';end
 catch err
  report.exact.status='ERROR';report.exact.identities=struct('ok',false,'identifier',err.identifier,'message',err.message);
 end
end
if o.Save
 if ~isfolder(out),mkdir(out);end
 if o.Plots
  try,report.graphics.files=plot_qstab_oneparameter(report,'Directory',out,'Visible',o.Visible);report.graphics.executed=true;
  catch err,report.graphics.executed=true;report.graphics.error=err.message;end
 end
end
report.elapsedSeconds=toc(started);
exactRequired=islogical(o.ExactCertificates)&&o.ExactCertificates;
report.ok=report.failed==0&&isempty(report.graphics.error)&&(~runExact||report.exact.ok)&&(~exactRequired||runExact);
report.articleComplete=report.ok&&runExact&&numel(C)==numel(catalogue)&& ...
 report.exact.radiusCertificates==9&&report.exact.rayCertificates==14&& ...
 report.exact.rayLiftingCertificates==13&&report.exact.identityChecks==42;
if report.articleComplete,report.status='ARTICLE_CHECKS_AND_EXACT_PROOFS_PASSED';
elseif report.ok,report.status='SELECTED_CHECKS_PASSED_ARTICLE_PROOFS_INCOMPLETE';
else,report.status='REVIEW_REQUIRED';end
if o.Save
 entries={'report.mat','report.txt','comparison.csv','radius_results.csv','ray_results.csv'};
 for k=1:numel(report.graphics.files),[~,n,e]=fileparts(report.graphics.files{k});entries{end+1}=[n e];end %#ok<AGROW>
 report.archiveEntries=entries;report.archivePolicy='Stage, unpack and compare every byte before replacing the fixed ZIP.';
 save(fullfile(out,'report.mat'),'report','-v7');writeText(report,fullfile(out,'report.txt'));writeCSV(report,fullfile(out,'comparison.csv'));
 writeResultsCSV(report,out);
 verifiedArchive(report.archiveFile,entries,out);
end
if o.Verbose
 fprintf('\nqstab one-parameter article: %s | %d groups matched, %d failed\n',report.status,report.passed,report.failed);
 fprintf('Fresh exact proof status: %s; complete article verification: %d\n',report.exact.status,report.articleComplete);
 if o.Save,fprintf('Results: %s\nSEND THIS ZIP: %s\n',out,report.archiveFile);end
end
end

function [D,K]=evaluateCase(C,o,exact)
D=struct('sourceAnchor',C.sourceAnchor,'exactRadiusCount',0,'exactRayCount',0,'exactLiftingCount',0,'unresolvedRays',0);
K=emptyChecks();d=C.data;
if isfield(d,'P0')
 [A,B,T]=normalized(d.P0,d.P1,d.frequencyScale);D.normalization=T;
 P0=poly(A,o);P1=poly(B,o);[tf,D.nominal]=qstab.isstable(P0);K(end+1)=ck('Nominal polynomial is Hurwitz','numerical',tf,1,0);
 D.problems=cell(1,numel(d.problems));
 for k=1:numel(d.problems)
  p=d.problems(k);side=p.side;w=p.frequency/T.frequencyScale;
  [rho,I]=qstab.freqmargin(P0,P1,w,'Side',side);
  [b,G]=qstab.stabilityradius(P0,P1,'Side',side,'AbsTol',o.AbsTol,'RelTol',o.RelTol,'MaxCells',o.MaxCells,'MaxTime',o.MaxTime);
  Q=struct('side',side,'reference',p,'frequency',w*T.frequencyScale,'frequencyRadius',rho, ...
   'frequencyInfo',I,'globalBounds',b,'globalInfo',G,'exact',struct('executed',false,'status','NOT_EXECUTED'));
  Q.svd=svdWitness(A,B,w,side);
  K(end+1)=ck([side ' fixed-frequency radius'],'pointwise',rho,p.radius,2e-8);
  K(end+1)=ck([side ' fixed-frequency witness valid'],'numerical',I.valid,true,0);
  K(end+1)=ck([side ' Hamilton and ordered SVD radii agree'],'pointwise',rho,Q.svd.radius,2e-8);
  K(end+1)=ck([side ' SVD original matrix residual'],'residual',Q.svd.originalResidual,0,1e-8);
  K(end+1)=ck([side ' global numerical bracket contains the reference'],'global-bracket', ...
   isfinite(G.lower)&&isfinite(G.upper)&&G.lower<=p.radius+1e-8*max(1,p.radius)&&G.upper>=p.radius-1e-8*max(1,p.radius),true,0);
  if exact
   if o.Verbose,fprintf('  %s: exact %s global-radius proof...\n',C.id,side);end
   rd=struct('P0',{rationalCells(d.rationalP0)},'P1',{rationalCells(d.rationalP1)});
   [eb,ei]=qstab.radiuscertificate(poly(d.P0,o),poly(d.P1,o),rationalCells(p.rationalBounds), ...
    'Side',side,'RationalData',rd,'UpperX',p.upperX);
   Q.exact=struct('executed',true,'bounds',eb,'info',ei,'status',ei.status);
   K(end+1)=ck([side ' exact global radius certificate'],'exact',ei.ok&&ei.certified,true,0);
   if ei.ok&&ei.certified,D.exactRadiusCount=D.exactRadiusCount+1;end
  end
  D.problems{k}=Q;
 end
end
D.rays=cell(1,numel(C.rays));
for k=1:numel(C.rays)
 rr=C.rays(k);tau=1;if strcmp(rr.articleCase,'E7'),tau=1e6;end
 [A,B,T]=normalized(rr.P0,rr.P1,tau);
 [ab,I]=qstab.raymargin(poly(A,o),poly(B,o),poly(rr.direction(:),o),'Side',rr.side);
 if rr.infinite,ref=Inf;else,ref=rr.reference;end
 Q=struct('id',rr.id,'name',rr.name,'normalization',T,'interval',ab,'info',I,'reference',ref, ...
  'directionComponents',rr.direction(:),'exact',struct('executed',false,'status','NOT_EXECUTED'));
 if strcmp(I.positiveStatus,'computed')
  K(end+1)=ck([rr.id ' first positive numerical contact'],'ray',ab(2),ref,2e-7);
 else
  D.unresolvedRays=D.unresolvedRays+1;
  K(end+1)=ck([rr.id ' unresolved contact remains NaN'],'diagnostic',isnan(ab(2)),true,0);
 end
 if exact
  if o.Verbose,fprintf('  %s/%s: exact first-contact and original-polynomial proof...\n',C.id,rr.id);end
  rd=struct('P0',{rationalCells(rr.rationalP0)},'P1',{rationalCells(rr.rationalP1)}, ...
   'Direction',{reshape(rationalCells(rr.rationalDirection),4,1)});
  [eb,ei]=qstab.raycertificate(poly(rr.P0,o),poly(rr.P1,o),poly(rr.direction(:),o), ...
   'Side',rr.side,'RationalData',rd,'AbsTol',1e-12,'MaxBisections',256,'VerifyLifting',true);
  Q.exact=struct('executed',true,'bounds',eb,'info',ei,'status',ei.status);
  K(end+1)=ck([rr.id ' exact first-contact proof'],'exact',ei.ok&&ei.certified,true,0);
  if rr.infinite,contains=all(isinf(eb)&eb>0);else,contains=eb(1)<=ref+2e-12*max(1,ref)&&eb(2)>=ref-2e-12*max(1,ref);end
  K(end+1)=ck([rr.id ' exact enclosure contains the source contact'],'exact',contains,true,0);
  if ei.ok&&ei.certified,D.exactRayCount=D.exactRayCount+1;end
  if ~rr.infinite
   liftingOK=isfield(ei,'lifting')&&ei.lifting.exact&& ...
    any(strcmp(ei.lifting.status,{'verified','verifiedImplicit'}));
   K(end+1)=ck([rr.id ' original-polynomial contact proved exactly'],'exact-lifting',liftingOK,true,0);
   if liftingOK,D.exactLiftingCount=D.exactLiftingCount+1;end
  end
 end
 D.rays{k}=Q;
end
switch C.id
 case 'E1'
  w=linspace(0,3,151);D.curves=struct('omega',w);
  D.curves.left=qstab.freqmargin(poly(d.P0,o),poly(d.P1,o),w,'Side','left');
  D.curves.right=qstab.freqmargin(poly(d.P0,o),poly(d.P1,o),w,'Side','right');
  D.curves.slice=zeros(size(w));
  for k=1:numel(w),D.curves.slice(k)=norm(evalq(d.P0,[0;w(k);0;0],'right'))/norm(evalq(d.P1,[0;w(k);0;0],'right'));end
  K(end+1)=ck('Left sweep identity','identity',D.curves.left,sqrt(2+w.^4./(1+w).^2),2e-9);
  K(end+1)=ck('Right sweep identity','identity',D.curves.right,sqrt((w.^4+2*w.^2-4*w+2)./(1+w).^2),2e-9);
  K(end+1)=ck('Fixed-slice quotient identity','identity',D.curves.slice,sqrt(1+w.^2+1./(1+w.^2)),2e-9);
  t=[0;1;0;1]/3;L=member(d.P0,d.P1,t,'left');R=member(d.P0,d.P1,t,'right');
  D.sampleParameter=t;D.sampleDecisions=[qstab.isstable(poly(L,o)),qstab.isstable(poly(R,o))];
  K(end+1)=ck('Same calibration parameter: actuator stable, sensor unstable','numerical',D.sampleDecisions,[1 0],0);
  D.rightHurwitzThird=minor(qstab.normpoly(poly(R,o)),3);
  K(end+1)=ck('Sensor error third Hurwitz minor','identity',D.rightHurwitzThird,-208/27,1e-11);
  [D.realGainInterval,D.realGainInfo]=qstab.raymargin(poly(d.P0,o),poly(d.P1,o),1);
  K(end+1)=ck('Pure real scale error has infinite ray margins','ray',D.realGainInterval,[-Inf Inf],0);
 case 'E2'
  P=complexPolynomial(d.P1,1i);D.perturbationSingularValues=svd(P);
  K(end+1)=ck('Perturbation matrix is rank one, not zero','numerical',sum(D.perturbationSingularValues>1e-10),1,0);
  D.criticalParameter=[-3;1;0;0]/4;Q=member(d.P0,d.P1,D.criticalParameter,'left');
  D.originalResidual=evalq(Q,[0;1;0;0],'right');
  K(end+1)=ck('Rank-one case original polynomial zero','residual',D.originalResidual,zeros(4,1),1e-12);
  D.forbiddenDenominator=evalq(d.P1,[0;-1;0;0],'right');
  K(end+1)=ck('Opposite point has a zero perturbation denominator','identity',D.forbiddenDenominator,zeros(4,1),1e-12);
 case 'E3'
  V=[eye(3),-eye(3),[1;2;3]/sqrt(14)];D.sphereParameters=[-.5*ones(1,size(V,2));-V/sqrt(2)];
  D.sphereResiduals=zeros(1,size(V,2));
  for k=1:size(V,2),P=member(d.P0,d.P1,D.sphereParameters(:,k),'left');D.sphereResiduals(k)=norm(evalq(P,[0;V(:,k)/sqrt(2)],'right'));end
  D.sphereCentre=[-.5;0;0;0];D.sphereRadius=1/sqrt(2);
  K(end+1)=ck('Nearest-gain sphere sample moduli','identity',sqrt(sum(D.sphereParameters.^2,1)),sqrt(3)/2*ones(1,size(V,2)),1e-12);
  K(end+1)=ck('Nearest-gain sphere original zeros','residual',D.sphereResiduals,zeros(size(D.sphereResiduals)),1e-12);
 case 'E4'
  w=d.problems(1).frequency;t=[w*w-1;0;w-1e-10;0];Q=D.problems{1};
  D.referenceParameter=t;D.HamiltonWitnessError=norm(Q.frequencyInfo.parameterComponents-t);
  D.svdWitnessError=norm(Q.svd.parameterComponents-t);D.singularValueGap=diff(flipud(Q.svd.singularValues));
  K(end+1)=ck('Hamilton orientation gives the analytic nearest gain','pointwise',D.HamiltonWitnessError,0,1e-10);
  % SVD witness error is intentionally recorded, not required to match a machine-specific value.
 case {'E5','E5V','E5C'}
  rr=C.rays(1);D.sampleGamma=[.5 1 1.25 1.5];D.sampleDecisions=NaN(size(D.sampleGamma));
  for k=1:numel(D.sampleGamma),P=member(rr.P0,rr.P1,D.sampleGamma(k)*rr.direction(:),rr.side);D.sampleDecisions(k)=qstab.isstable(poly(P,o));end
  P=member(rr.P0,rr.P1,rr.direction(:),rr.side);
  if any(strcmp(C.id,{'E5','E5V'}))
   K(end+1)=ck('Tangency has an original zero at zero frequency','residual',P(:,1),zeros(4,1),1e-12);
   K(end+1)=ck('Stability returns after the tangential first contact','numerical',D.sampleDecisions(3),1,0);
  end
  if strcmp(C.id,'E5V'),K(end+1)=ck('Noncommuting variant is stable at gamma=3/2','numerical',D.sampleDecisions(4),1,0);end
 case 'E6'
  D.referenceFrequency=1e-4;D.frequencyComparisonPolicy='Only radius accuracy is asserted at the nearly flat optimum; optimizer frequency is a diagnostic.';
 case 'E7'
  D.frequencyScale=d.frequencyScale;
  % Verify the original large coefficients at the actual frequency separately.
  D.originalFrequencyRadii=zeros(1,numel(d.problems));
  for k=1:numel(d.problems)
   p=d.problems(k);D.originalFrequencyRadii(k)=qstab.freqmargin(poly(d.P0,o),poly(d.P1,o),p.frequency,'Side',p.side,'CheckNominal',false);
   K(end+1)=ck([p.side ' combined scale invariance'],'pointwise',D.originalFrequencyRadii(k),p.radius,2e-8);
  end
 case 'A01'
  P0=[1 2 1;0 0 0;1 0 0;0 0 0];P1=[0 0;0 1;0 0;1 0];u=[0;1;0;1];
  [D.zeroInterval,D.zeroRay]=qstab.raymargin(poly(P0,o),0,1);[D.zeroRadius,D.zeroRadiusInfo]=qstab.stabilityradius(poly(P0,o),0);
  K(end+1)=ck('Zero perturbation ray','ray',D.zeroInterval,[-Inf Inf],0);
  K(end+1)=ck('Zero perturbation radius','global-bracket',D.zeroRadius,[Inf Inf],0);
  [D.baseInterval,D.baseInfo]=qstab.raymargin(poly(P0,o),poly(P1,o),poly(u,o));
  [D.nonunitInterval,D.nonunitInfo]=qstab.raymargin(poly(P0,o),poly(P1,o),poly(2*u,o));
  K(end+1)=ck('Direction doubling halves alpha, not the actual contact','ray',D.nonunitInterval,D.baseInterval/2,2e-8);
  a=[2;-1;3;1];A=zeros(size(P0));B=zeros(size(P1));
  for k=1:size(A,2),A(:,k)=mul(a,P0(:,k));end
  for k=1:size(B,2),B(:,k)=mul(a,P1(:,k));end
  % For a real ray the common left multiplier preserves both perturbation orders.
  [D.realBase,D.realBaseInfo]=qstab.raymargin(poly(P0,o),poly(P1,o),1);
  [D.nonmonicInterval,D.nonmonicInfo]=qstab.raymargin(poly(A,o),poly(B,o),1);
  K(end+1)=ck('Common nonreal leading factor on a real ray','ray',D.nonmonicInterval,D.realBase,2e-8);
 otherwise
end
end

function S=svdWitness(A,B,w,side)
P0=complexPolynomial(A,1i*w);P1=complexPolynomial(B,1i*w);
if strcmp(side,'left'),G=P1/P0;else,G=P0\P1;end
[U,D,V]=svd(G);sv=diag(D);v=V(:,1);u=U(:,1);
kv=[real(v(1));imag(v(1));-real(v(2));imag(v(2))];
ku=[real(u(1));imag(u(1));-real(u(2));imag(u(2))];
t=-mul(kv,conjq(ku))/sv(1);
if strcmp(side,'left'),Q=P0+chi(t)*P1;x=P0\v;else,Q=P0+P1*chi(t);x=u;end
S=struct('radius',1/sv(1),'singularValues',sv,'parameterComponents',t, ...
 'originalResidual',norm(Q*x)/max(realmin,(norm(P0)+norm(t)*norm(P1))*norm(x)), ...
 'assurance','Floating-point pointwise ordered SVD; no global norm proof.');
end
function P=complexPolynomial(C,z),P=zeros(2);for k=size(C,2):-1:1,P=P*z+chi(C(:,k));end;end
function A=chi(q),A=[q(1)+1i*q(2),q(3)+1i*q(4);-q(3)+1i*q(4),q(1)-1i*q(2)];end
function [A,B,T]=normalized(A,B,tau)
for k=1:size(A,2),A(:,k)=A(:,k)*tau^(k-1);end
for k=1:size(B,2),B(:,k)=B(:,k)*tau^(k-1);end
scale=max(abs([A(:);B(:)]));A=A/scale;B=B/scale;
T=struct('frequencyScale',tau,'commonAmplitudeScale',scale,'identity','P_normalized(s)=P_original(frequencyScale*s)/commonAmplitudeScale');
end
function P=poly(C,o)
if isvector(C)&&numel(C)==4,C=reshape(C,4,1);end
a=reshape(C(1,:),1,1,[]);b=reshape(C(2,:),1,1,[]);c=reshape(C(3,:),1,1,[]);d=reshape(C(4,:),1,1,[]);
if o.UseQuaternion,P=quaternion(a,b,c,d);else,P=struct('a',a,'b',b,'c',c,'d',d);end
end
function P=member(A,B,t,side)
P=A;for k=1:size(B,2),if strcmp(side,'left'),P(:,k)=P(:,k)+mul(t,B(:,k));else,P(:,k)=P(:,k)+mul(B(:,k),t);end;end
end
function y=evalq(P,z,side)
y=zeros(4,1);for k=size(P,2):-1:1,if strcmp(side,'right'),y=mul(y,z)+P(:,k);else,y=mul(z,y)+P(:,k);end;end
end
function z=mul(a,b),z=[a(1)*b(1)-a(2:4)'*b(2:4);a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];end
function q=conjq(q),q(2:4)=-q(2:4);end
function y=minor(N,k)
n=numel(N)-1;H=zeros(k);for i=1:k,for j=1:k,a=n-2*j+i;if a>=0&&a<=n,H(i,j)=N(a+1);end;end;end;y=det(H);
end
function C=rationalCells(x)
if isstring(x),C=cellstr(x);elseif ischar(x),C={x};
elseif iscell(x)&&~isempty(x)&&all(cellfun(@iscell,x(:)))
 rows=cell(size(x));for k=1:numel(x),rows{k}=reshape(rationalCells(x{k}),1,[]);end;C=vertcat(rows{:});
else,C=x;end
end
function K=emptyChecks(),K=struct('name',{},'kind',{},'pass',{},'actual',{},'reference',{},'tolerance',{});end
function K=ck(name,kind,a,b,tol)
same=isequal(size(a),size(b));good=false;
if same
 good=all((isinf(a(:))&isinf(b(:))&sign(a(:))==sign(b(:)))| ...
  (isfinite(a(:))&isfinite(b(:))&abs(a(:)-b(:))<=tol*max(1,abs(b(:)))));
end
K=struct('name',name,'kind',kind,'pass',logical(good),'actual',a,'reference',b,'tolerance',tol);
end
function [run,state]=exactAvailable(request)
run=false;state='NOT_REQUESTED';if islogical(request)&&~request,return;end
try,sym(1);run=true;state='READY';catch,state='SYMBOLIC_UNAVAILABLE';end
end
function H=sourceHashes(root)
names={'inputs.json','source_provenance.json'};H=struct('file',{},'sha256',{});
for k=1:numel(names)
 value='unavailable';try
  fid=fopen(fullfile(root,names{k}),'rb');cleanup=onCleanup(@()fclose(fid));bytes=fread(fid,Inf,'*uint8');clear cleanup
  if exist('OCTAVE_VERSION','builtin'),value=hash('sha256',char(bytes.'));
  else,md=java.security.MessageDigest.getInstance('SHA-256');md.update(bytes);value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));end
 catch,end
 H(end+1)=struct('file',names{k},'sha256',value); %#ok<AGROW>
end
end
function R=resolvedFunctions()
names={'run_qstab_oneparameter','qstab.freqmargin','qstab.raymargin','qstab.stabilityradius', ...
 'qstab.radiuscertificate','qstab.raycertificate','verify_qstab_oneparameter_identities'};
R=struct('name',{},'path',{});for k=1:numel(names),R(end+1)=struct('name',names{k},'path',which(names{k}));end
end
function y=flag(x),y=islogical(x)&&isscalar(x);end
function y=positive(x),y=isnumeric(x)&&isreal(x)&&isscalar(x)&&~isnan(x)&&x>0;end
function y=positiveFinite(x),y=positive(x)&&isfinite(x);end
function y=exactOption(x),y=flag(x)||(ischar(x)&&strcmpi(x,'auto'))||(isstring(x)&&isscalar(x)&&strcmpi(x,'auto'));end
function y=absolute(x),y=~isempty(x)&&(x(1)=='/'||x(1)=='\'||(numel(x)>1&&x(2)==':'));end
function writeText(R,file)
f=fopen(file,'w');if f<0,error('qstab:oneparameter:Write','Cannot write %s.',file);end;c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'qstab scalar one-parameter article\nVersion: %s\nRun: %s\nRuntime: %s | %s\nStatus: %s\nArticle complete: %d\nExact: %s\n%s\n\n', ...
 R.toolbox.version,R.runId,R.runtime,R.computer,R.status,R.articleComplete,R.exact.status,R.assurance);
for k=1:numel(R.cases)
 A=R.cases(k);fprintf(f,'%s %s | %s\n',A.id,A.status,A.title);
 if isfield(A.output,'problems')
  for j=1:numel(A.output.problems)
   P=A.output.problems{j};G=P.globalInfo;[lo,hi,state]=exactSummary(P.exact);
   fprintf(f,'  %s numerical global radius: [%.17g, %.17g] (%s)\n',P.side,G.lower,G.upper,G.status);
   fprintf(f,'    reference radius %.17g; value at reference frequency %.17g\n',P.reference.radius,P.frequencyRadius);
   fprintf(f,'    exact rational radius bounds: [%s, %s] (%s)\n',lo,hi,state);
  end
 end
 if isfield(A.output,'rays')
  for j=1:numel(A.output.rays)
   P=A.output.rays{j};[lo,hi,state]=exactSummary(P.exact);
   fprintf(f,'  %s direction %s; numerical signed interval %s (%s)\n',P.id,mat2str(P.directionComponents.',16),mat2str(P.interval,17),P.info.positiveStatus);
   fprintf(f,'    reference positive contact %.17g; exact positive enclosure [%s, %s] (%s)\n',P.reference,lo,hi,state);
   if P.exact.executed&&isfield(P.exact.info,'lifting'),fprintf(f,'    original-polynomial lifting: %s\n',P.exact.info.lifting.status);end
  end
 end
 if isfield(A.output,'svdWitnessError')
  fprintf(f,'  Nearest-gain errors: Hamilton %.17g; SVD %.17g; singular-value gap %.17g\n', ...
   A.output.HamiltonWitnessError,A.output.svdWitnessError,A.output.singularValueGap);
 end
 for j=1:numel(A.checks),B=A.checks(j);fprintf(f,'  pass=%d [%s] %s\n    actual: %s\n    reference: %s\n',B.pass,B.kind,B.name,mat2str(B.actual,16),mat2str(B.reference,16));end
 if strcmp(A.status,'ERROR'),fprintf(f,'  ERROR %s: %s\n',A.error.identifier,A.error.message);end
 fprintf(f,'\n');
end
fprintf(f,'Fresh radius certificates: %d\nFresh ray certificates: %d\nFresh original-polynomial ray-contact proofs: %d\nFresh exact identities: %d\nElapsed: %.6f s\nArchive: %s\n', ...
 R.exact.radiusCertificates,R.exact.rayCertificates,R.exact.rayLiftingCertificates,R.exact.identityChecks,R.elapsedSeconds,R.archiveFile);
end
function writeCSV(R,file)
f=fopen(file,'w');if f<0,error('qstab:oneparameter:Write','Cannot write %s.',file);end;c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'case,kind,check,passed,computed,reference,tolerance\n');
for k=1:numel(R.cases),A=R.cases(k);for j=1:numel(A.checks),B=A.checks(j);fprintf(f,'%s,%s,"%s",%d,"%s","%s",%.17g\n', ...
 A.id,B.kind,strrep(B.name,'"','""'),B.pass,mat2str(B.actual,16),mat2str(B.reference,16),B.tolerance);end;end
end
function writeResultsCSV(R,folder)
f=fopen(fullfile(folder,'radius_results.csv'),'w');if f<0,error('qstab:oneparameter:Write','Cannot write the radius table.');end
c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'case,side,numerical_lower,numerical_upper,numerical_status,estimated_frequency,reference_radius,reference_frequency,radius_at_reference_frequency,exact_lower,exact_upper,exact_status\n');
for k=1:numel(R.cases)
 A=R.cases(k);if ~isfield(A.output,'problems'),continue;end
 for j=1:numel(A.output.problems)
  P=A.output.problems{j};G=P.globalInfo;[lo,hi,state]=exactSummary(P.exact);
  fprintf(f,'%s,%s,%.17g,%.17g,%s,%.17g,%.17g,%.17g,%.17g,"%s","%s",%s\n', ...
   A.id,P.side,G.lower,G.upper,G.status,G.frequency*A.output.normalization.frequencyScale, ...
   P.reference.radius,P.reference.frequency,P.frequencyRadius,lo,hi,state);
 end
end
clear c
f=fopen(fullfile(folder,'ray_results.csv'),'w');if f<0,error('qstab:oneparameter:Write','Cannot write the ray table.');end
c=onCleanup(@()fclose(f)); %#ok<NASGU>
fprintf(f,'case,ray,side,direction_real,direction_i,direction_j,direction_k,numerical_negative_exit,numerical_positive_exit,numerical_positive_status,reference_positive_exit,exact_positive_lower,exact_positive_upper,exact_status,original_polynomial_lifting\n');
for k=1:numel(R.cases)
 A=R.cases(k);if ~isfield(A.output,'rays'),continue;end
 for j=1:numel(A.output.rays)
  P=A.output.rays{j};[lo,hi,state]=exactSummary(P.exact);lift='NOT_EXECUTED';
  if P.exact.executed&&isfield(P.exact.info,'lifting'),lift=P.exact.info.lifting.status;end
  u=P.directionComponents;
  fprintf(f,'%s,%s,%s,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%s,%.17g,"%s","%s",%s,%s\n', ...
   A.id,P.id,P.info.side,u(1),u(2),u(3),u(4),P.interval(1),P.interval(2),P.info.positiveStatus,P.reference,lo,hi,state,lift);
 end
end
end
function [lo,hi,state]=exactSummary(E)
lo='';hi='';state=E.status;
if E.executed&&isfield(E.info,'exactBounds')&&numel(E.info.exactBounds)==2
 B=E.info.exactBounds;lo=char(B{1});hi=char(B{2});
end
end
function verifiedArchive(target,entries,folder)
staged=[tempname(folder) '.zip'];checkdir=tempname(folder);
cleanup=onCleanup(@()archiveCleanup(staged,checkdir)); %#ok<NASGU>
zip(staged,entries,folder);mkdir(checkdir);unpacked=unzip(staged,checkdir);
if numel(unpacked)~=numel(entries),error('qstab:oneparameter:Archive','The staged ZIP has an unexpected entry count.');end
for k=1:numel(entries)
 source=fullfile(folder,entries{k});copy=fullfile(checkdir,entries{k});
 if ~isfile(copy)||~isequal(readBytes(source),readBytes(copy))
  error('qstab:oneparameter:Archive','The staged ZIP failed byte comparison for %s.',entries{k});
 end
end
[ok,msg]=movefile(staged,target,'f');
if ~ok,error('qstab:oneparameter:Archive','Verified ZIP could not replace the result archive: %s',msg);end
end
function bytes=readBytes(file)
f=fopen(file,'rb');if f<0,error('qstab:oneparameter:Read','Cannot read %s.',file);end
c=onCleanup(@()fclose(f));bytes=fread(f,Inf,'*uint8'); %#ok<NASGU>
end
function archiveCleanup(file,folder)
if isfile(file),delete(file);end
if isfolder(folder),rmdir(folder,'s');end
end
