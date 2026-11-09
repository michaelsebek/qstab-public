function E=run_qstab_paper_extensions(varargin)
%RUN_QSTAB_PAPER_EXTENSIONS Reproduce SELECTED source-linked qB01/qB02/qB03 examples.
% E=run_qstab_paper_extensions returns ordinary data only: catalogue/reference inputs,
% solver outputs, source checks and execution metadata. No figure, directory,
% path, RNG, archive or third-party files are modified. It is not a benchmark.
%
% Options: Articles={'qB01','qB02','qB03'}, Cases={} (all selected articles),
% UseQuaternion=true, Verbose=true. These finite checks have no subdivision budget.
% An empty article selection, duplicate IDs and unknown names are errors;
% Cases={} selects all cases of the requested articles. Selection is validated
% BEFORE the first solver call. Output order follows
% the catalogue, not the selection list. BOUNDED/UNRESOLVED are kept explicit.
% E.ok means no reference contradictions/errors, NOT complete reproduction.
% E.fullyMatched additionally requires all selected reference checks matched.
% To include E in the ONE standard report archive:
%   report=test_qstab('PaperExtensions',true);
% All original numbered figures/examples are NOT reproduced by this suite.
ip=inputParser;ip.FunctionName='run_qstab_paper_extensions';ip.PartialMatching=false;
addParameter(ip,'Articles',{'qB01','qB02','qB03'});
addParameter(ip,'Cases',{});
addParameter(ip,'UseQuaternion',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;
articles=textList(o.Articles,false);ids=textList(o.Cases,true);C=qstab_paper_extension_cases();
if any(~ismember(articles,{'qB01','qB02','qB03'})),error('qstab:paper:Article','Unknown article.');end
if any(~ismember(ids,{C.id})),error('qstab:paper:Selection','Unknown case ID.');end
keep=ismember({C.article},articles);
if ~isempty(ids)
 if any(~ismember(ids,{C(keep).id})),error('qstab:paper:Selection','Case does not belong to selected articles.');end
 keep=keep&ismember({C.id},ids);
end
C=C(keep);
if isempty(C),error('qstab:paper:Selection','At least one case must be selected.');end
if o.UseQuaternion
 try,q=quaternion(0,0,0,0);[a,b,c,d]=parts(q); %#ok<ASGLU>
 catch ME,error('qstab:paper:Quaternion','Native quaternion unavailable: %s',ME.message);end
end
E=struct('schema','qstab.paper-extensions.1','sourceVersion',qstab.about(), ...
 'runtime',version,'isOctave',exist('OCTAVE_VERSION','builtin')~=0, ...
 'nativeQuaternionInputs',o.UseQuaternion,'options',o,'catalogue',{C}, ...
 'status','RUNNING','ok',true,'fullyMatched',false,'completePaperReproduction',false, ...
 'matched',0,'bounded',0,'unresolved',0,'failed',0,'errors',0,'elapsedSeconds',0, ...
 'cases',struct('id',{},'article',{},'status',{},'output',{},'checks',{}, ...
   'elapsedSeconds',{},'message',{},'identifier',{}));
clock=tic;
for j=1:numel(C)
 row=struct('id',C(j).id,'article',C(j).article,'status','ERROR','output',struct(), ...
    'checks',{{}},'elapsedSeconds',0,'message','','identifier','');
 one=tic;
 try
  [row.output,row.checks]=qstab_paper_extension_evaluate(C(j),o);
  states=cellfun(@(x)x.status,row.checks,'UniformOutput',false);
  if isempty(states)||any(strcmp(states,'MISMATCH')),row.status='MISMATCH';E.failed=E.failed+1;
  elseif any(strcmp(states,'UNRESOLVED')),row.status='UNRESOLVED';E.unresolved=E.unresolved+1;
  elseif any(strcmp(states,'BOUNDED')),row.status='BOUNDED';E.bounded=E.bounded+1;
  else,row.status='MATCHED';E.matched=E.matched+1;end
  if ~strcmp(row.status,'MATCHED'),row.message=checkMessages(row.checks);end
 catch ME
  row.status='ERROR';row.message=ME.message;row.identifier=ME.identifier;E.errors=E.errors+1;
 end
 row.elapsedSeconds=toc(one);E.cases(end+1)=row;
 if o.Verbose
  fprintf('  paper %-30s %-12s %.3f s\n',row.id,row.status,row.elapsedSeconds);
  if ~isempty(row.message),fprintf('    %s\n',row.message);end
 end
end
E.elapsedSeconds=toc(clock);E.ok=E.failed==0&&E.errors==0;
E.fullyMatched=E.ok&&E.bounded==0&&E.unresolved==0;
if ~E.ok,E.status='REFERENCE_ISSUES';
elseif E.fullyMatched,E.status='COMPLETED';
else,E.status='COMPLETED_WITH_UNRESOLVED';end
E.note=['Selected source checks only; source references are separate from calculated outputs. ', ...
 'A bracket is not an exact radius; a pointwise witness is not a new global proof.'];
if o.Verbose
 fprintf('qstab additional paper examples: %s | %d matched, %d bounded, %d unresolved, %d mismatches, %d errors\n', ...
  E.status,E.matched,E.bounded,E.unresolved,E.failed,E.errors);
end
end
function c=textList(x,allowEmpty)
if ischar(x)&&isrow(x)&&~isempty(x),c={x};
elseif isstring(x)&&isvector(x)&&~isempty(x),c=cellstr(x(:).');
elseif iscell(x)&&(isvector(x)||isempty(x)),c=x(:).';
else,error('qstab:paper:Selection','Selections must be text or a cell vector of text.');end
if isempty(c)&&~allowEmpty,error('qstab:paper:Selection','Article selection must not be empty.');end
for j=1:numel(c)
 if isstring(c{j})&&isscalar(c{j}),c{j}=char(c{j});end
 if ~ischar(c{j})||~isrow(c{j})||isempty(strtrim(c{j})),error('qstab:paper:Selection','Invalid empty/numeric name.');end
end
if numel(unique(c))~=numel(c),error('qstab:paper:Selection','Duplicate selection names are not allowed.');end
end

function message=checkMessages(checks)
% Show which reference check is unresolved/contradictory, without changing it.
notes={};
for j=1:numel(checks)
 A=checks{j};
 if strcmp(A.status,'MATCHED'),continue;end
 note=sprintf('Check %d (%s): %s. %s',j,A.kind,A.status,A.message);
 if strcmp(A.kind,'status')&&ischar(A.actual)&&ischar(A.reference)
  note=[note,' actual=',A.actual,', reference=',A.reference,'.'];
 end
 notes{end+1}=note; %#ok<AGROW>
end
message=strjoin(notes,' | ');
end
