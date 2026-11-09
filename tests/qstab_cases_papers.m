function T=qstab_cases_papers(writer)
%QSTAB_CASES_PAPERS Reproducibility metadata and source-derived fixtures.
T={ ...
 'Selected paper catalogue has 20 unique source-anchored ordinary-data cases',@catalogue; ...
 'Internal qB01 Revision 2 is not silently renamed to Revision 3',@revision; ...
 'Source radius references remain separate from calculated outputs',@referenceSeparation; ...
 'Unknown status is not promoted to a matched source property',@unknownStatus; ...
 'Opposite decided status is a reference mismatch',@oppositeStatus; ...
 'A finite nonconverged bracket is containment rather than exact reproduction',@bounded; ...
 'A one-sided bound is retained as unresolved',@oneSided; ...
 'A source radius outside the returned bracket is rejected',@badBracket; ...
 'Malformed or NaN reference bounds cannot be called matched',@invalidBracket; ...
 'Numeric reference comparisons reject shape mismatches and nonfinite results',@numeric; ...
 'Unknown selections fail before a solver is called',@unknownSelection; ...
 'Duplicate and empty article selections are rejected',@duplicateSelection; ...
 'Case selection respects article constraints and deterministic catalogue order',@selection; ...
 'The paper runner leaves directory, path and random state unchanged',@sideEffects; ...
 'qB01 source side-separated real-ray contacts reproduce without normalization',@raySides; ...
 'qB01 cubic first contact is not replaced by renewed stability',@reentry; ...
 'qB02 stable endpoints and actual unstable interior retain the exact norm curve',@segment; ...
 'qB02 epsilon box checks all 256 vertices and its unstable interior separately',@vertices; ...
 'qB02 rational edge example checks all 256 vertices without a vertex theorem',@edgeBox; ...
 'qB02 hidden-edge point is a witness, not a new global minimization claim',@edgePoint; ...
 'qB02 eight-dimensional source-margin point has frequency two',@fullPoint; ...
 'qB03 parabola retains mixed controls and its non-certificate sample status',@parabola; ...
 'qB03 false alarm is separated from the complete-family decision',@bernstein; ...
 'qB03 master factorization and boundary quaternion are evaluated directly',@master; ...
 'qB03 tetrahedral face distances and unstable centre keep their meanings',@tetra; ...
 'qB03 Schur tetrahedron preserves its varying-leading Cayley coefficients',@schur; ...
 'qB03 degree loss is not a finite unstable root claim',@degree; ...
 'qB03 shared parameter is not replaced by independent coefficient choices',@correlation; ...
 'Paper payloads contain no function handles, classes or figure objects',@ordinary; ...
 'Paper payload round-trips through the unchanged two-entry result archive',@()roundtrip(writer)};
end
function catalogue()
C=qstab_paper_cases();assert(numel(C)==20&&numel(unique({C.id}))==20);
assert(sum(strcmp({C.article},'qB01'))==4&&sum(strcmp({C.article},'qB02'))==9&&sum(strcmp({C.article},'qB03'))==7);
assert(all(~cellfun(@isempty,{C.sourceAnchor}))&&plain(C));
end
function revision()
C=qstab_paper_cases();a=C(strcmp({C.article},'qB01'));
assert(all(cellfun(@(s)~isempty(strfind(s,'Revision 2')), {a.sourceRevision})));
end
function referenceSeparation()
E=run('qb03_parabola');assert(isfield(E.catalogue,'reference')&&isfield(E.cases,'output'));
assert(~E.completePaperReproduction&&E.fullyMatched);
end
function unknownStatus()
A=qstab_paper_check('status','undetermined','stable');assert(A.ok&&~A.matched&&strcmp(A.status,'UNRESOLVED'));
end
function oppositeStatus()
A=qstab_paper_check('status','notStable','stable');assert(~A.ok&&strcmp(A.status,'MISMATCH'));
end
function bounded()
A=qstab_paper_check('bracket',struct('lower',.5,'upper',1.5,'status','bounded'),1);
assert(A.ok&&~A.matched&&strcmp(A.status,'BOUNDED'));
end
function oneSided()
A=qstab_paper_check('bracket',struct('lower',0,'upper',Inf,'status','undetermined'),1);
assert(A.ok&&~A.matched&&strcmp(A.status,'UNRESOLVED'));
end
function badBracket()
for b={[1.1 2],[0 .9]}
 A=qstab_paper_check('bracket',struct('lower',b{1}(1),'upper',b{1}(2),'status','converged'),1);assert(~A.ok);
end
end
function invalidBracket()
for b={[NaN Inf],[-1 2],[2 1]}
 A=qstab_paper_check('bracket',struct('lower',b{1}(1),'upper',b{1}(2),'status','bounded'),1);assert(~A.ok);
end
fails(@()qstab_paper_check('bracket',struct(),1));
end
function numeric()
for a={NaN,Inf,[1;1]}
 A=qstab_paper_check('numeric',a{1},[1 1]);assert(~A.ok);
end
A=qstab_paper_check('numeric',1+2i,1+2i);assert(A.matched);
end
function unknownSelection()
fails(@()run_qstab_papers('Cases',{'not-a-case'},'UseQuaternion',false));
fails(@()run_qstab_papers('Articles','qB04','UseQuaternion',false));
end
function duplicateSelection()
fails(@()run_qstab_papers('Articles',{},'UseQuaternion',false));
fails(@()run_qstab_papers('Cases',{'qb03_master','qb03_master'},'UseQuaternion',false));
fails(@()run_qstab_papers('Cases',{1},'UseQuaternion',false));
fails(@()run_qstab_papers('MaxCells',0,'UseQuaternion',false));
end
function selection()
fails(@()run_qstab_papers('Articles','qB01','Cases','qb03_master','UseQuaternion',false));
E=run_qstab_papers('Cases',{'qb03_master','qb03_parabola'},'UseQuaternion',false,'Verbose',false);
assert(isequal({E.cases.id},{'qb03_parabola','qb03_master'}));assert(E.ok);
end
function sideEffects()
p=pwd;pth=path;rs=rng;before=dir(p);E=run('qb03_master');after=dir(p);
assert(E.ok&&strcmp(p,pwd)&&strcmp(pth,path)&&isequal(rs,rng)&&isequal({before.name},{after.name}));
end
function raySides(),matched('qb01_ray_sides');end
function reentry(),matched('qb01_reentry');end
function segment(),matched('qb02_segment');end
function vertices(),E=matched('qb02_fullbox_vertices');assert(E.cases.output.vertexCount==256);end
function edgeBox(),E=matched('qb02_fullbox_edge');assert(E.cases.output.vertexCount==256);end
function edgePoint()
E=matched('qb02_edge_witness');assert(strcmp(E.catalogue.scope,'pointwiseWitnessNotGlobalProof'));
end
function fullPoint()
E=matched('qb02_fullbox_margin_point');assert(abs(norm(E.catalogue.data.point)-2)<1e-12);
assert(E.catalogue.data.point(1)==0);
end
function parabola(),matched('qb03_parabola');end
function bernstein(),matched('qb03_bernstein');end
function master(),matched('qb03_master');end
function tetra(),matched('qb03_tetrahedron');end
function schur(),matched('qb03_schur_tetrahedron');end
function degree(),matched('qb03_degree_loss');end
function correlation(),matched('qb03_correlation');end
function ordinary()
E=run('qb03_bernstein');assert(plain(E));
end
function E=run(id)
E=run_qstab_papers('Cases',id,'UseQuaternion',false,'Verbose',false,'MaxTime',Inf);
end
function E=matched(id)
E=run(id);assert(E.ok&&E.fullyMatched, ['Paper case failed: ',id]);
end
function roundtrip(writer)
base=tempname;mkdir(base);cleanup=onCleanup(@()rmdir(base,'s')); %#ok<NASGU>
item=struct('group','fixture','name','Artificial paper-data persistence fixture', ...
 'status','PASS','seconds',0,'message','','identifier','');
r=struct('toolbox',struct('version','fixture'), ...
 'environment',struct('version','artificial','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','20000101T000000','runId','paper-persistence-fixture', ...
 'nativeMATLABRun',false,'nativeQuaternionOK',false,'ok',true,'status','OVERALL_OK', ...
 'passed',1,'failed',0,'skipped',0,'elapsedSeconds',0,'cases',item);
r.paperExamples=run('qb03_parabola');r=writer(r,base,'qstab_test_results');
folder=fullfile(base,'unpacked');mkdir(folder);unzip(r.archiveFile,folder);
f=dir(folder);f=f(~ismember({f.name},{'.','..'}));
assert(numel(f)==2&&isequal(sort({f.name}),sort({'report.mat','report.txt'})));
s=load(fullfile(folder,'report.mat'),'report');assert(isequaln(s.report.paperExamples,r.paperExamples));
assert(strcmp(r.archiveFile,fullfile(base,'qstab_test_results.zip')));
end
function tf=plain(x)
tf=isnumeric(x)||islogical(x)||ischar(x);
if isstruct(x)
 tf=true;names=fieldnames(x);
 for a=1:numel(x),for b=1:numel(names),if ~plain(x(a).(names{b})),tf=false;return;end,end,end
elseif iscell(x)
 tf=true;for a=1:numel(x),if ~plain(x{a}),tf=false;return;end,end
end
end
function fails(f)
caught=false;try,f();catch,caught=true;end;assert(caught,'Expected invalid input rejection.');
end
