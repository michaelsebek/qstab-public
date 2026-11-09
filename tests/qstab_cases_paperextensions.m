function T=qstab_cases_paperextensions(writer)
%QSTAB_CASES_PAPEREXTENSIONS Additional sources, plot data and safe persistence.
T={ ...
 'Extended catalogue has 16 distinct source-linked cases, separate from the original 20',@catalogue; ...
 'Extension selection retains catalogue order and rejects incompatible article selection',@selection; ...
 'Unknown, duplicate and malformed extension options fail before evaluation',@invalidSelection; ...
 'Source scope distinguishes finite identities from a general nonlinear-family certificate',@scope; ...
 'Additional checks do not change directory, path, random state or create files',@sideEffects; ...
 'qB01 several nearest-sphere parameters have the reference radius and actual zeros',@()matched('qb01_nearest_sphere'); ...
 'qB01 fixed real degree lifting preserves both frequency quotient values',@()matched('qb01_degree_lift'); ...
 'qB01 coefficient conjugation interchanges left and right frequency quotients',@()matched('qb01_conjugate_sides'); ...
 'qB02 larger real damping produces the source loss and reentry',@()matched('qb02_damping'); ...
 'qB02 unstable disk keeps sampled members and analytic edge distance distinct',@()matched('qb02_unstable_disk'); ...
 'qB02 cubic real constant direction retains the norm identity and unstable middle',@()matched('qb02_cubic_constant'); ...
 'qB02 fixed-slice omission and false projection are different phenomena',@()matched('qb02_fixed_slice'); ...
 'qB02 parabolic bands retain the zero-frequency and explicit boundary cases',@()matched('qb02_bands'); ...
 'qB02 zero in a convex hull is not a zero of the attainable curve',@()matched('qb02_nonconvex_hull'); ...
 'qB02 nonlinear rotating coefficient checks do not invoke a generic nonlinear solver',@()matched('qb02_nonlinear'); ...
 'qB02 exact marginals cannot be independently paired across different parameters',@()matched('qb02_correlated_marginals'); ...
 'qB02 ordered rotating factors preserve norm polynomials and zero classes',@()matched('qb02_rotating_factors'); ...
 'qB02 isolated contact is missed by a stable grid but verified explicitly',@()matched('qb02_isolated_touch'); ...
 'qB03 ordered polynomial multiplication differs from pointwise evaluation',@()matched('qb03_ordered_product'); ...
 'qB03 correlated square has unstable centre despite the source boundary distance',@()matched('qb03_square'); ...
 'qB03 reversed feedback multiplication removes the actual uncertain term',@()matched('qb03_feedback_order'); ...
 'Figure-data catalogue has eight source-labelled items without figure handles',@figureCatalogue; ...
 'Computed sweep and norm curves agree with separately stored source expressions',@figureIdentities; ...
 'Zero-frequency bands are real segments, not division by zero',@zeroFrequency; ...
 'Refined Bernstein picture retains two local hulls and mixed controls',@localHulls; ...
 'Picture metadata never claims full paper or continuum certification',@figureScope; ...
 'Figure-data generation leaves directory, path, random state and graphics unchanged',@figureSideEffects; ...
 'Figure choices and sample counts are validated before numerical work',@invalidFigures; ...
 'CSV export preserves every saved coordinate at round-trip precision',@csvRoundTrip; ...
 'Export overwrite refusal leaves existing outputs and unrelated files unchanged',@overwrite; ...
 'Invalid figure data and format requests create neither folders nor figures',@invalidExport; ...
 'Absolute figure output is independent of the current working directory',@exportDirectory; ...
 'Figure export manifest retains source anchors and distinguishes data from export runtime',@manifest; ...
 'New reproduction and figure payloads round-trip through the unchanged two-entry ZIP',@()roundTrip(writer)};
end
function catalogue()
C=qstab_paper_extension_cases();B=qstab_paper_cases();assert(numel(C)==16&&numel(B)==20);
assert(numel(unique({C.id}))==16&&isempty(intersect({C.id},{B.id})));
assert(sum(strcmp({C.article},'qB01'))==3&&sum(strcmp({C.article},'qB02'))==10&&sum(strcmp({C.article},'qB03'))==3);
assert(plain(C)&&all(~cellfun(@isempty,{C.sourceAnchor})));
end
function selection()
E=run_qstab_paper_extensions('Cases',{'qb03_ordered_product','qb01_nearest_sphere'},'UseQuaternion',false,'Verbose',false);
assert(isequal({E.cases.id},{'qb01_nearest_sphere','qb03_ordered_product'})&&E.ok);
fails(@()run_qstab_paper_extensions('Articles','qB01','Cases','qb03_square','UseQuaternion',false));
end
function invalidSelection()
fails(@()run_qstab_paper_extensions('Cases',{'x'},'UseQuaternion',false));
fails(@()run_qstab_paper_extensions('Articles',{},'UseQuaternion',false));
fails(@()run_qstab_paper_extensions('Cases',{'qb03_square','qb03_square'},'UseQuaternion',false));
fails(@()run_qstab_paper_extensions('Cases',{1},'UseQuaternion',false));
fails(@()run_qstab_paper_extensions('UseQuaternion',3));
end
function scope()
E=matched('qb02_nonlinear');assert(~E.completePaperReproduction&&~E.cases.output.isContinuumCertificate);
assert(~isempty(strfind(E.catalogue.scope,'NotGenericNonlinearSolver')));
end
function sideEffects()
p=pwd;pth=path;r=rng;f=dir(p);E=matched('qb03_ordered_product');g=dir(p);
assert(E.ok&&strcmp(p,pwd)&&strcmp(pth,path)&&isequal(r,rng)&&isequal({f.name},{g.name}));
end
function E=matched(id)
E=run_qstab_paper_extensions('Cases',id,'UseQuaternion',false,'Verbose',false);
assert(E.ok&&E.fullyMatched,['Extension failed: ',id,' ',E.cases.message]);
end
function F=figdata(ids)
F=qstab_paper_figures('Figures',ids,'UseQuaternion',false,'Samples',21);assert(F.ok);
end
function figureCatalogue()
F=figdata({});assert(numel(F.items)==8&&numel(unique({F.items.id}))==8&&plain(F));
for j=1:numel(F.items),qstab_validate_paper_figure(F.items(j));end
end
function figureIdentities()
F=figdata({});for j=1:numel(F.items),K=F.items(j).checks;assert(~isempty(K)&&all(cellfun(@(c)c.matched,K)));end
end
function zeroFrequency()
F=figdata('qb02_bands');c=F.items.curves(1);
assert(isequal(c.x,[1;7])&&isequal(c.y,[0;0])&&strcmp(c.role,'analyticBoundary'));
end
function localHulls()
F=figdata('qb03_parabola');D=F.items.raw;
assert(numel(D.refined.regions)==2&&numel(D.coarse.regions{1}.controls)==3);
assert(sum(strcmp({F.items.curves.role},'outerBoundary'))==3);
end
function figureScope()
F=figdata({});assert(~F.isCertificate&&~F.completePaperReproduction);
assert(all(~[F.items.isCertificate])&&all(~cellfun(@isempty,{F.items.sourceRevision})));
end
function figureSideEffects()
p=pwd;pth=path;r=rng;f=dir(p);old=findall(0,'Type','figure');F=figdata('qb02_counter');g=dir(p);
assert(F.ok&&strcmp(p,pwd)&&strcmp(pth,path)&&isequal(r,rng)&&isequal({f.name},{g.name})&&isequal(old,findall(0,'Type','figure')));
end
function invalidFigures()
fails(@()qstab_paper_figures('Figures','bad','UseQuaternion',false));
fails(@()qstab_paper_figures('Figures',{'qb02_hull','qb02_hull'},'UseQuaternion',false));
fails(@()qstab_paper_figures('Samples',100,'UseQuaternion',false));
fails(@()qstab_paper_figures('Samples',3,'UseQuaternion',false));
end
function csvRoundTrip()
F=figdata('qb02_hull');p=tempname;mkdir(p);c=onCleanup(@()rmdir(p,'s')); %#ok<NASGU>
X=export_qstab_paper_figures(F,'OutputDir',p,'Formats',{'csv'});assert(X.ok&&X.solverCalls==0&&numel(X.files)==1);
A=dlmread(X.files.file,',',1,0);B=[];
for j=1:numel(F.items.curves)
 v=F.items.curves(j);n=numel(v.x);B=[B;[j*ones(n,1),(1:n).',v.x,v.y,zeros(n,1)]]; %#ok<AGROW>
end
assert(isequal(A,B));
end
function overwrite()
F=figdata('qb02_counter');p=tempname;mkdir(p);c=onCleanup(@()rmdir(p,'s')); %#ok<NASGU>
X=export_qstab_paper_figures(F,'OutputDir',p,'Formats','csv');before=readBytes(X.files.file);
fid=fopen(fullfile(p,'keep.txt'),'w');fprintf(fid,'do not touch');fclose(fid);
fails(@()export_qstab_paper_figures(F,'OutputDir',p,'Formats','csv','Overwrite',false));
assert(isequal(before,readBytes(X.files.file))&&strcmp(fileread(fullfile(p,'keep.txt')),'do not touch'));
X=export_qstab_paper_figures(F,'OutputDir',p,'Formats','csv');assert(X.ok&&isequal(before,readBytes(X.files.file)));
end
function invalidExport()
F=figdata('qb02_hull');p=tempname;old=findall(0,'Type','figure');
G=F;G.items.id='../not-safe';fails(@()export_qstab_paper_figures(G,'OutputDir',p,'Formats','csv'));
G=F;G.items.curves(1).x(1)=NaN;fails(@()export_qstab_paper_figures(G,'OutputDir',p,'Formats','csv'));
fails(@()export_qstab_paper_figures(F,'OutputDir',p,'Formats',{'csv','csv'}));
fails(@()export_qstab_paper_figures(F,'OutputDir',p,'Formats','exe'));
assert(exist(p,'dir')~=7&&isequal(old,findall(0,'Type','figure')));
end
function exportDirectory()
F=figdata('qb02_counter');base=tempname;mkdir(base);p=fullfile(base,'figs');cwd=pwd;
c=onCleanup(@()restoreAndRemove(cwd,base)); %#ok<NASGU>
cd(base);X=export_qstab_paper_figures(F,'OutputDir',p,'Formats','csv');
assert(strcmp(pwd,base)&&strcmp(X.outputDirectory,p)&&exist(X.files.file,'file')==2);
assert(exist(fullfile(base,'qstab_test_results.zip'),'file')==0);
end
function manifest()
F=figdata('qb03_tetrahedron');p=tempname;mkdir(p);c=onCleanup(@()rmdir(p,'s')); %#ok<NASGU>
X=export_qstab_paper_figures(F,'OutputDir',p,'Formats','csv');s=fileread(X.manifestFile);
assert(~isempty(strfind(s,'eq:facegap'))&&~isempty(strfind(s,F.toolbox.version))&&~isempty(strfind(s,'Data runtime:'))&&plain(X));
end
function roundTrip(writer)
base=tempname;mkdir(base);c=onCleanup(@()rmdir(base,'s')); %#ok<NASGU>
item=struct('group','fixture','name','Artificial extension persistence fixture','status','PASS','seconds',0,'message','','identifier','');
r=struct('toolbox',struct('version','fixture'), ...
 'environment',struct('version','artificial','computer','fixture','quaternionPath','','qstabPath',''), ...
 'timestamp','20000101T000000','runId','paper-extension-fixture','nativeMATLABRun',false, ...
 'nativeQuaternionOK',false,'ok',true,'status','OVERALL_OK','passed',1,'failed',0,'skipped',0,'elapsedSeconds',0,'cases',item);
r.paperExtensions=matched('qb03_ordered_product');r.paperFigures=figdata('qb03_parabola');assert(plain(r));
r=writer(r,base,'qstab_test_results');p=fullfile(base,'unpacked');mkdir(p);unzip(r.archiveFile,p);
f=dir(p);f=f(~ismember({f.name},{'.','..'}));assert(isequal(sort({f.name}),sort({'report.mat','report.txt'})));
s=load(fullfile(p,'report.mat'),'report');assert(isequaln(s.report.paperExtensions,r.paperExtensions)&&isequaln(s.report.paperFigures,r.paperFigures));
end
function b=readBytes(p),f=fopen(p,'rb');c=onCleanup(@()fclose(f));b=fread(f,Inf,'*uint8');end %#ok<NASGU>
function restoreAndRemove(cwd,p),cd(cwd);rmdir(p,'s');end
function tf=plain(x)
tf=isnumeric(x)||islogical(x)||ischar(x);
if isstruct(x)
 tf=true;f=fieldnames(x);for j=1:numel(x),for k=1:numel(f),if ~plain(x(j).(f{k})),tf=false;return;end,end,end
elseif iscell(x),tf=true;for j=1:numel(x),if ~plain(x{j}),tf=false;return;end,end
end
end
function fails(f),caught=false;try,f();catch,caught=true;end;assert(caught,'Expected invalid input rejection.');end
