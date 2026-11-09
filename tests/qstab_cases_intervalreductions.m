function cases = qstab_cases_intervalreductions()
%QSTAB_CASES_INTERVALREDUCTIONS Whole qB02 cube/face subfamilies.
cases={ ...
 'Kharitonov cube reproduces the exact even-odd box at fixed frequencies',@identity; ...
 'Kharitonov cube uses ordinary affine data and preserves its input',@ordinary; ...
 'The whole Kharitonov cube certifies a stable higher-degree interval',@cubeStable; ...
 'Kharitonov conversion never applies its CT identity to a unit circle',@wrongVariable; ...
 'Eight active cube coordinates give exactly 1024 pair-free three-faces',@fullCount; ...
 'A paired even-odd component is reduced to four complete edges',@paired; ...
 'Three independent components retain the complete three-dimensional face',@three; ...
 'Monic quadratic reduction has at most 192 complete two-faces',@quadraticCount; ...
 'Quadratic reduction fixes the minimal real constant and retains imaginary vertices',@quadraticStable; ...
 'Stable quadratic corners do not hide an unstable two-face centre',@quadraticUnstable; ...
 'Nonmonic, wrong-degree and discrete quadratic reductions remain unsupported',@unsupported; ...
 'Face MaxCells and MaxSubfamilies are global budgets with exact pending faces',@budgets; ...
 'Face certificates retain source parameter bounds and face identifiers',@certificates; ...
 'An original interval degree loss is caught before any Kharitonov face',@degreeLoss; ...
 'robuststable forwards all three interval reduction names without changing results',@dispatch; ...
 'All 256 active-coordinate masks have unique pair-free face descriptors',@allMasks; ...
 'An entirely degenerate interval is a single continuous subfamily',@degenerate; ...
 'Noninterval families are not silently reinterpreted as independent boxes',@noninterval};
end
function [L,U]=smallBox()
L=poly([1.9,2.9,1;0,0,0;0,0,0;0,0,0]);
U=poly([2.1,3.1,1;.02,.02,0;.02,.02,0;0,0,0]);
end
function identity()
L=[1,2,3,4,1;-.2,-.3,-.4,-.5,0;-.4,-.1,-.2,-.3,0;-.1,-.2,-.3,-.4,0];
U=L+[.3,.2,.4,.1,0;.1,.2,.3,.4,0;.2,.3,.1,.4,0;.4,.1,.2,.3,0];
[K,I]=qstab.kharitonovcube({poly(L),poly(U)});F=qstab.familycomponents(K);
assert(I.parameterCount==8&&I.sameFixedFrequencyValueSets&&~I.vertexOnlyTest);
t=(1:8)'/9;
for w=[0,.5,2]
 C=F.base;for j=1:8,C=C+t(j)*F.generators(:,:,j);end
 y=C*(1i*w).^(0:4).';lo=zeros(8,1);hi=lo;
 for k=0:4
  f=(-1)^floor(k/2)*w^k;idx=(1:4)+4*mod(k,2);
  lo(idx)=lo(idx)+min(f*L(:,k+1),f*U(:,k+1));
  hi(idx)=hi(idx)+max(f*L(:,k+1),f*U(:,k+1));
 end
 target=lo+t.*(hi-lo);
 assert(norm(real(y)-target(1:4))<1e-10&&norm(imag(y)-target(5:8))<1e-10);
 assert(all(C(:)>=L(:)-1e-14)&&all(C(:)<=U(:)+1e-14));
end
end
function ordinary()
[L,U]=smallBox();oldL=L;oldU=U;[K,I]=qstab.kharitonovcube({L,U});
assert(isstruct(K)&&iscell(K.generators)&&numel(K.generators)==8&&I.degree==2);
assert(isequaln(L,oldL)&&isequaln(U,oldU));
F=qstab.familycomponents(K);C=F.base+.5*sum(F.generators,3);
assert(norm(C-.5*(qstab.polycomponents(L)+qstab.polycomponents(U)),'fro')<1e-14);
end
function cubeStable()
L=poly([5.98,10.98,5.98,1;-.01,-.01,-.01,0;-.01,-.01,-.01,0;-.01,-.01,-.01,0]);
U=poly([6.02,11.02,6.02,1;.01,.01,.01,0;.01,.01,.01,0;.01,.01,.01,0]);
R=run({L,U},'kharitonov-cube');stable(R);
assert(R.degree==3);
assert(R.reduction.totalSubfamilies==1&&R.reduction.complete);
end
function wrongVariable()
caught=false;try,qstab.kharitonovcube({[1,1],[2,1]},'Variable','z');catch,caught=true;end
assert(caught);
[~,R] = qstab.zeroexclude({[1,1],[2,1]},'Uncertainty','interval', ...
 'Reduction','kharitonov-faces','Variable','zinv');unknown(R);
end
function [L,U]=fullBox()
L=poly([2,3,1;-.01,-.01,0;-.01,-.01,0;-.01,-.01,0]);
U=poly([2.01,3.01,1;.01,.01,0;.01,.01,0;.01,.01,0]);
end
function fullCount()
[L,U]=fullBox();R=run({L,U},'kharitonov-faces','MaxSubfamilies',1);
unknown(R);assert(R.reduction.totalSubfamilies==1024&&R.reduction.processedSubfamilies==1);
assert(R.reduction.unvisitedSubfamilies==1023);
f=R.reduction.subfamilies{1}.freeParameters;
assert(numel(f)==3&&numel(unique(mod(f-1,4)))==3);
end
function paired()
R=run({[1,1],[2,2]},'kharitonov-faces');stable(R);
assert(R.reduction.totalSubfamilies==4&&R.reduction.parameterDimensionBound==1);
end
function three()
L=poly([2,1;-.1,0;-.1,0;-.1,0]);U=poly([2,1;.1,0;.1,0;.1,0]);
R=run({L,U},'kharitonov-faces');stable(R);
assert(R.reduction.totalSubfamilies==1&&R.reduction.parameterDimensionBound==3);
end
function quadraticCount()
[L,U]=fullBox();R=run({L,U},'quadratic-faces','MaxTime',0);unknown(R);
assert(R.reduction.totalSubfamilies==192&&R.reduction.processedSubfamilies==0);
assert(R.reduction.unvisitedSubfamilies==192&&~isempty(R.unresolvedRegions));
end
function quadraticStable()
L=poly([2,3,1;-.01,0,0;0,-.02,0;0,-.02,0]);
U=poly([4,3,1;.01,0,0;0,.02,0;0,.02,0]);
R=run({L,U},'quadratic-faces');stable(R);assert(R.reduction.totalSubfamilies==2);
for j=1:numel(R.reduction.subfamilies)
 B=R.reduction.subfamilies{j}.parameterBounds;
 F=qstab.familycomponents({L,U},'Uncertainty','interval');
 realB=find(F.parameterIndex(:,1)==1&F.parameterIndex(:,2)==0);
 assert(B(realB,1)==2&&B(realB,2)==2);
end
end
function quadraticUnstable()
L=poly([0,1,1;2,1,0;0,-2,0;0,-2,0]);
U=poly([0,1,1;2,1,0;0,2,0;0,2,0]);
for x=[-2,2]
 for y=[-2,2]
  P=poly([0,1,1;2,1,0;0,x,0;0,y,0]);[~,S] = qstab.isstable(P);assert(strcmp(S.status,'stable'));
 end
end
R=run({L,U},'quadratic-faces');
assert(strcmp(R.status,'notStable')&&~R.isStable&&R.reduction.totalSubfamilies==1);
assert(R.witness.originalComponents(3,2)==0&&R.witness.originalComponents(4,2)==0);
end
function unsupported()
R=run({[2,3,2],[3,4,2]},'quadratic-faces');unknown(R);
R=run({[1,1],[2,1]},'quadratic-faces');unknown(R);
[L,U]=smallBox();R=run({L,U},'quadratic-faces','Variable','z');unknown(R);
end
function budgets()
[L,U]=fullBox();A=run({L,U},'kharitonov-faces','MaxCells',1);unknown(A);
assert(A.processedCells==1&&A.reduction.processedSubfamilies==1&&A.reduction.unvisitedSubfamilies==1023);
B=run({L,U},'quadratic-faces','MaxSubfamilies',1);unknown(B);
assert(B.reduction.processedSubfamilies==1&&B.reduction.unvisitedSubfamilies==191);
end
function certificates()
R=run({[1,1],[2,2]},'kharitonov-faces','StoreCertificates',true);stable(R);
assert(numel(R.certificates)==R.certifiedCells);
for j=1:numel(R.certificates)
 c=R.certificates{j};assert(c.faceIndex>=1&&c.faceIndex<=4);
 assert(isequal(size(c.cell.parameterBounds),[8,2])&&c.diagnostic.gap>0);
end
end
function degreeLoss()
for mode={'kharitonov-cube','kharitonov-faces'}
 R=run({[1,-1],[2,1]},mode{1});
 assert(strcmp(R.status,'assumptionFailed')&&R.reduction.processedSubfamilies==0);
 assert(all(R.witness.originalComponents(:,end)==0));
end
end
function dispatch()
[L,U]=smallBox();V={L,U};
for mode={'kharitonov-cube','kharitonov-faces','quadratic-faces'}
 A=run(V,mode{1},'MaxCells',1);
 [~,B] = qstab.robuststable(V,'Uncertainty','interval','Method','zeroexclude', ...
  'Reduction',mode{1},'MaxCells',1,'MaxTime',Inf);
 assert(strcmp(A.status,B.status)&&A.processedCells==B.processedCells);
 assert(isequaln(A.reduction,B.reduction));
end
end
function allMasks()
for mask=0:255
 L=zeros(4,3);L(1,:)=[2,3,1];U=L;
 for j=1:8
  if bitget(mask,j)
   nu=mod(j-1,4)+1;k=1+(j>4);U(nu,k)=U(nu,k)+.01;
  end
 end
 R=run({poly(L),poly(U)},'kharitonov-faces','MaxTime',0);
 assert(R.reduction.totalSubfamilies>=1&&R.reduction.totalSubfamilies<=1024);
 faces=R.unresolvedRegions{end}.faces;keys=zeros(numel(faces),16);
 for j=1:numel(faces)
  t=faces{j};f=t.freeParameters;
  assert(numel(f)<=3&&numel(unique(mod(f-1,4)))==numel(f));
  keys(j,:)=t.bounds(:).';
 end
 assert(size(unique(keys,'rows'),1)==numel(faces));
end
end
function degenerate()
for mode={'kharitonov-cube','kharitonov-faces','quadratic-faces'}
 R=run({[2,3,1],[2,3,1]},mode{1});stable(R);
 assert(R.reduction.totalSubfamilies==1);
end
end
function noninterval()
F=struct('center',[2,1],'generators',{{[.1,0]}},'bounds',[-1,1]);
[~,R] = qstab.zeroexclude(F,'Reduction','kharitonov-cube');unknown(R);
[~,R] = qstab.zeroexclude({[1,1],[2,1]},'Reduction','quadratic-faces');unknown(R);
end
function R=run(V,mode,varargin)
[~,R] = qstab.zeroexclude(V,'Uncertainty','interval','Reduction',mode, ...
 'MaxTime',Inf,varargin{:});
end
function stable(R)
assert(strcmp(R.status,'stable')&&R.certificateComplete&&R.reduction.complete&&~R.validated,R.message);
end
function unknown(R)
assert(strcmp(R.status,'undetermined')&&isempty(R.isStable)&&~R.certificateComplete,R.message);
end
function P=poly(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
 'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
