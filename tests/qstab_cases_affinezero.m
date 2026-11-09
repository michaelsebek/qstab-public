function cases = qstab_cases_affinezero()
%QSTAB_CASES_AFFINEZERO Joint parameters on the axis and unit circle.
cases={ ...
 'A correlated Hurwitz affine box has a whole-cell norm certificate',@hurwitz; ...
 'Shared leading and constant coefficients are not independently recombined',@strongCorrelation; ...
 'Discrete Cartesian quaternion intervals use the joint parameter image',@discrete; ...
 'Delay reversal preserves the affine grade and forward zero mode',@delay; ...
 'Noncommuting one-quaternion sides reconstruct distinct actual witnesses',@sides; ...
 'An unstable discrete interval returns admissible original coefficients',@unstable; ...
 'An exact unit-circle midpoint disproves strict stability; nearby ambiguity stays unresolved',@boundary; ...
 'Nonzero correlated constant families need no frequency sweep',@constant; ...
 'Affine MaxCells preserves covering unresolved parameter regions',@cellBudget; ...
 'Affine MaxTime and MaxDepth never create a false complete certificate',@otherBudgets; ...
 'Stored affine leaves cover the complete parameter-frequency domain',@coverage; ...
 'Zero-exclusion dispatch preserves affine options and the original input',@dispatch};
end
function F=example()
F=struct('center',[2,1],'generators',{{poly([.1,0;0,0;.2,0;0,0])}},'bounds',[-1,1]);
end
function hurwitz()
[~,R] = qstab.zeroexclude(example(),'MaxTime',Inf);stable(R);
assert(strcmp(R.method,'affineNormValueZeroExclusion')&&strcmp(R.uncertainty,'affine'));
end
function strongCorrelation()
F=struct('center',poly([0,0;1,1;0,0;0,0]),'generators',{{[1,1]}},'bounds',[-2,2]);
[~,R] = qstab.zeroexclude(F,'MaxTime',Inf,'MaxCells',12000);stable(R);
% Independent mixing theta0=-2 and theta1=2 would instead be unstable.
[~,Q] = qstab.isstable(poly([-2,2;1,1;0,0;0,0]));assert(strcmp(Q.status,'notStable'));
end
function [L,U]=dtBox()
L=poly([-.3,1;-.05,0;-.1,0;-.05,0]);U=poly([-.1,1;.05,0;.1,0;.05,0]);
end
function discrete()
[L,U]=dtBox();[~,R] = qstab.zeroexclude({L,U},'Uncertainty','interval','Variable','z','MaxTime',Inf);
stable(R);assert(abs(R.frequencyLimit-pi)<1e-14);
end
function delay()
[L,U]=dtBox();a=qstab.polycomponents(L);b=qstab.polycomponents(U);
a=[zeros(4,1),a];b=[zeros(4,1),b];
[~,R] = qstab.zeroexclude({poly(a(:,end:-1:1)),poly(b(:,end:-1:1))}, ...
 'Uncertainty','interval','Variable','zinv','MaxTime',Inf);
stable(R);assert(R.degree==2&&strcmp(R.testVariable,'z'));
end
function sides()
F=struct('p0',[-1,1],'p1',poly([0;0;1;0]), ...
 'side','left','bounds',[0,0;1,1;0,0;0,0]);
[~,L] = qstab.zeroexclude(F,'MaxTime',Inf);F.side='right';[~,R] = qstab.zeroexclude(F,'MaxTime',Inf);
assert(strcmp(L.status,'notStable')&&strcmp(R.status,'notStable'));
assert(L.witness.originalComponents(4,1)==1&&R.witness.originalComponents(4,1)==-1);
assert(isequal(L.witness.parameter,[0;1;0;0]));
end
function unstable()
[~,R] = qstab.zeroexclude({[1.4,1],[1.6,1]},'Uncertainty','interval','Variable','z');
assert(strcmp(R.status,'notStable')&&~R.isStable);
c=R.witness.originalComponents(1,1);assert(c>=1.4&&c<=1.6);
end
function boundary()
[tf,R] = qstab.zeroexclude({[.9,1],[1.1,1]},'Uncertainty','interval','Variable','z');
assert(tf==0 && strcmp(R.status,'notStable') && isequal(R.witness.originalComponents(1,:),[1 1]));
A=[1+1e-12 1];[u,S]=qstab.zeroexclude({A,A},'Uncertainty','interval','Variable','z');
assert(isnan(u) && strcmp(S.status,'undetermined') && ~S.certificateComplete);
end
function constant()
F=struct('center',2,'generators',{{.2}},'bounds',[-1,1]);
[~,R] = qstab.zeroexclude(F);stable(R);assert(R.processedCells==0);
end
function cellBudget()
[~,R] = qstab.zeroexclude(example(),'MaxCells',1,'MaxTime',Inf);
assert(strcmp(R.status,'undetermined')&&~R.certificateComplete&&R.processedCells==1);
assert(~isempty(R.unresolvedRegions));
end
function otherBudgets()
[~,A] = qstab.zeroexclude(example(),'MaxTime',0);[~,B] = qstab.zeroexclude(example(),'MaxDepth',0,'MaxTime',Inf);
assert(strcmp(A.status,'undetermined')&&A.processedCells==0&&A.unresolvedCells>0);
assert(strcmp(B.status,'undetermined')&&B.unresolvedCells>0);
end
function coverage()
F=example();[~,R] = qstab.zeroexclude(F,'MaxTime',Inf,'StoreCertificates',true);stable(R);
assert(numel(R.certificates)==R.certifiedCells);
for t=linspace(-1,1,7)
 for f=linspace(0,R.frequencyLimit,11)
  covered=false;
  for k=1:numel(R.certificates)
   c=R.certificates{k};B=c.cell.parameterBounds;W=c.cell.frequencyBounds;
   covered=covered||(t>=B(1,1)&&t<=B(1,2)&&f>=W(1)&&f<=W(2));
   assert(c.diagnostic.gap>0);
  end
  assert(covered);
 end
end
end
function dispatch()
F=example();old=F;args={'MaxCells',1,'MaxTime',Inf,'StoreCertificates',true};
[~,A] = qstab.zeroexclude(F,args{:});[~,B] = qstab.robuststable(F,'Method','zeroexclusion',args{:});
assert(isequaln(F,old)&&strcmp(A.status,B.status)&&A.processedCells==B.processedCells);
assert(isequaln(A.unresolvedRegions,B.unresolvedRegions)&&B.options.MaxCells==1);
end
function P=poly(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
 'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
function stable(R)
assert(strcmp(R.status,'stable')&&R.certificateComplete&&~R.validated&&R.unresolvedCells==0,R.message);
end
