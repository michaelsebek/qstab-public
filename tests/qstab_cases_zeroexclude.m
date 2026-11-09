function cases = qstab_cases_zeroexclude()
%QSTAB_CASES_ZEROEXCLUDE qB02/qB03 norm-value exclusion regressions.
% Each case is executable with numeric arrays/component structs only.
% Stable means every frequency/parameter cell has a sufficient exclusion;
% no sample-only/vertex-only claim is used.
cases={ ...
    'qB03 Bernstein polygon certifies a real Hurwitz segment',@localRealSegment; ...
    'qB03 false Bernstein alarm is removed by parameter subdivision',@localBernsteinFalseAlarm; ...
    'qB03 sharp tetrahedron: stable faces but unstable centre',@localTetrahedron; ...
    'qB03 quaternion Schur polytope has a full-circle certificate',@localSchur; ...
    'qB03 delay convention reverses the same polynomial polytope',@localDelay; ...
    'qB02 exact even-odd Cartesian boxes certify Hurwitz interval',@localInterval; ...
    'qB02 permits a varying nonzero leading coefficient',@localVariableLeadingInterval; ...
    'qB02 zero leading coefficient is an assumption witness',@localIntervalDegreeLoss; ...
    'An unstable interval midpoint is a real admissible witness',@localUnstableInterval; ...
    'A varying-leading polytope uses a full-hull degree certificate',@localVariableLeadingPoly; ...
    'Zero-exclusion budgets preserve unresolved cells, not false stability',@localBudget; ...
    'The robuststable Method selector preserves default and new paths',@localSelector; ...
    'A discrete interval with an exact unit-root midpoint is not strictly stable',@localUnsupported; ...
    'Retained Bernstein certificates cover all accepted leaves',@localRetained; ...
    'Omitted and explicitly empty Uncertainty have identical inference',@localEmptyUncertainty; ...
    'Method aliases and explicit interval dispatch preserve options and budgets',@localDispatchOptions; ...
    'Empty Uncertainty does not admit numeric, cell or matrix options',@localBadUncertainty; ...
};
end

function localRealSegment()
[~,R] = qstab.zeroexclude({[1 1],[3 1]},'Uncertainty','polytope','Variable','s','MaxTime',Inf);
localStable(R);
assert(strcmp(R.method,'normValueZeroExclusion') && R.degree==1);
assert(R.processedCells>0 && R.frequencyLimit>=4);
end

function localBernsteinFalseAlarm()
% qB03: p_t=s+1+2*t*j, -1<=t<=1. Coarse control at w=0
% includes -3 and 5, while the actual norm set is [1,5].
V={localQ([1,1],[0,0],[-2,0],[0,0]), ...
   localQ([1,1],[0,0],[2,0],[0,0])};
[~,R] = qstab.zeroexclude(V,'Variable','s','MaxTime',Inf,'MaxCells',1024);
localStable(R);
assert(R.processedCells>1 && R.maximumDepth>0, ...
    'The coarse control hull cannot already certify this family.');
end

function localTetrahedron()
% qB03 master family p_y=s^2+(1+3*y)*s+(y-1), pure y.
% With four tetrahedral vertices y=2*(+,+,+),2*(+,-,-),
% 2*(-,+,-),2*(-,-,+), every face has |y|^2 >= 4/3 > 1/2
% and is Hurwitz, but the centre y=0 has an unstable real root.
sgn=[1,1,1;1,-1,-1;-1,1,-1;-1,-1,1];
V=cell(1,4);
for j=1:4
    y=2*sgn(j,:);
    V{j}=localQ([-1,1,1],[y(1),3*y(1),0], ...
       [y(2),3*y(2),0],[y(3),3*y(3),0]);
end
for skip=1:4
    idx=setdiff(1:4,skip);
    F=qstab.familycomponents(V(idx),'Uncertainty','polytope');
    C=mean(F.vertices,3);
    nominal=localComponent(C);
    [~,S] = qstab.isstable(nominal,'Variable','s');
    assert(strcmp(S.status,'stable'),'A face-centre should be stable.');
end
[~,R] = qstab.zeroexclude(V,'Variable','s','MaxCells',64,'MaxTime',Inf);
assert(strcmp(R.status,'notStable') && isequal(R.isStable,false));
assert(strcmp(R.witness.kind,'actualAdmissiblePolynomial'));
assert(max(abs(R.witness.originalWeights-1/4))<1e-14);
end

function localSchur()
V={localQ([-.3,1],[0,0],[-.2,0],[0,0]), ...
   localQ([-.3,1],[0,0],[ .2,0],[0,0])};
[~,R] = qstab.zeroexclude(V,'Variable','z','MaxTime',Inf);
localStable(R);
assert(abs(R.frequencyLimit-pi)<1e-14);
end

function localDelay()
V={localQ([1,-.3],[0,0],[0,-.2],[0,0]), ...
   localQ([1,-.3],[0,0],[0, .2],[0,0])};
[~,R] = qstab.zeroexclude(V,'Variable','zinv','MaxTime',Inf);
localStable(R);
assert(strcmp(R.variable,'zinv') && strcmp(R.testVariable,'z'));
end

function localInterval()
L=localQ([1,1],[-.05,0],[-.05,0],[-.05,0]);
U=localQ([1.2,1],[.05,0],[.05,0],[.05,0]);
[~,R] = qstab.zeroexclude({L,U},'Uncertainty','interval', ...
  'Variable','s','MaxTime',Inf);
localStable(R);
assert(R.frequencyLimit>0 && R.processedCells>0);
end

function localVariableLeadingInterval()
[~,R] = qstab.zeroexclude({[2,1],[2,2]},'Uncertainty','interval', ...
  'Variable','s','MaxTime',Inf);
localStable(R);
end

function localIntervalDegreeLoss()
[~,R] = qstab.zeroexclude({[2,-1],[2,1]},'Uncertainty','interval','Variable','s');
assert(strcmp(R.status,'assumptionFailed') && isempty(R.isStable));
assert(all(R.witness.components(:,end)==0));
end

function localUnstableInterval()
[~,R] = qstab.zeroexclude({[-1.2,1],[-.8,1]},'Uncertainty','interval','Variable','s');
assert(strcmp(R.status,'notStable') && ~R.isStable);
assert(R.witness.components(1,1)<0 && R.witness.components(1,2)==1);
end

function localVariableLeadingPoly()
[~,R] = qstab.zeroexclude({[1,1],[2,2]},'Variable','s');
localStable(R);
assert(R.leading.isSeparated&&R.leading.lowerBound>0);
end

function localBudget()
V={localQ([1,1],[0,0],[-2,0],[0,0]), ...
   localQ([1,1],[0,0],[ 2,0],[0,0])};
[~,R] = qstab.zeroexclude(V,'Variable','s','MaxCells',1,'MaxTime',Inf);
assert(strcmp(R.status,'undetermined') && ~R.certificateComplete);
assert(R.unresolvedCells>0 && ~isempty(R.unresolvedRegions));
end

function localSelector()
V={localQ([2,1],[0,0],[0,0],[0,0]), ...
   localQ([2,1],[0,0],[.1,0],[0,0])};
[~,R] = qstab.robuststable(V,'Variable','s','Method','zeroexclusion', ...
  'MaxTime',Inf,'StoreCertificates',true);
localStable(R);
assert(numel(R.certificates)==R.certifiedCells);
[~,Q] = qstab.robuststable(V,'Variable','s','MaxTime',Inf,'MaxCells',128);
localStable(Q);
assert(strcmp(Q.method,'subdivisionLyapunovPerturbation'));
end

function localUnsupported()
[tf,R] = qstab.zeroexclude({[.9,1],[1.1,1]},'Uncertainty','interval','Variable','z');
assert(tf==0 && strcmp(R.status,'notStable') && R.witness.nominal.diagnostics.boundaryWitness.exactForSuppliedComponents);
end

function localRetained()
[~,R] = qstab.zeroexclude({[1,1],[3,1]},'Variable','s', ...
 'StoreCertificates',true,'MaxTime',Inf);
localStable(R);
assert(numel(R.certificates)==R.certifiedCells);
for j=1:numel(R.certificates)
    c=R.certificates{j};
    assert(strcmp(c.cell.kind,'polytope'));
    assert(c.cell.frequencyBounds(1)<=c.cell.frequencyBounds(2));
end
end

function localStable(R)
assert(strcmp(R.status,'stable') && isequal(R.isStable,true) && ...
   R.certificateComplete && ~R.validated && R.unresolvedCells==0, ...
   ['No complete sufficient zero-exclusion certificate: ',R.message]);
end
function P=localQ(a,b,c,d)
P=struct('a',a,'b',b,'c',c,'d',d);
end
function S=localComponent(C)
S=struct('a',reshape(C(1,:),1,1,[]), ...
'b',reshape(C(2,:),1,1,[]),'c',reshape(C(3,:),1,1,[]), ...
'd',reshape(C(4,:),1,1,[]));
end


function localEmptyUncertainty()
V={[1,1],[3,1]};
args={'Variable','s','MaxTime',Inf,'MaxCells',2048,'MaxDepth',30, ...
      'Tolerance',1e-10,'StoreCertificates',true};
[~,A] = qstab.zeroexclude(V,args{:});
[~,B] = qstab.zeroexclude(V,'Uncertainty','',args{:});
[~,C] = qstab.zeroexclude(V,'Uncertainty','polytope',args{:});
[~,D] = qstab.robuststable(V,'Method','zeroexclusion',args{:});
[~,E] = qstab.robuststable(V,'Method','zeroexclusion','Uncertainty','',args{:});
localStable(A);
localSameCertificate(A,B);localSameCertificate(A,C);
localSameCertificate(A,D);localSameCertificate(A,E);
% MATLAB string scalar support is checked in the native run as well.
if exist('OCTAVE_VERSION','builtin')==0
    [~,Q] = qstab.zeroexclude(V,'Uncertainty',string(''),args{:});
    localSameCertificate(A,Q);
end
end

function localDispatchOptions()
V={localQ([1,1],[0,0],[-2,0],[0,0]), ...
   localQ([1,1],[0,0],[ 2,0],[0,0])};
args={'Variable','s','Uncertainty','polytope','MaxTime',Inf, ...
      'MaxCells',1,'MaxDepth',3,'Tolerance',1e-9,'StoreCertificates',true};
[~,A] = qstab.zeroexclude(V,args{:});
[~,B] = qstab.robuststable(V,'Method','zeroexclude',args{:});
[~,C] = qstab.robuststable(V,'Method','zeroexclusion',args{:});
assert(strcmp(A.status,'undetermined')&&~A.certificateComplete);
localSameCertificate(A,B);localSameCertificate(A,C);
assert(B.options.MaxCells==1&&B.options.MaxDepth==3&& ...
       B.options.Tolerance==1e-9&&B.options.StoreCertificates);
I={[1,1],[1.2,1]};
a={'Uncertainty','interval','Variable','s','MaxTime',Inf, ...
   'MaxCells',2048,'MaxDepth',30,'StoreCertificates',true};
[~,A] = qstab.zeroexclude(I,a{:});
[~,B] = qstab.robuststable(I,'Method','zeroexclusion',a{:});
localStable(A);localSameCertificate(A,B);
end

function localBadUncertainty()
V={[1,1],[3,1]};
bad={[],0,{'polytope'},char('aa','bb')};
for k=1:numel(bad)
    didFail=false;
    try
        qstab.zeroexclude(V,'Uncertainty',bad{k});
    catch
        didFail=true;
    end
    assert(didFail,'Malformed Uncertainty must be rejected.');
end
% Empty Variable is still invalid: the fix is specific to Uncertainty.
didFail=false;
try
    qstab.zeroexclude(V,'Variable','');
catch
    didFail=true;
end
assert(didFail,'Empty Variable must not acquire an implicit meaning.');
end

function localSameCertificate(A,B)
fields={'status','isStable','method','uncertainty','variable','testVariable', ...
    'degree','certificateComplete','processedCells','certifiedCells', ...
    'unresolvedCells','maximumDepth','frequencyLimit','witness', ...
    'certificates','unresolvedRegions'};
for k=1:numel(fields)
    f=fields{k};
    assert(isequaln(A.(f),B.(f)),['Dispatch changed result field: ',f]);
end
end
