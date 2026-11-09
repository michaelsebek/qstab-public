function cases=qstab_cases_fourgenerator()
%QSTAB_CASES_FOURGENERATOR Complete-subfamily reduction, not vertex tests.
cases={ ...
 'One to four generators retain the entire original family',@localSmall; ...
 'Five generators visit exactly five four-generator subfamilies',@localFive; ...
 'Six generators use all 15 lexicographic subsets without duplication',@localSix; ...
 'Noncommuting quadratic coefficients use full subfamily certificates',@localQuadratic; ...
 'A tetrahedron is not replaced by its stable triangular boundary',@localTetrahedron; ...
 'A subfamily counterexample reconstructs original five-generator weights',@localWitness; ...
 'Schur reduction tests all complete four-generator subfamilies',@localSchur; ...
 'Delay reversal agrees with the forward Schur reduction',@localDelay; ...
 'Declared delay grade is retained in every subfamily',@localGrade; ...
 'Delay witnesses distinguish input and tested coefficient variables',@localDelayWitness; ...
 'A common nonreal leading coefficient is preserved without normalization',@localNonrealLead; ...
 'Redundant and duplicate generators retain original indices',@localDuplicates; ...
 'MaxSubfamilies retains a compact exact unvisited suffix',@localSubfamilyBudget; ...
 'MaxCells is a global budget and unresolved weights map to the full input',@localCellBudget; ...
 'MaxDepth remains local and cannot produce a false complete certificate',@localDepthBudget; ...
 'MaxTime zero visits no subfamily and preserves the whole obligation',@localTimeBudget; ...
 'A 512-generator family is enumerated lazily under a one-subfamily budget',@localLazy; ...
 'Full leading-hull checks separate varying-leading and interior degree-loss families',@localVariableLead; ...
 'A zero leading generator is an admissible assumption witness',@localZeroLead; ...
 'Intervals and correlated affine families are not silently convexified',@localUnsupported; ...
 'Omitted and empty Uncertainty preserve four-generator inference',@localInference; ...
 'Method aliases preserve the reduction and budgets',@localDispatch; ...
 'Reduction none leaves the default zero-exclusion route unchanged',@localDefault; ...
 'Stored leaf certificates use original weights and cover every subfamily',@localCertificates; ...
 'Later stable subfamilies do not erase earlier boundary ambiguity',@localPersistentUnknown; ...
 'Malformed reduction names and budgets fail instead of falling back',@localInvalid; ...
 'Input generator cells and component arrays are not modified',@localNoMutation};
end

function R=localRun(V,varargin)
[~,R] = qstab.zeroexclude(V,'Reduction','four-generator','MaxTime',Inf,varargin{:});
end

function localSmall()
for m=1:4
    V=cell(1,m);
    for j=1:m, V{j}=[1+.1*j,1]; end
    R=localRun(V);
    localStable(R);
    assert(R.reduction.totalSubfamilies==1 && R.reduction.processedSubfamilies==1);
    assert(isequal(R.reduction.subfamilies{1}.generatorIndices,1:m));
    assert(R.reduction.parameterDimensionBound==m-1);
end
end

function localFive()
R=localRun(localVertices(5)); localStable(R);
assert(R.reduction.totalSubfamilies==5 && R.reduction.stableSubfamilies==5);
assert(R.reduction.fullVertexCount==5 && R.reduction.subsetSize==4);
localEnumeration(R,5);
end

function localSix()
R=localRun(localVertices(6)); localStable(R);
assert(R.reduction.totalSubfamilies==15 && R.reduction.stableSubfamilies==15);
localEnumeration(R,6);
cells=0; accepted=0;
for j=1:numel(R.reduction.subfamilies)
    S=R.reduction.subfamilies{j};
    cells=cells+S.processedCells; accepted=accepted+S.certifiedCells;
end
assert(R.processedCells==cells && R.certifiedCells==accepted);
end

function localQuadratic()
V=cell(1,5);
for j=1:5
    V{j}=localQ([2,3,1],[.01*j,0,0], ...
        [0,.01*(j-3),0],[.01*mod(j,2),0,0]);
end
R=localRun(V,'MaxCells',16384); localStable(R);
assert(R.degree==2 && R.reduction.stableSubfamilies==5);
end

function localTetrahedron()
V=localTetra();
% For p_y=s^2+(1+3y)s+y-1 the qB03 calculation gives stability iff
% |y|^2>1/2. All tetrahedral faces have |y|^2>=4/3, but y=0 is unstable.
assert(4/3>1/2);
for j=1:4
    [~,S] = qstab.isstable(V{j}); assert(strcmp(S.status,'stable'));
end
R=localRun(V);
assert(strcmp(R.status,'notStable') && isequal(R.isStable,false));
assert(R.reduction.subsetSize==4 && R.reduction.processedSubfamilies==1);
assert(isequal(R.witness.generatorIndices,1:4));
assert(max(abs(R.witness.originalWeights-.25))<1e-14);
end

function localWitness()
V=localTetra(); V{5}=[2,3,1];
R=localRun(V); assert(strcmp(R.status,'notStable') && ~R.certificateComplete);
assert(numel(R.witness.originalWeights)==5 && R.witness.originalWeights(5)==0);
assert(R.reduction.unvisitedSubfamilies==4);
localReconstruct(R,V);
[~,S] = qstab.isstable(localComponent(R.witness.originalComponents));
assert(strcmp(S.status,'notStable'));
end

function localSchur()
R=localRun(localDiscrete(5),'Variable','z'); localStable(R);
assert(strcmp(R.testVariable,'z') && abs(R.frequencyLimit-pi)<1e-14);
assert(R.reduction.stableSubfamilies==5);
end

function localDelay()
V=localDiscrete(5); D=localReverse(V);
A=localRun(V,'Variable','z'); B=localRun(D,'Variable','zinv');
localStable(A);localStable(B);
assert(A.processedCells==B.processedCells && A.certifiedCells==B.certifiedCells);
assert(strcmp(B.variable,'zinv') && strcmp(B.testVariable,'z'));
assert(isequal(A.reduction.subfamilies,B.reduction.subfamilies));
end

function localGrade()
V=cell(1,5);
for j=1:5, V{j}=[1,-(.1+.01*j),0]; end
R=localRun(V,'Variable','zinv'); localStable(R);
assert(R.degree==2 && R.nominal.grade==2);
F=qstab.familycomponents(V);
assert(size(F.vertices,2)==3);
end

function localDelayWitness()
V=repmat({[1,2,0]},1,5);
R=localRun(V,'Variable','zinv');
assert(strcmp(R.status,'notStable'));
assert(isequal(R.witness.originalComponents(1,:),[1,2,0]));
assert(isequal(R.witness.components(1,:),[0,2,1]));
assert(strcmp(R.witness.originalCoefficientVariable,'zinv'));
assert(strcmp(R.witness.coefficientVariable,'z'));
localReconstruct(R,V);
end

function localNonrealLead()
V=cell(1,5);
for j=1:5
    a=1+.1*j;
    V{j}=localQ([a,1],[a,1],[0,0],[0,0]);
end
R=localRun(V); localStable(R);
assert(R.reduction.stableSubfamilies==5);
end

function localDuplicates()
V=repmat({[1,1]},1,6);
R=localRun(V); localStable(R);
assert(R.reduction.fullVertexCount==6 && R.reduction.processedSubfamilies==15);
localEnumeration(R,6);
end

function localSubfamilyBudget()
R=localRun(localVertices(6),'MaxSubfamilies',1);
localUnknown(R);
assert(R.reduction.processedSubfamilies==1 && R.reduction.stableSubfamilies==1);
assert(R.reduction.unvisitedSubfamilies==14);
assert(isequal(R.reduction.nextGeneratorIndices,[1,2,3,5]));
S=R.unresolvedRegions{end};
assert(strcmp(S.kind,'unvisitedGeneratorSubfamilies') && S.count==14 && S.countExact);
assert(isequal(S.nextGeneratorIndices,[1,2,3,5]));
assert(~R.reduction.allVisited);
end

function localCellBudget()
R=localRun(localVertices(6),'MaxCells',1,'StoreCertificates',true);
localUnknown(R);
assert(R.processedCells<=1 && R.reduction.processedSubfamilies<=1);
assert(R.reduction.unvisitedSubfamilies>=14);
for j=1:numel(R.unresolvedRegions)
    U=R.unresolvedRegions{j};
    if isfield(U,'vertexWeights') && ~isempty(U.vertexWeights)
        assert(size(U.vertexWeights,1)==6);
        assert(max(abs(sum(U.vertexWeights,1)-1))<1e-12);
    end
end
end

function localDepthBudget()
V={localQ([1,1],[0,0],[-2,0],[0,0]), ...
   localQ([1,1],[0,0],[2,0],[0,0])};
V={V{1},V{2},V{1},V{2},V{1}};
R=localRun(V,'MaxDepth',0);
localUnknown(R);
assert(R.maximumDepth==0 && R.reduction.allVisited);
assert(R.reduction.unresolvedSubfamilies>0);
end

function localTimeBudget()
[~,R] = qstab.zeroexclude(localVertices(6),'Reduction','four-generator','MaxTime',0);
localUnknown(R);
assert(R.processedCells==0 && R.reduction.processedSubfamilies==0);
assert(R.reduction.unvisitedSubfamilies==15);
assert(isequal(R.reduction.nextGeneratorIndices,1:4));
end

function localLazy()
V=repmat({1},1,512);
R=localRun(V,'MaxSubfamilies',1);
localUnknown(R);
assert(R.reduction.totalSubfamilies==512*511*510*509/24);
assert(R.reduction.totalCountExact && R.reduction.processedSubfamilies==1);
assert(R.reduction.stableSubfamilies==1 && numel(R.reduction.subfamilies)==1);
assert(numel(R.unresolvedRegions)==1);
end

function localVariableLead()
A=localRun({[1,1],[2,2]}); localStable(A);
assert(A.reduction.processedSubfamilies==1&&A.leading.isSeparated);
% Five affinely independent leading quaternions enclose zero. No <=4
% generator face contains zero. A face-only degree check would be unsound.
Q=[eye(4),-ones(4,1)]; V=cell(1,5);
for j=1:5, V{j}=localComponent([Q(:,j),Q(:,j)]); end
F=qstab.familycomponents(V);
assert(norm(sum(F.vertices,3),'fro')==0);
R=localRun(V);assert(strcmp(R.status,'assumptionFailed'));
assert(R.reduction.processedSubfamilies==0);
assert(all(R.witness.originalComponents(:,end)==0));
end

function localZeroLead()
R=localRun({[1,1],[2,0]});
assert(strcmp(R.status,'assumptionFailed') && isempty(R.isStable));
assert(strcmp(R.witness.kind,'zeroLeadingCoefficient'));
assert(isequal(R.witness.originalWeights,[0;1]));
assert(all(R.witness.components(:,end)==0));
R=localRun({[1,1],[0,1]},'Variable','zinv');
assert(strcmp(R.status,'assumptionFailed') && all(R.witness.originalComponents(:,1)==0));
end

function localUnsupported()
[~,A] = qstab.zeroexclude({[1,1],[2,1]},'Uncertainty','interval', ...
    'Reduction','four-generator');localUnknown(A);
F=struct('center',[2,1],'generators',{{[.1,0]}},'bounds',[-1,1]);
[~,A] = qstab.zeroexclude(F,'Reduction','four-generator');localUnknown(A);
assert(A.reduction.processedSubfamilies==0 && strcmp(A.uncertainty,'affine'));
end

function localInference()
V=localVertices(5);
A=localRun(V); B=localRun(V,'Uncertainty','');
C=localRun(V,'Uncertainty','polytope');
localSame(A,B);localSame(A,C);
end

function localDispatch()
V=localVertices(5);
a={'Reduction','four-generator','MaxTime',Inf, ...
    'MaxCells',2048,'MaxDepth',30,'MaxSubfamilies',2,'StoreCertificates',true};
[~,A] = qstab.zeroexclude(V,a{:});
[~,B] = qstab.robuststable(V,'Method','zeroexclusion',a{:});
[~,C] = qstab.robuststable(V,'Method','zeroexclude',a{:});
localSame(A,B);localSame(A,C);
assert(B.options.MaxSubfamilies==2);
rejected=false;
try
    qstab.robuststable(V,'Reduction','four-generator');
catch ME
    rejected=strcmp(ME.identifier,'qstab:robuststable:ReductionMethod');
end
assert(rejected,'Do not silently apply or ignore reduction on the Lyapunov path.');
end

function localDefault()
V=localVertices(5);
[~,A] = qstab.zeroexclude(V,'MaxTime',Inf);
[~,B] = qstab.zeroexclude(V,'Reduction','none','MaxTime',Inf);
f={'status','isStable','method','witness','frequencyLimit', ...
   'processedCells','certifiedCells','maximumDepth','certificates','unresolvedRegions'};
for j=1:numel(f),assert(isequaln(A.(f{j}),B.(f{j})));end
assert(strcmp(A.method,'normValueZeroExclusion') && ~isfield(A,'reduction'));
end

function localCertificates()
V=repmat({[1,1]},1,5);
R=localRun(V,'StoreCertificates',true);localStable(R);
assert(numel(R.certificates)==R.certifiedCells && ~isempty(R.certificates));
for t=1:numel(R.reduction.subfamilies)
    indices=R.reduction.subfamilies{t}.generatorIndices;
    intervals=[]; count=0;
    for j=1:numel(R.certificates)
        C=R.certificates{j};
        if ~isequal(C.generatorIndices,indices),continue;end
        W=C.cell.vertexWeights;
        assert(size(W,1)==5 && all(W(:)>=0));
        assert(max(abs(sum(W,1)-1))<1e-12);
        absent=setdiff(1:5,indices);
        assert(all(all(W(absent,:)==0)));
        intervals(end+1,:)=C.cell.frequencyBounds; %#ok<AGROW>
        count=count+1;
    end
    % Identical generators do not require parameter bisection, so the
    % frequency leaves alone form a partition for each generator subset.
    intervals=sortrows(intervals,1);
    assert(count==R.reduction.subfamilies{t}.certifiedCells);
    assert(intervals(1,1)==0);
    assert(abs(intervals(end,2)-R.reduction.subfamilies{t}.frequencyLimit)<1e-12);
    assert(all(abs(intervals(1:end-1,2)-intervals(2:end,1))<1e-12));
end
end

function localPersistentUnknown()
V={[0,1],[2,1],[2,1],[2,1],[2,1]};
R=localRun(V,'MaxDepth',8,'MaxCells',8192);
localUnknown(R);
assert(R.reduction.allVisited && R.reduction.unresolvedSubfamilies>0);
assert(R.reduction.stableSubfamilies>0);
assert(strcmp(R.reduction.subfamilies{end}.status,'stable'));
end

function localInvalid()
V=localVertices(5);
bad={'',[],0,{'four-generator'},'edges','vertices',char('aa','bb')};
for j=1:numel(bad)
    rejected=false;
    try,qstab.zeroexclude(V,'Reduction',bad{j});catch,rejected=true;end
    assert(rejected,'Malformed reduction must be rejected.');
end
bad={0,-1,.5,NaN,Inf,[1,2]};
for j=1:numel(bad)
    rejected=false;
    try,localRun(V,'MaxSubfamilies',bad{j});catch,rejected=true;end
    assert(rejected,'MaxSubfamilies must be a positive finite integer.');
end
end

function localNoMutation()
V=localVertices(5); saved=V;
R=localRun(V,'MaxSubfamilies',1); %#ok<NASGU>
assert(isequaln(saved,V));
end

function V=localVertices(m)
V=cell(1,m);
for j=1:m
    V{j}=localQ([1+.1*j,1],[.02*j,0], ...
        [.03*(-1)^j,0],[.01*mod(j,3),0]);
end
end

function V=localDiscrete(m)
V=cell(1,m);
for j=1:m
    V{j}=localQ([-.2-.01*j,1],[.01*j,0], ...
        [.02*(-1)^j,0],[.01*mod(j,3),0]);
end
end

function V=localTetra()
signRows=[1,1,1;1,-1,-1;-1,1,-1;-1,-1,1];
V=cell(1,4);
for j=1:4
    y=2*signRows(j,:);
    V{j}=localQ([-1,1,1],[y(1),3*y(1),0], ...
        [y(2),3*y(2),0],[y(3),3*y(3),0]);
end
end

function D=localReverse(V)
D=cell(size(V));
for j=1:numel(V)
    C=qstab.polycomponents(V{j});D{j}=localComponent(C(:,end:-1:1));
end
end

function localReconstruct(R,V)
F=qstab.familycomponents(V); C=zeros(4,F.degree+1);
w=R.witness.originalWeights;
assert(all(w>=0) && abs(sum(w)-1)<1e-12);
for j=1:numel(w),C=C+w(j)*F.vertices(:,:,j);end
assert(norm(C-R.witness.originalComponents,'fro')<1e-12);
Ct=C;
if strcmp(R.variable,'zinv'),Ct=C(:,end:-1:1);end
assert(norm(Ct-R.witness.components,'fro')<1e-12);
end

function localEnumeration(R,m)
expected=nchoosek(1:m,4); actual=zeros(numel(R.reduction.subfamilies),4);
for j=1:size(actual,1),actual(j,:)=R.reduction.subfamilies{j}.generatorIndices;end
assert(isequal(expected,actual));
assert(R.reduction.totalCountExact && R.reduction.allVisited);
end

function localSame(A,B)
f={'status','isStable','method','variable','testVariable','degree', ...
 'certificateComplete','processedCells','certifiedCells','unresolvedCells', ...
 'maximumDepth','frequencyLimit','witness','certificates','unresolvedRegions','reduction'};
for j=1:numel(f)
    assert(isequaln(A.(f{j}),B.(f{j})),['Dispatch changed ',f{j}]);
end
end

function localStable(R)
assert(strcmp(R.status,'stable') && isequal(R.isStable,true),R.message);
assert(R.certificateComplete && R.reduction.complete && R.reduction.allVisited);
assert(R.reduction.unresolvedSubfamilies==0 && R.reduction.unvisitedSubfamilies==0);
assert(R.unresolvedCells==0 && isempty(R.unresolvedRegions) && ~R.validated);
end

function localUnknown(R)
assert(strcmp(R.status,'undetermined') && isempty(R.isStable),R.message);
assert(~R.certificateComplete && ~R.reduction.complete && ~R.validated);
assert(~isempty(R.unresolvedRegions));
end

function S=localComponent(C)
S=localQ(reshape(C(1,:),1,1,[]),reshape(C(2,:),1,1,[]), ...
    reshape(C(3,:),1,1,[]),reshape(C(4,:),1,1,[]));
end

function S=localQ(a,b,c,d)
S=struct('a',a,'b',b,'c',c,'d',d);
end
