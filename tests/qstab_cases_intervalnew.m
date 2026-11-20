function cases=qstab_cases_intervalnew()
%QSTAB_CASES_INTERVALNEW Degree/parity selection from the interval article.
cases={ ...
 'Interval-face output consists of reusable ordinary affine families',@publicFamilies; ...
 'Full degree/parity face counts are 256, 1536, 1536, 896 and 1024',@fullCounts; ...
 'All 256 active masks preserve pair-free coverage and remove duplicates',@allMasks; ...
 'Inactive odd directions retain required EO and OO two-faces',@degenerateQuartic; ...
 'Degree-one vertex reduction decides a complete positive interval',@linearStable; ...
 'Stable quadratic boundary cannot hide its unstable centre from two-face testing',@quadraticPocket; ...
 'Global degree loss is checked before the new face reduction',@degreeLoss; ...
 'The Hurwitz face selection is never applied to a Schur interval',@schurRejected; ...
 'Incomplete new reductions retain complete pending face descriptors',@budgets; ...
 'robuststable forwards the new opt-in selection without changing decisions',@dispatch; ...
 'Legacy pair-free and monic quadratic selections remain unchanged',@legacy; ...
 'Constant intervals retain their full family and report no enumeration certificate',@constant; ...
 'Public interval-face calls preserve their input and reject conflicting metadata',@inputRules};
end

function publicFamilies()
F=qstab.family('interval',[1.9,2.9,1],[2.1,3.1,1]);
[faces,I]=qstab.intervalfaces(F);
assert(iscell(faces)&&iscolumn(faces)&&numel(faces)==4);
assert(~I.stabilityTestPerformed&&~I.certificateComplete&&~I.validated);
original=qstab.familycomponents(F);
for j=1:numel(faces)
    assert(isstruct(faces{j})&&strcmp(faces{j}.kind,'affine'));
    G=qstab.familycomponents(faces{j});T=I.descriptors{j};
    assert(isequal(G.bounds,T.bounds)&&G.degree==2);
    theta=mean(G.bounds,2);C=G.base;
    for k=1:numel(theta),C=C+theta(k)*G.generators(:,:,k);end
    assert(all(C(:)>=original.intervalLower(:)-1e-14));
    assert(all(C(:)<=original.intervalUpper(:)+1e-14));
end
[tf,R]=qstab.zeroexclude(faces{1},'MaxTime',Inf);
assert(tf==1&&R.certificateComplete);
[raw,J]=qstab.intervalfaces({F.lower,F.upper});
assert(isequaln(faces,raw)&&isequaln(I,J));
end

function fullCounts()
expected=[256,1536,1536,896,1024];
for n=1:5
    F=fullBox(n);[~,I]=qstab.intervalfaces(F);
    assert(I.totalSubfamilies==expected(n));
    if n==4
        assert(I.faceCountsByDimension(3)==384&&I.faceCountsByDimension(4)==512);
    end
end
end

function allMasks()
% Independent check: enumerate all ternary cube faces, then require the
% mathematically selected full two/three skeleton, including degeneracies.
ternary=zeros(3^8,8);integers=(0:3^8-1)';
for j=1:8,ternary(:,j)=mod(floor(integers/3^(j-1)),3);end
pairFree=~any(ternary(:,1:4)==2 & ternary(:,5:8)==2,2);
dimension=sum(ternary==2,2);
for mask=0:255
    active=find(bitget(mask,1:8));inactive=setdiff(1:8,active);
    eligible=pairFree & all(ternary(:,inactive)==0,2);
    m=numel(unique(mod(active-1,4)));
    for n=[1,2,4,5]
        [~,tasks,I]=qstab.internal.intervalfaceplan(qstab.familycomponents(maskBox(n,mask)));
        states=taskStates(tasks);
        assert(size(unique(states,'rows'),1)==size(states,1));
        assert(all(all(states(:,inactive)==0)));
        assert(all(~any(states(:,1:4)==2 & states(:,5:8)==2,2)));
        if n==1,k=0;elseif n==2,k=min(2,m);else,k=min(3,m);end
        expected=ternary(eligible & dimension==k,:);
        if n==4 && k==3
            keep=sum(expected(:,5:8)==2,2)>=2;
            triples=expected(keep,:);
            twos=ternary(eligible & dimension==2,:);
            covered=false(size(twos,1),1);
            for j=1:size(triples,1)
                fixed=find(triples(j,:)~=2);
                covered=covered | all(bsxfun(@eq,twos(:,fixed),triples(j,fixed)),2);
            end
            expected=[triples;twos(~covered,:)];
        end
        assert(isequal(sortrows(states),sortrows(expected)));
        assert(I.totalSubfamilies==size(expected,1));
    end
end
end

function degenerateQuartic()
% E_real, O_i, E_j: no EOO/OOO triple can absorb the EO boundaries.
mask=sum(2.^([1,6,3]-1));[~,I]=qstab.intervalfaces(maskBox(4,mask));
assert(I.totalSubfamilies==6&&I.faceCountsByDimension(3)==6);
types=cellfun(@(x)x.parityType,I.descriptors,'UniformOutput',false);
assert(sum(strcmp(types,'EO'))==4&&sum(strcmp(types,'EE'))==2);
% With only two odd components, the entire two-face must remain.
mask=sum(2.^([5,6]-1));[~,I]=qstab.intervalfaces(maskBox(4,mask));
assert(I.totalSubfamilies==1&&strcmp(I.descriptors{1}.parityType,'OO'));
end

function linearStable()
F=qstab.family('interval',[1,1],[2,2]);
[tf,R]=qstab.zeroexclude(F,'Reduction','interval-faces','MaxTime',Inf);
assert(tf==1&&R.certificateComplete&&R.reduction.complete);
assert(R.reduction.vertexOnlyTest&&R.reduction.parameterDimensionBound==0);
assert(R.reduction.totalSubfamilies==4);
end

function quadraticPocket()
L=[1,1,1;2,0,0;0,-2,0;0,-2,0];U=L;U(3:4,2)=2;
F=qstab.family('interval',parts(L),parts(U));
for y=[-2,2]
    for z=[-2,2]
        C=L;C(3,2)=y;C(4,2)=z;
        assert(qstab.isstable(parts(C))==1);
    end
end
[tf,R]=qstab.zeroexclude(F,'Reduction','interval-faces','MaxTime',Inf);
assert(tf==0&&strcmp(R.status,'notStable'));
assert(R.reduction.totalSubfamilies==1&&R.reduction.parameterDimensionBound==2);
assert(R.witness.originalComponents(3,2)==0&&R.witness.originalComponents(4,2)==0);
end

function degreeLoss()
F=qstab.family('interval',[1,-1],[2,1]);
[~,I]=qstab.intervalfaces(F);assert(~I.cube.degreeInvariant);
[tf,R]=qstab.zeroexclude(F,'Reduction','interval-faces','MaxTime',Inf);
assert(isnan(tf)&&strcmp(R.status,'assumptionFailed'));
assert(R.reduction.processedSubfamilies==0&&all(R.witness.originalComponents(:,end)==0));
end

function schurRejected()
F=qstab.family('interval',[-.5,1],[.5,1],'Variable','z');
try
    qstab.intervalfaces(F);error('test:ExpectedError','Expected a CT-only error.');
catch E
    assert(strcmp(E.identifier,'qstab:intervalfaces:Variable'));
end
[tf,R]=qstab.zeroexclude(F,'Reduction','interval-faces','MaxTime',Inf);
assert(isnan(tf)&&strcmp(R.status,'undetermined')&&~R.certificateComplete);
assert(R.reduction.processedSubfamilies==0);
end

function budgets()
F=fullBox(2);
[tf,R]=qstab.zeroexclude(F,'Reduction','interval-faces','MaxTime',0);
assert(isnan(tf)&&~R.certificateComplete&&R.reduction.totalSubfamilies==1536);
assert(R.reduction.unvisitedSubfamilies==1536);
pending=R.unresolvedRegions{end};assert(numel(pending.faces)==1536);
[tf,R]=qstab.zeroexclude(F,'Reduction','interval-faces', ...
    'MaxSubfamilies',1,'MaxTime',Inf,'MaxCells',100);
assert(isnan(tf)&&~R.certificateComplete&&R.reduction.processedSubfamilies==1);
assert(R.reduction.unvisitedSubfamilies==1535);
assert(R.processedCells<=100);
end

function dispatch()
F=qstab.family('interval',[1,1],[2,2]);
[a,A]=qstab.zeroexclude(F,'Reduction','interval-faces','MaxTime',0);
[b,B]=qstab.robuststable(F,'Method','zeroexclusion', ...
    'Reduction','interval-faces','MaxTime',0);
assert(isequaln(a,b)&&strcmp(A.status,B.status));
assert(isequaln(A.reduction,B.reduction));
end

function legacy()
F=fullBox(2);
[~,A]=qstab.zeroexclude(F,'Reduction','kharitonov-faces','MaxTime',0);
[~,B]=qstab.zeroexclude(F,'Reduction','quadratic-faces','MaxTime',0);
[~,C]=qstab.zeroexclude(F,'Reduction','interval-faces','MaxTime',0);
assert(A.reduction.totalSubfamilies==1024&&B.reduction.totalSubfamilies==192);
assert(C.reduction.totalSubfamilies==1536);
end

function constant()
F=qstab.family('interval',parts([1;-.1;-.1;-.1]),parts([2;.1;.1;.1]));
[faces,I]=qstab.intervalfaces(F);
assert(numel(faces)==1&&strcmp(I.selection,'constant-family'));
assert(~I.stabilityTestPerformed&&~I.certificateComplete);
[tf,R]=qstab.zeroexclude(F,'Reduction','interval-faces','MaxTime',Inf);
assert(tf==1&&R.certificateComplete);
end

function inputRules()
F=fullBox(3);old=F;
folder=pwd;random=rng;
qstab.intervalfaces(F);
assert(isequaln(F,old)&&strcmp(folder,pwd)&&isequaln(random,rng));
try
    qstab.intervalfaces(F,'Variable','z');error('test:ExpectedError','Expected a metadata conflict.');
catch E
    assert(strcmp(E.identifier,'qstab:family:OptionConflict'));
end
try
    qstab.intervalfaces(qstab.family('polytope',{[1,1],[2,1]}));
    error('test:ExpectedError','Expected a family kind conflict.');
catch E
    assert(strcmp(E.identifier,'qstab:family:OptionConflict'));
end
end

function F=fullBox(n)
F=maskBox(n,255);
end
function F=maskBox(n,mask)
C=zeros(4,n+1);C(1,:)=arrayfun(@(k)nchoosek(n,k),0:n);L=C;U=C;
for j=1:8
    if bitget(mask,j)
        nu=mod(j-1,4)+1;k=1+(j>4);
        L(nu,k)=L(nu,k)-.001;U(nu,k)=U(nu,k)+.001;
    end
end
F=qstab.family('interval',parts(L),parts(U));
end
function P=parts(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
function states=taskStates(tasks)
states=zeros(numel(tasks),8);
for j=1:numel(tasks)
    B=tasks{j}.bounds;row=B(:,1)';row(B(:,1)~=B(:,2))=2;states(j,:)=row;
end
end
