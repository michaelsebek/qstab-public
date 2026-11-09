function R=qstab_guide_example(id)
%QSTAB_GUIDE_EXAMPLE Run the small, checked examples printed in the public guide.
%   qstab_guide_example() or qstab_guide_example('list') returns the catalogue.
%   R=qstab_guide_example('G06') runs one example. R.ok means the example's
%   assertions passed; a robust stability verdict, if present, is separate.
%   A=qstab_guide_example('all') returns a cell array of completed results or
%   raises the first failed assertion. No figure, file, path change, download or
%   caller variable is created. A fixed local seed is used in the sampling lab.
%   Each case below is the exact source of its corresponding guide listing.
names={'G01','Algebraic class versus point evaluation';'G02','Forward and delay variables';
 'G03','Right spectra of a nonnormal matrix';'G04','First exit and re-entry';
 'G05','Independent components and a shared affine parameter';
 'G06','Two local hulls exclude the origin';'G07','Six-generator cubic';
 'G08','A certificate and a separate admissible witness';'G09','Two perturbation geometries';
 'G10','Explore a random interval';'G11','A physical window and the all-direction radius'};
if nargin==0, R=names; return; end
if ~(ischar(id)&&isrow(id)) && ~(isstring(id)&&isscalar(id))
    error('qstab:guide:ID','Supply one example ID, list, or all.'); end
id=upper(strtrim(char(id)));
if strcmp(id,'LIST'), R=names; return; end
if strcmp(id,'ALL')
    R=cell(size(names,1),1);
    for n=1:numel(R), R{n}=qstab_guide_example(names{n,1}); end
    return
end
R=struct('id',id,'ok',false);
switch id
case 'G01'
%<G01>
P = struct('a',[0 1],'b',[0 0],'c',[-1 0],'d',[0 0]); % s-j
[N,~] = qstab.normpoly(P);
[lambda,S] = qstab.roots(P);
assert(S.isValid && abs(lambda-1i)<1e-8)
assert(isequal(N(:).',[1 0 1]))
%</G01>
R.normCoefficients=N; R.lambda=lambda;
case 'G02'
%<G02>
zP = [-0.5 1];   % z-1/2
dP = [1 -0.5];   % 1-d/2, d=z^-1
tfz = qstab.isstable(zP,'Variable','z');
tfd = qstab.isstable(dP,'Variable','zinv','Grade',2);
[delayRoots,I] = qstab.roots(dP,'Variable','zinv','Grade',2);
assert(tfz==1 && tfd==1 && I.isValid)
assert(any(isinf(delayRoots)) && any(abs(delayRoots-2)<1e-8))
%</G02>
R.decisions=[tfz tfd]; R.delayRoots=delayRoots;
case 'G03'
%<G03>
A = struct('a',[-1 4;0 -2],'b',[1 0;0 0], ...
           'c',[0 1;0 1],'d',zeros(2));
[tf,I] = qstab.isstable(A,'Type','matrix');
[V,D,E] = qstab.eigpairs(A,'Output','components');
res = norm(qstab.chi(A)*qstab.chi(V)-qstab.chi(V)*qstab.chi(D),'fro');
assert(tf==1 && E.eigenpairsValid && res<1e-8)
%</G03>
R.decision=tf; R.status=I.status; R.residual=res;
case 'G04'
%<G04>
P0 = [.5 1 1 1]; P1 = [4 1 1];
[ab,I] = qstab.raymargin(P0,P1,1,'Side','left');
alpha = [0 .5 2]; tf = zeros(size(alpha));
for n=1:numel(alpha)
    tf(n) = qstab.isstable(P0+alpha(n)*[P1 0],'Method','spectral');
end
assert(norm(ab-[-1/8,1-1/sqrt(2)],Inf)<1e-8)
assert(isequal(tf,[1 0 1]))
%</G04>
R.interval=ab; R.decisions=tf; R.status=I.status;
case 'G05'
%<G05>
L = [1 2 1]; U = [2 3 1];
F = qstab.family('interval',L,U);
line = qstab.family('affine',L,{U-L},[0 1]);
tfBox = qstab.robuststable(F);
tfLine = qstab.robuststable(line);
[K,~] = qstab.kharitonovcube({L,U});
tfCube = qstab.zeroexclude(K,'Uncertainty','affine');
assert(isequal([tfBox tfLine tfCube],[1 1 1]))
%</G05>
R.decisions=[tfBox tfLine tfCube];
case 'G06'
%<G06>
Pm = struct('a',[1 1],'b',[0 0],'c',[-2 0],'d',[0 0]);
Pp = Pm; Pp.c(1)=2;
vertices = {Pm,Pp}; beta = linspace(0,1,41).';
[points,S] = qstab.valueset(vertices,0,'Weights',[1-beta beta]);
[coarse,C] = qstab.valueenclosure(vertices,0);
left=[1;0]; mid=[.5;.5]; right=[0;1];
[local,E] = qstab.valueenclosure(vertices,0, ...
                'VertexCells',{[left mid],[mid right]});
ranges = [coarse{1}.valueBounds(1,:); local{1}.valueBounds(1,:);
          local{2}.valueBounds(1,:)];
assert(isequal(ranges,[-3 5;1 5;1 5]))
assert(~C.zeroExcludedNumerically && E.zeroExcludedNumerically)
%</G06>
R.ranges=ranges; R.points=points; R.samples=S; R.enclosure=E;
case 'G07'
%<G07>
vertices = cubicGenerators();  % six ordered quaternion polynomials
W = [eye(6); .25 .25 .25 .25 0 0];
[points,S] = qstab.valueset(vertices,0,'Weights',W);
[~,E] = qstab.valueenclosure(vertices,0);
C=zeros(4,6);
for n=1:6
    B=qstab.polycomponents(vertices{n}); C(:,n)=B(:,1);
end
G=C.'*C; hullZero=(5*G(1,1)+31*G(1,2))/36;
assert(abs(hullZero)<1e-9 && ~E.zeroExcludedNumerically)
assert(norm(points(:)-[620;620;620;620;628.2;620;80],Inf)<1e-8)
%</G07>
R.points=points; R.samples=S; R.hullZero=hullZero; R.enclosure=E;
case 'G08'
%<G08>
good = {[1 1],[2 1]}; bad = {[-1 1],[-2 1]};
[tfGood,I] = qstab.zeroexclude(good,'StoreCertificates',true);
[rows,C] = qstab.certificateData(I);
[tfBad,J] = qstab.robuststable(bad);
[P,W] = qstab.witnessData(bad,J);
assert(tfGood==1 && ~isempty(rows))
assert(tfBad==0 && W.isAdmissible && ~isempty(P))
%</G08>
R.decisions=[tfGood tfBad]; R.certificate=C; R.witness=W;
case 'G09'
%<G09>
P = [2 1]; gamma=[1 0];
W = zeros(4,2); W(1,1)=1;    % only the real constant may change
[ball,B] = qstab.ballvalue(P,gamma,1i);
[box,C] = qstab.boxvalue(P,W,0);
assert(B.valid && C.valid && abs(ball-sqrt(5))<1e-8 && abs(box-2)<1e-8)
%</G09>
R.ball=ball; R.box=box;
case 'G10'
%<G10>
F = qstab.randfamily(3,'interval','Seed',17);
[~,Z] = qstab.plotFamilyZeros(F,9,'Plot',false);
[~,V] = qstab.plotFamilyValues(F,1,9,'Projection','class','Plot',false);
[~,G] = qstab.plotFamilyValues(F,1,9,'Projection','normpoly','Plot',false);
dims = [size(V.sample.parameters); size(V.sample.weights)];
assert(isequal(dims,[9 12;9 0]) && isequal(size(Z.lambda),[3 9]))
assert(isequal(V.coefficients,G.coefficients))
%</G10>
R.dimensions=dims; R.family=F; R.values=V; R.zeros=Z;
case 'G11'
%<G11>
B = [-3 -1;0 0;0 0;0 0]; % real part in [-3,-1], imaginary parts fixed
F = qstab.family('oneparameter',[2 1],1,'Bounds',B,'Side','left');
[ab,I] = qstab.raymargin(F,-2);       % t=-2*alpha
admissible = I.parameterDomain.admissibleAlpha;
rho = qstab.freqmargin([2 1],1,[0 1]);
[bracket,J] = qstab.stabilityradius([2 1],1);
assert(ab(1)==-Inf && abs(ab(2)-1)<1e-8 && isequal(admissible,[.5 1.5]))
assert(norm(rho-[2 sqrt(5)],Inf)<1e-8)
assert(bracket(1)<=2+1e-9 && bracket(2)>=2-1e-9)
%</G11>
R.interval=ab; R.window=admissible; R.frequencyRadius=rho; R.bracket=bracket; R.status=J.status;
otherwise
    error('qstab:guide:ID','Unknown example ID. Use qstab_guide_example(''list'').');
end
R.ok=true;
end

function V=cubicGenerators()
% Coefficient columns are real,i,j,k. Construct A*C*B in this exact order.
one=[1;0;0;0]; A=[[1;1;0;0],one]; B=[[3;0;0;1],one];
U=[1 1 1;1 -1 -1;-1 1 -1;-1 -1 1]; X=cell(1,6);
for h=1:4
    C=[[2;3*U(h,:).'],one]; X{h}=polyprod(polyprod(A,C),B);
end
AB=polyprod(A,B); X{5}=X{1}; X{5}(:,1:3)=X{5}(:,1:3)+AB/10;
X{6}=X{1}; X{6}(1,3)=X{6}(1,3)+.1;
V=cell(1,6);
for h=1:6
    C=X{h}; V{h}=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
end
function P=polyprod(A,B)
% Ordered scalar Hamilton products; powers of the indeterminate commute.
P=zeros(4,size(A,2)+size(B,2)-1);
for i=1:size(A,2)
    for j=1:size(B,2)
        a=A(:,i); b=B(:,j);
        P(:,i+j-1)=P(:,i+j-1)+[a(1)*b(1)-a(2:4).'*b(2:4);
            a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
    end
end
end
