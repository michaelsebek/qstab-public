function cases=test_qstab_families_public()
%TEST_QSTAB_FAMILIES_PUBLIC Tests for definitions, generators and bounded rays.
cases={
 'Interval parameters and correlated convex weights describe different sets', @independence;
 'Seeded monic integer and delay families preserve draws, order and caller RNG', @generator;
 'One-parameter models preserve left and right Hamilton order', @sides;
 'Fixed-ray first exit can precede re-entry and is independent of parameter bounds', @cubic;
 'Closed parameter windows can miss the ray, cross a contact, or be singletons', @windows;
 'Unbounded laws are accepted only by appropriate searches; boxed radii need opt-in', @unbounded;
 'Robust calls inherit forward/delay variable and retain declared grade', @modelVariable;
 'Invalid bounds, conflicting metadata and mismatched saved reports are rejected', @invalid;
 'Constants have no degree loss and default calls generate a usable finite family', @constants;
};
end
function independence()
L=struct('a',[1 1],'b',[-1 0],'c',[0 0],'d',[0 0]);U=L;U.a(1)=3;U.b(1)=1;
F=qstab.family('interval',L,U);H=qstab.familycomponents(F);
assert(strcmp(H.kind,'interval')&&isequal(H.parameterIndex,[1 0;2 0])&&isequal(H.bounds,[1 3;-1 1]));
G=qstab.family('polytope',{L,U});p=qstab.valueset(G,0,'Weights',[.5 .5]);
assert(abs(p-4)<1e-12);
[G0,I]=qstab.valueenclosure(G,0);assert(~isempty(G0)&&strcmp(I.uncertainty,'polytope'));
% Two constant generators stay two generators, not one polynomial.
V=qstab.family('polytope',{2,3}); assert(isequal(qstab.valueset(V,0),[4;9]));
A=qstab.family('affine',1,{[1 1]},[-2 3]);E=qstab.familycomponents(A);assert(E.degree==1&&isequal(E.generators(1,:,1),[1 1])&&isequal(E.bounds,[-2 3]));
end
function generator()
before=rng; F=qstab.randfamily(3,'interval','Seed',17);F2=qstab.randfamily('interval',3,'Seed',17);after=rng;
assert(isequaln(before,after)&&isequaln(F,F2));
H=qstab.familycomponents(F);assert(size(H.bounds,1)==12&&H.degree==3&&isequal(H.intervalLower(:,4),[1;0;0;0]));
V=qstab.randfamily(2,'polytope',4,'z','int','Seed',2);
D=qstab.randfamily(2,'polytope',4,'zinv','int','SamplingPeriod',.1,'Seed',2);
A=qstab.familycomponents(V);B=qstab.familycomponents(D);
assert(isequal(B.vertices,A.vertices(:,end:-1:1,:))&&isequal(B.vertices,fix(B.vertices))&&B.samplingPeriod==.1);
assert(qstab.isstable(V.vertices{1},'Variable','z')==1);
end
function sides()
L=qstab.randfamily(2,'oneparameter','Side','left','Seed',17);
R=qstab.randfamily(2,'oneparameter','Side','right','Seed',17);
A=qstab.familycomponents(L);B=qstab.familycomponents(R);
assert(isequal(L.p0,R.p0)&&isequal(L.p1,R.p1)&&~isequal(A.generators,B.generators));
% Exact orientation: t=j, P1=i -> j*i=-k and i*j=k.
L=qstab.family('oneparameter',[2 1],1i,'Bounds',repmat([-1 1],4,1)); R=L;R.side='right';
A=qstab.familycomponents(L);B=qstab.familycomponents(R);
assert(isequal(A.generators(:,1,3),[0;0;0;-1])&&isequal(B.generators(:,1,3),[0;0;0;1]));
end
function cubic()
p=[.5 1 1 1];d=[4 1 1];F=qstab.family('oneparameter',p,d,'Bounds',[-.1 .2;zeros(3,2)]);
[ab,info]=qstab.raymargin(F,1);assert(norm(ab-[-1/8,1-1/sqrt(2)],inf)<1e-8);
assert(isequal(info.parameterDomain.admissibleAlpha,[-.1 .2])&&info.parameterDomain.nominalAdmissible);
a=[0 .5 2];tf=arrayfun(@(x)qstab.isstable(p+x*[d 0],'Method','spectral'),a);
assert(isequal(tf,[1 0 1]));
end
function windows()
F=qstab.family('oneparameter',[2 1],1,'Bounds',[-1 .5;zeros(3,2)]);
[ab,I]=qstab.raymargin(F,-2);assert(ab(1)==-Inf&&abs(ab(2)-1)<1e-8&&isequal(I.parameterDomain.admissibleAlpha,[-.25,.5]));
F.bounds(1,:)=[-3 -1];[a,I]=qstab.raymargin(F,-2); assert(isequaln(a,ab)&&isequal(I.parameterDomain.admissibleAlpha,[.5 1.5])&&~I.parameterDomain.nominalAdmissible);
F.bounds(2,:)=[.5 1];[~,I]=qstab.raymargin(F,-2);assert(isempty(I.parameterDomain.admissibleAlpha)&&strcmp(I.parameterDomain.alphaRangeStatus,'empty'));
F.bounds=[-1 -1;zeros(3,2)];[~,I]=qstab.raymargin(F,-2);assert(isequal(I.parameterDomain.admissibleAlpha,[.5 .5]));
end
function unbounded()
F=qstab.family('oneparameter',[2 1],1);[~,I]=qstab.raymargin(F,-2);
assert(isequal(I.parameterDomain.admissibleAlpha,[-Inf Inf])&&~I.parameterDomain.isBounded);
mustFail(@()qstab.robuststable(F),'qstab:family:Unbounded');
[r,I]=qstab.freqmargin(F,0);assert(abs(r-2)<1e-10&&~I.parameterDomain.explicitExtension);
F.bounds=repmat([-.1 .1],4,1);mustFail(@()qstab.freqmargin(F,0),'qstab:family:BoundedMargin');
mustFail(@()qstab.stabilityradius(F),'qstab:family:BoundedMargin');
[r,I]=qstab.freqmargin(F,0,'ExtendParameter',true);assert(abs(r-2)<1e-10&&I.parameterDomain.explicitExtension);
[b,I]=qstab.stabilityradius(F,'ExtendParameter',true,'MaxCells',2);
assert(b(1)<=2&&b(2)>=2&&I.parameterDomain.isBounded&&I.parameterDomain.explicitExtension);
end
function modelVariable()
% Delay input 1-0.5d has z=0.5; padding to declared grade 2 adds a zero mode.
F=qstab.family('polytope',{[1 -.5]},'Variable','zinv','Grade',2,'SamplingPeriod',.25);
[tf,I]=qstab.robuststable(F);[tz,Z]=qstab.zeroexclude(F);
assert(tf==1&&tz==1&&strcmp(I.variable,'zinv')&&strcmp(Z.variable,'zinv')&&I.degree==2);
[lo,J]=qstab.leadingbound(F);assert(lo>0&&strcmp(J.variable,'zinv'));
F=qstab.family('oneparameter',[1 -.5],[0 1],'Variable','zinv','Grade',2);
[a,R]=qstab.raymargin(F,1);b=qstab.raymargin([1 -.5],[0 1],1,'Variable','zinv','Grade',2);assert(isequaln(a,b)&&R.parameterDomain.grade==2);
end
function invalid()
mustFail(@()qstab.family('interval',[2 1],[1 1]),'qstab:family:IntervalOrder');
mustFail(@()qstab.family('affine',1,{1},[3 1]),'qstab:family:Bounds');
mustFail(@()qstab.family('oneparameter',1,1,'Bounds',[],'Bounds',[]),'qstab:family:Option');
F=qstab.family('oneparameter',[2 1],1,'Side','right');mustFail(@()qstab.raymargin(F,1,'Side','left'),'qstab:family:OptionConflict');
F=qstab.family('interval',[.5 1],[.6 1],'Variable','z');
mustFail(@()qstab.robuststable(F,'Variable','s'),'qstab:family:OptionConflict');
mustFail(@()qstab.robuststable(F,'Uncertainty','polytope'),'qstab:family:OptionConflict');
mustFail(@()qstab.raymargin(F,1),'qstab:family:MarginKind');
mustFail(@()qstab.raymargin(F,1,'Grade',2),'qstab:family:MarginKind');
mustFail(@()qstab.kharitonovcube(F),'qstab:kharitonovcube:Variable');
[~,R]=qstab.robuststable(F); F.variable='s';F.samplingPeriod=0;
mustFail(@()qstab.witnessData(F,R),'qstab:family:ReportConflict');
mustFail(@()qstab.randfamily(3,'oneparameter','Bounds',NaN(4,2)),'qstab:randfamily:Bounds');
mustFail(@()qstab.randfamily(2,'interval','Side','left'),'qstab:randfamily:Side');
end
function constants()
% Restore the RNG for this test of the intentionally unseeded, interactive call.
before=rng;clean=onCleanup(@()rng(before)); %#ok<NASGU>
F=qstab.randfamily();assert(ismember(F.grade,1:3)&&strcmp(F.kind,'interval')&&~isequaln(before,rng));
for kind={'interval','affine','polytope'}
 F=qstab.randfamily(0,kind{1},'int','Seed',13);[tf,S]=qstab.robuststable(F);assert(tf==1&&S.degree==0);
end
mustFail(@()qstab.randfamily(0,'oneparameter'),'qstab:randfamily:Degree');
end
function mustFail(fn,id)
try,fn();catch E,assert(strcmp(E.identifier,id),'Expected %s; got %s: %s',id,E.identifier,E.message);return;end
error('qstab:test:ExpectedError','Expected %s.',id);
end
