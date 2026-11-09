function cases=test_qstab_exploration_public(group)
%TEST_QSTAB_EXPLORATION_PUBLIC Sampled-geometry regression cases for test_qstab.
if nargin==0,group='numeric';end
switch group
case 'numeric'
cases={
 'Interval rectangle of complex algebraic zero representatives',@zeroRectangle;
 'Shared weights on a polynomial polytope',@polytope;
 'One shared affine coordinate at every degree',@affine;
 'Noncommuting quaternion side survives sampling',@quaternionSide;
 'Finite sampling window versus unbounded parameter law',@windows;
 'Physical value, two projections and central norm are distinct',@orientations;
 'Frequency sphere retains direction and member indices',@sphericalPlot;
 'Forward/delay variables and opposite frequency signs',@discrete;
 'Zero polynomials and infinity at preserved grade',@singularCases;
 'Deterministic sampling respects total member budget',@budget;
 'Repeatable seeded draws preserve caller RNG',@sharedRandom;
 'Malformed canonical families, grades and projection options are rejected',@rejectBad;
 'Tagged families carry variable, grade, parameter domain and sampling period',@taggedModels;
 'Raw interval at declared grade reproduces four checks on one cubic',@cubicInterval};
case 'graphics'
cases={'Existing axes, hold state and graphics-free mode',@parentAxes};
case 'native'
cases={'Native quaternion family and coordinate/axis agreement',@nativeInputs};
otherwise
error('qstab:tests:UnknownGroup','Unknown exploration test group.');
end
end
function q=coeff(a,b,c,d)
q=struct('a',a,'b',b,'c',c,'d',d);
end

function equalNear(a,b,tol)
if nargin<3,tol=8e-10;end
assert(isequal(size(a),size(b)) && all(isfinite(a(:))) && all(isfinite(b(:))) && ...
    all(abs(a(:)-b(:))<=tol*(1+abs(b(:)))),'Values or shapes differ.');
end

function zeroRectangle()
% P(s)=s+a+b*j. Its single algebraic class has representative -a+|b|i.
L=coeff([-1 1],[0 0],[0 0],[0 0]); U=L; U.a(1)=1; U.c(1)=2;
F=qstab.familycomponents({L,U},'Uncertainty','interval');
[ax,d]=qstab.plotFamilyZeros(F,9,'Plot',false);
assert(isempty(ax)&&all(d.valid)&&isequal(size(d.lambda),[1 9]));
equalNear(d.lambda(:),-d.sample.parameters(:,1)+1i*d.sample.parameters(:,2));
assert(size(unique(d.sample.parameters,'rows'),1)==9);
end
function polytope()
F={[1 1],[3 2],[4 4]};
[~,d]=qstab.plotFamilyZeros(F,10,'Plot',false);
W=d.sample.weights; equalNear(sum(W,2),ones(10,1));assert(all(W(:)>=0));
assert(size(unique(W,'rows'),1)==10 && all(ismember(eye(3),W,'rows')));
equalNear(d.lambda(:),-(W*[1;3;4])./(W*[1;2;4]));
end
function affine()
F=struct('center',[1 2],'generators',{{[3 -1]}},'bounds',[-.5 1]);
[~,d]=qstab.plotFamilyValues(F,2,7,'Plot',false);t=d.sample.parameters;
equalNear(d.valueComponents,[1+3*t  2*(2-t) zeros(7,2)]);
end
function quaternionSide()
% P1=k; j*k=i while k*j=-i.
p1=coeff(0,0,0,1);B=[0 0;0 0;1 1;0 0];
for s={'left','right'}
    F=struct('p0',[2 1],'p1',p1,'bounds',B,'side',s{1});
    [~,d]=qstab.plotFamilyValues(F,0,3,'Plot',false);
    sg=2*strcmp(s{1},'left')-1;
    equalNear(d.valueComponents,repmat([2 sg 0 0],3,1));
end
end
function windows()
% s+2+alpha has its first real contact at alpha=-2 irrespective of this B.
F=struct('p0',[2 1],'p1',1,'bounds',[-.5 .25;0 0;0 0;0 0],'side','left');
ab=qstab.raymargin(F.p0,F.p1,1,'Side',F.side);
assert(abs(ab(1)+2)<1e-8 && isinf(ab(2)) && ab(2)>0);
Fu=struct('p0',[2 1],'p1',1,'bounds',[],'side','left');
expectError(@()qstab.plotFamilyZeros(Fu,3,'Plot',false));
narrow=[-.25 .125;0 0;0 0;0 0];
[~,d]=qstab.plotFamilyZeros(Fu,3,'ParameterBox',narrow,'Plot',false);
equalNear(d.lambda(:),-2-d.sample.parameters(:,1));assert(isempty(Fu.bounds));
expectError(@()qstab.plotFamilyZeros(F,3,'ParameterBox',[-.6 .25;0 0;0 0;0 0],'Plot',false));
end
function orientations()
% P=s-j. At i its norm polynomial vanishes, but its right evaluation is i-j.
p=coeff([0 1],[0 0],[-1 0],[0 0]); F={p};
[~,a]=qstab.plotFamilyValues(F,1,1,'Plot',false);
[~,b]=qstab.plotFamilyValues(F,1,1,'AxisUnit','j','Plot',false);
[~,n]=qstab.plotFamilyValues(F,1,1,'Projection','normpoly','Plot',false);
[~,s]=qstab.plotFamilyValues(F,0,1,'Projection','slice','Plot',false);
equalNear(a.valueComponents,[0 1 -1 0]); equalNear(a.projected,[0 sqrt(2)]);
equalNear(b.valueComponents,zeros(1,4)); equalNear(n.normValues,0);
equalNear(s.projected,[0 0]); equalNear(s.offSlice,1);
% Left coefficient j times right argument i = -k. Reversing order gives +k.
G={coeff([0 0],[0 0],[0 1],[0 0])};
[~,g]=qstab.plotFamilyValues(G,1,1,'Plot',false);equalNear(g.valueComponents,[0 0 0 -1]);
% A genuinely complex norm value must not be replaced by |P|^2.
H={[1 1]};
[~,h]=qstab.plotFamilyValues(H,1,1,'Projection','normpoly','Plot',false);
equalNear(h.normValues,2i);equalNear(h.projected,[0 2]);
end
function sphericalPlot()
p=coeff([0 1],[0 0],[-1 0],[0 0]); F={p,p};
[~,d]=qstab.plotFamilyValues(F,1,3,'AxisUnit','sphere','DirectionCount',8,'Plot',false);
assert(isequal(size(d.valueComponents),[24,4])&&d.sample.count==3&&size(d.axisUnit,1)==8);
ix=d.directionIndex(:);u=d.axisUnit(ix,:);
equalNear(d.valueComponents,[zeros(24,1),u(:,1),u(:,2)-1,u(:,3)]);
assert(any(d.valueNorm(:)==0)&&all(accumarray(d.memberIndex(:),1,[3,1])==8));
end
function discrete()
Fz={[-.5 1]};
Fd={[1 -.5]};
[~,zf]=qstab.plotFamilyZeros(Fz,1,'Variable','z','Plot',false);
[~,df]=qstab.plotFamilyZeros(Fd,1,'Variable','zinv','Plot',false);equalNear(zf.lambda,.5);equalNear(df.lambda,2);
[~,a]=qstab.plotFamilyValues(Fz,pi/2,1,'Variable','z','Plot',false);[~,b]=qstab.plotFamilyValues(Fd,pi/2,1,'Variable','zinv','Plot',false);
equalNear(a.argumentComponents,[0 1 0 0]);equalNear(b.argumentComponents,[0 -1 0 0]);
equalNear(a.valueComponents,[-.5 1 0 0]);equalNear(b.valueComponents,[1 .5 0 0]);
end
function singularCases()
F={[1 1],1};
[~,d]=qstab.plotFamilyZeros(F,2,'Plot',false);
assert(all(d.valid)&&nnz(isinf(d.lambda))==1&&sum(d.infiniteCount)==1);
G={0,1};
[~,z]=qstab.plotFamilyZeros(G,2,'Plot',false);
assert(nnz(z.zeroPolynomial)==1&&nnz(z.valid)==1&&size(z.lambda,1)==0);
end
function budget()
M=22; G=repmat({[1]},1,M); F=struct('center',[2 1],'generators',{G},'bounds',repmat([-1 1],M,1));
[~,d]=qstab.plotFamilyZeros(F,11,'Sampling','grid','Plot',false);
T=d.sample.parameters;assert(isequal(size(T),[11 M])&&all(abs(T(:))<=1));
equalNear(d.lambda(:),-2-sum(T,2));
% Singleton domain: N counts evaluations, even if all describe one point.
S=qstab.familycomponents({[1 1],[1 1]},'Uncertainty','interval');
[~,s]=qstab.plotFamilyZeros(S,5,'Plot',false);equalNear(s.lambda,-ones(1,5));
end
function sharedRandom()
F={[1 1],[2 1],[3 1]};
before=rng;
[~,z]=qstab.plotFamilyZeros(F,13,'Sampling','random','Seed',902,'Plot',false);
[~,v]=qstab.plotFamilyValues(F,.4,13,'Sampling','random','Seed',902,'Plot',false);
after=rng;assert(isequaln(before,after)&&isequaln(z.coefficients,v.coefficients));
assert(all(z.sample.weights(:)>=0));equalNear(sum(z.sample.weights,2),ones(13,1));
assert(size(unique(z.sample.weights,'rows'),1)==13 && ~v.isCertificate);
end
function parentAxes()
F={[1 1],[2 1]};
oldFigure=get(0,'CurrentFigure');fig=figure('Visible','off');
cleanup=onCleanup(@()restoreFigure(fig,oldFigure)); %#ok<NASGU>
p=axes('Parent',fig);original=plot(p,1,2);hold(p,'off');
[ax,d]=qstab.plotFamilyZeros(F,4,'Parent',p);
assert(isequal(ax,p)&&ishandle(original)&&~ishold(p)&&all(d.valid));
hold(p,'on');[ax,d]=qstab.plotFamilyValues(F,1,4,'Parent',p,'Projection','slice');
assert(isequal(ax,p)&&ishandle(original)&&ishold(p)&&size(d.projected,1)==4);
n=numel(findall(0,'Type','figure'));
[ax,~]=qstab.plotFamilyValues(F,1,2,'Plot',false);
assert(isempty(ax)&&numel(findall(0,'Type','figure'))==n);
end
function restoreFigure(fig,old)
if ishandle(fig),close(fig);end
if ~isempty(old)&&ishandle(old),set(0,'CurrentFigure',old);end
end
function expectError(f)
try, f(); catch ME
    assert(strncmp(ME.identifier,'qstab:',6),'Unexpected diagnostic: %s',ME.identifier); return
end
error('testqstabexplore:MissingError','An invalid request was accepted.');
end

function rejectBad()
F=qstab.familycomponents({[1 1],[2 1]},'Uncertainty','interval');F.base(1,1)=77;
expectError(@()qstab.plotFamilyValues(F,1,'Plot',false));
expectError(@()qstab.plotFamilyValues({[1 2]},1,'Grade',0,'Plot',false));
expectError(@()qstab.plotFamilyValues({[1 2]},1,'AxisUnit','sphere','Projection','slice','Plot',false));
expectError(@()qstab.plotFamilyZeros({1},3,'AxisUnit','j','Plot',false));
end
function cubicInterval()
lo=coeff([-1 0 1 1],[-2 -1 0 0],[-1 1 2 0],[-1 0 0 0]);
hi=coeff([2 3 4 1],[1 2 3 0],[1 2 3 0],[1 1 1 0]);
raw={lo,hi}; args={'Uncertainty','interval','Plot',false};
[~,v]=qstab.plotFamilyValues(raw,1,9,args{:}); [~,s]=qstab.plotFamilyValues(raw,1,9,args{:},'Projection','slice');
[~,n]=qstab.plotFamilyValues(raw,1,9,args{:},'Projection','normpoly');
assert(isequal(size(v.sample.parameters),[9 12]) && isequal(size(v.sample.weights),[9 0]));
row=5; C=v.coefficients(:,:,row); E=C(:,1)-C(:,3);O=C(:,2)-C(:,4);w=E+[-O(2);O(1);O(4);-O(3)];
equalNear(v.valueComponents(row,:).',w);equalNear(v.projected(row,:),[w(1) norm(w(2:4))]);
equalNear([s.projected(row,:) s.offSlice(row)],[w(1) w(2) norm(w(3:4))]);equalNear(n.normValues(row),sum((E+1i*O).^2));
% A constant at declared grade 2 has two infinite representatives.
[~,r]=qstab.plotFamilyZeros({1},1,'Grade',2,'Variable','zinv','Plot',false);
assert(r.valid && r.infiniteCount==2 && isequal(size(r.lambda),[2 1]));
end
function nativeInputs()
p=quaternion([0 1],[0 0],[-1 0],[0 0]); unit=quaternion(0,0,1,0);
[~,a]=qstab.plotFamilyValues({p},1,1,'AxisUnit',unit,'Plot',false);
[~,z]=qstab.plotFamilyZeros({p},1,'Plot',false);
equalNear(a.valueComponents,zeros(1,4));equalNear(z.lambda,1i);assert(z.valid);
end

function taggedModels()
F=qstab.family('interval',[.5 1],[.75 1],'Variable','zinv','SamplingPeriod',.125,'Grade',3);
[~,D]=qstab.plotFamilyValues(F,pi/2,3,'Plot',false);
assert(strcmp(D.sample.variable,'zinv') && D.sample.samplingPeriod==.125 && D.sample.grade==3);
equalNear(D.valueComponents,[linspace(.5,.75,3).' -ones(3,1) zeros(3,2)]);
C=qstab.familycomponents(F);[~,J]=qstab.plotFamilyValues(C,pi/2,3,'Plot',false);
equalNear(J.valueComponents,D.valueComponents);
expectError(@()qstab.plotFamilyValues(F,1,2,'Variable','z','Plot',false));
expectError(@()qstab.plotFamilyZeros(F,2,'Uncertainty','polytope','Plot',false));
Fu=qstab.family('oneparameter',[2 1],1,'Side','left');
box=[-1 1;0 0;0 0;0 0];[~,z]=qstab.plotFamilyZeros(Fu,3,'ParameterBox',box,'Plot',false);
assert(isempty(z.sample.originalBounds) && isequal(z.sample.bounds,box));
equalNear(z.lambda,[-1 -2 -3]);
end
