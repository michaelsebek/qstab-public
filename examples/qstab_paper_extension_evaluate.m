function [D,K]=qstab_paper_extension_evaluate(C,o)
%QSTAB_PAPER_EXTENSION_EVALUATE Additional finite identities and witnesses.
% Uses ordinary stored data; explicit spectral checks are labelled as such.
D=struct();K={};d=C.data;ref=C.reference;qi=[0;1;0;0];qj=[0;0;1;0];qk=[0;0;0;1];
switch C.id
 case 'qb01_nearest_sphere'
  n=size(d.directions,2);D.parameters=zeros(4,n);D.points=zeros(4,n);D.residuals=zeros(1,n);D.normAtClass=zeros(1,n);
  for j=1:n
   u=[0;d.directions(:,j)];t=[-.5;zeros(3,1)]-u/sqrt(2);lambda=u/sqrt(2);
   P=d.P;P(:,1)=P(:,1)+t;D.parameters(:,j)=t;D.points(:,j)=lambda;
   D.residuals(j)=norm(evaluate(P,lambda));N=qstab.normpoly(poly(P,o));D.normAtClass(j)=abs(polyval(fliplr(N),1i/sqrt(2)));
  end
  K={ck('numeric',sqrt(sum(D.parameters.^2,1)),ref.radius*ones(1,n)), ...
     ck('numeric',D.residuals,zeros(1,n)),ck('numeric',D.normAtClass,zeros(1,n))};
 case 'qb01_degree_lift'
  D.P0=lift(d.P0,ref.factor);D.P1=lift(d.P1,ref.factor);D.omega=d.omega;
  D.norm=qstab.normpoly(poly(D.P0,o));D.referenceNorm=conv(qstab.normpoly(poly(d.P0,o)),conv(ref.factor,ref.factor));
  K={ck('numeric',D.norm,D.referenceNorm)};
  for side={'left','right'}
   s=side{1};D.base.(s)=qstab.freqmargin(poly(d.P0,o),poly(d.P1,o),d.omega,'Side',s);
   D.lifted.(s)=qstab.freqmargin(poly(D.P0,o),poly(D.P1,o),d.omega,'Side',s);
   K{end+1}=ck('numeric',D.lifted.(s),D.base.(s));
  end
 case 'qb01_conjugate_sides'
  A=d.P0;A(2:4,:)=-A(2:4,:);B=d.P1;B(2:4,:)=-B(2:4,:);D.omega=d.omega;
  D.originalLeft=qstab.freqmargin(poly(d.P0,o),poly(d.P1,o),d.omega,'Side','left');
  D.originalRight=qstab.freqmargin(poly(d.P0,o),poly(d.P1,o),d.omega,'Side','right');
  D.conjugateLeft=qstab.freqmargin(poly(A,o),poly(B,o),d.omega,'Side','left');
  D.conjugateRight=qstab.freqmargin(poly(A,o),poly(B,o),d.omega,'Side','right');
  K={ck('numeric',D.originalLeft,D.conjugateRight),ck('numeric',D.originalRight,D.conjugateLeft)};
 case 'qb02_damping'
  D.alpha=d.alpha;D.members=cell(1,numel(d.alpha));D.decisions=cell(size(D.members));
  for j=1:numel(d.alpha)
   P=d.base;P(1,2)=d.alpha(j);D.members{j}=P;D.decisions{j}=spectral(P,o);
   K{end+1}=ck('status',D.decisions{j}.status,ref.statuses{j});
  end
  P=d.base;P(1,2)=.5;D.boundaryHalf=qstab.normpoly(poly(P,o));P(1,2)=1;D.boundaryOne=qstab.normpoly(poly(P,o));
  K{end+1}=ck('numeric',D.boundaryHalf,conv([4 0 1],[21/4 1 1]));
  K{end+1}=ck('numeric',D.boundaryOne,conv([3 0 1],[7 2 1]));
 case 'qb02_unstable_disk'
  n=size(d.yz,1);D.yz=d.yz;D.criteria=zeros(n,1);D.norms=zeros(n,5);D.decisions=cell(n,1);
  for j=1:n
   y=d.yz(j,1);z=d.yz(j,2);P=d.base;P(3,2)=y;P(4,2)=z;N=qstab.normpoly(poly(P,o));
   D.norms(j,:)=N;D.criteria(j)=H3(N);D.decisions{j}=spectral(P,o);
   K{end+1}=ck('numeric',N,[5 2 3+y*y+z*z 2 1]);
   if y*y+z*z>3,s='stable';else,s='notStable';end
   K{end+1}=ck('status',D.decisions{j}.status,s);
  end
  K{end+1}=ck('numeric',D.criteria,4*(sum(d.yz.^2,2)-3));
  D.edgeDistanceSquared=4;D.claim='The continuum edge inequality is a source-derived geometric reference, not inferred from these samples.';
 case 'qb02_cubic_constant'
  D.x=d.x;D.decisions=cell(1,numel(d.x));D.norms=zeros(numel(d.x),7);
  for j=1:numel(d.x)
   x=d.x(j);P=d.base;P(1,1)=x;N=qstab.normpoly(poly(P,o));D.norms(j,:)=N;D.decisions{j}=spectral(P,o);
   K{end+1}=ck('numeric',N,[x*x+2 4*x+4 4*x+25 2*x+10 17 4 1]);
   K{end+1}=ck('status',D.decisions{j}.status,ref.statuses{j});
  end
  P=d.base;P(1,1)=2;D.boundaryNorm=qstab.normpoly(poly(P,o));
  K{end+1}=ck('numeric',D.boundaryNorm,conv([2 0 1],[3 6 15 4 1]));
 case 'qb02_fixed_slice'
  D.alpha=d.alpha;D.omega=d.omega;D.fixedSliceModulus=zeros(numel(d.alpha),numel(d.omega));D.normValues=complex(D.fixedSliceModulus);
  for j=1:numel(d.alpha)
   a=d.alpha(j);P=[a 1;0 0;1 0;0 0];[~,S] = qstab.valueset(poly(P,o),1i*d.omega);D.normValues(j,:)=S.points;
   for h=1:numel(d.omega),D.fixedSliceModulus(j,h)=norm(evaluate(P,d.omega(h)*qi));end
   K{end+1}=ck('numeric',S.points,a*a+1-d.omega.^2+2i*a*d.omega);
  end
  P=[0 1;0 0;1 0;0 0];D.actualZeroResidual=evaluate(P,-qj);D.projectedValue=[0;0];
  D.projectionPolynomial=[qj,qj];D.projectionDecision=spectral(D.projectionPolynomial,o);
  K=[K,{ck('numeric',D.actualZeroResidual,zeros(4,1)), ...
    ck('truth',all(D.fixedSliceModulus(:)>=ref.minimumFixedSliceModulus),true),ck('status',D.projectionDecision.status,'stable')}];
 case 'qb02_bands'
  D.widths=d.widths;D.omega=d.omega;D.members=cell(1,numel(d.widths));D.decisions=cell(size(D.members));D.values=complex(zeros(numel(d.widths),numel(d.omega)));
  for j=1:numel(d.widths)
   r=d.widths(j);P=[1 2 1;r 0 0;r 0 0;r 0 0];D.members{j}=P;D.decisions{j}=spectral(P,o);
   [~,S] = qstab.valueset(poly(P,o),1i*d.omega);D.values(j,:)=S.points;w=d.omega;
   K{end+1}=ck('numeric',S.points,(1-w.^2).^2+3*r*r-4*w.^2+4i*w.*(1-w.^2));
  end
  D.boundaryPoint=[0;-.5*ref.criticalWidth*ones(3,1)];D.boundaryResidual=evaluate(D.members{2},D.boundaryPoint);
  K=[K,{ck('status',D.decisions{1}.status,'stable'),ck('status',D.decisions{3}.status,'notStable'),ck('numeric',D.boundaryResidual,zeros(4,1))}];
 case 'qb02_nonconvex_hull'
  D.b=d.b;D.values=complex(zeros(size(d.b)));
  for j=1:numel(d.b),[~,S] = qstab.valueset(poly([d.b(j) .5 1;zeros(3,3)],o),1i);D.values(j)=S.points;end
  D.referenceCurve=(d.b-1).^2-1/4+1i*(d.b-1);D.chordEndpoints=[-.5i .5i];D.chordCentre=sum(D.chordEndpoints)/2;
  K={ck('numeric',D.values,D.referenceCurve),ck('numeric',min(abs(D.values)),ref.minimumCurveModulus),ck('numeric',D.chordCentre,0)};
 case 'qb02_nonlinear'
  D.theta=d.theta;D.norms=zeros(numel(d.theta),5);D.decisions=cell(1,numel(d.theta));
  for j=1:numel(d.theta)
   t=d.theta(j);u=[0;cos(t);sin(t);0];P=[1 1 1;zeros(3,3)];P(:,1)=P(:,1)+2*(1-t*t)*u;
   N=qstab.normpoly(poly(P,o));D.norms(j,:)=N;D.decisions{j}=spectral(P,o);Nref=conv([1 1 1],[1 1 1]);Nref(1)=Nref(1)+4*(1-t*t)^2;
   K{end+1}=ck('numeric',N,Nref);K{end+1}=ck('status',D.decisions{j}.status,ref.statuses{j});
  end
  D.boundaryResiduals=zeros(4,2);
  for j=1:2
   t=(2*j-3)/sqrt(2);u=[0;cos(t);sin(t);0];P=[1 1 1;zeros(3,3)];P(:,1)=P(:,1)+2*(1-t*t)*u;
   D.boundaryResiduals(:,j)=evaluate(P,-u);
  end
  K{end+1}=ck('numeric',D.boundaryResiduals,zeros(4,2));
 case 'qb02_correlated_marginals'
  D.theta=d.theta;D.norms=zeros(numel(d.theta),3);
  for j=1:numel(d.theta)
   c=[0;cos(d.theta(j));sin(d.theta(j));0];D.norms(j,:)=qstab.normpoly(poly([c c],o));
  end
  D.falseEven=qi;D.falseOdd=qj;D.falseDirection=-qk;D.falsePairResidual=qi+mulQ(qj,-qk);
  D.admissibilityGap=norm(D.falseEven-D.falseOdd);
  K={ck('numeric',D.norms,repmat(ref.norm,numel(d.theta),1)),ck('numeric',D.falsePairResidual,zeros(4,1)),ck('truth',D.admissibilityGap>1,true)};
 case 'qb02_rotating_factors'
  D.theta=d.theta;D.norms=zeros(numel(d.theta),5);D.representatives=complex(zeros(numel(d.theta),2));
  for j=1:numel(d.theta)
   u=[0;cos(d.theta(j));sin(d.theta(j));0];a=[1;0;0;0]+u;b=[2;0;0;1];P=[mulQ(a,b),a+b,[1;0;0;0]];
   D.norms(j,:)=qstab.normpoly(poly(P,o));[lam,rootInfo]=qstab.roots(poly(P,o));K{end+1}=ck('truth',rootInfo.isValid,true);
   [~,order]=sort(real(lam));D.representatives(j,:)=reshape(lam(order),1,[]);
  end
  K=[K,{ck('numeric',D.norms,repmat(ref.norm,numel(d.theta),1)),ck('numeric',D.representatives,repmat([-2+1i,-1+1i],numel(d.theta),1))}];
 case 'qb02_isolated_touch'
  D.grid=d.grid;D.thetaStar=ref.thetaStar;D.decisions=cell(size(d.grid));D.values=complex(zeros(size(d.grid)));
  for j=1:numel(d.grid)
   a=(d.grid(j)-ref.thetaStar)^2;P=[a 1;0 0;1 0;0 0];D.decisions{j}=spectral(P,o);[~,S] = qstab.valueset(poly(P,o),1i);D.values(j)=S.points;
   K{end+1}=ck('status',D.decisions{j}.status,'stable');
  end
  D.boundaryPolynomial=[0 1;0 0;1 0;0 0];D.boundaryResidual=evaluate(D.boundaryPolynomial,-qj);D.boundaryNorm=qstab.normpoly(poly(D.boundaryPolynomial,o));
  K=[K,{ck('truth',all(d.grid~=ref.thetaStar),true),ck('numeric',D.values,(d.grid-ref.thetaStar).^4+2i*(d.grid-ref.thetaStar).^2), ...
   ck('numeric',D.boundaryResidual,zeros(4,1)),ck('numeric',D.boundaryNorm,[1 0 1])}];
 case 'qb03_ordered_product'
  D.P=[mulQ(qi,qj),qi+qj,[1;0;0;0]];D.atMinusJ=evaluate(D.P,-qj);D.atMinusI=evaluate(D.P,-qi);D.norm=qstab.normpoly(poly(D.P,o));
  K={ck('numeric',D.atMinusJ,zeros(4,1)),ck('numeric',D.atMinusI,2*qk),ck('numeric',D.norm,ref.norm)};
 case 'qb03_square'
  D.xz=d.xz;D.criteria=zeros(size(d.xz,1),1);D.decisions=cell(size(D.criteria));
  for j=1:size(d.xz,1)
   y=[0;d.xz(j,1);d.xz(j,2);0];P=[[-1;0;0;0]+y,[1;0;0;0]+3*y,[1;0;0;0]];
   N=qstab.normpoly(poly(P,o));D.criteria(j)=H3(N);D.decisions{j}=spectral(P,o);rr=sum(d.xz(j,:).^2);
   K{end+1}=ck('numeric',D.criteria(j),4*(2*rr-1)*(9*rr+1));
   if rr>.5,s='stable';else,s='notStable';end
   K{end+1}=ck('status',D.decisions{j}.status,s);
  end
  D.edgeDistanceSquared=1;D.claim='Analytic edge distance is source data, not a conclusion from sampled members.';
 case 'qb03_feedback_order'
  D.t=d.t;D.correct=zeros(4,2,numel(d.t));D.reversed=D.correct;
  for j=1:numel(d.t)
   t=d.t(j);d0=[2;0;0;0]+t*qj;n0=[1;0;0;0]+t*qk;
   D.correct(:,:,j)=[d0+mulQ(n0,qi),[1;0;0;0]];D.reversed(:,:,j)=[d0+mulQ(qi,n0),[1;0;0;0]];
   K{end+1}=ck('numeric',D.correct(:,1,j)-D.reversed(:,1,j),t*ref.differenceDirection);
   K{end+1}=ck('numeric',D.correct(:,:,j),[2 1;1 0;2*t 0;0 0]);
  end
  D.claim='Ordered denominator identity only; no internal-stability theorem is asserted.';
 otherwise,error('qstab:paperExtension:Case','Unknown extension case: %s.',C.id);
end
D.evidenceScope=C.scope;D.isContinuumCertificate=false;
end
function A=ck(kind,a,b),A=qstab_paper_check(kind,a,b);end
function P=poly(C,o)
if o.UseQuaternion
 P=quaternion(reshape(C(1,:),1,1,[]),reshape(C(2,:),1,1,[]),reshape(C(3,:),1,1,[]),reshape(C(4,:),1,1,[]));
else,P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));end
end
function R=spectral(P,o),[~,R] = qstab.isstable(poly(P,o),'Method','spectral');end
function h=H3(N),h=N(4)*N(3)*N(2)-N(2)^2-N(4)^2*N(1);end
function Q=lift(P,h),Q=zeros(4,size(P,2)+numel(h)-1);for j=1:4,Q(j,:)=conv(P(j,:),h);end,end
function y=evaluate(P,q),y=zeros(4,1);for j=size(P,2):-1:1,y=mulQ(y,q)+P(:,j);end,end
function c=mulQ(a,b),c=[a(1)*b(1)-dot(a(2:4),b(2:4));a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];end
