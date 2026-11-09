function [D,K]=qstab_paper_evaluate(C,o)
%QSTAB_PAPER_EVALUATE Execute one catalogue item without writes or plotting.
% K contains explicit reference checks, D retains unmodified solver outputs.
D=struct();K={};d=C.data;ref=C.reference;
common={'MaxTime',o.MaxTime,'MaxCells',o.MaxCells};
switch C.id
 case 'qb01_ray_sides'
  P=poly(d.P0,o);Q=poly(d.P1,o);u=poly(d.u,o);
  [~,D.left] = qstab.raymargin(P,Q,u,'Side','left');[~,D.right] = qstab.raymargin(P,Q,u,'Side','right');
  K={check('numeric',[D.left.positive,D.right.positive],ref.positive), ...
     check('truth',isinf(D.left.negative)&&D.left.negative<0&&isinf(D.right.negative)&&D.right.negative<0,true)};
 case {'qb01_radius_left','qb01_radius_right'}
  if strcmp(C.id,'qb01_radius_left'),side='left';else,side='right';end
  [~,D.result] = qstab.stabilityradius(poly(d.P0,o),poly(d.P1,o),'Side',side, ...
      'AbsTol',1e-4,'RelTol',5e-3,common{:});
  K={check('bracket',D.result,ref.radius)};
 case 'qb01_reentry'
  [~,D.result] = qstab.raymargin(poly(d.P0,o),poly(d.P1,o),1);
  D.alpha=[0 .5 1];D.members=cell(1,3);D.decisions=cell(1,3);expected={'stable','notStable','stable'};
  for j=1:3
   P=d.P0;P(:,1:3)=P(:,1:3)+D.alpha(j)*d.P1;
   D.members{j}=P;[~,D.decisions{j}] = qstab.isstable(poly(P,o));
   K{end+1}=check('status',D.decisions{j}.status,expected{j});
  end
  K{end+1}=check('numeric',D.result.interval,ref.interval);
 case 'qb02_segment'
  D.members=cell(1,3);D.decisions=cell(1,3);D.H3=zeros(1,3);
  for j=1:3
   P=d.base;P(2,2)=ref.x(j);D.members{j}=P;
   [~,D.decisions{j}] = qstab.isstable(poly(P,o));N=qstab.normpoly(poly(P,o));h=quartic(N);D.H3(j)=h(5);
  end
  K={check('status',D.decisions{1}.status,'stable'),check('status',D.decisions{2}.status,'stable'), ...
     check('status',D.decisions{3}.status,'notStable'),check('numeric',D.H3,ref.H3)};
  D.x=linspace(-1,1,41).';V={d.base,d.base};V{1}(2,2)=-1;V{2}(2,2)=1;
  [~,D.values] = qstab.valueset(polys(V,o),1i*sqrt(5),'Weights',[(1-D.x)/2,(1+D.x)/2]);
  D.expectedCurve=-5*D.x.^2+2i*sqrt(5)*D.x;
  K{end+1}=check('numeric',D.values.points,D.expectedCurve);
  K{end+1}=check('truth',~D.values.isCertificate,true);
 case {'qb02_fullbox_vertices','qb02_fullbox_edge'}
  D=enumerateBox(d,o);K={check('numeric',D.minima,ref.minimum), ...
      check('truth',D.vertexCount==256&&all(strcmp(D.statuses,'stable')),true), ...
      check('truth',all(d.interior(:)>=d.L(:))&&all(d.interior(:)<=d.U(:)),true), ...
      check('status',D.interiorResult.status,'notStable')};
  if isfield(ref,'edgeH3'),K{end+1}=check('numeric',D.interiorCriteria(5),ref.edgeH3);end
 case {'qb02_edge_witness','qb02_fullbox_margin_point'}
  [D.radius,D.point]=qstab.boxvalue(poly(d.P,o),d.W,poly(d.point,o),'Witness',true);
  K={check('numeric',D.radius,ref.radius),check('truth',all(D.point.valid(:)),true)};
  w=D.point.witnesses{1};
  if ~isempty(w)
   K{end+1}=check('truth',w.valid,true);
   K{end+1}=check('numeric',w.evaluationResidual,0);
   K{end+1}=check('truth',isequal(w.coefficientComponents(d.W==0),d.P(d.W==0)),true);
  else,K{end+1}=check('truth',false,true);end
  D.claim='Pointwise attaining witness only; the source global proof is not re-run by this example.';
 case {'qb02_ball_constant','qb02_ball_joint2','qb02_ball_joint1','qb02_degree_loss'}
  [~,D.result] = qstab.ballradius(poly(d.P,o),d.weights,'Exponent',d.exponent, ...
      'AbsTol',1e-4,'RelTol',5e-3,common{:});
  K={check('bracket',D.result,ref.radius)};
  if ~isempty(D.result.witness),K{end+1}=check('truth',D.result.witness.valid,true);end
  if strcmp(C.id,'qb02_degree_loss')
   K{end+1}=check('status',D.result.limitingEvent,ref.event);
   D.claim='Regularity loss, not a finite imaginary-axis zero.';
  end
 case 'qb03_parabola'
  V=polys(d.vertices,o);D.t=linspace(0,1,41).';D.a=1+2*D.t;
  [~,D.values] = qstab.valueset(V,1i,'Weights',[1-D.t,D.t]);
  [~,D.enclosure] = qstab.valueenclosure(V,1i);
  D.referenceCurve=D.a.^2-1+2i*D.a;
  K={check('numeric',D.values.points,D.referenceCurve), ...
     check('numeric',D.enclosure.regions{1}.controls,ref.controls), ...
     check('truth',~D.values.isCertificate&&~D.enclosure.isCertificate,true)};
 case 'qb03_bernstein'
  V=polys(d.vertices,o);[~,D.coarse] = qstab.valueenclosure(V,0);
  [~,D.refined] = qstab.valueenclosure(V,0,'VertexCells',{[1 .5;0 .5],[.5 0;.5 1]});
  [~,D.result] = qstab.zeroexclude(V,'StoreCertificates',true,common{:});
  a=real(D.coarse.regions{1}.controls);
  K={check('numeric',[min(a),max(a)],ref.coarseBounds), ...
     check('truth',~D.coarse.zeroExcludedNumerically&&D.refined.zeroExcludedNumerically,true), ...
     check('status',D.result.status,'stable')};
 case 'qb03_master'
  D.norm=qstab.normpoly(poly(d.P,o));D.value=evaluate(d.P,d.point);
  [~,D.nominal] = qstab.isstable(poly(d.P,o));
  K={check('numeric',D.norm,ref.norm),check('numeric',D.value,zeros(4,1)), ...
     check('status',D.nominal.status,'undetermined')};
 case 'qb03_tetrahedron'
  [~,D.result] = qstab.zeroexclude(polys(d.vertices,o),common{:});
  D.faceCentres=zeros(4,3);D.faceDistanceSquared=zeros(1,4);D.faceOrthogonality=zeros(4,2);
  for j=1:4
   face=d.parameterVertices(setdiff(1:4,j),:);c=mean(face,1);
   D.faceCentres(j,:)=c;D.faceDistanceSquared(j)=dot(c,c);
   D.faceOrthogonality(j,:)=[dot(c,face(2,:)-face(1,:)),dot(c,face(3,:)-face(1,:))];
  end
  K={check('status',D.result.status,'notStable'), ...
     check('numeric',D.faceDistanceSquared,repmat(ref.faceDistanceSquared,1,4)), ...
     check('numeric',D.faceOrthogonality,zeros(4,2))};
  D.claim='Face-distance identities plus the source analytic criterion; not a finite sample proof of every face.';
 case 'qb03_schur_tetrahedron'
  % The original coefficient-only Schur test may be undecided even at a
  % genuinely unstable centre. Keep its complete zero-exclusion result;
  % verify the actual admissible member with an EXPLICIT spectral request.
  D.centre=mean(cat(3,d.vertices{:}),3);
  D.centreWeights=ones(numel(d.vertices),1)/numel(d.vertices);
  D.centreFromWeights=zeros(size(D.centre));
  for j=1:numel(d.vertices)
   D.centreFromWeights=D.centreFromWeights+D.centreWeights(j)*d.vertices{j};
  end
  [~,D.result] = qstab.zeroexclude(polys(d.vertices,o),'Variable','z',common{:});
  [~,D.centreResult] = qstab.isstable(poly(D.centre,o),'Variable','z','Method','spectral');
  D.decisionBasis='Explicit spectral test of the admissible barycentric centre.';
  D.comparisonRole='result retains the original zero-exclusion output, including undetermined.';
  D.claim='One admissible unstable member disproves robustness; no sampled proof of stable faces.';
  K={check('numeric',D.centre,ref.centre), ...
     check('truth',all(D.centreWeights>=0)&&sum(D.centreWeights)==1,true), ...
     check('numeric',D.centreFromWeights,D.centre), ...
     check('status',D.centreResult.status,'notStable'), ...
     check('truth',any(strcmp(D.result.status,{'notStable','undetermined'})),true)};
 case 'qb03_degree_loss'
  [~,D.result] = qstab.zeroexclude(polys(d.vertices,o),common{:});
  K={check('status',D.result.status,ref.status)};
 case 'qb03_correlation'
  % A coarse Bernstein hull can touch zero at an isolated frequency even
  % for this strictly stable segment. The generic refinement can exhaust
  % its depth budget. Preserve that result and use the existing whole-cell
  % Lyapunov method explicitly on the SAME correlated input, not on a box.
  V=polys(d.vertices,o);
  [~,D.correlated] = qstab.zeroexclude(V,common{:});
  [~,D.correlatedLyapunov] = qstab.robuststable(V,'Method','lyapunov', ...
      'Variable','s','StoreCertificates',true,common{:});
  [~,D.independent] = qstab.isstable(poly(d.independent,o));
  D.independentNorm=qstab.normpoly(poly(d.independent,o));
  D.decisionBasis='Complete-family Lyapunov subdivision of the original correlated segment.';
  D.comparisonRole='correlated retains the original zero-exclusion output, including undetermined.';
  D.claim='The independent member belongs to the larger component box, not to the correlated segment.';
  stableComplete=~strcmp(D.correlatedLyapunov.status,'stable')|| ...
      (D.correlatedLyapunov.certificateComplete&&D.correlatedLyapunov.unresolvedCells==0&& ...
       D.correlatedLyapunov.remainingCells==0&& ...
       numel(D.correlatedLyapunov.certificates)==D.correlatedLyapunov.certifiedCells);
  K={check('status',D.correlatedLyapunov.status,'stable'), ...
     check('truth',stableComplete,true), ...
     check('truth',any(strcmp(D.correlated.status,{'stable','undetermined'})),true), ...
     check('status',D.independent.status,'notStable'), ...
     check('numeric',D.independentNorm,ref.independentNorm)};
 otherwise,error('qstab:paper:UnknownCase','Unknown paper case %s.',C.id);
end
end
function A=check(kind,x,ref)
A=qstab_paper_check(kind,x,ref);
% Store concise assessment. Full raw data are retained once in D, not copied
% into every check (especially full radius-cell arrays).
if strcmp(kind,'bracket'),A.actual=struct('status',x.status,'lower',x.lower,'upper',x.upper);end
end
function P=poly(C,o)
if o.UseQuaternion
 P=quaternion(reshape(C(1,:),1,1,[]),reshape(C(2,:),1,1,[]), ...
   reshape(C(3,:),1,1,[]),reshape(C(4,:),1,1,[]));
else,P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));end
end
function V=polys(C,o)
V=cell(size(C));for j=1:numel(C),V{j}=poly(C{j},o);end
end
function h=quartic(N)
% Ascending real monic norm coefficients, no root computation.
g0=N(1);g1=N(2);g2=N(3);g3=N(4);
h=[g3,g1,g0,g3*g2-g1,g3*g2*g1-g1^2-g3^2*g0];
end
function D=enumerateBox(d,o)
active=find(d.L~=d.U);n=2^numel(active);
D=struct('vertexCount',n,'activeComponentIndices',active,'vertices',zeros(4,3,n), ...
 'normCoefficients',zeros(n,5),'criteria',zeros(n,5),'statuses',{cell(n,1)});
for j=1:n
 C=d.L;
 for k=1:numel(active),if bitget(uint32(j-1),k),C(active(k))=d.U(active(k));end,end
 D.vertices(:,:,j)=C;P=poly(C,o);N=qstab.normpoly(P);
 D.normCoefficients(j,:)=N;D.criteria(j,:)=quartic(N);
 [~,S] = qstab.isstable(P);D.statuses{j}=S.status;
end
D.minima=min(D.criteria,[],1);
[~,D.interiorResult] = qstab.isstable(poly(d.interior,o));
D.interiorNorm=qstab.normpoly(poly(d.interior,o));D.interiorCriteria=quartic(D.interiorNorm);
D.claim='Exhaustive enumeration of finite coefficient vertices only, not a vertex criterion for the whole family.';
end
function y=evaluate(C,q)
y=zeros(4,1);for j=size(C,2):-1:1,y=mulQ(y,q)+C(:,j);end
end
function c=mulQ(a,b)
c=[a(1)*b(1)-dot(a(2:4),b(2:4));a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
end
