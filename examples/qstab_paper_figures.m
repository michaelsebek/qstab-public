function F=qstab_paper_figures(varargin)
%QSTAB_PAPER_FIGURES Ordinary data for eight selected article illustrations.
% No graphics, exports, RNG/path/directory changes or robust solver calls.
% Numerical curves and source-formula references are stored separately.
% Options: Figures={} (all eight), UseQuaternion=true, Samples=101 (odd).
% These are regenerated mathematical illustrations, NOT identical typography
% or a claim to reproduce all numbered figures of qB01/qB02/qB03.
ip=inputParser;ip.FunctionName='qstab_paper_figures';ip.PartialMatching=false;
addParameter(ip,'Figures',{});
addParameter(ip,'UseQuaternion',true,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Samples',101,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=5&&x<=2001&&x==fix(x)&&mod(x,2)==1);
parse(ip,varargin{:});o=ip.Results;
allIds={'qb01_sweep','qb02_counter','qb02_hull','qb02_bands','qb02_width','qb02_nonlinear','qb03_parabola','qb03_tetrahedron'};
ids=selection(o.Figures,allIds);
if o.UseQuaternion
 try,quaternion(0,0,0,0);catch ME,error('qstab:paperFigures:Quaternion','Native quaternion unavailable: %s',ME.message);end
end
F=struct('schema','qstab.paper-figures.1','toolbox',qstab.about(),'runtime',version, ...
 'nativeQuaternionInputs',o.UseQuaternion,'options',o,'isCertificate',false, ...
 'completePaperReproduction',false,'ok',true,'status','COMPLETED','items',struct([]));
for j=1:numel(ids)
 item=makeItem(ids{j},o);qstab_validate_paper_figure(item);
 if isempty(F.items),F.items=item;else,F.items(end+1)=item;end
 F.ok=F.ok&&all(cellfun(@(c)c.ok&&c.matched,item.checks));
end
if ~F.ok,F.status='REFERENCE_ISSUES';end
F.note='Polylines sample the stated analytic curves; plots alone certify no continuum property. Original paper typography is not replicated.';
end
function I=makeItem(id,o)
n=o.Samples;I=template(id);K={};
switch id
 case 'qb01_sweep'
  I=meta(I,'qB01','fig:sweep; eq:sweep-example','Left and right one-parameter frequency margins','frequency omega','radius');
  w=linspace(0,3,n).';P=[1 2 1;0 0 0;1 0 0;0 0 0];Q=[0 0;0 1;0 0;1 0];
  L=qstab.freqmargin(poly(P,o),poly(Q,o),w,'Side','left');R=qstab.freqmargin(poly(P,o),poly(Q,o),w,'Side','right');
  rL=sqrt(2+w.^4./(1+w).^2);rR=sqrt(w.^4+2*w.^2-4*w+2)./(w+1);
  I.curves=[curve(w,L,[],'left parameter','computedSamples','-','none'),curve(w,R,[],'right parameter','computedSamples','--','none'), ...
   curve([0;sqrt(3)-1],[sqrt(2);sqrt(2)*(2-sqrt(3))],[],'source critical values','reference','none','o')];
  I.raw=struct('omega',w,'left',L,'right',R);I.reference=struct('left',rL,'right',rR);
  K={ck(L(:),rL),ck(R(:),rR)};
 case 'qb02_counter'
  I=meta(I,'qB02','fig:counter; eq:counterG','Stable endpoints do not exclude an interior boundary zero','Re N_p(i omega)','Im N_p(i omega)');
  x=linspace(-1,1,n).';A=[2 1 1;1 -1 0;-2 -1 0;1 1 0];B=A;B(2,2)=1;
  [~,S] = qstab.valueset({poly(A,o),poly(B,o)},1i*sqrt(5),'Weights',[(1-x)/2,(1+x)/2]);v=S.points;
  I.curves=[zcurve(v,'actual sampled image','computedSamples','-','none'),zcurve(v([1,end]),'coefficient endpoints','point','none','o'),zeroCurve()];
  I.raw=S;I.reference=-5*x.^2+2i*sqrt(5)*x;K={ck(v,I.reference)};I.equalAxes=true;
 case 'qb02_hull'
  I=meta(I,'qB02','fig:hull; eq:hullcurve','Actual norm curve and its convex hull are different','Re N_p(i)','Im N_p(i)');
  y=linspace(-.75,.75,n).';V={poly([.25 .5 1;zeros(3,3)],o),poly([1.75 .5 1;zeros(3,3)],o)};
  [~,S] = qstab.valueset(V,1i,'Weights',[(.75-y)/1.5,(.75+y)/1.5]);v=S.points;
  I.curves=[zcurve(v,'actual sampled curve','computedSamples','-','none'), ...
    zcurve(v([1,end]),'closing hull chord (not the curve)','outerBoundary',':','none'), ...
    zcurve([-.5i;.5i],'chord through zero (not the curve)','guide','--','o'),zeroCurve()];
  I.raw=S;I.reference=y.^2-.25+1i*y;K={ck(v,I.reference)};I.equalAxes=true;
 case {'qb02_bands','qb02_width'}
  if strcmp(id,'qb02_bands')
   I=meta(I,'qB02','fig:bands; eq:bandregion','Parabolic value-set bands at selected frequencies','Re N_p(i omega)','Im N_p(i omega)');
   widths=ones(1,5);omegas=[0 .5 1 sqrt(2) 2];
  else
   I=meta(I,'qB02','fig:width; eq:rcritical','Interval width: separation, contact, inclusion','Re N_p(i)','Im N_p(i)');
   widths=[1 2/sqrt(3) 1.3];omegas=ones(1,3);
  end
  b=linspace(1,2,n).';styles={'-','--',':','-.','-'};I.raw=struct('b',b,'widths',widths,'omega',omegas,'left',{{}},'right',{{}});I.reference=struct('left',{{}},'right',{{}});
  for j=1:numel(widths)
   r=widths(j);w=omegas(j);P=[1 2 1;zeros(3,3)];Q=P;Q(1,1)=2;
   [~,L] = qstab.valueset({poly(P,o),poly(Q,o)},1i*w,'Weights',[2-b,b-1]);
   P(2:4,1)=r;Q(2:4,1)=r;[~,R] = qstab.valueset({poly(P,o),poly(Q,o)},1i*w,'Weights',[2-b,b-1]);
   l=(b-w*w).^2-4*w*w+4i*w*(b-w*w);u=l+3*r*r;
   I.raw.left{j}=L.points;I.raw.right{j}=R.points;I.reference.left{j}=l;I.reference.right{j}=u;
   K=[K,{ck(L.points,l),ck(R.points,u)}];
   if w==0,edge=[1;4+3*r*r];else,edge=[L.points;flipud(R.points);L.points(1)];end
   if strcmp(id,'qb02_bands'),label=sprintf('omega=%.4g, r=1',w);else,label=sprintf('r=%.6g, omega=1',r);end
   I.curves(end+1)=zcurve(edge,label,'analyticBoundary',styles{j},'none');
  end
  I.curves(end+1)=zeroCurve();I.subtitle='Sampled analytic band boundaries; interpolation is not a zero-exclusion proof.';
 case 'qb02_nonlinear'
  I=meta(I,'qB02','fig:nonlinear; eq:nonlinearV','Nonlinear rotating uncertainty: exact segment images','Re N_p(i omega)','Im N_p(i omega)');
  ws=[0 .5 1 1.5];I.raw=complex(zeros(2,4));I.reference=I.raw;
  for j=1:4
   P=[1 1 1;zeros(3,3)];Q=P;Q(2,1)=2;
   [~,S] = qstab.valueset({poly(P,o),poly(Q,o)},1i*ws(j),'Weights',eye(2));I.raw(:,j)=S.points;
   w=ws(j);left=(1-w*w)^2-w*w+2i*w*(1-w*w);I.reference(:,j)=[left;left+4];
   I.curves(end+1)=zcurve(S.points,sprintf('omega=%.4g',w),'analyticBoundary','-','o');
  end
  I.curves(end+1)=zeroCurve();K={ck(I.raw,I.reference)};
  I.subtitle='Segment completeness follows from the source identity, not from the two computed endpoints.';
 case 'qb03_parabola'
  I=meta(I,'qB03','fig:parabola; ex:parabola','Bernstein outer polygons and the actual norm parabola','Re N_p(i)','Im N_p(i)');
  V={poly([1 1;zeros(3,2)],o),poly([3 1;zeros(3,2)],o)};a=linspace(1,3,n).';
  [~,S] = qstab.valueset(V,1i,'Weights',[(3-a)/2,(a-1)/2]);[~,E] = qstab.valueenclosure(V,1i);
  [~,Es] = qstab.valueenclosure(V,1i,'VertexCells',{[1 .5;0 .5],[.5 0;.5 1]});
  I.curves=[zcurve(S.points,'actual sampled parabola','computedSamples','-','none'), ...
   zcurve(closePolygon(E.regions{1}.polygon),'coarse Bernstein outer polygon','outerBoundary','--','none')];
  for j=1:2,I.curves(end+1)=zcurve(closePolygon(Es.regions{j}.polygon),sprintf('local outer polygon %d',j),'outerBoundary',':','none');end
  I.curves(end+1)=zcurve(E.regions{1}.controls,'Bernstein controls (not all attainable)','control','none','o');
  I.raw=struct('samples',S,'coarse',E,'refined',Es);I.reference=a.^2-1+2i*a;K={ck(S.points,I.reference)};I.equalAxes=true;
  I.subtitle='Two local outer polygons remain separate; mixed controls are not vertex-value samples.';
 case 'qb03_tetrahedron'
  I=meta(I,'qB03','fig:tetra; eq:tetra; eq:facegap','Stable entire boundary, unstable tetrahedral centre','y_i','y_j');
  I.zlabel='y_k';I.dimension=3;I.equalAxes=true;Y=2*[1 1 1;1 -1 -1;-1 1 -1;-1 -1 1];edges=nchoosek(1:4,2);
  for j=1:size(edges,1)
   p=Y(edges(j,:),:);c=curve(p(:,1),p(:,2),p(:,3),'tetrahedron edge','guide','-','none');
   c.showLegend=j==1;I.curves(end+1)=c;
  end
  I.curves(end+1)=curve(Y(:,1),Y(:,2),Y(:,3),'stable vertices','point','none','o');
  I.curves(end+1)=curve(0,0,0,'unstable centre','criticalPoint','none','x');
  D2=zeros(4,1);for j=1:4,b=Y(j,:)/2;D2(j)=4/sum(b.^2);end
  I.raw=struct('vertices',Y,'edges',edges,'faceDistanceSquared',D2);I.reference=4/3*ones(4,1);K={ck(D2,I.reference)};
  I.subtitle='3D display is illustrative; whole-face stability follows from the analytic distance bound.';
end
I.checks=K;
end
function I=template(id)
I=struct('schema','qstab.paper-figure.1','id',id,'article','','sourceRevision','', ...
 'sourceAnchor','','title','','subtitle','Sampled mathematical reconstruction; not an interval-arithmetic certificate.', ...
 'xlabel','','ylabel','','zlabel','','dimension',2,'equalAxes',false,'isCertificate',false, ...
 'curves',struct('x',{},'y',{},'z',{},'label',{},'role',{},'lineStyle',{},'marker',{},'showLegend',{}), ...
 'raw',struct(),'reference',[],'checks',{{}});
end
function I=meta(I,a,anchor,title,x,y)
I.article=a;I.sourceAnchor=anchor;I.title=title;I.xlabel=x;I.ylabel=y;
switch a
 case 'qB01',I.sourceRevision='Internal corrected Revision 2, 7 October 2026';
 case 'qB02',I.sourceRevision='Revision 3, 13 October 2026';
 case 'qB03',I.sourceRevision='Version 1, 16 October 2026';
end
end
function C=curve(x,y,z,label,role,style,marker)
C=struct('x',x(:),'y',y(:),'z',z(:),'label',label,'role',role,'lineStyle',style,'marker',marker,'showLegend',true);
end
function C=zcurve(v,label,role,style,marker),C=curve(real(v),imag(v),[],label,role,style,marker);end
function C=zeroCurve(),C=curve(0,0,[],'zero','reference','none','+');end
function p=closePolygon(p),p=p(:);if numel(p)>1&&p(1)~=p(end),p=[p;p(1)];end,end
function A=ck(a,b),A=qstab_paper_check('numeric',a,b);end
function P=poly(C,o)
if o.UseQuaternion,P=quaternion(reshape(C(1,:),1,1,[]),reshape(C(2,:),1,1,[]),reshape(C(3,:),1,1,[]),reshape(C(4,:),1,1,[]));
else,P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));end
end
function ids=selection(x,allIds)
if ischar(x)&&isrow(x)&&~isempty(x),x={x};elseif isstring(x)&&isvector(x)&&~isempty(x),x=cellstr(x(:).');end
if ~iscell(x)||(~isvector(x)&&~isempty(x)),error('qstab:paperFigures:Selection','Figures must be IDs in a cell vector.');end
for j=1:numel(x)
 if isstring(x{j})&&isscalar(x{j}),x{j}=char(x{j});end
 if ~ischar(x{j})||~isrow(x{j})||isempty(x{j}),error('qstab:paperFigures:Selection','Invalid figure ID.');end
end
if isempty(x),ids=allIds;return;end
if numel(unique(x))~=numel(x)||any(~ismember(x,allIds)),error('qstab:paperFigures:Selection','Duplicate or unknown figure ID.');end
ids=allIds(ismember(allIds,x));
end
