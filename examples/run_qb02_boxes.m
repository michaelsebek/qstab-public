function E=run_qb02_boxes(varargin)
%RUN_QB02_BOXES Component boxes, not coefficient balls or vertex-only tests.
ip=inputParser;addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));parse(ip,varargin{:});o=ip.Results;
E=struct('schema','qstab.component-box-examples.1','ok',true,'sourceVersion',qstab.about());
W=zeros(4,2);W(1,1)=1;
[~,E.ct] = qstab.boxradius([2 1],W,'RelTol',.02,'MaxTime',Inf,'MaxCells',2048);
[~,E.schur] = qstab.boxradius([-.5 1],W,'Variable','z','RelTol',.02,'MaxTime',Inf,'MaxCells',2048);
Wd=zeros(4,3);Wd(1,2)=1;
[~,E.delay] = qstab.boxradius([1 -.5 0],Wd,'Variable','zinv','RelTol',.02,'MaxTime',Inf,'MaxCells',2048);
C=[1 1 1;2 0 0;0 2.5 0;0 0 0];We=zeros(size(C));We(3:4,2)=1;
q=parts([0;-.5;0;-sqrt(3)/2]);
[edgeRadius,edgeData]=qstab.boxvalue(parts(C),We,q,'Witness',true);
E.edgePoint=struct('radius',edgeRadius,'data',edgeData);
[~,E.edge] = qstab.boxradius(parts(C),We,'BoundaryPoints',q,'RelTol',.03,'AbsTol',.005, ...
    'MaxIterations',12,'MaxCells',4096,'MaxTime',Inf);
Wl=zeros(4,2);Wl(1,2)=1;
[~,E.leading] = qstab.boxradius([1 1],Wl,'RelTol',.02,'MaxCells',2048,'MaxTime',Inf);
[~,E.budget] = qstab.boxradius([2 1],W,'MaxIterations',0);
E.references=struct('ct',2,'schur',.5,'delay',.5,'edge',2.5-sqrt(3),'leading',1);
keys=fieldnames(E.references);
for j=1:numel(keys)
 k=keys{j};R=E.(k);ref=E.references.(k);
 E.ok=E.ok&&isfinite(R.lower)&&R.lower<=ref+1e-8&&R.upper>=ref-1e-8&&~isempty(R.witness)&&R.witness.valid;
 if o.Verbose,fprintf('  box %-10s %-10s [%g, %g] %s\n',k,R.status,R.lower,R.upper,R.limitingEvent);end
end
E.ok=E.ok&&E.edgePoint.data.valid&&abs(E.edgePoint.radius-E.references.edge)<1e-8;
E.ok=E.ok&&E.budget.lower==0&&E.budget.upper>=2;
if o.Verbose
 fprintf('qB02 component boxes: ok=%d\n',E.ok);
 fprintf('  Whole-family tests supply lower bounds; actual-point gauges supply only upper candidates.\n');
 fprintf('  Finite budgets may leave bounded brackets; no direction-grid completeness is claimed.\n');
end
end
function P=parts(C)
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end
