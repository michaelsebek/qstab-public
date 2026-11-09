function D = run_qb01_discrete(varargin)
%RUN_QB01_DISCRETE Selected analytic discrete one-quaternion margin examples.
% Return ordinary data only. test_qstab('DiscreteExamples',true) stores D in
% report.discreteMarginExamples before writing the unchanged two-entry ZIP.
ip=inputParser;addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));parse(ip,varargin{:});
D=struct('schema','qstab.discrete-examples.1','ok',false,'sourceVersion',qstab.about());
[~,D.rayForward] = qstab.raymargin([-.5 1],1,1,'Variable','z');
[~,D.rayDelay] = qstab.raymargin([1 -.5 0],[0 1 0],1,'Variable','zinv');
[~,D.rayLeading] = qstab.raymargin([-.5 1],[0 1],1,'Variable','z');
[~,D.rayProportional] = qstab.raymargin([-.5 1],[-.5 1],1,'Variable','z');
A=struct('a',.2,'b',0,'c',.3,'d',0);u=struct('a',0,'b',0,'c',0,'d',1);
[~,D.rayMatrix] = qstab.raymargin(A,1,u,'Type','matrix','Domain','discrete');
opts={'AbsTol',1e-5,'RelTol',.005,'MaxCells',2047};
[~,D.radiusForward] = qstab.stabilityradius([-.5 1],1,'Variable','z',opts{:});
[~,D.radiusDelay] = qstab.stabilityradius([1 -.5],[0 1],'Variable','zinv',opts{:});
P=struct('a',[-.2 1],'b',[0 0],'c',[-.3 0],'d',[0 0]);
[~,D.radiusNonreal] = qstab.stabilityradius(P,1,'Variable','z',opts{:});
[~,D.radiusCorrelated] = qstab.stabilityradius([-.25 1],[.5 1],'Variable','z',opts{:});
[~,D.radiusDegreeLoss] = qstab.stabilityradius([-.5 1],[-.5 1],'Variable','z',opts{:});
[~,D.budget] = qstab.stabilityradius(P,1,'Variable','z','MaxCells',1,'AbsTol',1e-12,'RelTol',1e-10);
D.references=struct('rayForward',[-.5 1.5],'rayDelay',[-.5 1.5], ...
    'rayLeading',[-.5 Inf],'rayProportional',[-1 Inf],'rayMatrix',[-sqrt(.87) sqrt(.87)], ...
    'radiusForward',.5,'radiusDelay',.5,'radiusNonreal',1-sqrt(.13),'radiusCorrelated',.5,'radiusDegreeLoss',1);
D.ok=true;
for name={'rayForward','rayDelay','rayLeading','rayProportional','rayMatrix'}
 R=D.(name{1});ref=D.references.(name{1});mask=isfinite(ref);
 D.ok=D.ok&&strcmp(R.status,'computed')&&all(abs(R.interval(mask)-ref(mask))<1e-8)&&all(isinf(R.interval(~mask)));
end
for name={'radiusForward','radiusDelay','radiusNonreal','radiusCorrelated','radiusDegreeLoss'}
 R=D.(name{1});ref=D.references.(name{1});
 D.ok=D.ok&&strcmp(R.status,'converged')&&R.lower<=ref+1e-8&&R.upper>=ref-1e-8&&R.witness.valid;
end
D.ok=D.ok&&strcmp(D.rayProportional.negativeEvent,'degreeLoss')&&strcmp(D.radiusDegreeLoss.limitingEvent,'degreeLoss')&& ...
    any(strcmp(D.budget.status,{'bounded','undetermined'}))&&D.budget.lower<=D.references.radiusNonreal;
if ip.Results.Verbose
 fprintf('qB01 discrete one-quaternion margins: ok=%d\n',D.ok);
 for name={'rayForward','rayDelay','rayLeading','rayProportional','rayMatrix'}
  R=D.(name{1});fprintf('  %-20s %-10s [% .9g, % .9g]\n',name{1},R.status,R.negative,R.positive);
 end
 for name={'radiusForward','radiusDelay','radiusNonreal','radiusCorrelated','radiusDegreeLoss'}
  R=D.(name{1});fprintf('  %-20s %-10s [%.9g, %.9g] %s\n',name{1},R.status,R.lower,R.upper,R.limitingEvent);
 end
 fprintf('  budget: %s (full angular bracket retained)\n',D.budget.status);
 fprintf('  Directions are not normalized; this is one shared quaternion parameter, not independent coefficient balls.\n');
end
end
