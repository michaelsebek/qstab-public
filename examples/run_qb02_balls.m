function out = run_qb02_balls(varargin)
%RUN_QB02_BALLS Selected weighted-coefficient-ball examples, ordinary data only.
% Retain via test_qstab('BallExamples',true); the existing report writer and
% qstab_test_results.zip location remain unchanged. These are selected qB02
% examples/analytic regressions, not an exhaustive reproduction of the paper.
ip=inputParser;addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});verbose=ip.Results.Verbose;
options={'AbsTol',2e-3,'RelTol',2e-3,'MaxCells',2047,'MaxTime',Inf};
out=struct('schema','qstab.ball-examples.1','ok',false,'sourceVersion',qstab.about(), ...
    'assurance','ordinary-double bounds; no validated interval arithmetic');
[~,out.ctConstant] = qstab.ballradius([2 1],[1 0],options{:});
[~,out.ctLeading] = qstab.ballradius([1 1],[0 1],options{:});
[~,out.ctQuadratic] = qstab.ballradius([1 1 1],[1 0 0],options{:});
[~,out.independentBalls] = qstab.ballradius([1 1],[1 1],'Exponent',Inf,options{:});
[~,out.schur] = qstab.ballradius([-.5 1],[1 0],'Variable','z',options{:});
[~,out.delay] = qstab.ballradius([1 -.5 0],[0 1 0],'Variable','zinv',options{:});
[~,out.budget] = qstab.ballradius([1 1 1],[1 0 0], ...
    'AbsTol',1e-12,'RelTol',1e-12,'MaxCells',1,'MaxTime',Inf);
P=struct('a',[2 .7 1],'b',[.1 .2 0],'c',[.4 -.3 0],'d',[-.2 .1 0]);
[r,I]=qstab.ballvalue(P,[.3 .5 .2],-.2+1.3i,'Exponent',2,'Witness',true);
out.noncommutingPoint=struct('center',P,'weights',[.3 .5 .2], ...
    'argument',-.2+1.3i,'radius',r,'details',I);
[xx,yy]=meshgrid(linspace(-3,0,17),linspace(1,5,17));
Q=struct('a',[2 1],'b',[0 0],'c',[3 0],'d',[0 0]);
[~,out.spectralSet] = qstab.ballset(Q,[1 0],complex(xx,yy),.5);
labels={'ctConstant','ctLeading','ctQuadratic','independentBalls','schur','delay'};
exact=[2,1,sqrt(3)/2,1/sqrt(2),.5,.5];okay=true;
for j=1:numel(labels)
    R=out.(labels{j});okay=okay&&strcmp(R.status,'converged')&& ...
        R.lower<=exact(j)+1e-10&&R.upper>=exact(j)-1e-10;
end
out.analyticMargins=struct('names',{labels},'values',exact);
out.ok=okay&&strcmp(out.ctLeading.limitingEvent,'degreeLoss')&& ...
    strcmp(out.budget.status,'bounded')&&I.valid&&I.witnesses{1}.valid&& ...
    out.spectralSet.isSampled&&~out.spectralSet.continuumCertified;
if verbose
    fprintf('qB02 weighted coefficient balls: ok=%d\n',out.ok);
    for j=1:numel(labels)
        R=out.(labels{j});fprintf('  %-19s %-10s [%.8g, %.8g] %s\n',labels{j},R.status,R.lower,R.upper,R.limitingEvent);
    end
    fprintf('  budget: %s (global bracket retained)\n',out.budget.status);
    fprintf('  noncommuting witness residual: %.3g\n',I.witnesses{1}.evaluationResidual);
    fprintf('  Spectral-set classifications apply only at stored query points.\n');
end
end
