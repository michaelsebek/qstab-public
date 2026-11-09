function cases = qstab_cases_robust()
%QSTAB_CASES_ROBUST Mathematical and contract tests of scalar robustness.
%   CASES is an N-by-2 cell array {description,@testFunction}. The caller
%   runs each zero-argument test and catches assertion/implementation errors.
%   These cases use only numeric arrays and native component structs, so
%   the quaternion product class is not required for this part of the suite.
%   The same tests may be run in MATLAB or as auxiliary checks in Octave.

    cases = { ...
        'Robust CT quaternion interval: real part gives an analytic certificate', @localCTInterval; ...
        'Robust CT quaternion polytope and retained numerical certificates', @localCTPolytope; ...
        'Unstable interval returns an admissible actual polynomial witness', @localUnstableInterval; ...
        'qB01 re-entry cubic: stable vertices and an unstable interior', @localReentry; ...
        'Interior degree loss is an assumption witness, not an unstable finite root', @localDegreeLoss; ...
        'Affine correlations survive subdivision and varying leading coefficients', @localAffineCorrelation; ...
        'One-quaternion parameter keeps the side of multiplication', @localOneParameterSides; ...
        'Schur stability of a full quaternion coefficient box', @localSchurInterval; ...
        'Delay reversal retains grade and handles zero/infinity endpoints', @localDelay; ...
        'Nonzero constant families and the zero-polynomial exception', @localConstants; ...
        'Cell, depth and time budgets preserve an undetermined result and its region', @localBudgets; ...
        'Malformed family bounds and representations are rejected', @localInputValidation ...
    };
end

function localCTInterval()
    % p=s+a0, Re(a0)>=1. The unique algebraic class has Re(lambda)<=-1.
    L = localPoly([1,1],[-0.05,0],[-0.05,0],[-0.05,0]);
    U = localPoly([1.2,1],[0.05,0],[0.05,0],[0.05,0]);
    [~,R] = qstab.robuststable({L,U},'Uncertainty','interval','Variable','s', ...
        'MaxTime',Inf,'MaxCells',64);
    localAssertStable(R);
    assert(R.degree==1,'The declared scalar degree was not retained.');
    assert(~R.validated,'A floating-point certificate must not claim interval validation.');
end

function localCTPolytope()
    % Every convex combination has Re(a0)=1, independently of its direction.
    V = {localPoly([1,1],[1,0],[0,0],[0,0]), ...
         localPoly([1,1],[0,0],[1,0],[0,0]), ...
         localPoly([1,1],[0,0],[0,0],[1,0])};
    [~,R] = qstab.robuststable(V,'Variable','s','StoreCertificates',true, ...
        'MaxTime',Inf,'MaxCells',128);
    localAssertStable(R);
    assert(numel(R.certificates)==R.certifiedCells, ...
        'Requested accepted-leaf certificates were not retained.');
    for j = 1:numel(R.certificates)
        Q = R.certificates{j};
        assert(Q.accepted && Q.slack>0,'An accepted certificate lacks a strict slack.');
        assert(~isempty(Q.X),'A nonconstant leaf certificate must retain its Lyapunov matrix.');
        assert(abs(sum(Q.vertexWeights)-1)<1e-12 && all(Q.vertexWeights>=-1e-14), ...
            'A certificate center is outside the original convex hull.');
    end
end

function localUnstableInterval()
    % Every member s+a0 has a positive real zero because a0 is negative.
    [~,R] = qstab.robuststable({[-1.2,1],[-0.8,1]},'Uncertainty','interval', ...
        'Variable','s','MaxTime',Inf,'MaxCells',16);
    assert(strcmp(R.status,'notStable') && isequal(R.isStable,false), ...
        'An actual unstable interval member was not reported as a counterexample.');
    C = R.witness.components;
    assert(C(1,1)>=-1.2 && C(1,1)<=-0.8 && C(1,2)==1, ...
        'The polynomial witness is not in the supplied interval.');
    assert(all(all(C(2:4,:)==0)) && -C(1,1)>0, ...
        'The supplied witness does not have its claimed unstable real zero.');
    assert(R.witness.theta(1)>=-1.2 && R.witness.theta(1)<=-0.8, ...
        'The original interval parameter was not retained.');
end

function localReentry()
    % qB01, real-ray re-entry example:
    % p(alpha)=s^3+(1+alpha)s^2+(1+alpha)s+13/16+3*alpha.
    % The remaining Routh inequality is (alpha-1/4)*(alpha-3/4)>0.
    p0 = [13/16,1,1,1];
    p1 = [61/16,2,2,1];
    [~,N0] = qstab.isstable(p0);
    [~,N1] = qstab.isstable(p1);
    assert(strcmp(N0.status,'stable') && strcmp(N1.status,'stable'), ...
        'The two analytically Hurwitz endpoints were not recognized.');
    [~,R] = qstab.robuststable({p0,p1},'Variable','s','MaxTime',Inf,'MaxCells',32);
    assert(strcmp(R.status,'notStable') && ~isempty(R.witness), ...
        'Stable vertices were incorrectly promoted to robust stability.');
    weights = R.witness.vertexWeights;
    assert(numel(weights)==2 && all(weights>=0) && abs(sum(weights)-1)<1e-12, ...
        'Counterexample barycentric weights are invalid.');
    a = R.witness.components(1,:);
    assert(norm(a-(weights(1)*p0+weights(2)*p1),Inf)<1e-12, ...
        'The polynomial witness disagrees with its original vertex weights.');
    alpha = weights(2);
    assert(alpha>1/4 && alpha<3/4 && a(2)*a(3)<a(1), ...
        'The witness is not in the analytically unstable interior interval.');
end

function localDegreeLoss()
    % Both endpoints are nonzero real multiples of s+1. Their midpoint is
    % the zero polynomial, so this is outside the fixed-degree assumptions.
    [~,R] = qstab.robuststable({[1,1],[-1,-1]},'Variable','s','MaxTime',Inf);
    assert(strcmp(R.status,'assumptionFailed') && isempty(R.isStable), ...
        'Interior degree loss must be distinguished from an unstable finite zero.');
    assert(strcmp(R.witness.kind,'zeroLeadingCoefficient'), ...
        'The violated leading-coefficient assumption is not identified.');
    assert(all(R.witness.components(:)==0), ...
        'The cancellation witness is not the zero polynomial.');
    assert(norm(R.witness.vertexWeights-[0.5;0.5],Inf)<1e-12, ...
        'The interior witness lost its original vertex coordinates.');
end

function localAffineCorrelation()
    % p(s,t)=(1+t*i)*(s+1), -2<=t<=2, has the fixed zero -1.
    % If the two imaginary coefficient components were made independent,
    % a0=1-2*i, a1=1+2*i would instead give Re(-a1\a0)=3/5>0.
    data = localCorrelatedFamily();
    [~,R] = qstab.robuststable(data,'Uncertainty','affine','Variable','s', ...
        'MaxTime',Inf,'MaxCells',256,'MaxDepth',20);
    localAssertStable(R);
end

function localOneParameterSides()
    % With p1=i and t=j, t*p1=-k whereas p1*t=k.
    p1 = localPoly(0,1,0,0);
    bounds = zeros(4,2);
    bounds(3,:) = [1,1];
    data = struct('p0',[2,1],'p1',p1,'side','left','bounds',bounds);
    FL = qstab.familycomponents(data,'Uncertainty','oneparameter');
    [~,RL] = qstab.robuststable(data,'Uncertainty','oneparameter','MaxTime',Inf);
    data.side = 'right';
    FR = qstab.familycomponents(data,'Uncertainty','oneparameter');
    [~,RR] = qstab.robuststable(data,'Uncertainty','oneparameter','MaxTime',Inf);
    assert(isequal(FL.generators(:,1,3),[0;0;0;-1]), ...
        'Left multiplication failed j*i=-k.');
    assert(isequal(FR.generators(:,1,3),[0;0;0;1]), ...
        'Right multiplication failed i*j=k.');
    assert(isequal(FL.generators(:,1,1),[0;1;0;0]), ...
        'The real component of the quaternion parameter has the wrong generator.');
    localAssertStable(RL);
    localAssertStable(RR);
end

function localSchurInterval()
    % p=z+a0 and |a0|<=sqrt(0.2^2+3*0.1^2)<1 for the entire box.
    L = localPoly([-0.2,1],[-0.1,0],[-0.1,0],[-0.1,0]);
    U = localPoly([0.2,1],[0.1,0],[0.1,0],[0.1,0]);
    [~,R] = qstab.robuststable({L,U},'Uncertainty','interval','Variable','z', ...
        'MaxTime',Inf,'MaxCells',64);
    localAssertStable(R);
    assert(strcmp(R.variable,'z'),'The reported discrete variable is wrong.');
end

function localDelay()
    [~,D] = qstab.robuststable({[1,-0.5]},'Variable','zinv','MaxTime',Inf);
    [~,Z] = qstab.robuststable({[-0.5,1]},'Variable','z','MaxTime',Inf);
    localAssertStable(D);
    localAssertStable(Z);
    assert(strcmp(D.variable,'zinv'),'Delay reporting silently changed the user variable.');

    % A missing highest delay coefficient corresponds to a stable z=0 mode.
    [~,D0] = qstab.robuststable({[1,0]},'Variable','zinv','MaxTime',Inf);
    localAssertStable(D0);
    assert(D0.degree==1,'The delay grade was trimmed and the zero mode was lost.');

    % A zero current-sample coefficient is the different causal endpoint.
    [~,bad] = qstab.robuststable({[0,1]},'Variable','zinv','MaxTime',Inf);
    assert(strcmp(bad.status,'assumptionFailed') && isempty(bad.isStable), ...
        'A zero delay constant coefficient escaped the causal leading check.');
end

function localConstants()
    % The segment 1+t*i, -1<=t<=1, never vanishes although its first ball
    % enclosure touches zero. Subdivision must retain its real-part constraint.
    V = {localPoly(1,1,0,0),localPoly(1,-1,0,0)};
    [~,R] = qstab.robuststable(V,'Variable','s','MaxTime',Inf,'MaxCells',64);
    localAssertStable(R);
    assert(R.degree==0,'A constant polynomial family acquired a spurious degree.');
    [~,zero] = qstab.robuststable({0},'Variable','s','MaxTime',Inf);
    assert(strcmp(zero.status,'assumptionFailed') && isempty(zero.isStable), ...
        'The zero polynomial was treated as a stable nonzero constant.');
end

function localBudgets()
    data = localCorrelatedFamily();
    [~,A] = qstab.robuststable(data,'Uncertainty','affine','MaxCells',1, ...
        'MaxDepth',20,'MaxTime',Inf);
    [~,B] = qstab.robuststable(data,'Uncertainty','affine','MaxCells',64, ...
        'MaxDepth',0,'MaxTime',Inf);
    [~,C] = qstab.robuststable(data,'Uncertainty','affine','MaxCells',64, ...
        'MaxDepth',20,'MaxTime',0);
    reports = {A,B,C};
    for j = 1:numel(reports)
        R = reports{j};
        assert(strcmp(R.status,'undetermined') && isempty(R.isStable), ...
            'An exhausted budget was silently interpreted as stability or instability.');
        assert(isempty(R.witness) && ~R.certificateComplete, ...
            'An exhausted budget fabricated a witness or complete certificate.');
        assert(R.unresolvedCells>0 && numel(R.unresolvedRegions)==R.unresolvedCells, ...
            'The unresolved parameter regions were not retained.');
        regions = R.unresolvedRegions;
        lower = zeros(numel(regions),1);
        upper = lower;
        for k = 1:numel(regions)
            bounds = regions{k}.parameterBounds;
            assert(isequal(size(bounds),[1,2]) && bounds(1)<=bounds(2), ...
                'An unresolved affine parameter interval is malformed.');
            assert(~isempty(regions{k}.reason),'An unresolved cell has no reason.');
            lower(k)=bounds(1); upper(k)=bounds(2);
        end
        [lower,index] = sort(lower);
        upper = upper(index);
        assert(lower(1)==-2 && upper(end)==2, ...
            'The unresolved cells lost an endpoint of the original interval.');
        for k = 1:numel(lower)-1
            assert(lower(k+1)<=upper(k), ...
                'The reported unresolved cover has a gap.');
        end
    end
    assert(A.processedCells==1 && B.maximumDepth==0 && C.processedCells==0, ...
        'A cell, depth, or immediate time budget was not respected.');
end

function localInputValidation()
    localMustError(@()qstab.robuststable({[2,1],[1,1]},'Uncertainty','interval'));
    localMustError(@()qstab.robuststable({[1,1],[1,1,1]},'Uncertainty','interval'));
    data = struct('center',[1,1],'generators',{{1}},'bounds',[0,NaN]);
    localMustError(@()qstab.familycomponents(data,'Uncertainty','affine'));
    data.bounds = [1,-1];
    localMustError(@()qstab.familycomponents(data,'Uncertainty','affine'));
    one = struct('p0',[1,1],'p1',1,'side','unspecified','bounds',zeros(4,2));
    localMustError(@()qstab.familycomponents(one,'Uncertainty','oneparameter'));
    localMustError(@()qstab.robuststable({eye(2)},'Uncertainty','polytope'));
    localMustError(@()qstab.robuststable({[1,1]},'Variable','unsupported'));
end

function data = localCorrelatedFamily()
    generator = localPoly([0,0],[1,1],[0,0],[0,0]);
    data = struct('center',[1,1],'generators',{{generator}},'bounds',[-2,2]);
end

function P = localPoly(a,b,c,d)
    P = struct('a',a,'b',b,'c',c,'d',d);
end

function localAssertStable(R)
    assert(strcmp(R.status,'stable') && isequal(R.isStable,true), ...
        'An analytically stable family was not certified by the allocated budget.');
    assert(R.certificateComplete && R.unresolvedCells==0 && isempty(R.unresolvedRegions), ...
        'A stable result does not cover the entire uncertainty family.');
    assert(R.certifiedCells>0 && R.certificate.minimumSlack>0, ...
        'A stable result lacks accepted cells with a strict sufficient inequality.');
end

function localMustError(f)
    caught = false;
    try
        f();
    catch
        caught = true;
    end
    assert(caught,'Malformed family data were not rejected.');
end
