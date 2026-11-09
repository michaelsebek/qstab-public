function A = qstab_benchmark_assess(fixture,R,tolerance)
%QSTAB_BENCHMARK_ASSESS Check semantics and reconstruct a returned witness.
% Executed OUTSIDE timing. A stable certificate is checked for report
% consistency, not re-proved with validated arithmetic. Counterexamples
% are reconstructed in the original family and checked by an independent
% nominal SPECTRAL route. Unknown is not instability and is never rankable.
A=struct('ok',true,'status','','resolved',false,'comparisonEligible',false, ...
    'witnessChecked',false,'witnessMember',false,'witnessNominalStatus','', ...
    'membershipResidual',NaN,'referenceCompatible',true,'message','');
if ~isstruct(R)||~isscalar(R)||~isfield(R,'status')
    A.ok=false;A.message='Solver did not return a scalar status report.';return
end
A.status=R.status;
if ~any(strcmp(R.status,{'stable','notStable','undetermined','assumptionFailed'}))
    A.ok=false;A.message='Unknown solver status.';return
end
if ~isfield(R,'isStable') || ...
        (strcmp(R.status,'stable')&&~isequal(R.isStable,true)) || ...
        (strcmp(R.status,'notStable')&&~isequal(R.isStable,false)) || ...
        (any(strcmp(R.status,{'undetermined','assumptionFailed'}))&&~isempty(R.isStable))
    A.ok=false;A.message='Status and isStable disagree.';return
end
ref=fixture.referenceStatus;
if strcmp(R.status,'stable')
    A.referenceCompatible=strcmp(ref,'stable');
elseif strcmp(R.status,'notStable')
    A.referenceCompatible=any(strcmp(ref,{'notStable','boundary'}));
elseif strcmp(R.status,'assumptionFailed')
    A.referenceCompatible=strcmp(ref,'assumptionFailed');
end
if ~A.referenceCompatible
    A.ok=false;A.message='Resolved result contradicts the analytic fixture reference.';return
end
if strcmp(R.status,'undetermined')
    A.message='Unresolved is retained; elapsed time is not time to a solution.';return
end
if strcmp(R.status,'stable')
    A.ok=isfield(R,'certificateComplete')&&isequal(R.certificateComplete,true);
    if isfield(R,'unresolvedCells'),A.ok=A.ok&&R.unresolvedCells==0;end
    A.resolved=A.ok;A.comparisonEligible=A.ok;
    if A.ok,A.message='Complete numerical certificate reported (not an interval proof).';
    else,A.message='A stable report is missing complete coverage.';end
    return
end
F=qstab.familycomponents(fixture.data,'Uncertainty',fixture.uncertainty);
if ~isfield(R,'witness')||~isstruct(R.witness)||~isscalar(R.witness)
    A.ok=false;A.message='A resolved negative/assumption report has no witness.';return
end
W=R.witness;
C=[];
if isfield(W,'originalComponents')
    C=W.originalComponents;
elseif isfield(W,'coefficients')
    C=qstab.polycomponents(W.coefficients);
elseif isfield(W,'components')
    C=W.components;
    % Legacy direct zero-exclusion witnesses are in the TEST variable;
    % the Lyapunov engine stores its ORIGINAL-variable coefficient member.
    if strcmp(fixture.variable,'zinv') && ...
            ~strcmp(R.method,'subdivisionLyapunovPerturbation')
        C=C(:,end:-1:1);
    end
end
if ~isnumeric(C)||~isreal(C)||~isequal(size(C),[4,F.degree+1])||any(~isfinite(C(:)))
    A.ok=false;A.message='Witness coefficients are absent or have the wrong original grade.';return
end
scale=max(1,max(abs(C(:))));tol=100*tolerance*scale;
A.witnessChecked=true;
switch F.kind
    case 'interval'
        A.membershipResidual=max([0;F.intervalLower(:)-C(:);C(:)-F.intervalUpper(:)]);
        A.witnessMember=A.membershipResidual<=tol;
    case 'polytope'
        w=[];
        if isfield(W,'originalWeights'),w=W.originalWeights;
        elseif isfield(W,'vertexWeights'),w=W.vertexWeights;end
        w=w(:);
        if numel(w)==size(F.vertices,3)&&isreal(w)&&all(isfinite(w))
            D=zeros(size(C));for k=1:numel(w),D=D+w(k)*F.vertices(:,:,k);end
            A.membershipResidual=max(abs(D(:)-C(:)));
            A.witnessMember=min(w)>=-100*tolerance && ...
                abs(sum(w)-1)<=100*tolerance && A.membershipResidual<=tol;
        end
    otherwise
        t=[];
        if isfield(W,'parameter'),t=W.parameter;
        elseif isfield(W,'theta'),t=W.theta;end
        t=t(:);
        if numel(t)==size(F.bounds,1)&&isreal(t)&&all(isfinite(t))
            D=F.base;for k=1:numel(t),D=D+t(k)*F.generators(:,:,k);end
            A.membershipResidual=max(abs(D(:)-C(:)));
            A.witnessMember=all(t>=F.bounds(:,1)-tol)&& ...
                all(t<=F.bounds(:,2)+tol)&&A.membershipResidual<=tol;
        end
end
if ~A.witnessMember
    A.ok=false;A.message='Witness did not reconstruct in the original correlated family.';return
end
if strcmp(R.status,'assumptionFailed')
    lead=F.degree+1;if strcmp(fixture.variable,'zinv'),lead=1;end
    A.ok=all(C(:,lead)==0);
    A.message='Degree loss is a separate diagnosis, excluded from speed comparisons.';
    if ~A.ok,A.message='An assumption witness did not have an exactly zero required coefficient.';end
    return
end
P=struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
[~,S] = qstab.isstable(P,'Variable',fixture.variable,'Method','spectral','Tolerance',tolerance);
A.witnessNominalStatus=S.status;
A.ok=strcmp(S.status,'notStable');A.resolved=A.ok;
A.comparisonEligible=A.ok&&strcmp(ref,'notStable');
if A.ok,A.message='Admissible original-variable counterexample agrees with spectral nominal check.';
else,A.message='Independent spectral check did not resolve the reported negative witness.';end
end
