function T = boxtrial(M,r,variable,method,maxCells,maxDepth,maxTime,tol,store)
%BOXTRIAL Whole-family decision at a scale, with checked original components.
% Two original-box diagonal corners are NEGATIVE-ONLY probes. A stable
% corner is never promoted to whole-family stability. In particular, the CT
% E/O zero-exclusion backend need not reconstruct an unstable noncentral
% member by itself. Positive decisions still come only from that selected
% whole-family backend. Time limits are shared with the corner checks.
started=tic;
T=struct('status','undetermined','radius',r,'processedCells',0,'result',[], ...
    'witness',[],'endpoints',[],'message','','roundingOverbound',false, ...
    'backendCalled',false,'memberChecks',{{}},'negativeScreeningOnly',true, ...
    'validated',false);
D=r*M.weights;
if ~isfinite(r)||any(~isfinite(D(:)))
    T.message='Scaled component widths overflowed.';return
end
if any(M.weights(:)>0 & D(:)==0) && r>0
    T.message='A positive component width underflowed.';return
end
% Fixed components stay BIT-IDENTICAL. Adjustable endpoints are conservatively
% padded in double. This safety padding is not a validated interval library.
amp=abs(M.C)+D;
if any(~isfinite(amp(:)))
    T.message='Endpoint magnitudes overflowed.';return
end
guard=zeros(size(D));mask=M.weights>0;
guard(mask)=16*eps(amp(mask));
L=M.C-D-guard;U=M.C+D+guard;
L(~mask)=M.C(~mask);U(~mask)=M.C(~mask);
if any(~isfinite([L(:);U(:)]))
    T.message='Padded endpoints are not finite.';return
end
Plo=localParts(L);Phi=localParts(U);T.endpoints={Plo,Phi};T.roundingOverbound=any(guard(:)>0);
% Candidate construction uses UNPADDED original bounds. The helper retracts
% rounding-outside coordinates inward and independently checks the resulting
% member. These are two corners, NOT an exhaustive vertex-only stability test.
previous=[];
for corner=1:2
    if toc(started)>=maxTime,break;end
    if corner==1,C=M.C-D;origin='originalLowerCorner';
    else,C=M.C+D;origin='originalUpperCorner';end
    if isequal(C,M.C)||isequal(C,previous),continue;end
    previous=C;
    [member,check]=qstab.internal.boxcandidate(M,C,r,variable,tol);
    check.origin=origin;T.memberChecks{end+1}=check; %#ok<AGROW>
    if ~isempty(member)&&member.valid&&member.insideRequestedBox
        member.source=origin;T.witness=member;T.status=check.status;
        T.message=check.message;return
    end
end
T.backendCalled=true;
Q=qstab.internal.robuststableReport(T.endpoints,'Uncertainty','interval','Variable',variable, ...
    'Method',method,'MaxCells',maxCells,'MaxDepth',maxDepth, ...
    'MaxTime',max(0,maxTime-toc(started)), ...
    'Tolerance',tol,'StoreCertificates',store);
T.result=Q;T.processedCells=Q.processedCells;T.message=Q.message;
if strcmp(Q.status,'stable')
    T.status='stable';return
end
if ~any(strcmp(Q.status,{'notStable','assumptionFailed'}))||isempty(Q.witness),return;end
% The robust solver worked on a slightly OUTER box. First confirm its actual
% member. A larger-radius negative witness may bound boxradius from above,
% but can never disprove stability of the original smaller requested box.
W=qstab.internal.witnessDataReport(T.endpoints,Q,'Uncertainty','interval','Tolerance',tol);
if ~W.isAdmissible||isempty(W.components),T.message='Backend witness could not be reconstructed.';return;end
C=W.components;delta=C-M.C;
if any(delta(~mask)~=0),T.message='A backend witness changed a fixed component.';return;end
ratios=zeros(size(delta));ratios(mask)=abs(delta(mask))./M.weights(mask);actual=max(ratios(:));
if ~isfinite(actual),T.message='Witness component norm overflows.';return;end
if strcmp(Q.status,'notStable')&&W.isConfirmedCounterexample
    kind='unstableMember';negativeStatus='notStable';
elseif strcmp(Q.status,'assumptionFailed')&&W.isDegreeLoss
    kind='degreeLoss';negativeStatus='assumptionFailed';
else
    T.message='Backend candidate lacks an independent confirmed negative check.';return
end
inside=actual<=r && all(abs(delta(mask))<=D(mask));
T.witness=struct('kind',kind,'valid',true,'variable',variable,'grade',M.grade, ...
    'coefficientComponents',C,'perturbationComponents',delta,'polynomial',localParts(C), ...
    'weightedNorm',actual,'weights',M.weights,'spectralCheck',W, ...
    'requestedRadius',r,'insideRequestedBox',inside, ...
    'fixedComponentsPreserved',true,'boundaryComponents',[], ...
    'retractedToOriginalBox',false,'source','paddedBackend','validated',false);
if inside
    T.status=negativeStatus;return
end
T.status='upperWitnessOnly';
T.message='Only a larger-radius negative member is confirmed; not a counterexample to this requested box.';
if toc(started)>=maxTime,return;end
% Clipping an unstable member can make it stable. Never reuse the old
% spectral decision: the helper checks the moved polynomial from scratch.
[member,check]=qstab.internal.boxcandidate(M,C,r,variable,tol);
check.origin='retractedBackendCandidate';T.memberChecks{end+1}=check;
if ~isempty(member)&&member.valid&&member.insideRequestedBox
    member.source='retractedBackendCandidate';T.witness=member;
    T.status=check.status;T.message=check.message;
end
end
function P=localParts(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
