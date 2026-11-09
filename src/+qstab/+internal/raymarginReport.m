function R = raymarginReport(X0, X1, u, varargin)
%RAYMARGINREPORT Internal implementation returning the full report.

    o=struct('Type','polynomial','Side','left','Variable','s',...
        'Domain','continuous','Tolerance',1e-9,'VerifyBoundary',true,'Grade',[], ...
        'MaxGuardCells',2047,'MaxGuardDepth',50,'MaxGuardTime',5);
    o=parseopts(o,varargin{:});
    typ=lower(char(o.Type)); side=lower(char(o.Side));
    if ~any(strcmp(typ,{'polynomial','matrix'})) || ...
            ~any(strcmp(side,{'left','right'}))
        error('qstab:raymargin:Option','Invalid Type or Side.');
    end
    if (strcmp(typ,'polynomial') && any(strcmpi(o.Variable,{'z','zinv'}))) || ...
            (strcmp(typ,'matrix') && strcmpi(o.Domain,'discrete'))
        R=qstab.internal.discreteray(X0,X1,u,o); return
    end
    if ~isempty(o.Grade)
        error('qstab:raymargin:Grade','Grade is currently a discrete-ray option only.');
    end
    if ~strcmpi(o.Variable,'s') || ~strcmpi(o.Domain,'continuous')
        error('qstab:raymargin:Domain','raymargin currently implements continuous-time first exit.');
    end
    if ~(isnumeric(o.Tolerance)&&isreal(o.Tolerance)&&isscalar(o.Tolerance)&&...
            isfinite(o.Tolerance)&&o.Tolerance>0&&o.Tolerance<1)
        error('qstab:raymargin:Tolerance','Tolerance must lie strictly between 0 and 1.');
    end
    if ~(isscalar(o.VerifyBoundary)&&(islogical(o.VerifyBoundary)||...
            (isnumeric(o.VerifyBoundary)&&any(o.VerifyBoundary==[0,1]))))
        error('qstab:raymargin:VerifyBoundary','VerifyBoundary must be a logical scalar.');
    end
    [ua,ub,uc,ud]=qstab.components(u);
    if numel(ua)~=1
        error('qstab:raymargin:Direction','u must be scalar.');
    end
    uv=[ua;ub;uc;ud]; un=norm(uv);
    if un==0 || ~isfinite(un)
        error('qstab:raymargin:Direction','u must be a finite nonzero quaternion.');
    end
    R=struct('status','undetermined','method','quaternion-hermitian-lyapunov-pencil',...
        'assurance','floating-point','type',typ,'side',side,'variable','s',...
        'directionComponents',uv,'directionNorm',un,'interval',[NaN NaN],...
        'candidateInterval',[NaN NaN],'negative',NaN,'positive',NaN,...
        'negativeStatus','undetermined','positiveStatus','undetermined',...
        'positiveParameterNorm',NaN,'message','','diagnostics',struct());
    if strcmp(typ,'polynomial')
        C0=qstab.polycomponents(X0); C1=qstab.polycomponents(X1);
        last0=find(any(C0~=0,1),1,'last'); last1=find(any(C1~=0,1),1,'last');
        if isempty(last0)||last0<2||(~isempty(last1)&&last1>=last0)
            R.status='assumptionFailed';
            R.message='Require degree(P0)>=1 and degree(P1)<degree(P0).'; return
        end
        C0=C0(:,1:last0); n=last0-1;
        Cp=zeros(4,n); if ~isempty(last1), Cp(:,1:last1)=C1(:,1:last1); end
        A=zeros(n,n,4); B=A;
        if n>1, A(:,:,1)=diag(ones(n-1,1),1); end
        lead=C0(:,end); leadScale=max(abs(lead)); lead=lead/leadScale;
        il=conjq(lead)/(lead'*lead);
        for k=1:n
            A(n,k,:)=reshape(-mulq(il,C0(:,k)/leadScale),1,1,4);
            if strcmp(side,'left'), b=mulq(uv,Cp(:,k)/leadScale); else, b=mulq(Cp(:,k)/leadScale,uv); end
            B(n,k,:)=reshape(-mulq(il,b),1,1,4);
        end
        nominal=qstab.internal.isstableReport(X0,'Variable','s','Tolerance',o.Tolerance);
    else
        [a,b,c,d]=qstab.components(X0);
        if ndims(a)>2, error('qstab:raymargin:MatrixShape','A0 must be a two-dimensional matrix.'); end
        A=cat(3,a,b,c,d);
        [a,b,c,d]=qstab.components(X1);
        if ndims(a)>2, error('qstab:raymargin:MatrixShape','A1 must be a two-dimensional matrix.'); end
        E=cat(3,a,b,c,d);
        if isempty(A)||ndims(A)>3||size(A,1)~=size(A,2)||~isequal(size(A),size(E))
            error('qstab:raymargin:MatrixShape','A0 and A1 must be square matrices of equal size.');
        end
        n=size(A,1); B=zeros(size(A));
        for i=1:n
            for j=1:n
                e=reshape(E(i,j,:),4,1);
                if strcmp(side,'left'), b=mulq(uv,e); else, b=mulq(e,uv); end
                B(i,j,:)=reshape(b,1,1,4);
            end
        end
        nominal=qstab.internal.isstableReport(X0,'Type','matrix','Domain','continuous','Tolerance',o.Tolerance);
    end
    R.nominal=nominal;
    if any(~isfinite(A(:)))||any(~isfinite(B(:)))
        R.message='The normalized matrix family exceeds floating-point range.'; return
    end
    if ~strcmp(nominal.status,'stable')
        if strcmp(nominal.status,'notStable'), R.status='assumptionFailed'; end
        R.message='A verified numerically Hurwitz nominal member is required.'; return
    end
    SA=asstruct(A); SB=asstruct(B);
    [L0,linfo]=qstab.lyapoperator(SA); L1=qstab.lyapoperator(SB);
    R.diagnostics.operatorDimension=linfo.dimension;
    if any(~isfinite(L0(:))) || any(~isfinite(L1(:)))
        R.message='The Lyapunov-operator pencil exceeds floating-point range.'; return
    end
    R.diagnostics.nominalOperatorRcond=rcond(L0);
    if rcond(L0)<eps*max(100,size(L0,1))
        R.message='Nominal Lyapunov operator is numerically singular.'; return
    end
    s0=norm(L0,2); s1=norm(L1,2);
    R.diagnostics.pencilScale=[s0,s1];
    R.diagnostics.pencilEigenvalueMeaning='finite alpha contacts, not their reciprocals';
    R.diagnostics.boundaryVerification=[];
    if all(L1(:)==0)
        % Test exact zero separately: a small but nonzero direction can
        % have a very distant finite first contact.
        lo=-Inf; hi=Inf; negativeOK=true; positiveOK=true;
        R.candidateInterval=[lo,hi];
        R.diagnostics.pencilEigenvalues=Inf(size(L0,1),1);
        R.diagnostics.scaledPencilEigenvalues=Inf(size(L0,1),1);
        R.diagnostics.infiniteEigenvalueCount=size(L0,1);
        R.diagnostics.eigenResidual=zeros(0,1);
        R.diagnostics.eigenUncertaintyEstimate=zeros(0,1);
        R.diagnostics.ambiguousRealEigenvalues=zeros(0,1);
    else
        if ~(isfinite(s0)&&isfinite(s1)&&s0>0&&s1>0)
            R.message='The separate pencil scales cannot be represented reliably.'; return
        end
        ratio=s0/s1;
        if ~isfinite(ratio)||ratio==0
            R.message='The parameter scale exceeds the available floating-point range.'; return
        end
        AA=-L0/s0; BB=L1/s1;
        try
            [V,D,W]=eig(AA,BB,'qz');
        catch ME
            R.message=['Generalized QZ evaluation failed: ',ME.message]; return
        end
        beta=diag(D);
        [eta,backwardResidual,pairing,roundoffBand]=localEigenChecks(AA,BB,V,W,beta);
        finite=isfinite(beta);
        exactlyReal=finite & imag(beta)==0;
        % A genuine nonreal pair remains nonreal even when it is closer to
        % the real axis than the user's requested parameter tolerance.
        % If floating-point uncertainty reaches the axis, it is unresolved.
        possible=finite & (exactlyReal | ...
            abs(imag(beta))<=max(roundoffBand,10*eta));
        accuracy=max(128*eps,o.Tolerance);
        reliable=exactlyReal & isfinite(eta) & ...
            eta<=accuracy*abs(beta) & real(beta)~=0;
        [hi,hiCandidate,positiveOK,positiveIndex]=localEndpoint(beta,eta, ...
            roundoffBand,reliable,possible,1,ratio);
        [lo,loCandidate,negativeOK,negativeIndex]=localEndpoint(beta,eta, ...
            roundoffBand,reliable,possible,-1,ratio);
        if any(isnan(beta))
            % A 0/0 generalized eigenvalue is incompatible with a regular
            % pencil with nonsingular L0 and cannot be silently ignored.
            positiveOK=false; negativeOK=false;
        end
        R.candidateInterval=[loCandidate,hiCandidate];
        R.diagnostics.pencilScaleRatio=ratio;
        R.diagnostics.pencilEigenvalues=ratio*beta;
        R.diagnostics.scaledPencilEigenvalues=beta;
        R.diagnostics.infiniteEigenvalueCount=sum(isinf(beta));
        R.diagnostics.eigenResidual=backwardResidual;
        R.diagnostics.eigenUncertaintyEstimate=ratio*eta;
        R.diagnostics.finiteEigenpairPairing=pairing;
        R.diagnostics.realEigenvalueMask=exactlyReal;
        R.diagnostics.reliableRealEigenvalueMask=reliable;
        R.diagnostics.ambiguousRealEigenvalues=ratio*beta(possible & ~reliable);
        R.diagnostics.selectedEndpointIndices=[negativeIndex,positiveIndex];
    end
    if o.VerifyBoundary
        H0=qstab.chi(SA); H1=qstab.chi(SB);
        endpoints=[lo,hi]; check=NaN(2,3);
        for k=1:2
            endpointOK=(k==1 && negativeOK)||(k==2 && positiveOK);
            if isfinite(endpoints(k)) && endpointOK
                H=H0+endpoints(k)*H1; ev=eig(H);
                scale=max(1,norm(H,2));
                check(k,:)=[min(abs(real(ev)))/scale,max(real(ev))/scale,scale];
                if check(k,1)>max(100*eps,o.Tolerance)||check(k,2)>max(100*eps,o.Tolerance)
                    if k==1, negativeOK=false; else, positiveOK=false; end
                end
            end
        end
        R.diagnostics.boundaryVerification=check;
    end
    if negativeOK
        R.negativeStatus='computed'; R.negative=lo; R.interval(1)=lo;
    end
    if positiveOK
        R.positiveStatus='computed'; R.positive=hi; R.interval(2)=hi;
        R.positiveParameterNorm=hi*un;
    end
    if positiveOK&&negativeOK
        R.status='computed';
        R.message='Connected open Hurwitz interval containing zero; floating-point pencil evaluation.';
    elseif positiveOK||negativeOK
        R.status='partial';
        R.message='One first-exit endpoint is resolved; the other is numerically ambiguous.';
    else
        R.message='Pencil contacts are ill-conditioned or ambiguous; inspect candidateInterval and diagnostics.';
    end
end

function [eta,residual,pairing,roundoffBand]=localEigenChecks(A,B,V,W,beta)
% Per-finite-eigenvalue first-order sensitivity estimates. They are numerical
% diagnostics, not rigorous perturbation bounds. Defective eigenvalues have
% a small left/right B pairing and therefore cannot pass the accuracy gate.
    m=numel(beta);
    eta=Inf(m,1); residual=NaN(m,1); pairing=NaN(m,1);
    roundoffBand=NaN(m,1);
    normA=norm(A,2); normB=norm(B,2);
    arithmetic=64*eps*max(1,size(A,1));
    for j=1:m
        if ~isfinite(beta(j)), continue; end
        v=V(:,j); w=W(:,j); nv=norm(v); nw=norm(w);
        if ~(isfinite(nv)&&isfinite(nw)&&nv>0&&nw>0), continue; end
        v=v/nv; w=w/nw;
        eigenScale=max(1,abs(beta(j)));
        K=A/eigenScale-(beta(j)/eigenScale)*B;
        backward=max(norm(K*v),norm(K'*w));
        coefficientScale=normA/eigenScale+(abs(beta(j))/eigenScale)*normB;
        denominator=abs(w'*B*v);
        pairing(j)=denominator;
        residual(j)=backward/max(coefficientScale,realmin);
        roundoffBand(j)=arithmetic*eigenScale;
        if denominator>0
            eta(j)=((backward+arithmetic*coefficientScale)/denominator)*eigenScale;
        end
    end
end

function [value,candidate,ok,index]=localEndpoint(beta,eta,roundoffBand,reliable,possible,direction,ratio)
% Work in the scaled parameter to avoid reciprocating tiny pencil values.
% An uncertain candidate can block only endpoints that it could precede.
    oriented=direction*real(beta);
    eligible=find(reliable & oriented>0);
    index=0; best=Inf; bestUncertainty=0;
    if ~isempty(eligible)
        [best,position]=min(oriented(eligible));
        index=eligible(position);
        bestUncertainty=10*eta(index)+roundoffBand(index);
    end
    candidates=oriented(possible & oriented>0);
    if isempty(candidates), approximate=Inf; else, approximate=min(candidates); end
    candidate=direction*(ratio*approximate);
    ok=true;
    uncertain=find(possible & ~reliable);
    for j=reshape(uncertain,1,[])
        uncertainty=max(roundoffBand(j),10*eta(j));
        if ~isfinite(uncertainty)
            ok=false;
            continue
        end
        left=oriented(j)-uncertainty;
        right=oriented(j)+uncertainty;
        if right>0 && left<=best+bestUncertainty
            ok=false;
        end
    end
    value=direction*(ratio*best);
    if isfinite(best) && (~isfinite(value)||value==0)
        ok=false;
    end
end

function z=mulq(x,y)
    z=[x(1)*y(1)-x(2:4)'*y(2:4);...
        x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
function z=conjq(x)
    z=[x(1);-x(2:4)];
end
function S=asstruct(A)
    S=struct('a',A(:,:,1),'b',A(:,:,2),'c',A(:,:,3),'d',A(:,:,4));
end
function o=parseopts(o,varargin)
    if mod(numel(varargin),2)~=0, error('qstab:raymargin:Options','Use name/value pairs.'); end
    keys=fieldnames(o);
    for k=1:2:numel(varargin)
        if ~(ischar(varargin{k})||(isstring(varargin{k})&&isscalar(varargin{k})))
            error('qstab:raymargin:Options','Option names must be text.');
        end
        j=find(strcmpi(char(varargin{k}),keys),1);
        if isempty(j), error('qstab:raymargin:Options','Unknown option %s.',char(varargin{k})); end
        o.(keys{j})=varargin{k+1};
    end
end
