function R = robuststableReport(data,varargin)
%ROBUSTSTABLEREPORT Internal implementation returning the full report.

    opts = localOptions(varargin{:});
    variable = localVariable(opts.Variable);
    method = lower(strtrim(char(opts.Method)));
    if strcmp(method,'zeroexclusion') || strcmp(method,'zeroexclude')
        R = qstab.internal.zeroexcludeReport(data,'Uncertainty',opts.Uncertainty, ...
            'Variable',opts.Variable,'MaxCells',opts.MaxCells, ...
            'MaxDepth',opts.MaxDepth,'MaxTime',opts.MaxTime, ...
            'Tolerance',opts.Tolerance,'StoreCertificates',opts.StoreCertificates, ...
            'Reduction',opts.Reduction,'MaxSubfamilies',opts.MaxSubfamilies);
        return
    elseif ~strcmp(method,'lyapunov')
        error('qstab:robuststable:Method', ...
            'Method must be lyapunov or zeroexclusion.');
    end
    if ~strcmp(lower(strtrim(char(opts.Reduction))),'none')
        error('qstab:robuststable:ReductionMethod', ...
            'Reduction requires Method=zeroexclusion; the Lyapunov path is unchanged.');
    end
    F = qstab.familycomponents(data,'Uncertainty',opts.Uncertainty);
    timer = tic;
    R = localInitialReport(F,variable,opts);
    testVariable = variable;
    if strcmp(variable,'zinv')
        testVariable = 'z';
    end

    initial = localInitialNode(F);
    stack = {initial};
    processed = 0;
    certified = 0;
    unresolved = 0;
    maxDepth = 0;
    minimumSlack = Inf;
    maximumResidual = 0;
    maximumPerturbation = 0;
    stored = {};
    lastReason = '';

    % An interval coefficient box can have a zero leading quaternion exactly
    % when each of its four leading component intervals contains zero.
    if strcmp(F.kind,'interval')
        leadIndex = F.degree+1;
        if strcmp(variable,'zinv'), leadIndex = 1; end
        if all(F.intervalLower(:,leadIndex)<=0 & F.intervalUpper(:,leadIndex)>=0)
            theta = 0.5*F.bounds(:,1)+0.5*F.bounds(:,2);
            affected = find(F.parameterIndex(:,2)==leadIndex-1);
            theta(affected) = 0;
            C = localAffineValue(F,theta);
            if all(C(:,leadIndex)==0)
                R.status = 'assumptionFailed';
                R.witness = localWitness(C,F,variable,theta,[],[]);
                R.witness.kind = 'zeroLeadingCoefficient';
                R.message = 'The family contains a zero leading coefficient in the tested variable.';
                R.assurance = 'floatingPointAssumptionWitness';
                R = localAppendPending(R,F,stack,'Not evaluated after an assumption witness');
                R = localFinish(R,timer,0,0,0,0,1,Inf,0,0,stored);
                return
            end
        end
    end

    while ~isempty(stack)
        if processed>=opts.MaxCells
            lastReason = 'MaxCells reached';
            break
        end
        if toc(timer)>=opts.MaxTime
            lastReason = 'MaxTime reached';
            break
        end
        node = stack{end};
        stack(end) = [];
        processed = processed+1;
        maxDepth = max(maxDepth,node.depth);

        [C,radius,theta,weights] = localCell(F,node);
        [zeroC,zeroTheta,zeroWeights] = localZeroLeadingWitness(F,node,C,theta,weights,variable);
        if ~isempty(zeroC)
            R.status = 'assumptionFailed';
            R.witness = localWitness(zeroC,F,variable,zeroTheta,zeroWeights,[]);
            R.witness.kind = 'zeroLeadingCoefficient';
            R.message = 'The family contains a zero leading coefficient in the tested variable.';
            R.assurance = 'floatingPointAssumptionWitness';
            R = localAppendPending(R,F,stack,'Not evaluated after an assumption witness');
            R = localFinish(R,timer,processed,certified,unresolved,maxDepth, ...
                numel(stack),minimumSlack,maximumResidual,maximumPerturbation,stored);
            return
        end

        Ct = C;
        rt = radius;
        if strcmp(variable,'zinv')
            Ct = C(:,end:-1:1);
            rt = radius(end:-1:1);
        end
        certificate = localCertificate(Ct,rt,testVariable,opts.Tolerance);
        if certificate.accepted
            certified = certified+1;
            minimumSlack = min(minimumSlack,certificate.slack);
            maximumResidual = max(maximumResidual,certificate.residualNorm);
            maximumPerturbation = max(maximumPerturbation,certificate.perturbationBound);
            if opts.StoreCertificates
                certificate.centerComponents = C;
                certificate.coefficientRadii = radius;
                certificate.parameterPoint = theta;
                certificate.vertexWeights = weights;
                certificate.cell = node;
                stored{end+1,1} = certificate; %#ok<AGROW>
            end
            continue
        end
        lastReason = certificate.reason;

        nominal = localNominal(C,variable,opts.Tolerance);
        if strcmp(nominal.status,'notStable')
            R.status = 'notStable';
            R.isStable = false;
            R.witness = localWitness(C,F,variable,theta,weights,nominal);
            R.witness.kind = 'notStrictlyStablePolynomial';
            R.message = 'A parameter point is not strictly stable according to the nominal test.';
            R.assurance = 'floatingPointCounterexample';
            R = localAppendPending(R,F,stack,'Not evaluated after a counterexample');
            R = localFinish(R,timer,processed,certified,unresolved,maxDepth, ...
                numel(stack),minimumSlack,maximumResidual,maximumPerturbation,stored);
            return
        end

        if node.depth>=opts.MaxDepth
            unresolved = unresolved+1;
            lastReason = 'MaxDepth reached in an uncertified cell';
            R.unresolvedRegions{end+1,1} = localRegion(F,node, ...
                [lastReason,': ',certificate.reason]);
            continue
        end
        [left,right,canSplit] = localSplit(F,node);
        if ~canSplit
            unresolved = unresolved+1;
            lastReason = ['Uncertified cell cannot be subdivided: ',certificate.reason];
            R.unresolvedRegions{end+1,1} = localRegion(F,node,lastReason);
            continue
        end
        stack{end+1} = right; %#ok<AGROW>
        stack{end+1} = left; %#ok<AGROW>
    end

    remaining = numel(stack);
    R = localAppendPending(R,F,stack,lastReason);
    if unresolved==0 && remaining==0
        R.status = 'stable';
        R.isStable = true;
        R.certificateComplete = true;
        R.message = 'Every cell in a covering subdivision satisfies a sufficient Lyapunov perturbation inequality.';
    else
        R.status = 'undetermined';
        R.isStable = [];
        R.message = ['No robust decision within the available certificate/budget: ',lastReason,'.'];
    end
    R = localFinish(R,timer,processed,certified,unresolved,maxDepth,remaining, ...
        minimumSlack,maximumResidual,maximumPerturbation,stored);
end

function opts = localOptions(varargin)
    ip = inputParser;
    ip.FunctionName = 'qstab.robuststable';
    addParameter(ip,'Uncertainty','',@localIsText);
    addParameter(ip,'Variable','s',@localIsText);
    addParameter(ip,'Type','polynomial',@localIsText);
    addParameter(ip,'Method','lyapunov',@localIsText);
    addParameter(ip,'Reduction','none',@localIsText);
    addParameter(ip,'MaxSubfamilies',1024,@(x)localInteger(x,1));
    addParameter(ip,'MaxCells',2048,@(x)localInteger(x,1));
    addParameter(ip,'MaxDepth',30,@(x)localInteger(x,0));
    addParameter(ip,'MaxTime',20,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&x>=0&&~isnan(x));
    addParameter(ip,'Tolerance',1e-10,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<0.1);
    addParameter(ip,'StoreCertificates',false,@(x)islogical(x)&&isscalar(x));
    parse(ip,varargin{:});
    opts = ip.Results;
    kind = lower(strtrim(char(opts.Type)));
    if ~ismember(kind,{'polynomial','scalarpolynomial','poly'})
        error('qstab:robuststable:UnsupportedType', ...
            'robuststable currently accepts scalar polynomial families. Use isstable for nominal matrices.');
    end
end

function variable = localVariable(value)
    variable = lower(strtrim(char(value)));
    if ismember(variable,{'continuous','continuous-time','hurwitz'}), variable='s'; end
    if ismember(variable,{'discrete','discrete-time','schur'}), variable='z'; end
    if ismember(variable,{'z^-1','z^{-1}','delay','d'}), variable='zinv'; end
    if ~ismember(variable,{'s','z','zinv'})
        error('qstab:robuststable:InvalidVariable','Variable must be s, z, or zinv.');
    end
end

function tf = localIsText(x)
    tf = (ischar(x)&&(isrow(x)||isempty(x))) || (isa(x,'string')&&isscalar(x));
end

function tf = localInteger(x,minimum)
    tf = isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=minimum&&x==floor(x);
end

function R = localInitialReport(F,variable,opts)
    R = struct('status','undetermined','isStable',[], ...
        'type','polynomial','variable',variable,'uncertainty',F.kind, ...
        'degree',F.degree,'method','subdivisionLyapunovPerturbation', ...
        'assurance','floatingPointSufficientCertificate','validated',false, ...
        'certificateComplete',false,'witness',[],'message','', ...
        'options',opts,'processedCells',0,'certifiedCells',0, ...
        'unresolvedCells',0,'remainingCells',0,'maximumDepth',0, ...
        'elapsedTime',0,'certificate',struct(),'certificates',{{}}, ...
        'unresolvedRegions',{{}},'parameterNames',{F.parameterNames}, ...
        'parameterIndex',F.parameterIndex);
    R.assumptions = { ...
        'Scalar polynomial with central variable and constant real uncertainty parameters.', ...
        'The declared degree is retained; leading coefficients must remain nonzero in the tested variable.', ...
        'Delay input is reversed at its declared degree; its original constant coefficient must remain nonzero.', ...
        'Reported certificates use double precision, not validated interval arithmetic.'};
end

function R = localFinish(R,timer,processed,certified,unresolved,maxDepth,remaining,minSlack,maxResidual,maxPerturbation,stored)
    R.processedCells = processed;
    R.certifiedCells = certified;
    R.unresolvedCells = unresolved+remaining;
    R.remainingCells = remaining;
    R.maximumDepth = maxDepth;
    R.elapsedTime = toc(timer);
    if isinf(minSlack), minSlack = []; end
    R.certificate = struct('minimumSlack',minSlack, ...
        'maximumLyapunovResidual',maxResidual,'maximumPerturbationBound',maxPerturbation, ...
        'description','Slack belongs to the Lyapunov inequality; it is not a spectral stability margin.');
    R.certificates = stored;
end

function node = localInitialNode(F)
    if strcmp(F.kind,'polytope')
        m = size(F.vertices,3);
        node = struct('vertices',F.vertices,'weights',eye(m),'depth',0);
    else
        node = struct('lower',F.bounds(:,1),'upper',F.bounds(:,2),'depth',0);
    end
end

function [C,radius,theta,weights] = localCell(F,node)
    count = F.degree+1;
    radius = zeros(1,count);
    theta = [];
    weights = [];
    if strcmp(F.kind,'polytope')
        m = size(node.vertices,3);
        C = sum(node.vertices/m,3);
        weights = sum(node.weights/m,2);
        for k = 1:count
            for j = 1:m
                radius(k) = max(radius(k),norm(node.vertices(:,k,j)-C(:,k)));
            end
        end
    else
        theta = 0.5*node.lower+0.5*node.upper;
        halfwidth = 0.5*node.upper-0.5*node.lower;
        C = localAffineValue(F,theta);
        componentRadii = zeros(4,count);
        triangleRadii = zeros(1,count);
        for j = 1:numel(theta)
            G = F.generators(:,:,j);
            componentRadii = componentRadii+halfwidth(j)*abs(G);
            for k = 1:count
                triangleRadii(k) = triangleRadii(k)+halfwidth(j)*norm(G(:,k));
            end
        end
        for k = 1:count
            radius(k) = min(norm(componentRadii(:,k)),triangleRadii(k));
        end
    end
end

function C = localAffineValue(F,theta)
    C = F.base;
    for j = 1:numel(theta)
        C = C+theta(j)*F.generators(:,:,j);
    end
end

function [C0,theta0,weights0] = localZeroLeadingWitness(F,node,C,theta,weights,variable)
    C0 = [];
    theta0 = [];
    weights0 = [];
    k = F.degree+1;
    if strcmp(variable,'zinv'), k=1; end
    if all(isfinite(C(:))) && all(C(:,k)==0)
        C0=C; theta0=theta; weights0=weights;
        return
    end
    if strcmp(F.kind,'polytope')
        for j = 1:size(node.vertices,3)
            if all(node.vertices(:,k,j)==0)
                C0=node.vertices(:,:,j); weights0=node.weights(:,j);
                return
            end
        end
    end
end

function cert = localCertificate(C,radius,variable,tolerance)
    cert = struct('accepted',false,'reason','','slack',NaN,'residualNorm',NaN, ...
        'perturbationBound',NaN,'delta',NaN,'leadingLowerBound',NaN, ...
        'positiveDefiniteTolerance',NaN,'X',[],'companion',[]);
    if any(~isfinite(C(:))) || any(~isfinite(radius(:)))
        cert.reason = 'Coefficient enclosure exceeds finite numerical range';
        return
    end
    scale = max([abs(C(:));radius(:)]);
    if scale==0
        cert.reason = 'Zero polynomial';
        return
    end
    C = C/scale;
    radius = radius/scale;
    lead = norm(C(:,end));
    leadRadius = radius(end);
    cert.leadingLowerBound = scale*(lead-leadRadius);
    if ~(lead>leadRadius)
        cert.reason = 'Nonzero leading coefficient not certified over this cell';
        return
    end
    n = size(C,2)-1;
    if n==0
        cert.slack = 1-leadRadius/lead;
        cert.residualNorm = 0;
        cert.perturbationBound = leadRadius/lead;
        cert.delta = 0;
        cert.accepted = cert.slack>tolerance;
        cert.reason = 'Constant-polynomial regularity bound';
        return
    end

    aunit = C(:,end)/lead;
    aunit(2:4) = -aunit(2:4);
    normalized = zeros(4,n);
    perturbations = zeros(1,n);
    relativeLeadRadius = leadRadius/lead;
    for k = 1:n
        normalized(:,k) = localQmul(aunit,C(:,k)/lead);
        perturbations(k) = (radius(k)/lead+relativeLeadRadius*norm(normalized(:,k))) ...
                          /(1-relativeLeadRadius);
    end
    delta = norm(perturbations);
    F = zeros(2*n,2*n);
    for k = 1:n-1
        F(2*k-1:2*k,2*k+1:2*k+2) = eye(2);
    end
    for k = 1:n
        F(end-1:end,2*k-1:2*k) = -localChi(normalized(:,k));
    end
    if any(~isfinite(F(:))) || ~isfinite(delta)
        cert.reason = 'Normalized companion exceeds finite numerical range';
        return
    end
    [X,solved] = localLyapunov(F,variable);
    if ~solved
        cert.reason = 'Lyapunov equation could not be solved reliably';
        return
    end
    normX = norm(X,2);
    positivityTolerance = tolerance*normX;
    [~,flag] = chol(X-positivityTolerance*eye(size(X)));
    if flag~=0 || ~(normX>0)
        cert.reason = 'A positive definite Lyapunov solution was not certified';
        return
    end
    XE = X(:,end-1:end);
    if strcmp(variable,'s')
        residual = F'*X+X*F+eye(size(F));
        perturbation = 2*norm(XE,2)*delta;
    else
        residual = F'*X*F-X+eye(size(F));
        perturbation = 2*norm(F'*XE,2)*delta+norm(X(end-1:end,end-1:end),2)*delta^2;
    end
    residualNorm = norm(residual,2);
    slack = 1-residualNorm-perturbation;
    safety = tolerance*(1+residualNorm+perturbation);
    cert.slack = slack;
    cert.residualNorm = residualNorm;
    cert.perturbationBound = perturbation;
    cert.delta = delta;
    cert.positiveDefiniteTolerance = positivityTolerance;
    cert.X = X;
    cert.companion = F;
    cert.accepted = isfinite(slack) && slack>safety;
    if cert.accepted
        cert.reason = 'Strict Lyapunov perturbation inequality';
    else
        cert.reason = 'Lyapunov perturbation bound is too large';
    end
end

function [X,ok] = localLyapunov(F,variable)
    X = [];
    ok = false;
    try
        [U,T] = schur(F,'complex');
        n = size(T,1);
        Y = zeros(n,n);
        scale = max(norm(T,1),realmin);
        for j = 1:n
            for i = 1:n
                if strcmp(variable,'s')
                    denominator = conj(T(i,i))+T(j,j);
                    if abs(denominator)<=16*eps*scale, return; end
                    numerator = -double(i==j);
                    if i>1
                        numerator = numerator-T(1:i-1,i)'*Y(1:i-1,j);
                    end
                    if j>1
                        numerator = numerator-Y(i,1:j-1)*T(1:j-1,j);
                    end
                else
                    product = conj(T(i,i))*T(j,j);
                    denominator = 1-product;
                    if abs(denominator)<=16*eps*max(1,abs(product)), return; end
                    numerator = double(i==j);
                    if i>1
                        numerator = numerator+T(1:i-1,i)'*Y(1:i-1,1:j)*T(1:j,j);
                    end
                    if j>1
                        numerator = numerator+conj(T(i,i))*Y(i,1:j-1)*T(1:j-1,j);
                    end
                end
                Y(i,j) = numerator/denominator;
            end
        end
        X = U*Y*U';
        X = 0.5*X+0.5*X';
        ok = all(isfinite(X(:)));
    catch
        X = [];
        ok = false;
    end
end

function nominal = localNominal(C,variable,tolerance)
    S = localComponentStruct(C);
    try
        nominal = qstab.internal.isstableReport(S,'Variable',variable,'Tolerance',tolerance);
    catch ME
        nominal = struct('status','undetermined','isStable',[], ...
            'message',['Nominal calculation failed: ',ME.message]);
    end
end

function R = localAppendPending(R,F,stack,reason)
    for j = 1:numel(stack)
        R.unresolvedRegions{end+1,1} = localRegion(F,stack{j},reason);
    end
end

function region = localRegion(F,node,reason)
    region = struct('uncertainty',F.kind,'depth',node.depth,'reason',reason, ...
        'parameterBounds',[],'vertexWeights',[]);
    if strcmp(F.kind,'polytope')
        region.vertexWeights = node.weights;
    else
        region.parameterBounds = [node.lower,node.upper];
    end
end

function witness = localWitness(C,F,variable,theta,weights,nominal)
    witness = struct('kind','','coefficients',localComponentStruct(C), ...
        'components',C,'variable',variable,'theta',theta, ...
        'vertexWeights',weights,'parameterNames',{F.parameterNames}, ...
        'nominalResult',nominal);
    if strcmp(F.kind,'oneparameter')
        witness.tComponents = theta;
        witness.side = F.side;
    end
    if strcmp(F.kind,'interval')
        witness.parameterIndex = F.parameterIndex;
    end
end

function S = localComponentStruct(C)
    S = struct('a',C(1,:),'b',C(2,:),'c',C(3,:),'d',C(4,:));
end

function [left,right,ok] = localSplit(F,node)
    left = node;
    right = node;
    ok = false;
    if strcmp(F.kind,'polytope')
        m = size(node.vertices,3);
        if m<2, return; end
        scale = max(abs(node.vertices(:)));
        if scale==0, return; end
        V = node.vertices/scale;
        best = 0;
        pair = [1,1];
        for i = 1:m-1
            for j = i+1:m
                distance = norm(V(:,:,i)-V(:,:,j),'fro');
                if distance>best
                    best=distance; pair=[i,j];
                end
            end
        end
        if best==0 || ~isfinite(best), return; end
        i = pair(1); j = pair(2);
        middle = 0.5*node.vertices(:,:,i)+0.5*node.vertices(:,:,j);
        if isequal(middle,node.vertices(:,:,i)) || isequal(middle,node.vertices(:,:,j))
            return
        end
        middleWeights = 0.5*node.weights(:,i)+0.5*node.weights(:,j);
        left.vertices(:,:,j) = middle;
        left.weights(:,j) = middleWeights;
        right.vertices(:,:,i) = middle;
        right.weights(:,i) = middleWeights;
    else
        halfwidth = 0.5*node.upper-0.5*node.lower;
        score = zeros(size(halfwidth));
        for j = 1:numel(score)
            score(j) = halfwidth(j)*norm(F.generators(:,:,j),'fro');
        end
        [best,index] = max(score);
        if isempty(best) || best==0 || ~isfinite(best), return; end
        middle = 0.5*node.lower(index)+0.5*node.upper(index);
        if ~(middle>node.lower(index) && middle<node.upper(index)), return; end
        left.upper(index) = middle;
        right.lower(index) = middle;
    end
    left.depth = node.depth+1;
    right.depth = node.depth+1;
    ok = true;
end

function c = localQmul(a,b)
    c = [a(1)*b(1)-a(2:4).'*b(2:4); ...
         a(1)*b(2:4)+b(1)*a(2:4)+cross(a(2:4),b(2:4))];
end

function Z = localChi(q)
    z = complex(q(1),q(2));
    w = complex(q(3),q(4));
    Z = [z,w;-conj(w),conj(z)];
end
