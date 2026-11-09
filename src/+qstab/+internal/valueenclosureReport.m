function E = valueenclosureReport(data,lambda,varargin)
%VALUEENCLOSUREREPORT Internal implementation returning the full report.
ip=inputParser;ip.FunctionName='qstab.valueenclosure';
addParameter(ip,'Uncertainty','',@localText);
addParameter(ip,'VertexCells',{},@(x)iscell(x));
addParameter(ip,'ParameterBoxes',{},@(x)iscell(x));
addParameter(ip,'Tolerance',1e-10,@localTolerance);
parse(ip,varargin{:});o=ip.Results;
if ~isnumeric(lambda)||~isscalar(lambda)||~isfinite(real(lambda))||~isfinite(imag(lambda))
    error('qstab:valueenclosure:Lambda','lambda must be one finite numeric complex argument.');
end
lambda=double(lambda);F=qstab.familycomponents(data,'Uncertainty',o.Uncertainty);
E=struct('schema','qstab.value-enclosure.1','lambda',lambda,'uncertainty',F.kind, ...
    'degree',F.degree,'regions',{{}},'scope','oneInputVariableArgument', ...
    'parameterCoverage','fullOriginalDomain','isCertificate',false,'validated',false, ...
    'assurance','ordinaryDoubleOuterEnclosure','tolerance',o.Tolerance, ...
    'zeroExcludedNumerically',false,'description','');
if strcmp(F.kind,'polytope')
    if ~isempty(o.ParameterBoxes),error('qstab:valueenclosure:ParameterBoxes','Use VertexCells for a polytope.');end
    cells=o.VertexCells;m=size(F.vertices,3);
    if isempty(cells),cells={eye(m)};else,E.parameterCoverage='userSubcellsCoverageNotChecked';end
    if ~isvector(cells),error('qstab:valueenclosure:VertexCells','VertexCells must be a cell vector.');end
    for j=1:numel(cells)
        W=cells{j};
        if ~isnumeric(W)||~isreal(W)||ndims(W)>2||size(W,1)~=m||size(W,2)<1|| ...
                any(~isfinite(W(:)))||any(W(:)<0)||any(abs(sum(W,1)-1)>64*eps*max(1,m))
            error('qstab:valueenclosure:Weights','Subcell columns must be original nonnegative convex weights.');
        end
        % Normalize only the accepted sum-roundoff; record exact used weights.
        W=bsxfun(@rdivide,double(W),sum(W,1));
        C=reshape(reshape(F.vertices,[],m)*W,4,F.degree+1,size(W,2));
        k=size(C,3);y=complex(zeros(4,k));
        for h=1:k,y(:,h)=localEval(C(:,:,h),lambda);end
        controls=complex(zeros(k*(k+1)/2,1));pairs=zeros(numel(controls),2);u=0;
        for a=1:k
            for b=a:k
                u=u+1;controls(u)=y(:,a).'*y(:,b);pairs(u,:)=[a b];
            end
        end
        localFinite(controls);polygon=qstab.internal.diagnosticGeometry('hull',controls);
        region=localRegion(polygon,controls,o.Tolerance,'BernsteinPolygon');
        region.vertexWeights=W;region.parameterBounds=[];region.controlPairs=pairs;
        region.valueBounds=[min(real(controls)),max(real(controls));min(imag(controls)),max(imag(controls))];
        region.componentValues=y;E.regions{end+1,1}=region;
    end
else
    if ~isempty(o.VertexCells),error('qstab:valueenclosure:VertexCells','Use ParameterBoxes for affine or interval families.');end
    cells=o.ParameterBoxes;
    if isempty(cells),cells={F.bounds};else,E.parameterCoverage='userSubcellsCoverageNotChecked';end
    if ~isvector(cells),error('qstab:valueenclosure:ParameterBoxes','ParameterBoxes must be a cell vector.');end
    for j=1:numel(cells)
        B=cells{j};
        if ~isnumeric(B)||~isreal(B)||~isequal(size(B),size(F.bounds))||any(~isfinite(B(:)))|| ...
                any(B(:,1)>B(:,2))||any(B(:,1)<F.bounds(:,1))||any(B(:,2)>F.bounds(:,2))
            error('qstab:valueenclosure:Bounds','Each parameter box must lie inside the original domain.');
        end
        B=double(B);t=.5*B(:,1)+.5*B(:,2);r=.5*B(:,2)-.5*B(:,1);
        C=F.base;H=bsxfun(@times,F.generators,reshape(r,1,1,[]));
        for a=1:numel(t),C=C+t(a)*F.generators(:,:,a);end
        y=localEval(C,lambda);dy=complex(zeros(4,numel(t)));
        for a=1:numel(t),dy(:,a)=localEval(H(:,:,a),lambda);end
        evenBox=[];oddBox=[];
        if strcmp(F.kind,'interval') && real(lambda)==0
            % On the imaginary axis every coefficient contributes to E OR O.
            er=sum(abs(real(dy)),2);orr=sum(abs(imag(dy)),2);
            evenBox=[real(y)-er,real(y)+er];oddBox=[imag(y)-orr,imag(y)+orr];
            [emin,emax]=localSquare(evenBox);[omin,omax]=localSquare(oddBox);
            lo=0;hi=0;
            for a=1:4
                prod=[evenBox(a,1)*oddBox(a,1),evenBox(a,1)*oddBox(a,2), ...
                      evenBox(a,2)*oddBox(a,1),evenBox(a,2)*oddBox(a,2)];
                lo=lo+min(prod);hi=hi+max(prod);
            end
            bounds=[emin-omax,emax-omin;2*lo,2*hi];kind='EvenOddBoundingRectangle';
        else
            q=y.'*y;lin=2*dy.'*y;Q=dy.'*dy;d=diag(Q);off=triu(Q,1);
            lo=[real(q)-sum(abs(real(lin)))+sum(min(0,real(d)))-2*sum(abs(real(off(:)))); ...
                imag(q)-sum(abs(imag(lin)))+sum(min(0,imag(d)))-2*sum(abs(imag(off(:))))];
            hi=[real(q)+sum(abs(real(lin)))+sum(max(0,real(d)))+2*sum(abs(real(off(:)))); ...
                imag(q)+sum(abs(imag(lin)))+sum(max(0,imag(d)))+2*sum(abs(imag(off(:))))];
            bounds=[lo hi];kind='JointQuadraticRectangle';
        end
        localFinite(bounds);polygon=qstab.internal.diagnosticGeometry('rectangle',bounds);
        region=localRegion(polygon,polygon,o.Tolerance,kind);
        region.parameterBounds=B;region.vertexWeights=[];region.controlPairs=[];
        region.valueBounds=bounds;region.evenBox=evenBox;region.oddBox=oddBox;
        region.centerComponents=C;region.componentValues=y;E.regions{end+1,1}=region;
    end
end
E.zeroExcludedNumerically=all(cellfun(@(x)x.zeroExcludedNumerically,E.regions));
E.description=['Local outer sets at one argument; retain their UNION. ', ...
    'Zero inside an outer set is inconclusive, not an actual zero. ', ...
    'No robust-stability conclusion follows from this function alone.'];
end
function r=localRegion(p,c,tol,kind)
d=qstab.internal.diagnosticGeometry('distance',p);s=max([1;abs(c(:))]);
r=struct('kind',kind,'polygon',p(:),'controls',c(:),'distanceToZero',d, ...
    'errorAllowance',tol*s,'zeroExcludedNumerically',d>tol*s,'inflationRadius',0);
end
function y=localEval(C,z)
y=complex(zeros(4,1));for k=size(C,2):-1:1,y=y*z+C(:,k);end
end
function [lo,hi]=localSquare(B)
near=max(0,max(B(:,1),-B(:,2)));lo=sum(near.^2);hi=sum(max(abs(B),[],2).^2);
end
function localFinite(x)
if any(~isfinite(real(x(:)))|~isfinite(imag(x(:))))
    error('qstab:valueenclosure:Overflow','The enclosure exceeds finite numerical range; no plot or verdict was fabricated.');
end
end
function tf=localText(x)
tf=(ischar(x)&&(isrow(x)||isempty(x)))||(isstring(x)&&isscalar(x));
end
function tf=localTolerance(x)
tf=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<.01;
end
