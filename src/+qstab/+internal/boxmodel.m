function M = boxmodel(P,weights,grade)
%BOXMODEL Scalar centre and NONNEGATIVE REAL COMPONENT widths, not ball radii.
if isempty(grade),[C,I]=qstab.polycomponents(P);
else,[C,I]=qstab.polycomponents(P,'Grade',grade);end
if isstruct(weights)
    W=qstab.polycomponents(weights);
elseif isnumeric(weights)&&isreal(weights)
    W=double(weights);
else
    error('qstab:box:Weights','Use a real 4-by-(grade+1) width array or a real nonnegative a,b,c,d struct.');
end
if ~isequal(size(W),size(C))||any(~isfinite(W(:)))||any(W(:)<0)
    error('qstab:box:Weights','Widths must be finite nonnegative real components, one for each coefficient at the declared grade.');
end
cs=max(abs(C(:)));ws=max(W(:));
M=struct('C',C,'weights',W,'grade',I.grade,'inputInfo',I, ...
    'coefficientScale',cs,'weightScale',ws,'Cscaled',C,'Wscaled',W, ...
    'hasPerturbation',ws>0,'logScaleRatio',0);
if cs>0,M.Cscaled=C/cs;end
if ws>0,M.Wscaled=W/ws;end
if any(C(:)~=0 & M.Cscaled(:)==0)||any(W(:)>0 & M.Wscaled(:)==0)
    error('qstab:box:DynamicRange','Normalization would erase a nonzero coefficient or component width.');
end
if cs>0&&ws>0,M.logScaleRatio=log(cs)-log(ws);end
end
