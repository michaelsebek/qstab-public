function M = ballmodel(P,weights,exponent,grade)
%BALLMODEL Validate scalar isotropic coefficient-ball data, without classes.
if isempty(grade),[C,I]=qstab.polycomponents(P);
else,[C,I]=qstab.polycomponents(P,'Grade',grade);end
if ~(isnumeric(weights)&&isreal(weights)&&isvector(weights)&& ...
        ~isempty(weights)&&numel(weights)==size(C,2)&& ...
        all(isfinite(weights(:)))&&all(weights(:)>=0))
    error('qstab:ball:Weights','Weights must be a finite nonnegative real vector, one per ascending coefficient at the declared grade.');
end
if ~(isnumeric(exponent)&&isreal(exponent)&&isscalar(exponent)&& ...
        ~isnan(exponent)&&exponent>=1)
    error('qstab:ball:Exponent','Exponent must be a real scalar in [1,Inf].');
end
exponent=double(exponent);weights=double(weights(:).');
if exponent==1,dual=Inf;
elseif isinf(exponent),dual=1;
else,dual=1+1/(exponent-1);end
cs=max(abs(C(:)));gs=max(weights);
M=struct('C',C,'weights',weights,'exponent',exponent,'dualExponent',dual, ...
    'grade',I.grade,'inputInfo',I,'coefficientScale',cs,'weightScale',gs, ...
    'Cscaled',C,'weightsScaled',weights,'hasPerturbation',gs>0,'logScaleRatio',NaN);
if cs>0,M.Cscaled=C/cs;end
if gs>0,M.weightsScaled=weights/gs;end
if any(C(:)~=0 & M.Cscaled(:)==0)||any(weights>0 & M.weightsScaled==0)
    error('qstab:ball:DynamicRange','Normalization would erase a nonzero coefficient or weight; rescale the input.');
end
if cs>0&&gs>0,M.logScaleRatio=log(cs)-log(gs);end
end
