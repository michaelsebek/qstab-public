function [radius,info] = ballvalue(P,weights,z,varargin)
%BALLVALUE Exact-formula pointwise coefficient-ball distance to a zero class.
% [R,INFO]=qstab.ballvalue(P,GAMMA,Z,'Exponent',ETA) evaluates the smallest
% || (|delta_k|/gamma_k) ||_eta that permits an algebraic zero in each class
% represented by complex Z. Zero weights fix their coefficients. ETA is in
% [1,Inf] (default 2). Coefficients/weights are ascending at the same grade.
% Z can be a finite numeric array of any shape; R preserves that shape.
%
% This is the qB02 weighted-coefficient-ball spectral-set formula, evaluated
% via the smaller singular value of a 2-by-2 complex polynomial value. It is
% NOT abs(p(i*omega)) on one arbitrarily chosen quaternion slice. Its small
% singular values use SVD, not a subtraction of two nearly equal squares.
% 'Witness',true reconstructs an admissible critical polynomial per query.
% 'Grade',G retains zero padding; weights must include all G+1 coefficients.
% Z is always in the INPUT polynomial variable. No delay reversal occurs.
% INFO.lower/upper are guarded floating-point pointwise estimates, not
% outward-rounded interval proofs. Inspect INFO.valid before using results.
%
%   Examples (ascending coefficients):
%   [r,I] = qstab.ballvalue([2 1],[1 0],0,'Witness',true);
%   disp(r); disp(I.valid);
%
%   Interpretation:
%   Z describes a complete similarity class. Full quaternion coefficient balls are isotropic. Witness polynomials may have nonreal coefficients even for a real nominal polynomial. Pointwise validity is not a global stability certificate.
%
%   Reference: docs/reference/ballvalue.html.
%   See also qstab.about.
try
ip=inputParser;
addParameter(ip,'Exponent',2);
addParameter(ip,'Grade',[]);
addParameter(ip,'Witness',false,@(x)islogical(x)&&isscalar(x));
addParameter(ip,'Tolerance',1e-10,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<1);
parse(ip,varargin{:});o=ip.Results;
M=qstab.internal.ballmodel(P,weights,o.Exponent,o.Grade);
if ~(isnumeric(z)&&~isempty(z)&&all(isfinite(z(:))))
    error('qstab:ballvalue:Arguments','Z must be a nonempty finite real or complex numeric array.');
end
radius=NaN(size(z));
info=struct('schema','qstab.ballvalue.1','exponent',M.exponent,'weights',M.weights, ...
    'grade',M.grade,'classRepresentatives',complex(real(double(z)),abs(imag(double(z)))), ...
    'valid',false(size(z)),'lower',zeros(size(z)),'upper',Inf(size(z)), ...
    'nearZero',false(size(z)),'scaledMinimumModulus',NaN(size(z)), ...
    'scaledEnvelope',NaN(size(z)),'coefficientScale',M.coefficientScale, ...
    'weightScale',M.weightScale,'powerScaleLog',NaN(size(z)), ...
    'witnesses',{cell(size(z))},'messages',{cell(size(z))}, ...
    'assurance','floating-point pointwise bounds','validated',false);
for j=1:numel(z)
    E=qstab.internal.ballpoint(M,double(z(j)),o.Witness,o.Tolerance);
    radius(j)=E.radius;info.valid(j)=E.valid;info.lower(j)=E.lower;info.upper(j)=E.upper;
    info.nearZero(j)=E.nearZero;info.scaledMinimumModulus(j)=E.scaledMinimumModulus;
    info.scaledEnvelope(j)=E.scaledEnvelope;info.powerScaleLog(j)=E.powerScaleLog;
    info.witnesses{j}=E.witness;info.messages{j}=E.message;
end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
