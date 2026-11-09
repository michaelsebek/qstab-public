function [radius,info] = boxvalue(P,widths,lambda,varargin)
%BOXVALUE Component-box distance to a zero at an ACTUAL quaternion point.
% [R,I]=qstab.boxvalue(P,W,LAMBDA) minimizes max(abs(delta_knu)/W_knu)
% subject to right evaluation of P+delta vanishing at LAMBDA. W is an
% explicit NONNEGATIVE REAL 4-by-(grade+1) array, or an a,b,c,d component
% struct. Zero entries fix the corresponding real components exactly.
% P has the usual ascending scalar coefficient representation.
%
% IMPORTANT: unlike ballvalue, this evaluates an ACTUAL quaternion, NOT a
% whole similarity class. Complex LAMBDA lies in the 1,i slice. No arbitrary
% direction grid certifies the class minimum. A native quaternion array or
% an ordinary a,b,c,d point struct can specify other directions explicitly.
%
% R is a numerical point estimate; I.lower/I.upper retain the bracket.
% A bounded pointData status does not assert that the optimum was resolved.
% The qB02 zonotope gauge is evaluated in its intrinsic real dimension.
% At imaginary points the generators are compressed to at most eight, also
% for high-degree polynomials. Primal/dual checks are numerical only.
% Options: Evaluation='right'|'left', Grade=[], Tolerance=1e-10,
% MaxFacets=4096 (total tuple budget per point, recursive reconstruction
% included), Witness=false. Query shapes and coefficient order are retained.
%
%   Examples (ascending coefficients):
%   W = zeros(4,2); W(1,1) = 1;
%   [r,I] = qstab.boxvalue([2 1],W,0,'Witness',true);
%   disp(r); disp(I.valid);
%
%   Interpretation:
%   LAMBDA is an ACTUAL quaternion point, not just a class representative. Component boxes are anisotropic. A point witness supplies only an upper candidate for a global radius. Zero width fixes the component exactly.
%
%   Reference: docs/reference/boxvalue.html.
%   See also qstab.about.
try
ip=inputParser;
addParameter(ip,'Evaluation','right',@(x)ischar(x)||(isstring(x)&&isscalar(x)));
addParameter(ip,'Grade',[]);
addParameter(ip,'Tolerance',1e-10,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<1);
addParameter(ip,'MaxFacets',4096,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x==fix(x));
addParameter(ip,'Witness',false,@(x)islogical(x)&&isscalar(x));
parse(ip,varargin{:});o=ip.Results;o.MaxFacets=double(o.MaxFacets);o.Tolerance=double(o.Tolerance);
side=lower(strtrim(char(o.Evaluation)));
if ~any(strcmp(side,{'right','left'})),error('qstab:boxvalue:Evaluation','Evaluation must be right or left.');end
M=qstab.internal.boxmodel(P,widths,o.Grade);[a,b,c,d]=qstab.components(lambda);
if isempty(a),error('qstab:boxvalue:Point','Point array must not be empty.');end
radius=NaN(size(a));
info=struct('schema','qstab.boxvalue.1','scope','actualQuaternionPoints', ...
    'evaluationSide',side,'grade',M.grade,'widths',M.weights, ...
    'valid',false(size(a)),'lower',zeros(size(a)),'upper',Inf(size(a)), ...
    'pointData',{cell(size(a))},'witnesses',{cell(size(a))}, ...
    'classMinimumClaim',false,'validated',false,'assurance','floating-point primal/dual checks');
for j=1:numel(a)
    E=qstab.internal.boxpoint(M,[a(j);b(j);c(j);d(j)],side,o.Tolerance,o.MaxFacets,o.Witness);
    radius(j)=E.radius;info.valid(j)=E.valid;info.lower(j)=E.lower;info.upper(j)=E.upper;
    info.witnesses{j}=E.witness;info.pointData{j}=E;
end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
