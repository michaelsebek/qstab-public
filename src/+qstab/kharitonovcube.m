function [K,info] = kharitonovcube(data,varargin)
%KHARITONOVCUBE Exact qB02 reduction of a CT quaternion interval polynomial.
%   A qstab.family definition supplies its variable, grade and parameter model.
% [K,INFO]=qstab.kharitonovcube({L,U}) returns an ordinary affine-family
% struct (center,generators,bounds), accepted directly by zeroexclude.
% Coefficients are ascending. Parameters 1:4 are even real,i,j,k and 5:8
% are odd real,i,j,k, all in [0,1]. Some generators can be identically zero.
% K is contained in the original interval family and has EXACTLY its norm
% value sets at each imaginary-axis frequency, in exact arithmetic.
% Stability of 256 cube VERTICES alone is NOT a robustness test.
% No corresponding independent-even/odd claim is made for z or zinv.
% Use zeroexclude(...,'Reduction','kharitonov-faces') for pair-free faces.
%
%   Examples (ascending coefficients):
%   [K,I] = qstab.kharitonovcube({[1 1],[2 1]});
%   [tf,R] = qstab.zeroexclude(K,'Uncertainty','affine');
%   disp(R.status);
%
%   Interpretation:
%   This is the CT eight-parameter continuous cube, not a list of 256 independently sufficient vertex polynomials. Its independent even/odd construction does not apply to a unit-circle boundary.
%
%   Reference: docs/reference/kharitonovcube.html.
%   See also qstab.about.
try
if qstab.internal.isFamilyDefinition(data)
    varargin=qstab.internal.familyOptions(data,varargin,{'Variable',data.variable,'Grade',data.grade},{'Variable'});
end
ip=inputParser;ip.FunctionName='qstab.kharitonovcube';
addParameter(ip,'Variable','s',@(x)(ischar(x)&&isrow(x))||(isstring(x)&&isscalar(x)));
parse(ip,varargin{:});v=lower(strtrim(char(ip.Results.Variable)));
if ~ismember(v,{'s','hurwitz','continuous'})
    error('qstab:kharitonovcube:Variable','The exact independent-even/odd reduction is for s only.');
end
F=qstab.familycomponents(data,'Uncertainty','interval');
[A,info]=qstab.internal.cubefamily(F);
if ~info.representable
    error('qstab:kharitonovcube:CoefficientRange','Endpoint differences exceed the representable range.');
end
generators=cell(1,8);
for j=1:8,generators{j}=localPolynomial(A.generators(:,:,j));end
K=struct('center',localPolynomial(A.base),'generators',{generators}, ...
    'bounds',A.bounds,'parameterNames',{A.parameterNames});
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
function P=localPolynomial(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
