function [n,info] = normpoly(P,varargin)
%NORMPOLY Real scalar norm polynomial, with ascending coefficients.
%   N = qstab.normpoly(P) forms N(v)=p(v)*conj(p)(v), where conjugation
%   acts on coefficients and v is central. If p=a+b*i+c*j+d*k, then
%   N=a^2+b^2+c^2+d^2, using polynomial multiplication, not pointwise
%   squaring. N has declared grade twice that of P; no zeros are trimmed.
%   [N,INFO] = ... returns degree/grade and input metadata.
%   ...,'Grade',G has the same meaning as in qstab.polycomponents.
%
%   Examples (ascending coefficients):
%   [n,I] = qstab.normpoly([2 1]);
%   disp(n);  % [4 4 1], ascending, not MATLAB roots order
%
%   Interpretation:
%   Returns a real ascending coefficient vector of declared grade 2*g. Multiplication takes place in the central polynomial ring. N(z) uses complex squares, not absolute squares; N(lambda) is not generally |p(lambda)|^2.
%
%   Reference: docs/reference/normpoly.html.
%   See also qstab.about.

[C,pinfo] = qstab.polycomponents(P,varargin{:});
n = zeros(1,2*pinfo.grade+1);
for k = 1:4
    n = n + conv(C(k,:),C(k,:));
end
if any(isnan(n) | isinf(n))
    error('qstab:ArithmeticOverflow','Norm-polynomial coefficients overflowed; rescale the input.');
end
if ~pinfo.isZero && n(2*pinfo.actualDegree+1)==0
    error('qstab:ArithmeticUnderflow','The leading norm-polynomial coefficient underflowed; rescale the input.');
end
info = struct('grade',2*pinfo.grade, ...
    'actualDegree',2*pinfo.actualDegree,'input',pinfo, ...
    'coefficientOrder','ascending','assurance','floating-point');
end
