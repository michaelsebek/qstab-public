function [a,b,c,d] = components(Q)
%COMPONENTS Real components of a numeric or native quaternion array.
%   [A,B,C,D] = qstab.components(Q) represents Q = A+B*i+C*j+D*k.
%   Q may be a MathWorks quaternion array, a real or complex numeric array,
%   or an ordinary scalar struct with real numeric fields a,b,c,d of the
%   same size. Complex numeric inputs use the 1,i subalgebra. Outputs are
%   real double arrays with the same shape as Q (or its component fields).
%   No normalization as a rotation is performed. NaN and Inf are rejected.
%
%   Examples (ascending coefficients):
%   Q = quaternion(1,2,3,4);
%   [a,b,c,d] = qstab.components(Q);
%   disp([a b c d]);
%
%   Interpretation:
%   The four outputs are real numeric component arrays of exactly the original shape. They are not normalized rotation parameters.
%
%   Reference: docs/reference/components.html.
%   See also qstab.about.

if isa(Q,'quaternion')
    [a,b,c,d] = parts(Q);
elseif isnumeric(Q) || islogical(Q)
    a = real(Q);
    b = imag(Q);
    c = zeros(size(Q));
    d = zeros(size(Q));
elseif isstruct(Q) && isscalar(Q) && all(isfield(Q,{'a','b','c','d'}))
    a = Q.a; b = Q.b; c = Q.c; d = Q.d;
    if ~(isnumeric(a) && isnumeric(b) && isnumeric(c) && isnumeric(d))
        error('qstab:InvalidComponents','Component fields must be numeric.');
    end
    if ~(isreal(a) && isreal(b) && isreal(c) && isreal(d))
        error('qstab:InvalidComponents','Component fields must be real.');
    end
    if ~isequal(size(a),size(b),size(c),size(d))
        error('qstab:ComponentSizeMismatch','Component fields must have identical sizes.');
    end
else
    error('qstab:UnsupportedInput', ...
        'Use a numeric array, a native quaternion array, or a struct with a,b,c,d.');
end
a = double(a); b = double(b); c = double(c); d = double(d);
if any(isnan(a(:)) | isinf(a(:))) || any(isnan(b(:)) | isinf(b(:))) || ...
        any(isnan(c(:)) | isinf(c(:))) || any(isnan(d(:)) | isinf(d(:)))
    error('qstab:NonfiniteInput','Quaternion components must be finite.');
end
end
