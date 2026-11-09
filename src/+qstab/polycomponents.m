function [C,info] = polycomponents(P,varargin)
%POLYCOMPONENTS Ascending coefficients of a scalar quaternion polynomial.
%   C = qstab.polycomponents(P) returns a real 4-by-(grade+1) matrix, with
%   columns containing the a,b,c,d parts of coefficients p_0,...,p_grade.
%   P is a numeric/native-quaternion coefficient vector, a 1-by-1-by-(g+1)
%   coefficient tensor, a vector cell array of scalar coefficients, or a
%   component struct a,b,c,d with the same coefficient layout.
%   [C,INFO] = ... also records grade, actualDegree and exact zero padding.
%   ...,'Grade',G pads to grade G; it never silently removes coefficients.
%   All coefficients are ascending. A raw numeric 4-by-N component matrix
%   is deliberately not inferred; use the explicit component struct.
%
%   Examples (ascending coefficients):
%   [C,I] = qstab.polycomponents([2 3 1]);
%   disp(C); disp(I.grade);
%
%   Interpretation:
%   The returned 4-by-(g+1) array is an explicitly documented real representation. A raw 4-by-N input is NOT interpreted as quaternion polynomial coefficients.
%
%   Reference: docs/reference/polycomponents.html.
%   See also qstab.about.

grade = [];
if mod(numel(varargin),2) ~= 0
    error('qstab:InvalidOptions','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = lower(char(varargin{k}));
    switch name
        case 'grade'
            grade = varargin{k+1};
        otherwise
            error('qstab:UnknownOption','Unknown option: %s.',name);
    end
end
if iscell(P)
    if isempty(P) || ~isvector(P)
        error('qstab:InvalidPolynomial','Scalar polynomial coefficient cells must form a nonempty vector.');
    end
    C = zeros(4,numel(P));
    for k = 1:numel(P)
        [a,b,c,d] = qstab.components(P{k});
        if ~isscalar(a)
            error('qstab:NotScalarPolynomial','Each scalar polynomial coefficient must be scalar.');
        end
        C(:,k) = [a;b;c;d];
    end
else
    [a,b,c,d] = qstab.components(P);
    if isempty(a) || ~(isvector(a) || isscalar(a) || ...
            (size(a,1)==1 && size(a,2)==1 && ndims(a)<=3))
        error('qstab:NotScalarPolynomial', ...
            'Use a coefficient vector or a 1-by-1-by-(grade+1) tensor.');
    end
    C = [a(:).';b(:).';c(:).';d(:).'];
end
inputGrade = size(C,2)-1;
if isempty(grade)
    grade = inputGrade;
end
if ~(isnumeric(grade) && isreal(grade) && isscalar(grade) && ...
        ~isnan(grade) && ~isinf(grade) && grade>=0 && grade==fix(grade))
    error('qstab:InvalidGrade','Grade must be a finite nonnegative integer.');
end
if grade<inputGrade
    error('qstab:GradeTooSmall','Grade cannot be smaller than the supplied coefficient layout.');
end
if grade>inputGrade
    C(:,end+1:grade+1) = 0;
end
nz = find(any(C~=0,1),1,'last');
if isempty(nz)
    actualDegree = -Inf;
else
    actualDegree = nz-1;
end
info = struct('grade',double(grade),'inputGrade',inputGrade, ...
    'actualDegree',actualDegree,'isZero',isempty(nz), ...
    'leadingIsZero',all(C(:,end)==0),'constantIsZero',all(C(:,1)==0), ...
    'coefficientOrder','ascending','componentOrder','a,b,c,d', ...
    'assurance','floating-point');
end
