function cases = qstab_cases_values()
%QSTAB_CASES_VALUES Scalar norm-value geometry and input-contract checks.
%   CASES is an N-by-2 cell array {name,@zeroArgumentTest}. These groups
%   use numeric arrays and ordinary component structs, without graphics.
%   Running them in Octave is an auxiliary check, not native MATLAB
%   quaternion-class validation.

cases = { ...
    'Norm zero of s-j identifies a class without asserting evaluation at i', @localClassVsEvaluation; ...
    'Complex norm values use squares rather than absolute squares', @localComplexSquares; ...
    'Shared convex weights preserve coefficient correlations before the norm map', @localSharedWeights; ...
    'A norm value of a convex combination can lie outside the vertex-value hull', @localNonconvexImage; ...
    'Malformed arguments, weight dimensions and inadmissible convex weights fail', @localValidation; ...
    'Ascending vectors, coefficient tensors and coefficient cells have explicit meanings', @localLayouts ...
};
end

function localClassVsEvaluation()
% p=s-j has N_p=s^2+1, but p_r(i)=i-j is nonzero.
P = localQ([0,1],[0,0],[-1,0],[0,0]);
[~,S] = qstab.valueset(P,1i);
assert(abs(S.points)<1e-14, 'The algebraic class represented by i was lost.');
evaluationComponents = [0;1;-1;0];
assert(norm(evaluationComponents)==sqrt(2), 'The reference evaluation must be nonzero.');
assert(strcmp(S.kind,'singleton') && ~S.isCertificate, ...
    'One fixed-polynomial value must not be described as a robustness certificate.');
end

function localComplexSquares()
% For the real polynomial p=1+s, N_p=(1+s)^2. At i it equals 2i,
% whereas abs(1+i)^2=2. The complex value cannot be replaced by a norm.
z = [0;1i;1+2i;-0.5+0.3i];
[~,S] = qstab.valueset([1,1],z);
expected = (1+z.').^2;
assert(isequal(size(S.points),[1,numel(z)]) && isequal(size(S.lambda),[1,numel(z)]));
assert(max(abs(S.points-expected))<1e-13);
assert(abs(S.points(2)-2i)<1e-14 && abs(S.points(2)-abs(1+1i)^2)>1, ...
    'Squared absolute values replaced the central norm polynomial.');
assert(all(S.finiteMask(:)) && S.sampleCount==1 && S.vertexCount==1);
end

function localSharedWeights()
% Vertices (1+i)*(1+s) and (1-i)*(1+s) share the SAME coefficient
% parameter. At weights (w,1-w), N=(1+(2*w-1)^2)*(1+s)^2.
Pplus = localQ([1,1],[1,1],[0,0],[0,0]);
Pminus = localQ([1,1],[-1,-1],[0,0],[0,0]);
W = [1,0;.75,.25;.5,.5;.25,.75;0,1];
z = [0,1i,1+1i];
[~,S] = qstab.valueset({Pplus,Pminus},z,'Weights',W);
t = 2*W(:,1)-1;
expected = (1+t.^2)*((1+z).^2);
assert(isequal(size(S.points),[5,3]) && max(abs(S.points(:)-expected(:)))<1e-12);
assert(isequal(S.parameters,W) && S.weightAdjustment==0, ...
    'The supplied shared parameter coordinates were not retained.');
assert(abs(S.points(3,2)-2i)<1e-14, ...
    'The midpoint must be formed before applying the quadratic norm map.');
assert(strcmp(S.kind,'samples') && ~S.isCertificate);
end

function localNonconvexImage()
% At lambda=0 the vertices s+1 and s-1 both have norm value 1.
% Their polynomial midpoint is s and has norm value 0. Therefore the
% image of the coefficient segment is not the convex hull {1} of its
% endpoint norm values. A vertex-value hull would miss a boundary zero.
V = {[1,1],[-1,1]};
[~,vertices] = qstab.valueset(V,0);
[~,middle] = qstab.valueset(V,0,'Weights',[.5,.5]);
assert(isequal(vertices.points,[1;1]) && middle.points==0);
assert(abs(middle.points-vertices.points(1))==1, ...
    'The polynomial midpoint was replaced by interpolation of vertex values.');
assert(~vertices.isCertificate && ~middle.isCertificate, ...
    'Finite samples cannot certify zero exclusion.');
assert(strcmp(vertices.kind,'samples') && vertices.sampleCount==2);
end

function localValidation()
V = {[1,1],[2,1]};
localMustFail(@() qstab.valueset([1,1],ones(2)));
localMustFail(@() qstab.valueset([1,1],NaN));
localMustFail(@() qstab.valueset([1,1],Inf));
localMustFail(@() qstab.valueset(V,1i,'Weights',[1,0,0]));
localMustFail(@() qstab.valueset(V,1i,'Weights',[-.1,1.1]));
localMustFail(@() qstab.valueset(V,1i,'Weights',[.4,.4]));
localMustFail(@() qstab.valueset(V,1i,'Weights',[.5,.5+1i]));
localMustFail(@() qstab.valueset(V,1i,'Weights',[NaN,1]));
localMustFail(@() qstab.valueset(eye(2),1i));
localMustFail(@() qstab.valueset([1,Inf],1i));
localMustFail(@() qstab.valueset([1,1],1i,'UnknownOption',1));

[~,S] = qstab.valueset(V,0,'Weights',[.25;.75]);
assert(isequal(S.parameters,[.25,.75]) && abs(S.points-1.75^2)<1e-14);
[~,S] = qstab.valueset(V,0,'Weights',[-eps,1+eps]);
assert(isequal(S.parameters,[0,1]) && S.points==4 && S.weightAdjustment>0, ...
    'Roundoff-sized simplex corrections must be explicit in the output.');
end

function localLayouts()
z = [0,1i,2+1i];
P = [1,2,1];
tensor = reshape(P,1,1,[]);
coefficients = {1,2,1};
[~,A] = qstab.valueset(P,z);
[~,B] = qstab.valueset(tensor,z);
[~,C] = qstab.valueset(coefficients,z);
assert(isequal(A.points,B.points) && isequal(A.points,C.points));
assert(strcmp(C.kind,'singleton') && C.grade==2 && C.vertexCount==1, ...
    'A scalar coefficient cell without Weights must mean one polynomial.');

V = {reshape([1,1],1,1,[]),reshape([2,0,1],1,1,[])};
[~,S] = qstab.valueset(V,z,'Weights',[.5,.5]);
expected = (1.5+.5*z+.5*z.^2).^2;
assert(max(abs(S.points-expected))<1e-12 && S.grade==2);
assert(isequal(S.vertexGrades,[1;2]), 'High-degree padding lost vertex metadata.');

[~,singlePolynomial] = qstab.valueset({1,2},0);
[~,constantFamily] = qstab.valueset({1,2},0,'Weights',[.5,.5]);
[~,constantVertices] = qstab.valueset({1,2},0,'Weights',[]);
assert(singlePolynomial.points==1 && strcmp(singlePolynomial.kind,'singleton'));
assert(constantFamily.points==2.25 && strcmp(constantFamily.kind,'samples'));
assert(isequal(constantVertices.points,[1;4]));
end

function P = localQ(a,b,c,d)
P = struct('a',a,'b',b,'c',c,'d',d);
end

function localMustFail(f)
failed = false;
try
    f();
catch ME
    failed = strncmp(ME.identifier,'qstab:',6);
end
assert(failed,'An invalid input did not raise a qstab input error.');
end
