function [V,D,info] = eigpairs(A,varargin)
%EIGPAIRS Right quaternion eigenvectors, diagonal representatives and checks.
%   [V,D,INFO] = qstab.eigpairs(A) computes RIGHT eigenpairs A*V=V*D.
%   A may be native quaternion, numeric, or an ordinary a,b,c,d component
%   struct. D is diagonal with canonical representatives a+1i*b, b>=0,
%   retaining algebraic multiplicities. Columns of V have unit Euclidean
%   quaternion norm. This is NOT a left-eigenvalue solver.
%
%   Output='like' (default) returns native quaternion V,D for a native
%   quaternion input, and ordinary a,b,c,d structs otherwise. Explicit
%   Output='components' always uses plain data; Output='quaternion' requests
%   native objects. Numeric/component inputs need only base MATLAB.
%
%   INFO.eigenpairsValid checks pairing, finite outputs and EACH eigenpair
%   residual. INFO.fullBasisObserved is independent: a defective matrix can
%   have valid eigenpairs but no full eigenbasis. INFO.status='rankDeficient'
%   is not a claim of exact defectiveness. Inspect the basis singular values,
%   rank threshold and reciprocal condition. Never infer V is invertible
%   from A*V-V*D alone. No generalized eigenvectors or Jordan chains are made.
%
%   The existing [lambda,INFO]=qstab.eig(A) has NOT changed. Use it when only
%   the classes and its established diagnostics are needed. This separate
%   function adds eigenvectors without overloading that second output.
%   In particular qstab.eigpairs(A) with only ONE output returns V, NOT the
%   eigenvalue vector. The standard representatives are diag(D), or read
%   INFO.standardEigenvalues when using component outputs.
%
%   The numerical core pairs full complex conjugates, then selects columns
%   using quaternionic independence (complex span of z and J*conj(z)).
%   Separate eigenvalue entries are not collapsed when close numerical
%   classes are grouped for vector selection. Spectrum/vector ordering is
%   ascending real part, then ascending nonnegative imaginary part.
%   No call to qeig, no global settings and no private toolbox are used.
%   All decisions are floating point, not interval-certified.
%
%   [V,D,I]=qstab.eigpairs(-eye(4),'Output','components');
%   assert(I.eigenpairsValid && I.fullBasisObserved)
%   X=qstab.chi(-eye(4));Y=qstab.chi(V);T=qstab.chi(D);
%   norm(X*Y-Y*T,'fro')
%
%   See also qstab.eig, qstab.chi, qstab.components.
%
%   Examples (ascending coefficients):
%   [V,D,I] = qstab.eigpairs(-eye(4),'Output','components');
%   assert(I.eigenpairsValid && I.fullBasisObserved);
%   norm(qstab.chi(-eye(4))*qstab.chi(V)-qstab.chi(V)*qstab.chi(D),'fro')
%
%   Interpretation:
%   A single output is V. D contains canonical standard representatives and V is oriented consistently with A*V=V*D. Small residuals do not imply a full basis; check eigenpairsValid and fullBasisObserved separately.
%
%   Reference: docs/reference/eigpairs.html.
%   See also qstab.about.
try
narginchk(1,Inf);nargoutchk(0,3);
ip=inputParser;addParameter(ip,'Output','like',@localText);parse(ip,varargin{:});
output=lower(char(ip.Results.Output));
if ~ismember(output,{'like','components','quaternion'})
    error('qstab:eigpairs:Output','Output must be like, components or quaternion.');
end
[a,b,c,d]=qstab.components(A);
if ~ismatrix(a)||size(a,1)~=size(a,2)
    error('qstab:eigpairs:NotSquare','Right eigenpairs require a two-dimensional square matrix.');
end
[a,b,c,d]=deal(full(a),full(b),full(c),full(d));
[Z,e,info]=qstab.internal.rightpairs(a,b,c,d);
n=size(a,1);u=Z(1:n,:);w=-conj(Z(n+1:2*n,:));
va=real(u);vb=imag(u);vc=real(w);vd=imag(w);
da=diag(real(e));db=diag(imag(e));dc=zeros(n);dd=zeros(n);
native=strcmp(output,'quaternion')||(strcmp(output,'like')&&isa(A,'quaternion'));
if native
    V=quaternion(va,vb,vc,vd);D=quaternion(da,db,dc,dd);
    info.outputType='quaternion';
else
    V=struct('a',va,'b',vb,'c',vc,'d',vd);
    D=struct('a',da,'b',db,'c',dc,'d',dd);
    info.outputType='components';
end
info.standardEigenvalues=e;
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'), throwAsCaller(problem); end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end
function tf=localText(x)
tf=(ischar(x)&&isrow(x)&&~isempty(x))||(isstring(x)&&isscalar(x)&&~ismissing(x));
end
