function cases = qstab_cases_nominal()
%QSTAB_CASES_NOMINAL Nominal algebra, spectral and root-free stability cases.
%   N-by-2 cell array {name,@testFunction}. These tests use numeric arrays
%   and ordinary component structs; native quaternion-class integration is
%   tested separately by the main runner in MATLAB. Octave execution is an
%   auxiliary check and must not be reported as native MATLAB validation.

cases={ ...
    'Numeric and component arrays preserve their shape and Hamilton parts',@localComponents; ...
    'Invalid component arrays and ambiguous raw component matrices are rejected',@localInvalidInputs; ...
    'Ascending scalar vector, tensor and coefficient-cell layouts agree',@localLayouts; ...
    'Complex adjoint preserves noncommuting matrix products and adjoints',@localEmbedding; ...
    'Real norm polynomial uses central polynomial products',@localNormPolynomial; ...
    'Nonmonic quaternion companion preserves the left coefficient convention',@localNonmonicCompanion; ...
    'Right eigenvalue representatives preserve real doubled multiplicities',@localEigenRealMultiplicity; ...
    'Genuinely quaternion triangular right spectrum and eigen residuals',@localEigenQuaternion; ...
    'The algebraic class of s-j is distinct from point evaluation at i',@localAlgebraicVsEvaluation; ...
    'Repeated real scalar zeros retain their algebraic multiplicity',@localRootMultiplicity; ...
    'Noncommuting polynomial factors retain their algebraic zero classes',@localNoncommutingFactors; ...
    'Declared grade preserves zero and infinity in algebraic zero lists',@localInfinity; ...
    'Routh resolves stable and unstable polynomials without computing roots',@localRouth; ...
    'Exact Hurwitz boundary witnesses are separated from the unresolved band',@localRouthBoundary; ...
    'Schur recursion resolves stable and unstable polynomials without roots',@localSchur; ...
    'Exact unit roots are distinguished from a near-unit-circle tolerance band',@localSchurBoundary; ...
    'Delay reversal and grade padding preserve forward zero modes',@localDelay; ...
    'Explicit spectral tests agree with root-free nominal decisions',@localSpectralAgreement; ...
    'Continuous matrices use a Lyapunov certificate for nonnormal dynamics',@localContinuousMatrix; ...
    'Discrete matrices use Stein certificates and preserve the boundary band',@localDiscreteMatrix; ...
    'Nonzero constants, zero polynomials and singular constants are distinct',@localConstants; ...
    'Experimental nominal diagonal polynomial matrix: continuous time',@localPolynomialMatrixCT; ...
    'Experimental nominal diagonal polynomial matrix: delay time',@localPolynomialMatrixDelay ...
};
end

function localComponents()
A=[1+2i,3-4i;5,6i]; [a,b,c,d]=qstab.components(A);
assert(isequal(a,real(A)) && isequal(b,imag(A)) && isequal(c,zeros(2)) && isequal(d,zeros(2)));
S=localQ(reshape(1:6,2,3),ones(2,3),2*ones(2,3),-ones(2,3));
[a,b,c,d]=qstab.components(S);
assert(isequal(size(a),[2,3]) && isequal(a,S.a) && isequal(b,S.b) && isequal(c,S.c) && isequal(d,S.d));
end

function localInvalidInputs()
localMustFail(@() qstab.components(struct('a',1,'b',[0,0],'c',0,'d',0)));
localMustFail(@() qstab.components(localQ(1,1i,0,0)));
localMustFail(@() qstab.components(NaN));
localMustFail(@() qstab.polycomponents(zeros(4,5)));
localMustFail(@() qstab.polycomponents([1,2,3],'Grade',1));
[~,R] = qstab.isstable(ones(2)); assert(strcmp(R.status,'assumptionFailed'));
[~,R] = qstab.isstable([1,1],'Type','matrix'); assert(strcmp(R.status,'assumptionFailed'));
[~,R] = qstab.isstable([1,1,0]); assert(strcmp(R.status,'assumptionFailed'));
end

function localLayouts()
S=localQ([1,2,3],[4,5,6],[7,8,9],[-1,-2,-3]);
T=localQ(reshape(S.a,1,1,[]),reshape(S.b,1,1,[]), ...
    reshape(S.c,1,1,[]),reshape(S.d,1,1,[]));
V=cell(1,3);
for k=1:3, V{k}=localQ(S.a(k),S.b(k),S.c(k),S.d(k)); end
[C,I]=qstab.polycomponents(S);
assert(isequal(C,qstab.polycomponents(T)) && isequal(C,qstab.polycomponents(V)));
assert(I.grade==2 && I.actualDegree==2 && strcmp(I.coefficientOrder,'ascending'));
[C,I]=qstab.polycomponents(S,'Grade',4);
assert(I.grade==4 && I.actualDegree==2 && isequal(C(:,4:5),zeros(4,2)));
end

function localEmbedding()
A=cat(3,[1,2;0,-1],[0,1;2,0],[2,0;1,3],[1,-1;0,2]);
B=cat(3,[2,0;1,1],[1,2;0,1],[0,1;-1,0],[2,0;1,-2]);
AB=localMatrixProduct(A,B); BA=localMatrixProduct(B,A);
XA=qstab.chi(localTensorQ(A)); XB=qstab.chi(localTensorQ(B));
assert(norm(XA*XB-qstab.chi(localTensorQ(AB)),'fro')<1e-12);
assert(norm(AB(:)-BA(:))>1,'The fixture must be noncommutative.');
AH=cat(3,A(:,:,1).',-A(:,:,2).',-A(:,:,3).',-A(:,:,4).');
assert(norm(XA'-qstab.chi(localTensorQ(AH)),'fro')<1e-12);
end

function localNormPolynomial()
C=[1,2,1;2,-1,0;-1,3,0;1,0,0];
P=localRows(C); [N,I]=qstab.normpoly(P);
expected=zeros(1,5);
for k=1:4, expected=expected+conv(C(k,:),C(k,:)); end
assert(isequal(N,expected) && I.grade==4 && I.actualDegree==4);
conjugate=C; conjugate(2:4,:)=-conjugate(2:4,:);
product=localPolynomialProduct(C,conjugate);
assert(norm(product(1,:)-N,Inf)<1e-12 && norm(product(2:4,:),'fro')<1e-12);
end

function localNonmonicCompanion()
% i+j*s has its right evaluation zero -k. Its algebraic class is [i].
P=localQ([0,0],[1,0],[0,1],[0,0]);
M=qstab.companion(P); [A,B,I]=qstab.companion(P);
expected=qstab.chi(localQ(0,0,0,-1));
assert(norm(M-expected,'fro')<1e-12 && norm(B\A-M,'fro')<1e-12);
assert(I.grade==1 && I.quaternionOrder==1);
[z,J]=qstab.roots(P); assert(J.isValid && abs(z-1i)<1e-10);
end

function localEigenRealMultiplicity()
[z,I]=qstab.eig(diag([-2,-2,1]));
assert(I.isValid && numel(z)==3 && numel(I.fullSpectrum)==6);
assert(norm(z-[-2;-2;1],Inf)<1e-12);
assert(isequal(I.classMultiplicities,[2;1]) && sum(I.classMultiplicities)==3);
assert(I.pairingValid && I.eigenResidual<1e-12 && all(I.pairingResidual<1e-12));
end

function localEigenQuaternion()
A=localQ([-1,0;0,-2],[1,0;0,0],[0,0;0,1],[0,2;0,0]);
[z,I]=qstab.eig(A);
assert(I.isValid && numel(z)==2 && all(imag(z)>=0));
localAssertSpectrum(z,[-2+1i;-1+1i],1e-10);
assert(I.eigenResidual<1e-12 && I.embeddingResidual==0);
end

function localAlgebraicVsEvaluation()
P=localQ([0,1],[0,0],[-1,0],[0,0]);
[z,I]=qstab.roots(P);
assert(I.isValid && numel(z)==1 && abs(z-1i)<1e-10);
evaluationAtRepresentative=[0;imag(z);-1;0];
assert(norm(evaluationAtRepresentative)>1,'A class representative was incorrectly treated as an evaluation zero.');
assert(I.eigenResidual<1e-12 && all(I.pairingResidual<1e-10));
end

function localRootMultiplicity()
[z,I]=qstab.roots([1,2,1],'Tolerance',1e-6);
assert(I.isValid && numel(z)==2 && numel(I.fullSpectrum)==4);
assert(all(abs(z+1)<1e-6) && sum(I.classMultiplicities)==2);
assert(numel(I.classMultiplicities)==1 && I.classMultiplicities(1)==2);
end

function localNoncommutingFactors()
u=[-1;1;0;0]; v=[-2;0;1;0];
C=localPolynomialProduct([-u,[1;0;0;0]],[-v,[1;0;0;0]]);
[z,I]=qstab.roots(localRows(C));
assert(I.isValid); localAssertSpectrum(z,[-2+1i;-1+1i],1e-9);
end

function localInfinity()
[z,I]=qstab.roots([1,0]);
assert(I.isValid && numel(z)==1 && isinf(z) && I.infiniteCount==1 && I.grade==1);
[z,I]=qstab.roots(1,'Grade',3);
assert(I.isValid && numel(z)==3 && all(isinf(z)) && I.infiniteCount==3);
[z,I]=qstab.roots([0,1,0]);
assert(I.isValid && numel(z)==2 && any(abs(z)==0) && sum(isinf(z))==1);
end

function localRouth()
[~,S] = qstab.isstable([6,11,6,1]); [~,U] = qstab.isstable([5,1,1,1]);
[~,Q] = qstab.isstable(localQ([1,1],[0,0],[2,0],[0,0]));
localStatus(S,'stable'); localStatus(U,'notStable'); localStatus(Q,'stable');
assert(~S.spectrumComputed && ~U.spectrumComputed && ~Q.spectrumComputed);
assert(strcmp(S.method,'routh') && all(S.diagnostics.test.firstColumn>0));
end

function localRouthBoundary()
% Test signed identities on ORIGINAL input, before numerical scaling.
% The scaled subnormal imaginary part disappears; the original has no
% exact zero constant. A near-boundary spectral residual is insufficient.
a=realmin*eps; P=localQ([0 realmax],[a 0],[0 0],[0 0]);
[tf,U]=qstab.isstable(localQ([a realmax],[0 0],[0 0],[0 0]));
assert(isnan(tf) && ~isfield(U.diagnostics,'boundaryWitness'));
[tf,R]=qstab.isstable(P);
assert(tf==0 && strcmp(R.diagnostics.boundaryWitness.kind,'orthogonalLinearCoefficients'));
[tf,R]=qstab.isstable([0 1 0]);assert(isnan(tf) && strcmp(R.status,'assumptionFailed'));
% A structural scalar diagonal block survives arbitrary upper right
% entries. Removing triangularity removes the proof.
P={ [0 3;0 2],eye(2) };
[tf,R]=qstab.isstable(P,'Type','polynomialMatrix');
assert(tf==0 && R.diagnostics.boundaryWitness.diagonalIndex==1);
P{1}(2,1)=-1;
[tf,R]=qstab.isstable(P,'Type','polynomialMatrix');
assert(tf==1 && ~isfield(R.diagnostics,'boundaryWitness'));
P={ [0,1], [1,0,1], ...
    localQ([1e-12,1],[1,0],[0,0],[0,0]), ...
    localQ([-1e-12,1],[1,0],[0,0],[0,0]) };
for k=1:numel(P)
    [tf,R] = qstab.isstable(P{k},'Tolerance',1e-10);
    if k<=2, localStatus(R,'notStable'); assert(tf==0 && R.diagnostics.boundaryWitness.exactForSuppliedComponents);
    else, localStatus(R,'undetermined'); assert(isnan(tf) && ~isfield(R.diagnostics,'boundaryWitness'));end
    assert(~R.spectrumComputed);
end
end

function localSchur()
[~,S] = qstab.isstable([0.12,-0.7,1],'Variable','z');
[~,U] = qstab.isstable([1.2,1],'Variable','z');
[~,Q] = qstab.isstable(localQ([0.2,1],[0.2,0],[0.2,0],[0.2,0]),'Variable','z');
localStatus(S,'stable'); localStatus(U,'notStable'); localStatus(Q,'stable');
assert(strcmp(S.method,'schur') && ~S.spectrumComputed && ~U.spectrumComputed);
end

function localSchurBoundary()
for a=[1,1-1e-12,1+1e-12]
    [tf,R] = qstab.isstable([-a,1],'Variable','z','Tolerance',1e-10);
    if a==1, localStatus(R,'notStable');assert(tf==0);
    else,localStatus(R,'undetermined');assert(isnan(tf));end
    assert(~R.spectrumComputed);
end
end

function localDelay()
[~,D] = qstab.isstable([1,-0.5],'Variable','zinv');
[~,Z] = qstab.isstable([-0.5,1],'Variable','z');
localStatus(D,'stable'); localStatus(Z,'stable');
assert(strcmp(D.variable,'zinv') && strcmp(D.diagnostics.testVariable,'z'));
[z,I]=qstab.roots([1,-0.5],'Variable','zinv');
assert(I.isValid && abs(z-2)<1e-12 && strcmp(I.variable,'zinv'));
[~,D] = qstab.isstable([1,0],'Variable','zinv'); localStatus(D,'stable'); assert(D.grade==1);
[~,D] = qstab.isstable(1,'Variable','zinv','Grade',3); localStatus(D,'stable'); assert(D.grade==3);
[~,D] = qstab.isstable([0,1],'Variable','zinv'); localStatus(D,'assumptionFailed');
end

function localSpectralAgreement()
inputs={ [2,3,1],'s'; [3,-1,1],'s'; [0.12,-0.7,1],'z'; [1.2,1],'z'; [1,-0.5],'zinv' };
for k=1:size(inputs,1)
    [~,A] = qstab.isstable(inputs{k,1},'Variable',inputs{k,2});
    [~,B] = qstab.isstable(inputs{k,1},'Variable',inputs{k,2},'Method','spectral');
    assert(strcmp(A.status,B.status) && ~A.spectrumComputed && B.spectrumComputed);
    assert(B.diagnostics.spectrum.isValid && B.diagnostics.spectrum.eigenResidual<1e-10);
end
end

function localContinuousMatrix()
[~,S] = qstab.isstable([-1,4;0,-1],'Type','matrix','Domain','continuous');
localStatus(S,'stable'); assert(~S.spectrumComputed);
assert(strcmp(S.diagnostics.certificate.solver,'kronecker'), ...
    'The nonnormal fixture must exercise the Lyapunov solve.');
assert(S.diagnostics.certificate.equationResidual<1e-10);
[~,U] = qstab.isstable(diag([-1,2]),'Type','matrix'); localStatus(U,'notStable');
assert(~isempty(U.diagnostics.certificate.negativeDirection));
[~,N] = qstab.isstable(-1e-12,'Type','matrix','Tolerance',1e-10); localStatus(N,'undetermined');
[~,N] = qstab.isstable(1e-12,'Type','matrix','Tolerance',1e-10); localStatus(N,'undetermined');
end

function localDiscreteMatrix()
[~,S] = qstab.isstable([0.5,1;0,0.5],'Type','matrix','Domain','discrete');
localStatus(S,'stable'); assert(~S.spectrumComputed);
assert(strcmp(S.diagnostics.certificate.solver,'kronecker'));
[~,U] = qstab.isstable(1.2,'Type','matrix','Domain','discrete'); localStatus(U,'notStable');
for a=[1,1-1e-12,1+1e-12]
    [tf,R] = qstab.isstable(a,'Type','matrix','Domain','discrete','Tolerance',1e-10);
    if a==1, localStatus(R,'notStable'); assert(tf==0 && R.diagnostics.boundaryWitness.diagonalIndex==1);
    else, localStatus(R,'undetermined'); assert(isnan(tf));end
end
end

function localConstants()
P=localQ(1,2,3,4); [~,R] = qstab.isstable(P); localStatus(R,'stable');
[z,I]=qstab.roots(P); assert(isempty(z) && I.isValid && ~I.spectrumComputed);
[~,R] = qstab.isstable(0); localStatus(R,'assumptionFailed');
[z,I]=qstab.roots(0); assert(isempty(z) && ~I.isValid && strcmp(I.status,'assumptionFailed'));
localMustFail(@() qstab.companion(0));
S={zeros(2)};
[~,R] = qstab.isstable(S,'Type','polynomialMatrix'); localStatus(R,'assumptionFailed');
[z,I]=qstab.roots(S); assert(isempty(z) && ~I.isValid);
end

function localPolynomialMatrixCT()
P={diag([1,2]),eye(2)};
[z,I]=qstab.roots(P); assert(I.isValid); localAssertSpectrum(z,[-2;-1],1e-12);
[~,S] = qstab.isstable(P,'Type','polynomialMatrix'); localStatus(S,'stable'); assert(~S.spectrumComputed);
[~,R] = qstab.isstable({diag([2,3])},'Type','polynomialMatrix'); localStatus(R,'stable');
[~,R] = qstab.isstable({diag([-1,2]),eye(2)},'Type','polynomialMatrix','Method','spectral');
localStatus(R,'notStable');
end

function localPolynomialMatrixDelay()
P={eye(2),-diag([0.2,0.5])};
[~,S] = qstab.isstable(P,'Type','polynomialMatrix','Variable','zinv'); localStatus(S,'stable');
[~,Z] = qstab.isstable({-diag([0.2,0.5]),eye(2)},'Type','polynomialMatrix','Variable','z');
localStatus(Z,'stable');
[~,R] = qstab.isstable(P,'Type','polynomialMatrix','Variable','zinv','Method','spectral');
localStatus(R,'stable'); localAssertSpectrum(R.diagnostics.representatives,[0.2;0.5],1e-12);
end

function localStatus(R,expected)
assert(strcmp(R.status,expected),['Expected ',expected,'; got ',R.status,': ',R.message]);
assert(strcmp(R.assurance,'floating-point'));
if strcmp(expected,'stable'), assert(isequal(R.isStable,true));
elseif strcmp(expected,'notStable'), assert(isequal(R.isStable,false));
else, assert(isempty(R.isStable)); end
end

function localAssertSpectrum(actual,expected,tol)
assert(numel(actual)==numel(expected));
remaining=actual(:);
for k=1:numel(expected)
    [distance,j]=min(abs(remaining-expected(k)));
    assert(distance<=tol*max(1,abs(expected(k)))); remaining(j)=[];
end
end

function localMustFail(f)
failed=false;
try, value=f(); %#ok<NASGU>
catch, failed=true;
end
assert(failed,'The invalid input was silently accepted.');
end

function S=localQ(a,b,c,d)
S=struct('a',a,'b',b,'c',c,'d',d);
end

function S=localRows(C)
S=localQ(C(1,:),C(2,:),C(3,:),C(4,:));
end

function S=localTensorQ(C)
S=localQ(C(:,:,1),C(:,:,2),C(:,:,3),C(:,:,4));
end

function C=localMatrixProduct(A,B)
C=zeros(size(A,1),size(B,2),4);
for i=1:size(A,1)
    for j=1:size(B,2)
        v=zeros(4,1);
        for k=1:size(A,2)
            v=v+localMultiply(reshape(A(i,k,:),4,1),reshape(B(k,j,:),4,1));
        end
        C(i,j,:)=reshape(v,1,1,4);
    end
end
end

function C=localPolynomialProduct(A,B)
C=zeros(4,size(A,2)+size(B,2)-1);
for i=1:size(A,2)
    for j=1:size(B,2)
        C(:,i+j-1)=C(:,i+j-1)+localMultiply(A(:,i),B(:,j));
    end
end
end

function z=localMultiply(x,y)
z=[x(1)*y(1)-dot(x(2:4),y(2:4)); ...
    x(1)*y(2:4)+y(1)*x(2:4)+cross(x(2:4),y(2:4))];
end
