function report=verify_qstab_oneparameter_identities(varargin)
%VERIFY_QSTAB_ONEPARAMETER_IDENTITIES Recompute 42 exact article identities.
%   R=verify_qstab_oneparameter_identities uses Symbolic Math Toolbox to
%   verify the E1 feedback model, calibration margins and the E5 variant.
%   R.ok is true only after all 42 checks were freshly executed.
%   'Verbose',false suppresses the final count. No Python is launched.
%   These identities supplement the generic radiuscertificate/raycertificate
%   proofs; they do not replace those global and first-contact certificates.
%   See also run_qstab_oneparameter.
ip=inputParser;ip.PartialMatching=false;addParameter(ip,'Verbose',true,@(x)islogical(x)&&isscalar(x));parse(ip,varargin{:});
try,sym(1);catch,error('qstab:exact:SymbolicRequired','These exact identities require Symbolic Math Toolbox.');end
s=sym('qstab_identity_s','real');g=sym('qstab_identity_gamma','real');
w=sym('qstab_identity_w','positive');Z=sym(0);O=sym(1);I=sym(1i);
p0=[s*s+2*s+1;Z;O;Z];K=[Z;s;Z;O];plant=[s*s+2*s+1;-s;O;-O];D=[Z;O;Z;O];
checks=struct('name',{},'passed',{});started=tic;
check('plant + controller equals E1 nominal polynomial',zero(plant+K-p0));
check('nominal plant Hurwitz',positiveAll(hurwitz(normp(plant),s)));
check('nominal closed loop Hurwitz',positiveAll(hurwitz(normp(p0),s)));
x=[sym('qstab_identity_x0','real');sym('qstab_identity_x1','real');sym('qstab_identity_x2','real');sym('qstab_identity_x3','real')];
one=[O;Z;Z;Z];a=one+x;
check('real calibration error = L(t)',zero(L(a)-eye(4)-L(x)));
check('induced norm identity',zero(L(x).'*L(x)-sum(x.^2)*eye(4)));
V=[0 0;-1 0;0 1;0 -1];
for k=1:2,M=L(sym(V(:,k)));check(sprintf('mechanical skew coupling %d',k-1),zero(M.'+M));end
for side={'L','R'}
 if strcmp(side{1},'L'),actual=plant+mul(a,K);intended=p0+mul(x,K);G=L(a)*L(K);
 else,actual=plant+mul(K,a);intended=p0+mul(K,x);G=L(K)*L(a);end
 check([side{1} ' hypercomplex feedback identity'],zero(actual-intended));
 check([side{1} ' real feedback matrix identity'],zero(L(plant)+G-L(intended)));
end
P0=subs(chi(p0),s,I*w);P1=subs(chi(K),s,I*w);
for side={'L','R'}
 if strcmp(side{1},'L')
  rr=2+w^4/(1+w)^2;A0=P0'*P0;A1=P1'*P1;tr=4*w^5/(w+1)^2;
 else
  rr=(w^4+2*w^2-4*w+2)/(w+1)^2;A0=P0*P0';A1=P1*P1';tr=4*w*(w*w+2)^2/(w+1)^2;
 end
 B=simplify(A0-rr*A1);
 check([side{1} ' frequency formula determinant'],zero(det(B)));
 check([side{1} ' frequency formula PSD trace'],zero(trace(B)-tr)&&truth(tr>0)&&zero(subs(tr,w,0)));
 check([side{1} ' positive radius expression'],zero(subs(rr,w,0)-2));
 if strcmp(side{1},'L')
  check('L global minimum certificate',zero(rr-2-w^4/(1+w)^2));
 else
  b=sqrt(sym(3))-1;rs=14-8*sqrt(sym(3));
  check('R exact nonnegative factor certificate',zero((1+w)^2*(rr-rs)-(w-b)^2*(w*w+2*b*w+2*sqrt(sym(3)))));
  check('R minimum attained',zero(subs(rr,w,b)-rs));
 end
end
minorsR=sym([]);
for side={'L','R'}
 if strcmp(side{1},'L'),p=p0+g*mul(D,K);first=O;expected=8*(g-3)*(g-1)*(g*g-2*g+2);lam=[Z;Z;Z;Z];
 else,p=p0+g*mul(K,D);first=2-sqrt(sym(3));expected=24*(g*g-4*g+1)*(g*g-2*g+2);lam=[Z;Z;1-sqrt(sym(3));Z];end
 H=hurwitz(normp(p),s);h=H(3);
 check([side{1} ' guardian factorization'],zero(h-expected));
 if strcmp(side{1},'L')
  firstProof=zero(h-8*(g-first)*(g-3)*((g-1)^2+1))&&truth(first>0)&&truth(first<3)&&truth((g-1)^2+1>0);
 else
  later=2+sqrt(sym(3));firstProof=zero(h-24*(g-first)*(g-later)*((g-1)^2+1))&&truth(first>0)&&truth(first<later)&&truth((g-1)^2+1>0);
 end
 check([side{1} ' first positive guardian root'],firstProof);
 pc=subs(p,g,first);if strcmp(side{1},'L'),evaluation='right';else,evaluation='left';end
 check([side{1} ' original polynomial contact'],zero(evalq(pc,s,lam,evaluation)));
 vals=subs(H,g,sym(1)/3);
 if strcmp(side{1},'L'),classification=positiveAll(vals);else,classification=truth(vals(3)<0)&&all(arrayfun(@(k)~zero(vals(k)),1:4));minorsR=vals;end
 check([side{1} ' at 1/3 correct stability classification'],classification);
end
check('chosen witness inside half-unit ball',truth(sym(2)/9<sym(1)/4));
check('actuator radius above half-unit ball',truth(sym(2)>sym(1)/4));
check('sensor radius below half-unit ball',truth(14-8*sqrt(sym(3))<sym(1)/4));
check('right third minor exact',zero(minorsR(3)+sym(208)/27));
realp=p0+g*K;H=hurwitz(normp(realp),s);
check('real gain norm exact',zero(normp(realp)-((s+1)^4+1+g*g*(s*s+1))));
check('real scale reduced Hurwitz determinant constant',zero(H(3)-48));
pert=[-s-1;s;-s-1;Z];p=p0+g*pert;H=hurwitz(normp(p),s);h=H(3);cubic=3*g^3-9*g^2+10*g-6;
check('E5 variant guardian',zero(h+8*(g-1)^2*cubic));
check('E5 variant cubic increasing',zero(diff(cubic,g)-(9*(g-1)^2+1))&&truth(9*(g-1)^2+1>0));
check('E5 variant next root > 1',truth(subs(cubic,g,1)<0));
check('E5 variant next root rational bracket',truth(subs(cubic,g,sym(1747)/1000)<0)&&truth(subs(cubic,g,sym(1748)/1000)>0));
check('E5 variant exact stable after contact',positiveAll(subs(H,g,sym(3)/2)));
check('E5 variant contact factor',zero(subs(p,g,1)-[s*s+s;s;-s;Z]));
check('E5 variant zero on original polynomial',zero(evalq(subs(p,g,1),s,[Z;Z;Z;Z],'right')));
q=[Z;Z;O;Z];r=[Z;O;-O;Z];
check('E5 variant coefficients have no common slice',~zero(mul(q,r)-mul(r,q)));
report=struct('schema','qstab.oneparameter.exact-identities.1','ok',numel(checks)==42&&all([checks.passed]), ...
 'status','EXECUTED','computedThisRun',true,'assurance','exact-symbolic', ...
 'passed',sum([checks.passed]),'failed',sum(~[checks.passed]),'checks',checks,'elapsedSeconds',toc(started));
if report.ok,report.status='PASS';else,report.status='FAIL';end
if ip.Results.Verbose,fprintf('One-parameter exact identities: %d passed, %d failed.\n',report.passed,report.failed);end
 function check(name,condition)
  checks(end+1)=struct('name',name,'passed',logical(condition)); %#ok<AGROW>
 end
end
function H=hurwitz(N,s)
% Every polynomial in this specific identity suite has nominal grade two;
% its real norm has degree four. No numerical degree inference is used.
C=sym(zeros(1,5));for k=0:4,C(k+1)=simplify(subs(diff(N,s,k),s,0)/factorial(k));end
M=sym(zeros(4));for i=1:4,for j=1:4,k=4-2*j+i;if k>=0&&k<=4,M(i,j)=C(k+1);end;end;end
H=sym(zeros(4,1));for k=1:4,H(k)=simplify(det(M(1:k,1:k)));end
end
function n=normp(p),n=expand(sum(p.^2));end
function z=mul(a,b)
z=expand([a(1)*b(1)-sum(a(2:4).*b(2:4));a(1)*b(2)+a(2)*b(1)+a(3)*b(4)-a(4)*b(3); ...
 a(1)*b(3)-a(2)*b(4)+a(3)*b(1)+a(4)*b(2);a(1)*b(4)+a(2)*b(3)-a(3)*b(2)+a(4)*b(1)]);
end
function M=L(q),M=sym(zeros(4));E=sym(eye(4));for k=1:4,M(:,k)=mul(q,E(:,k));end;end
function M=chi(q),I=sym(1i);M=[q(1)+I*q(2),q(3)+I*q(4);-q(3)+I*q(4),q(1)-I*q(2)];end
function y=evalq(p,s,z,side)
y=sym(zeros(4,1));for k=2:-1:0,c=subs(diff(p,s,k),s,0)/factorial(k);if strcmp(side,'right'),y=mul(y,z)+c;else,y=mul(z,y)+c;end;end
end
function y=zero(x),x=simplify(x);y=isequal(x,sym(zeros(size(x))));end
function y=truth(x),if islogical(x),y=all(x(:));else,y=all(isAlways(x));end;y=logical(y);end
function y=positiveAll(x),y=true;for k=1:numel(x),y=y&&truth(x(k)>0);end;end
