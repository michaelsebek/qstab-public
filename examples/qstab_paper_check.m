function A=qstab_paper_check(kind,actual,reference,varargin)
%QSTAB_PAPER_CHECK Assess a source reference without turning unknown into true.
% Ordinary data only. BOUNDED is containment, not computed equality. Reference
% tolerances are comparison tolerances; they never change any solver tolerance.
ip=inputParser;addParameter(ip,'Tolerance',1e-8,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0);
parse(ip,varargin{:});tol=ip.Results.Tolerance;
A=struct('kind',kind,'status','MISMATCH','ok',false,'matched',false, ...
 'actual',actual,'reference',reference,'comparisonTolerance',tol,'message','');
switch kind
 case 'numeric'
  valid=isnumeric(actual)&&isnumeric(reference)&&isequal(size(actual),size(reference))&& ...
   all(isfinite(actual(:)))&&all(isfinite(reference(:)));
  if valid&&all(abs(actual(:)-reference(:))<=tol*(1+abs(reference(:))))
   A=finish(A,'MATCHED','Source formula reproduced within the stated comparison tolerance.');
  else,A.message='Nonfinite, mismatched shape or inconsistent numeric result.';end
 case 'status'
  if ~ischar(actual)||~ischar(reference),error('qstab:paper:Check','Status comparisons require character vectors.');end
  if strcmp(actual,reference),A=finish(A,'MATCHED','Expected solver status.');
  elseif any(strcmp(actual,{'undetermined','bounded','partial'}))
   A=finish(A,'UNRESOLVED','Solver did not decide the source property.');
  else,A.message='Resolved result does not agree with the source property.';end
 case 'truth'
  if isequal(actual,true)&&isequal(reference,true),A=finish(A,'MATCHED','Explicit finite check passed.');
  else,A.message='Explicit finite check failed.';end
 case 'bracket'
  if ~isstruct(actual)||~isscalar(actual)||~all(isfield(actual,{'lower','upper','status'}))|| ...
    ~isnumeric(reference)||~isscalar(reference)||~isreal(reference)||~isfinite(reference)||reference<0
   error('qstab:paper:Check','A radius bracket and one finite source radius are required.');
  end
  lo=actual.lower;hi=actual.upper;
  if ~isnumeric(lo)||~isnumeric(hi)||~isscalar(lo)||~isscalar(hi)||~isreal(lo)||~isreal(hi)|| ...
       isnan(lo)||isnan(hi)||lo<0||hi<lo
   A.message='Missing or inconsistent radius bounds.';return
  end
  if lo>reference+tol*(1+abs(reference))||hi<reference-tol*(1+abs(reference))
   A.message='Source radius is outside the numerical bracket.';return
  end
  if strcmp(actual.status,'converged')&&isfinite(hi)
   A=finish(A,'MATCHED','Converged numerical bracket contains source radius; this is not exact equality.');
  elseif isfinite(hi)&&any(strcmp(actual.status,{'bounded','undetermined','partial'}))
   A=finish(A,'BOUNDED','Finite bracket contains source radius but is not converged.');
  elseif any(strcmp(actual.status,{'bounded','undetermined','partial'}))
   A=finish(A,'UNRESOLVED','Only a one-sided or empty bound is available.');
  else,A.message='Unexpected radius status for a finite positive source margin.';end
 otherwise,error('qstab:paper:Check','Unknown comparison kind.');
end
end
function A=finish(A,s,msg)
A.status=s;A.ok=~strcmp(s,'MISMATCH');A.matched=strcmp(s,'MATCHED');A.message=msg;
end
