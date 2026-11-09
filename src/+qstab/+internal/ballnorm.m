function y = ballnorm(x,p)
%BALLNORM Scaled finite-dimensional vector p-norm (internal).
a=abs(x(:));
if isempty(a),y=0;return;end
s=max(a);
if s==0||isinf(s),y=s;return;end
if isinf(p),y=s;else,y=s*sum((a/s).^p)^(1/p);end
end
