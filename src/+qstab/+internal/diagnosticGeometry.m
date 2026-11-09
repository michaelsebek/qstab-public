function varargout = diagnosticGeometry(action,varargin)
%DIAGNOSTICGEOMETRY Small ordinary-double geometry for diagnostic displays.
% No function here calls a stability solver or supplies a stability verdict.
switch action
    case 'hull'
        varargout{1}=localHull(varargin{1});
    case 'distance'
        v=varargin{1};v=v(:);d=0;
        if isempty(v) || any(~isfinite(real(v))|~isfinite(imag(v)))
            varargout{1}=NaN;return
        end
        scale=max(abs(v));
        if scale==0,varargout{1}=0;return;end
        w=v/scale;
        angles=sort(mod(angle(w),2*pi));
        gaps=diff([angles;angles(1)+2*pi]);
        if any(w==0) || max(gaps)<=pi+64*eps*pi
            varargout{1}=0;return
        end
        h=localHull(w);d=min(abs(h));
        if numel(h)>1
            for k=1:numel(h)
                a=h(k);b=h(mod(k,numel(h))+1);u=b-a;
                if abs(u)>0
                    t=max(0,min(1,-real(conj(u)*a)/abs(u)^2));
                    d=min(d,abs(a+t*u));
                end
            end
        end
        varargout{1}=scale*d;
    case 'inflate'
        v=localHull(varargin{1});r=varargin{2};sides=64;
        if numel(varargin)>2,sides=varargin{3};end
        if ~isnumeric(r)||~isreal(r)||~isscalar(r)||~isfinite(r)||r<0
            error('qstab:diagnostics:Radius','Inflation radius must be finite and nonnegative.');
        end
        if ~isnumeric(sides)||~isscalar(sides)||~isreal(sides)||~isfinite(sides)||sides<8||sides~=fix(sides)
            error('qstab:diagnostics:Sides','Sides must be an integer of at least 8.');
        end
        if r==0,varargout{1}=v;return;end
        % CIRCUMSCRIBED polygon, not an inscribed rendering of a disk.
        disk=(r/cos(pi/sides))*exp(2i*pi*(0:sides-1)/sides);
        cloud=bsxfun(@plus,v(:),disk);
        varargout{1}=localHull(cloud(:));
    case 'rectangle'
        b=varargin{1};
        if ~isnumeric(b)||~isreal(b)||~isequal(size(b),[2 2])||any(~isfinite(b(:)))||any(b(:,1)>b(:,2))
            error('qstab:diagnostics:Rectangle','Bounds must be finite [realLo realHi; imagLo imagHi].');
        end
        varargout{1}=localHull([b(1,1)+1i*b(2,1);b(1,2)+1i*b(2,1); ...
            b(1,2)+1i*b(2,2);b(1,1)+1i*b(2,2)]);
    otherwise
        error('qstab:diagnostics:GeometryAction','Unknown diagnostic geometry action.');
end
end
function h=localHull(v)
if ~isnumeric(v)||isempty(v)||~isvector(v)||any(~isfinite(real(v(:)))|~isfinite(imag(v(:))))
    error('qstab:diagnostics:Hull','Supply a nonempty finite numeric vector.');
end
v=double(v(:));xy=unique([real(v),imag(v)],'rows');
if size(xy,1)<=2,h=complex(xy(:,1),xy(:,2));return;end
scale=max(abs(xy(:)));
if scale==0,h=0;return;end
q=xy/scale;n=size(q,1);lowerHull=zeros(n,1);upperHull=zeros(n,1);a=0;b=0;
for k=1:n
    while a>=2 && localCross(q(lowerHull(a-1),:),q(lowerHull(a),:),q(k,:))<=0,a=a-1;end
    a=a+1;lowerHull(a)=k;
end
for k=n:-1:1
    while b>=2 && localCross(q(upperHull(b-1),:),q(upperHull(b),:),q(k,:))<=0,b=b-1;end
    b=b+1;upperHull(b)=k;
end
idx=[lowerHull(1:a-1);upperHull(1:b-1)];h=complex(xy(idx,1),xy(idx,2));
end
function z=localCross(a,b,c)
u=b-a;v=c-a;z=u(1)*v(2)-u(2)*v(1);
end
