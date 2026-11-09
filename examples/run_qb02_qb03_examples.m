function out = run_qb02_qb03_examples(varargin)
%RUN_QB02_QB03_EXAMPLES Compact reproducible cases for qB02/qB03 methods.
%   OUT=run_qb02_qb03_examples calls the specialized zero-exclusion solver
%   on Hurwitz coefficient boxes, affine polynomial polytopes, and Schur
%   polytopes. This runner uses native MATLAB quaternion coefficients and
%   is not executed automatically by test_qstab.
%
%   OUT=run_qb02_qb03_examples('MaxTime',T) sets the per-case time limit.
%   See also qstab.zeroexclude, qstab.robuststable, test_qstab.

ip=inputParser;
addParameter(ip,'MaxTime',30,@(x)isnumeric(x)&&isscalar(x)&&isreal(x)&&x>0);
parse(ip,varargin{:});
T=ip.Results.MaxTime;

% qB02: full quaternion coefficient box with Re(q0) >= 1, p(s)=s+q0.
L=localQ([1,1],[-.05,0],[-.05,0],[-.05,0]);
U=localQ([1.2,1],[.05,0],[.05,0],[.05,0]);
[~,out.interval] = qstab.zeroexclude({L,U},'Uncertainty','interval', ...
  'Variable','s','MaxTime',T);

% qB03: coarse Bernstein control hull includes zero at w=0 but the
% underlying segment is Hurwitz: p_t=s+1+2*t*j, -1 <= t <= 1.
V={localQ([1,1],[0,0],[-2,0],[0,0]), ...
   localQ([1,1],[0,0],[ 2,0],[0,0])};
[~,out.bernstein] = qstab.zeroexclude(V,'Variable','s','MaxTime',T, ...
    'StoreCertificates',true);

% qB03: four vertices of a tetrahedron. Its boundary is Hurwitz by the
% exact |y|^2 > 1/2 criterion, while the centre y=0 is unstable.
sgn=[1,1,1;1,-1,-1;-1,1,-1;-1,-1,1];
W=cell(1,4);
for j=1:4
    y=2*sgn(j,:);
    W{j}=localQ([-1,1,1],[y(1),3*y(1),0], ...
          [y(2),3*y(2),0],[y(3),3*y(3),0]);
end
[~,out.tetrahedron] = qstab.zeroexclude(W,'Variable','s','MaxTime',T);

% Schur polynomial segment, plus the same family written in z^-1.
SZ={localQ([-.3,1],[0,0],[-.2,0],[0,0]), ...
    localQ([-.3,1],[0,0],[ .2,0],[0,0])};
SD={localQ([1,-.3],[0,0],[0,-.2],[0,0]), ...
    localQ([1,-.3],[0,0],[0, .2],[0,0])};
[~,out.schur] = qstab.zeroexclude(SZ,'Variable','z','MaxTime',T);
[~,out.delay] = qstab.zeroexclude(SD,'Variable','zinv','MaxTime',T);

names=fieldnames(out);
fprintf('qB02/qB03 specialized zero-exclusion examples\n');
for k=1:numel(names)
    r=out.(names{k});
    fprintf('  %-14s : %-15s (%d cells, %.3f s)\n', ...
      names{k},r.status,r.processedCells,r.elapsedTime);
end
end

function P=localQ(a,b,c,d)
% Native quaternion values, in ascending coefficient order.
P=quaternion(reshape(a,1,1,[]),reshape(b,1,1,[]), ...
    reshape(c,1,1,[]),reshape(d,1,1,[]));
end
