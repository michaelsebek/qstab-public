function out = run_qb02_extensions(varargin)
%RUN_QB02_EXTENSIONS qB02 reductions and correlated CT/DT zero-exclusion.
% Uses only ordinary component structs, so examples do not depend on a
% private quaternion class. This is an example runner, not a benchmark or
% a proof obtained by sampling. It does not alter the result archive.
% OUT=run_qb02_extensions('MaxTime',30) sets the budget for each example.
ip=inputParser;
addParameter(ip,'MaxTime',30,@(x)isnumeric(x)&&isreal(x)&&isscalar(x)&&x>0&&~isnan(x));
parse(ip,varargin{:});T=ip.Results.MaxTime;
args={'MaxTime',T,'MaxCells',12000,'MaxDepth',40};

% Full leading-coefficient hull is separated, but varies with the member.
V={localP([1,1;1,1;0,0;0,0]),localP([2,2;-1,-1;0,0;0,0])};
[~,out.varyingLead] = qstab.zeroexclude(V,args{:});

% p_theta=(i+theta)*(s+1), theta in [-2,2]. Both coefficients share theta.
% Independent coefficient recombination would introduce unstable members.
F=struct('center',localP([0,0;1,1;0,0;0,0]), ...
    'generators',{{[1,1]}},'bounds',[-2,2]);
[~,out.correlated] = qstab.zeroexclude(F,args{:});

% Full quaternion first-order Schur box and its delay-variable reversal.
L=[-.3,1;-.05,0;-.1,0;-.05,0];U=[-.1,1;.05,0;.1,0;.05,0];
[~,out.schurInterval] = qstab.zeroexclude({localP(L),localP(U)}, ...
    'Uncertainty','interval','Variable','z',args{:});
[~,out.delayInterval] = qstab.zeroexclude({localP(L(:,end:-1:1)),localP(U(:,end:-1:1))}, ...
    'Uncertainty','interval','Variable','zinv',args{:});

% A monic cubic coefficient box, reduced to an eight-parameter cube.
L=[5.98,10.98,5.98,1;-.01,-.01,-.01,0;-.01,-.01,-.01,0;-.01,-.01,-.01,0];
U=[6.02,11.02,6.02,1;.01,.01,.01,0;.01,.01,.01,0;.01,.01,.01,0];
[K,cubeInfo]=qstab.kharitonovcube({localP(L),localP(U)});
[~,out.cube] = qstab.zeroexclude({localP(L),localP(U)}, ...
    'Uncertainty','interval','Reduction','kharitonov-cube',args{:});
% K and cubeInfo are directly inspectable; do not replace K by its vertices.
assert(cubeInfo.parameterCount==8 && isstruct(K));

% Both even and odd real parameters are active: four complete edges, not
% four endpoint polynomials, suffice by the rectangle-opening reduction.
[~,out.pairfree] = qstab.zeroexclude({[1,1],[1.2,1.2]}, ...
    'Uncertainty','interval','Reduction','kharitonov-faces',args{:});

% Monic quadratic: at most two free components of the coefficient of s.
L=[1.9,2.9,1;-.01,-.01,0;0,0,0;0,0,0];U=[2.1,3.1,1;.01,.01,0;0,0,0;0,0,0];
[~,out.quadratic] = qstab.zeroexclude({localP(L),localP(U)}, ...
    'Uncertainty','interval','Reduction','quadratic-faces',args{:});

% Constructed regression: s^2+(1+i+x*j+y*k)*s+2*i. All square edges
% x/y=+-2 are Hurwitz, but x=y=0 is unstable. A whole two-face is required.
L=[0,1,1;2,1,0;0,-2,0;0,-2,0];U=[0,1,1;2,1,0;0,2,0;0,2,0];
[~,out.quadraticInterior] = qstab.zeroexclude({localP(L),localP(U)}, ...
    'Uncertainty','interval','Reduction','quadratic-faces',args{:});

% Deliberately incomplete global subfamily budget.
[~,out.budget] = qstab.zeroexclude({[1,1],[1.2,1.2]}, ...
    'Uncertainty','interval','Reduction','kharitonov-faces', ...
    'MaxSubfamilies',1,args{:});

fprintf('qB02 extensions: complete parameter families, not vertex-only tests\n');
names=fieldnames(out);
for j=1:numel(names)
    R=out.(names{j});fprintf('  %-20s : %-14s %d cells',names{j},R.status,R.processedCells);
    if isfield(R,'reduction')
        fprintf(', %d/%d subfamilies',R.reduction.processedSubfamilies,R.reduction.totalSubfamilies);
    end
    fprintf('\n');
end
fprintf('quadraticInterior=notStable and budget=undetermined are expected.\n');
end
function P=localP(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
