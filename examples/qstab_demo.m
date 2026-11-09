function out = qstab_demo(varargin)
%QSTAB_DEMO Small examples of quaternion nominal and robust stability.
%   OUT = qstab_demo uses the native MATLAB quaternion class and opens
%   spectral/value-set figures. The scalar examples use ASCENDING
%   coefficient tensors P(:,:,k+1)=P_k; matrices are ordinary quaternion
%   arrays. No custom quaternion or polynomial class is needed.
%
%   OUT = qstab_demo('Plots',false) runs without graphics.
%   OUT = qstab_demo('UseQuaternion',false,'Plots',false) uses ordinary
%   a,b,c,d component structs. This path can be used for an auxiliary
%   Octave check; it does not validate MATLAB's native quaternion class.
%
%   Demonstrations:
%     p(s)=s^2+3*s+2+j: nominal root-free Hurwitz test and algebraic zeros;
%     an upper-triangular quaternion matrix: Lyapunov test and right spectrum;
%     z-(0.4+0.2*j) and 1-(0.4+0.2*j)*z^-1: reciprocal conventions;
%     scalar first-degree quaternion intervals and a three-vertex polytope;
%     optional exact 2-D/3-D similarity-class views and sampled value sets.
%
%   The robustness examples have simple independent analytic explanations:
%   every member is s+a0 with Re(a0)>=1, hence its zero class is Hurwitz.
%   Numerical robustness outputs still carry their ordinary floating-point
%   assurance; the plotted value samples are not a robustness certificate.

opts = struct('UseQuaternion',true,'Plots',true);
if mod(numel(varargin),2) ~= 0
    error('qstab:demo:Options','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = varargin{k};
    if isstring(name) && isscalar(name), name = char(name); end
    if ~ischar(name) || size(name,1) ~= 1
        error('qstab:demo:Options','Option names must be character rows or scalar strings.');
    end
    switch lower(name)
        case 'usequaternion', opts.UseQuaternion = localLogical(varargin{k+1});
        case 'plots', opts.Plots = localLogical(varargin{k+1});
        otherwise, error('qstab:demo:Options','Unknown option "%s".',name);
    end
end
if opts.UseQuaternion
    try
        quaternion(0,0,0,0);
    catch
        error('qstab:demo:QuaternionUnavailable', ...
            ['The native MATLAB quaternion constructor is unavailable. ', ...
             'For a component-struct demonstration, set UseQuaternion to false.']);
    end
end

P = localPolynomial([2,3,1],[0,0,0],[1,0,0],[0,0,0],opts);
A = localArray([-1,.2;0,-2],[1,0;0,0],[0,0;0,1],[0,.3;0,0],opts);
Pz = localPolynomial([-.4,1],[0,0],[-.2,0],[0,0],opts);
Pd = localPolynomial([1,-.4],[0,0],[0,-.2],[0,0],opts);

out = struct();
out.options = opts;
out.polynomial = P;
out.matrix = A;
out.forwardPolynomial = Pz;
out.delayPolynomial = Pd;
[~,out.nominalPolynomial] = qstab.isstable(P,'Variable','s');
[~,out.nominalMatrix] = qstab.isstable(A,'Type','matrix','Domain','continuous');
[~,out.forward] = qstab.isstable(Pz,'Variable','z');
[~,out.delay] = qstab.isstable(Pd,'Variable','zinv');
[out.polynomialZeros,out.polynomialZeroInfo] = qstab.roots(P);
[out.matrixEigenvalues,out.matrixEigenInfo] = qstab.eig(A);
[out.forwardZeros,out.forwardZeroInfo] = qstab.roots(Pz,'Variable','z');
[out.delayZeros,out.delayZeroInfo] = qstab.roots(Pd,'Variable','zinv');

L = localPolynomial([1,1],[-.05,0],[-.05,0],[-.05,0],opts);
U = localPolynomial([1.2,1],[.05,0],[.05,0],[.05,0],opts);
[~,out.robustInterval] = qstab.robuststable({L,U},'Uncertainty','interval', ...
    'Variable','s','MaxCells',64,'MaxTime',10);
V = {localPolynomial([1,1],[.3,0],[0,0],[0,0],opts), ...
     localPolynomial([1,1],[0,0],[.3,0],[0,0],opts), ...
     localPolynomial([1,1],[0,0],[0,0],[.3,0],opts)};
[~,out.robustPolytope] = qstab.robuststable(V,'Uncertainty','polytope', ...
    'Variable','s','MaxCells',64,'MaxTime',10);

% The weights select the three vertices, their center, and edge midpoints.
% These are samples of the SAME family just passed to robuststable.
W = [eye(3);ones(1,3)/3;.5,.5,0;0,.5,.5;.5,0,.5];
[~,out.valueSamples] = qstab.valueset(V,1i*[0,.5,1,1.5,2],'Weights',W);
out.graphics = struct();
if opts.Plots
    [~,out.graphics.spectrum2D] = qstab.plotSpectrum(out.polynomialZeros, ...
        'Variable','s','Dimension',2);
    [~,out.graphics.spectrum3D] = qstab.plotSpectrum(out.polynomialZeros, ...
        'Variable','s','Dimension',3);
    [~,out.graphics.values] = qstab.plotValueSet(out.valueSamples);
end

if opts.UseQuaternion, kind = 'native MATLAB quaternion'; else, kind = 'ordinary component structs'; end
fprintf('qstab demonstration (%s)\n',kind);
localPrint('p(s)=s^2+3s+2+j',out.nominalPolynomial);
localPrint('Quaternion state matrix',out.nominalMatrix);
localPrint('Forward variable z',out.forward);
localPrint('Delay variable z^-1',out.delay);
localPrint('Quaternion coefficient interval',out.robustInterval);
localPrint('Quaternion polynomial polytope',out.robustPolytope);
fprintf('The forward and delay root representatives stay in their input variables.\n');
fprintf('Value-set points are samples; they do not prove zero exclusion.\n');
end

function P = localPolynomial(a,b,c,d,opts)
P = localArray(reshape(a,1,1,[]),reshape(b,1,1,[]), ...
    reshape(c,1,1,[]),reshape(d,1,1,[]),opts);
end

function A = localArray(a,b,c,d,opts)
if opts.UseQuaternion
    A = quaternion(a,b,c,d);
else
    A = struct('a',a,'b',b,'c',c,'d',d);
end
end

function x = localLogical(x)
if ~isscalar(x) || ~(islogical(x) || isnumeric(x)) || ~isreal(x) || ...
        ~isfinite(x) || ~ismember(x,[0,1])
    error('qstab:demo:Options','UseQuaternion and Plots must be logical scalars.');
end
x = logical(x);
end

function localPrint(label,R)
fprintf('  %-34s : %s\n',label,R.status);
end
