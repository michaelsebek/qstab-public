function [families,Info] = intervalfaces(data,varargin)
%INTERVALFACES Continuous test families for a Hurwitz polynomial interval.
%   F = qstab.family('interval',L,U);
%   families = qstab.intervalfaces(F)
%   [families,Info] = qstab.intervalfaces(F)
%   [families,Info] = qstab.intervalfaces({L,U},'Variable','s')
%
%   F describes independent intervals for all four real components of
%   ascending scalar polynomial coefficients. L and U may be MATLAB
%   quaternion arrays, numeric arrays, or a,b,c,d component structs.
%   families is a column CELL of ordinary qstab.family affine definitions.
%   Each families{k} is directly accepted by qstab.zeroexclude and
%   qstab.robuststable. The bounds in a face are continuous intervals;
%   fixed coordinates have equal endpoints. No new class is introduced.
%
%   The degree/parity selection uses the Hurwitz Kharitonov cube:
%     degree 1   vertices (at most 256);
%     degree 2/3 pair-free two-faces (at most 1536);
%     degree 4   EOO/OOO three-faces and uncovered two-faces
%                (896 for all eight active coordinates: 512 + 384);
%     degree >=5 pair-free three-faces (at most 1024).
%   A pair-free face frees at most one parity in each quaternion component.
%   Zero-width directions are removed from the enumeration. Required
%   lower-dimensional families are kept when coordinates are inactive.
%   A constant interval returns its complete family.
%
%   Info.descriptors{k} gives bounds, freeParameters, fixedParameters,
%   dimension and parityType for families{k}. Info.faceCountsByDimension
%   is indexed by dimension+1. Info.cube gives the endpoint construction,
%   active coordinates and the original leading-coefficient condition.
%
%   ENUMERATION IS NOT A STABILITY TEST. The face theorem requires the
%   original family to preserve degree, and every selected continuous
%   family must be tested. Use the combined call below to enforce the
%   degree condition and global work budgets:
%     [tf,R] = qstab.zeroexclude(F,'Reduction','interval-faces', ...
%                               'MaxSubfamilies',1536);
%   Alternatively use qstab.robuststable(F,'Method','zeroexclusion', ...
%                               'Reduction','interval-faces', ...
%                               'MaxSubfamilies',1536).
%   The default global MaxSubfamilies=1024 can leave a full degree-2/3
%   selection unfinished. Incomplete testing returns NaN, never stable.
%   The legacy kharitonov-faces and specialized monic quadratic-faces
%   selections remain available and unchanged.
%
%   Variable is 's' only (default). The family metadata is checked for
%   conflicts. The independent even/odd cube does not apply to Schur or
%   delay intervals; z and zinv are rejected by this enumeration function.
%
%   Example:
%     F = qstab.family('interval',[1.9 2.9 1],[2.1 3.1 1]);
%     [families,I] = qstab.intervalfaces(F);
%     I.totalSubfamilies
%     [tf,R] = qstab.zeroexclude(families{1});
%
%   See also qstab.family, qstab.kharitonovcube, qstab.zeroexclude,
%            qstab.robuststable.
try
    if qstab.internal.isFamilyDefinition(data)
        varargin=qstab.internal.familyOptions(data,varargin, ...
            {'Variable',data.variable,'Grade',data.grade},{'Variable'});
    end
    ip=inputParser;ip.FunctionName='qstab.intervalfaces';
    addParameter(ip,'Variable','s',@(x)(ischar(x)&&isrow(x)) || ...
        (isstring(x)&&isscalar(x)));
    parse(ip,varargin{:});v=lower(strtrim(char(ip.Results.Variable)));
    if ~ismember(v,{'s','hurwitz','continuous'})
        error('qstab:intervalfaces:Variable', ...
            'The degree/parity Kharitonov face selection applies to s only.');
    end
    F=qstab.familycomponents(data,'Uncertainty','interval');
    [A,tasks,Info]=qstab.internal.intervalfaceplan(F);
    if ~Info.cube.representable
        error('qstab:intervalfaces:CoefficientRange', ...
            'Endpoint differences exceed the representable range.');
    end
    generators=cell(1,8);
    for j=1:8,generators{j}=localPolynomial(A.generators(:,:,j));end
    template=qstab.family('affine',localPolynomial(A.base),generators, ...
        A.bounds,'Variable','s','Grade',F.degree);
    families=cell(numel(tasks),1);
    for j=1:numel(tasks)
        families{j}=template;families{j}.bounds=tasks{j}.bounds;
    end
catch problem
    if qstab.internal.isUserInputError(problem)
        if isa(problem,'MException'),throwAsCaller(problem);end
        problem.stack=dbstack(1,'-completenames');
    end
    rethrow(problem)
end
end

function P=localPolynomial(C)
P=struct('a',reshape(C(1,:),1,1,[]),'b',reshape(C(2,:),1,1,[]), ...
    'c',reshape(C(3,:),1,1,[]),'d',reshape(C(4,:),1,1,[]));
end
