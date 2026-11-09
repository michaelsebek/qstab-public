function cases=qstab_cases_friendly(group)
%QSTAB_CASES_FRIENDLY Distinguish invalid arguments from implementation errors.
% Returns {name,handle} rows for test_qstab; no global settings/files changed.
if nargin==0,group='numeric';end
switch group
case 'numeric'
    cases={'A reversed interval reports the public caller and unchanged identifier',@ordered; ...
           'Bad plot and random-family options point to the user',@plotsAndRandom; ...
           'An arithmetic range failure preserves its actual source and call stack',@numericFault; ...
           'Unknown, foreign, and caused exceptions are not trimmed',@foreignFaults};
case 'native'
    cases={'A genuine MATLAB inputParser error identifies the public caller',@matlabParser; ...
           'The public classifier preserves native MException causes',@nativeCause};
otherwise
    error('qstab:test:Group','Unknown friendly-errors test group.');
end
end

function ordered()
e=capture(@userDefinition);compact(e,'qstab:family:IntervalOrder');
% Test the raw, non-generator path independently of qstab.family metadata.
e=capture(@()qstab.robuststable({2,1},'Uncertainty','interval'));
compact(e,'qstab:familycomponents:IntervalOrder');
[~,s]=qstab.robuststable({1,2},'Uncertainty','interval');
assert(isstruct(s));
end

function userDefinition()
F=qstab.family('interval',2,1); %#ok<NASGU>
end

function plotsAndRandom()
e=capture(@()qstab.plotFamilyValues({[1 1]},1,'AxisUnit','sphere','Projection','slice','Plot',false));
compact(e,'qstab:sample:Projection');
e=capture(@()qstab.randfamily(-1));compact(e,'qstab:randfamily:Degree');
end

function numericFault()
e=capture(@()qstab.plotFamilyValues({1e200},0,1,'Projection','normpoly','Plot',false));
assert(strcmp(e.identifier,'qstab:sample:ArithmeticRange'));
assert(strcmpPath(e.stack(1).file,which('qstab.plotFamilyValues')));
assert(numel(e.stack)>=2,'The original numerical stack disappeared.');
end

function foreignFaults()
% A foreign caller may use the SAME identifier. Neither that frame nor
% subsequent frames in our implementation may disguise the true origin.
e=capture(@foreign); assert(~qstab.internal.isUserInputError(e));
s=struct('identifier','qstab:family:IntervalOrder','stack',e.stack);
s.stack(2)=s.stack(1);s.stack(2).file=which('qstab.family');
assert(~qstab.internal.isUserInputError(s));
s.stack(1).file=which('qstab.family');assert(qstab.internal.isUserInputError(s));
s.cause={'deliberate witness cause'};assert(~qstab.internal.isUserInputError(s));
s.cause={};s.identifier='qstab:DecisionStatus';assert(~qstab.internal.isUserInputError(s));
end

function foreign()
error('qstab:family:IntervalOrder','Independent user error.');
end

function nativeCause()
% Real native MException carrying a cause, with a known identifier.
e=capture(@raiseCaused);assert(numel(e.cause)==1&&strcmp(e.cause{1}.identifier,'qstab:test:Cause'));
assert(~qstab.internal.isUserInputError(e));
end

function raiseCaused()
e=MException('qstab:family:IntervalOrder','Known wording, but an underlying fault.');
cause=MException('qstab:test:Cause','Original diagnostic.');throw(addCause(e,cause));
end

function matlabParser()
parser=inputParser;parser.FunctionName='qstab.familycomponents';addParameter(parser,'Uncertainty','');
reference=capture(@()parse(parser,'BogusExampleOption',true));
assert(~isempty(reference.identifier));
e=capture(@()qstab.familycomponents({1,2},'BogusExampleOption',true));
compact(e,reference.identifier);
end

function compact(e,id)
assert(strcmp(id,e.identifier),'Expected %s, received %s: %s',id,e.identifier,e.message);
assert(~isempty(e.stack)&&strcmpPath(e.stack(1).file,which('qstab_cases_friendly')), ...
    'The first frame is not the public caller.');
end

function tf=strcmpPath(a,b)
a=strrep(a,'\','/');b=strrep(b,'\','/');if ispc,tf=strcmpi(a,b);else,tf=strcmp(a,b);end
end

function e=capture(call)
try,call();catch e,return;end
error('qstab:test:MissingError','Expected an error.');
end
