function [A,tasks,info] = intervalfaceplan(F)
%INTERVALFACEPLAN Degree/parity selection for a Hurwitz coefficient interval.
% Enumeration only. The caller must check uniform degree before using the
% face theorem as a stability test. Legacy face selections are unchanged.
[A,cube]=qstab.internal.cubefamily(F);
n=F.degree;active=cube.activeParameters;
components=unique(cube.parameterComponents(active));
if n==0
    % Constants need only the nonzero-family condition. Retain the complete
    % family for inspection, rather than suggesting that its vertices prove it.
    states=zeros(1,8);states(active)=2;
    selection='constant-family';dimension=numel(active);
    theorem='Nonzero constant interval; the full family is retained';
elseif n==1
    states=localStates(active,cube.parameterComponents,0);
    selection='vertices';dimension=0;
    theorem='Hurwitz degree-one interval theorem: coefficient-box vertices';
elseif n<=3
    dimension=min(2,numel(components));
    states=localStates(active,cube.parameterComponents,dimension);
    selection='pair-free-two-faces';
    theorem='Hurwitz interval degree two/three theorem: pair-free two-faces';
else
    dimension=min(3,numel(components));
    states=localStates(active,cube.parameterComponents,dimension);
    if n==4 && dimension==3
        retained=zeros(0,8);descended=zeros(0,8);
        for j=1:size(states,1)
            row=states(j,:);free=find(row==2);
            if sum(free>4)>=2
                retained(end+1,:)=row; %#ok<AGROW>
            else
                descended=[descended;localBoundary(row)]; %#ok<AGROW>
            end
        end
        descended=unique(descended,'rows');
        % A retained three-face already tests all six of its boundary
        % two-faces. No EO/OO face is discarded unless this containment is
        % actually established; inactive odd directions may prevent it.
        covered=zeros(0,8);
        for j=1:size(retained,1)
            covered=[covered;localBoundary(retained(j,:))]; %#ok<AGROW>
        end
        if ~isempty(covered) && ~isempty(descended)
            descended=descended(~ismember(descended,covered,'rows'),:);
        end
        states=[retained;descended];
    end
    if n==4
        selection='quartic-parity-faces';
        theorem='Hurwitz quartic interval theorem: EOO/OOO three-faces and uncovered two-faces';
    else
        selection='pair-free-three-faces';
        theorem='Hurwitz interval theorem: pair-free faces of dimension at most three';
    end
end
states=unique(states,'rows');
tasks=cell(size(states,1),1);counts=zeros(1,9);
for j=1:size(states,1)
    row=states(j,:);free=find(row==2);fixed=setdiff(1:8,free);
    B=[row(:),row(:)];B(free,1)=0;B(free,2)=1;
    parity=[repmat('E',1,sum(free<=4)),repmat('O',1,sum(free>4))];
    tasks{j}=struct('bounds',B,'freeParameters',free, ...
        'fixedParameters',fixed,'dimension',numel(free), ...
        'parityType',parity,'faceIndex',j);
    counts(numel(free)+1)=counts(numel(free)+1)+1;
end
if isempty(tasks),actualDimension=0;else,actualDimension=find(counts,1,'last')-1;end
info=struct('schema','qstab.intervalfaces.1','variable','s', ...
    'degree',n,'selection',selection,'theorem',theorem, ...
    'totalSubfamilies',numel(tasks),'parameterDimensionBound',actualDimension, ...
    'theoremDimensionBound',dimension,'faceCountsByDimension',counts, ...
    'activeParameters',active,'activeParameterCount',numel(active), ...
    'parameterNames',{cube.parameterNames},'cube',cube, ...
    'descriptors',{tasks},'vertexOnlyTest',n==1, ...
    'requiresDegreePreservation',true,'stabilityTestPerformed',false, ...
    'certificateComplete',false,'validated',false, ...
    'assurance','Combinatorial face enumeration; no stability decision');
end

function rows=localStates(active,labels,k)
% State 0/1 fixes an endpoint; state 2 frees that coordinate. Inactive
% coordinates are always fixed at zero, avoiding duplicate zero generators.
components=unique(labels(active));
if k==0,sets=zeros(1,0);else,sets=nchoosek(components,k);end
rows=zeros(0,8);
for row=1:size(sets,1)
    choices={zeros(1,0)};
    for component=sets(row,:)
        candidate=active(labels(active)==component);next={};
        for j=1:numel(choices)
            for index=candidate
                next{end+1}=[choices{j},index]; %#ok<AGROW>
            end
        end
        choices=next;
    end
    for j=1:numel(choices)
        free=choices{j};fixed=setdiff(active,free);
        block=zeros(2^numel(fixed),8);block(:,free)=2;
        for ell=1:numel(fixed)
            block(:,fixed(ell))=bitget((0:size(block,1)-1)',ell);
        end
        rows=[rows;block]; %#ok<AGROW>
    end
end
end

function rows=localBoundary(row)
free=find(row==2);rows=repmat(row,2*numel(free),1);
for j=1:numel(free)
    rows(2*j-1,free(j))=0;rows(2*j,free(j))=1;
end
end
