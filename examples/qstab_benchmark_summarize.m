function [summary,pairs] = qstab_benchmark_summarize(plan,rows,options)
%QSTAB_BENCHMARK_SUMMARIZE Preserve unresolved and failed observations.
% Ratios are supplied only for equal, resolved verdicts in ALL repetitions,
% with all samples long enough. No imputed time limit, zero, or dropped
% failure substitutes for the time of a solved problem.
summary=struct('caseId',{},'method',{},'applicable',{},'status',{}, ...
    'ok',{},'decided',{},'timingAdequate',{},'comparisonEligible',{}, ...
    'medianSeconds',{},'minimumSeconds',{},'maximumSeconds',{}, ...
    'medianCells',{},'medianSubfamilies',{},'measuredCalls',{},'message',{});
for j=1:numel(plan)
    p=plan(j);row=rows{j};S=struct('caseId',p.caseId,'method',p.method, ...
        'applicable',p.applicable,'status','notApplicable','ok',true, ...
        'decided',false,'timingAdequate',false,'comparisonEligible',false, ...
        'medianSeconds',NaN,'minimumSeconds',NaN,'maximumSeconds',NaN, ...
        'medianCells',NaN,'medianSubfamilies',NaN,'measuredCalls',0,'message',p.reason);
    if ~p.applicable,summary(end+1)=S;continue;end %#ok<AGROW>
    samples=row.samples;observations=row.warmups;
    for k=1:numel(samples),observations=[observations;samples{k}.observations(:)];end %#ok<AGROW>
    states=cellfun(@(x)x.status,observations,'UniformOutput',false);
    good=cellfun(@(x)x.ok,observations);
    eligible=cellfun(@(x)x.comparisonEligible,observations);
    complete=numel(samples)==options.Repetitions;
    S.ok=complete&&~isempty(good)&&all(good);
    if isempty(states),S.status='notRun';
    elseif all(strcmp(states,states{1})),S.status=states{1};
    else,S.status='mixed';end
    % Undetermined/stable mixtures can be ordinary budget variation, not
    % false correctness. They are still excluded from the timing ranking.
    if any(strcmp(states,'stable'))&&any(strcmp(states,'notStable')),S.ok=false;end
    if ~isempty(samples)
        times=cellfun(@(x)x.secondsPerCall,samples);
        S.medianSeconds=median(times);S.minimumSeconds=min(times);S.maximumSeconds=max(times);
        counts=cellfun(@(x)x.calls,samples);S.measuredCalls=sum(counts);
        S.timingAdequate=all(cellfun(@(x)x.elapsedSeconds>=options.MinSampleTime,samples));
        cells=[];subs=[];
        for k=1:numel(samples)
            obs=samples{k}.observations;
            for z=1:numel(obs),cells(end+1)=obs{z}.processedCells;subs(end+1)=obs{z}.processedSubfamilies;end %#ok<AGROW>
        end
        finite=cells(isfinite(cells));if ~isempty(finite),S.medianCells=median(finite);end
        finite=subs(isfinite(subs));if ~isempty(finite),S.medianSubfamilies=median(finite);end
    end
    S.decided=S.ok&&any(strcmp(S.status,{'stable','notStable'}))&&all(eligible);
    S.comparisonEligible=S.decided&&S.timingAdequate&&isfinite(S.medianSeconds)&&S.medianSeconds>0;
    if ~S.ok,S.message='Runtime, reference, witness, or completeness issue; inspect raw observations.';
    elseif S.decided&&~S.timingAdequate,S.message='Resolved, but samples are too short for a reported ratio.';
    elseif S.comparisonEligible,S.message='Comparable resolved verdict; timing remains environment-specific.';
    else,S.message='No solved-case comparison: unresolved, boundary, degree loss, or mixed verdicts.';end
    summary(end+1)=S; %#ok<AGROW>
end
pairs=struct('caseId',{},'methodA',{},'methodB',{},'eligible',{},'timeAoverB',{},'reason',{});
for i=1:numel(summary)
    for j=i+1:numel(summary)
        a=summary(i);b=summary(j);
        if ~strcmp(a.caseId,b.caseId)||~a.applicable||~b.applicable,continue;end
        ok=a.comparisonEligible&&b.comparisonEligible&&strcmp(a.status,b.status);
        q=struct('caseId',a.caseId,'methodA',a.method,'methodB',b.method, ...
            'eligible',ok,'timeAoverB',NaN,'reason','Different/unresolved verdict or insufficient timing coverage.');
        if ok,q.timeAoverB=a.medianSeconds/b.medianSeconds;q.reason='Median wall seconds per call, A divided by B.';end
        pairs(end+1)=q; %#ok<AGROW>
    end
end
end
