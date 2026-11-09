function plan = qstab_benchmark_plan(cases,options)
%QSTAB_BENCHMARK_PLAN Explicit applicable/nonapplicable comparison rows.
% No solver is executed. Do not silently try an unsupported reduction.
methods={'lyapunov','direct','four-generator','kharitonov-cube', ...
    'kharitonov-faces','quadratic-faces'};
if ~isempty(options.Methods),methods=options.Methods;end
plan=struct('caseIndex',{},'caseId',{},'method',{},'applicable',{}, ...
    'reason',{},'solverOptions',{},'comparisonGroup',{});
for c=1:numel(cases)
    F=qstab.familycomponents(cases(c).data,'Uncertainty',cases(c).uncertainty);
    for j=1:numel(methods)
        m=methods{j};ok=true;reason='';engine='zeroexclusion';reduction='none';
        switch m
            case 'lyapunov',engine='lyapunov';
            case 'direct'
            case 'four-generator'
                reduction=m;ok=strcmp(F.kind,'polytope');
                reason='Four-generator reduction requires a polynomial polytope.';
            case {'kharitonov-cube','kharitonov-faces','quadratic-faces'}
                reduction=m;ok=strcmp(F.kind,'interval')&&strcmp(cases(c).variable,'s');
                reason='This interval reduction is continuous-time and requires Cartesian bounds.';
                if ok&&strcmp(m,'quadratic-faces')
                    ok=F.degree==2 && isequal(F.intervalLower(:,end),[1;0;0;0]) && ...
                        isequal(F.intervalUpper(:,end),[1;0;0;0]);
                    reason='Quadratic faces require exactly monic degree-two input.';
                end
            otherwise
                error('qstab:benchmark:Method','Unknown benchmark method %s.',m);
        end
        if ok,reason='Applicable';end
        opts=struct('Method',engine,'Reduction',reduction, ...
            'Uncertainty',cases(c).uncertainty,'Variable',cases(c).variable, ...
            'Tolerance',options.Tolerance,'MaxCells',options.MaxCells, ...
            'MaxDepth',options.MaxDepth,'MaxTime',options.MaxTime, ...
            'MaxSubfamilies',options.MaxSubfamilies,'StoreCertificates',false);
        row=struct('caseIndex',c,'caseId',cases(c).id,'method',m, ...
            'applicable',ok,'reason',reason,'solverOptions',opts, ...
            'comparisonGroup',cases(c).id);
        plan(end+1)=row; %#ok<AGROW>
    end
end
end
