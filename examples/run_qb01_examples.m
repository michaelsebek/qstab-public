function results = run_qb01_examples(varargin)
%RUN_QB01_EXAMPLES Selected exact examples from the corrected qB01 notes.
%   RESULTS=run_qb01_examples uses native MATLAB quaternion inputs. It
%   reports the example separating left/right radii and the real-ray cubic
%   with loss and re-entry. It is not a runner for every example in qB01.
%   Options: UseQuaternion (true), Plots (true), OutputDir ('qb01_results').
%   The output MAT file contains only ordinary structs and numeric arrays.
%   Numerical estimates and analytical reference values are labelled apart.

    ip=inputParser;
    addParameter(ip,'UseQuaternion',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'Plots',true,@(x)islogical(x)&&isscalar(x));
    addParameter(ip,'OutputDir',fullfile(pwd,'qb01_results'),...
        @(x)ischar(x)||(isstring(x)&&isscalar(x)));
    parse(ip,varargin{:}); o=ip.Results;
    p0=makepoly([1 2 1;0 0 0;1 0 0;0 0 0],o.UseQuaternion);
    p1=makepoly([0 0;0 1;0 0;1 0],o.UseQuaternion);
    u=makepoly([0;1;0;1],o.UseQuaternion);
    results=struct('source','qB01 corrected internal Revision 2, 7 October 2026',...
        'runtime',version,'isOctave',exist('OCTAVE_VERSION','builtin')~=0,...
        'nativeQuaternionInputs',o.UseQuaternion,'timestamp',datestr(now,30));
    results.exact=struct('leftRadius',sqrt(2),'rightRadius',sqrt(2)*(2-sqrt(3)),...
        'rightFrequency',sqrt(3)-1,'leftRay',1,'rightRay',2-sqrt(3),...
        'reentryInterval',[-13/48,1/4]);
    [~,results.rayLeft] = qstab.raymargin(p0,p1,u,'Side','left');
    [~,results.rayRight] = qstab.raymargin(p0,p1,u,'Side','right');
    [~,results.radiusLeft] = qstab.stabilityradius(p0,p1,'Side','left',...
        'AbsTol',5e-4,'RelTol',5e-4,'MaxTime',30,'MaxCells',4096);
    [~,results.radiusRight] = qstab.stabilityradius(p0,p1,'Side','right',...
        'AbsTol',5e-4,'RelTol',5e-4,'MaxTime',30,'MaxCells',4096);
    results.omega=linspace(0,3,121);
    [results.frequencyLeft,results.frequencyLeftInfo]=qstab.freqmargin(p0,p1,results.omega,'Side','left');
    [results.frequencyRight,results.frequencyRightInfo]=qstab.freqmargin(p0,p1,results.omega,'Side','right');
    r0=[13/16 1 1 1]; r1=[3 1 1];
    [~,results.reentryRay] = qstab.raymargin(r0,r1,1);
    results.reentryVertices={r0,r0+[r1 0]};
    [~,results.reentryRobust] = qstab.robuststable(results.reentryVertices,'Uncertainty','polytope');
    fprintf('qB01 selected examples (%s)\n',version);
    fprintf('  Left global radius:  [%.10g, %.10g], analytical %.10g (%s)\n',...
        results.radiusLeft.lower,results.radiusLeft.upper,results.exact.leftRadius,results.radiusLeft.status);
    fprintf('  Right global radius: [%.10g, %.10g], analytical %.10g (%s)\n',...
        results.radiusRight.lower,results.radiusRight.upper,results.exact.rightRadius,results.radiusRight.status);
    fprintf('  Direction u=i+k: left alpha=%.10g, right alpha=%.10g\n',...
        results.rayLeft.positive,results.rayRight.positive);
    fprintf('  Re-entry cubic first interval: (%.10g, %.10g)\n',results.reentryRay.interval);
    fprintf('  Segment between two stable endpoints: %s\n',results.reentryRobust.status);
    output=char(o.OutputDir);
    if ~exist(output,'dir'), mkdir(output); end
    save(fullfile(output,'qb01_results.mat'),'results','-v7');
    if o.Plots
        f=figure('Color','w','Name','qB01 frequency margins');
        plot(results.omega,results.frequencyLeft,'LineWidth',1.5); hold on;
        plot(results.omega,results.frequencyRight,'LineWidth',1.5);
        grid on; xlabel('\omega'); ylabel('\rho(\omega)');
        legend({'p_0+t p_1','p_0+p_1 t'},'Location','best');
        title('Minimum parameter norm at each frequency');
        saveas(f,fullfile(output,'frequency_margins.fig'));
        print(f,fullfile(output,'frequency_margins.png'),'-dpng','-r160');
    end
end
function P=makepoly(C,useQuaternion)
    a=reshape(C(1,:),1,1,[]); b=reshape(C(2,:),1,1,[]);
    c=reshape(C(3,:),1,1,[]); d=reshape(C(4,:),1,1,[]);
    if useQuaternion
        P=quaternion(a,b,c,d);
    else
        P=struct('a',a,'b',b,'c',c,'d',d);
    end
end
