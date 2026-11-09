function tests = qstab_cases_graphics()
%QSTAB_CASES_GRAPHICS MATLAB hidden-figure integration and geometry checks.
    tests={...
        'Two-dimensional spectrum and separate infinite representative',@twoD;...
        'Three-dimensional class slice preserves norm and real part',@threeD;...
        'Plot functions preserve an existing axes and hold state',@existingAxes};
end
function twoD()
    [~,H] = qstab.plotSpectrum([-1+2i;.3+.1i;Inf],'Visible','off','Variable','s');
    cleanup=onCleanup(@()close(H.figure)); %#ok<NASGU>
    assert(H.infinityCount==1&&~H.isCertificate);
    assert(numel(H.canonicalLambda)==2);
    assert(max(abs(sort(imag(H.canonicalLambda))-[.1;2]))<1e-12);
    assert(numel(H.markers)==3&&all(isgraphics(H.markers,'line')));
    assert(isequal(get(H.markers(1),'XData'),-1)&&isequal(get(H.markers(2),'XData'),.3));
    assert(isempty(get(H.markers(3),'XData'))&&isempty(get(H.markers(3),'YData')), ...
        'An absent boundary group must retain one line handle with empty data.');
    % A row input with only boundary classes leaves both other groups empty.
    [~,B] = qstab.plotSpectrum([1i,2i],'Parent',H.axes,'Variable','s');
    assert(numel(B.markers)==3&&all(isgraphics(B.markers,'line')));
    assert(isempty(get(B.markers(1),'XData'))&&isempty(get(B.markers(2),'XData')));
    assert(isequal(get(B.markers(3),'YData'),[1,2]));
    % Constant/grade bookkeeping can produce an entirely infinite spectrum.
    [~,E] = qstab.plotSpectrum(Inf,'Parent',H.axes,'Variable','s');
    assert(isempty(E.canonicalLambda)&&E.infinityCount==1&&isgraphics(E.annotation));
    assert(numel(E.markers)==3&&all(isgraphics(E.markers,'line')));
    for k=1:3, assert(isempty(get(E.markers(k),'XData'))); end
end
function threeD()
    [~,H] = qstab.plotSpectrum(-.5+2i,'Visible','off','Variable','z','Dimension',3);
    cleanup=onCleanup(@()close(H.figure)); %#ok<NASGU>
    assert(numel(H.classes)==1);
    x=get(H.classes(1),'XData'); y=get(H.classes(1),'YData'); z=get(H.classes(1),'ZData');
    assert(max(abs(x+.5))<1e-12&&max(abs(x.^2+y.^2+z.^2-4.25))<1e-12);
    assert(strcmp(H.classification{1},'unstable'));
    assert(numel(H.markers)==3&&all(isgraphics(H.markers,'line')));
    assert(isequal(get(H.markers(2),'XData'),-.5)&&isequal(get(H.markers(2),'YData'),2));
    assert(isequal(get(H.markers(2),'ZData'),0));
    for k=[1,3]
        assert(isempty(get(H.markers(k),'XData'))&&isempty(get(H.markers(k),'YData'))&& ...
            isempty(get(H.markers(k),'ZData')), ...
            'An absent 3-D group must retain a scalar line handle with empty data.');
    end
end
function existingAxes()
    f=figure('Visible','off'); cleanup=onCleanup(@()close(f)); %#ok<NASGU>
    ax=axes('Parent',f); old=plot(ax,[0 1],[0 1]); set(ax,'NextPlot','replace');
    [returnedAx,P] = qstab.plotSpectrum(-1+1i,'Parent',ax);
    assert(isequal(returnedAx,ax)&&isequal(P.axes,returnedAx));
    assert(isgraphics(old)&&strcmp(get(ax,'NextPlot'),'replace'));
    assert(numel(P.markers)==3&&all(isgraphics(P.markers,'line')));
    assert(isempty(get(P.markers(2),'XData'))&&isempty(get(P.markers(3),'XData')));
    [~,S] = qstab.valueset({[1 1],[2 1]},1i,'Weights',[.5 .5]);
    [returnedAx,H] = qstab.plotValueSet(S,'Parent',ax);
    assert(isequal(returnedAx,ax)&&isequal(H.axes,returnedAx));
    assert(isgraphics(old)&&strcmp(get(ax,'NextPlot'),'replace')&&~H.isCertificate);
end
