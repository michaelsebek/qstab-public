function [fig,ax,restore] = diagnosticAxes(parent,visible,name)
%DIAGNOSTICAXES Explicit axes with exact NextPlot restoration.
if isstring(visible)&&isscalar(visible),visible=char(visible);end
if ~ischar(visible)||~any(strcmpi(visible,{'on','off'}))
    error('qstab:diagnostics:Visible','Visible must be on or off.');
end
if isempty(parent)
    fig=figure('Visible',lower(visible),'Name',name);ax=axes('Parent',fig);
else
    if ~isscalar(parent)||~isgraphics(parent,'axes')
        error('qstab:diagnostics:Parent','Parent must be one Cartesian axes.');
    end
    ax=parent;fig=ancestor(ax,'figure');
end
old=get(ax,'NextPlot');restore=onCleanup(@()localRestore(ax,old));set(ax,'NextPlot','add');
end
function localRestore(ax,old)
if isgraphics(ax,'axes'),set(ax,'NextPlot',old);end
end
