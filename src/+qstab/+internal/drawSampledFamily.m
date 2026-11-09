function ax=drawSampledFamily(mode,D,o)
% Graphics only. Numerical data are assembled before a figure is created.
if isempty(o.parent)
    fig=figure('Color','w','Name','qstab: sampled polynomial family','Visible',o.visible);
    ax=axes('Parent',fig);
else
    ax=o.parent;
end
np=get(ax,'NextPlot'); fresh=isempty(get(ax,'Children'));
restore=onCleanup(@()restoreHold(ax,np)); %#ok<NASGU>
set(ax,'NextPlot','add');
if strcmp(mode,'zeros')
    z=D.lambda(:,D.valid); z=z(isfinite(z));
    xp=real(z(:)); yp=imag(z(:));
    if strcmp(D.sample.variable,'s')
        extent=max([1;imag(z(:))]);
        line(ax,[0 0],[0 extent],'LineStyle','--','Color',[.34 .36 .38],'DisplayName','Re = 0');
        where='s: stable to the left';
        xp=[xp;0]; yp=[yp;0];
    else
        theta=linspace(0,pi,200);line(ax,cos(theta),sin(theta),'LineStyle','--','Color',[.34 .36 .38],'DisplayName','unit boundary');
        xp=[xp;-1;1]; yp=[yp;0;1];
        if strcmp(D.sample.variable,'zinv'),where='delay: stable outside circle';else,where='z: stable inside circle';end
    end
    if strcmp(D.sample.variable,'s')
        dist=real(z);
    elseif strcmp(D.sample.variable,'zinv')
        dist=1-abs(z);
    else
        dist=abs(z)-1;
    end
    % A narrow floating-point band affects marker styling only.
    band=100*eps*max(1,abs(z));
    sets={dist < -band, dist > band, abs(dist)<=band};
    palette=[.07 .45 .27;.72 .14 .19;.73 .43 .04];
    symbols={'o','x','d'}; captions={'stable side','unstable side','near boundary'};
    hh={};labels={};
    for group=1:3
        zg=z(sets{group});
        if ~isempty(zg)
            hh{end+1}=line('Parent',ax,'XData',real(zg(:)).','YData',imag(zg(:)).', ...
                'Marker',symbols{group},'MarkerSize',5,'LineWidth',1.0, ...
                'LineStyle','none','Color',palette(group,:),'DisplayName',captions{group}); %#ok<AGROW>
            labels{end+1}=captions{group}; %#ok<AGROW>
        end
    end
    if fresh
        xlabel(ax,'Real part');ylabel(ax,'Imaginary radius');
        title(ax,sprintf('Zeros of %d sampled polynomials  |  %s',D.sample.count,where));
        if ~isempty(hh),legend(ax,[hh{:}],labels,'Location','northeast');end
    end
    nbad=nnz(~D.valid);ninf=sum(D.infiniteCount(D.valid));
    if nbad>0 || ninf>0 || isempty(z)
        note=sprintf('Finite reps: %d, infinite reps: %d; unresolved/undefined members: %d',numel(z),ninf,nbad);
        text('Parent',ax,'Units','normalized','Position',[.015 .02 0],'String',note,'Interpreter','none','VerticalAlignment','bottom');
    end
else
    xy=D.projected; xp=xy(:,1); yp=xy(:,2);
    line('Parent',ax,'XData',0,'YData',0,'Marker','+','MarkerSize',12,'LineWidth',1.4, ...
        'LineStyle','none','Color',[.18 .20 .22],'DisplayName','origin');
    if strcmp(D.projection,'slice')
        scatter(ax,xy(:,1),xy(:,2),26,D.offSlice,'filled','DisplayName','values (color: perpendicular length)');
        if fresh
            cb=colorbar(ax);
            if isprop(cb,'Label'), set(get(cb,'Label'),'String','Perpendicular length');else,ylabel(cb,'Perpendicular length');end
        end
    else
        line('Parent',ax,'XData',xy(:,1).','YData',xy(:,2).','Marker','o','MarkerSize',5, ...
            'LineStyle','none','Color',[.04 .44 .49],'DisplayName','sampled values');
    end
    if fresh
        if strcmp(D.projection,'class')
            name='Quaternion values';xl='Real part';yl='Imaginary radius';
        elseif strcmp(D.projection,'slice')
            name='Signed projection';xl='Real part';yl='Imaginary component along u';
        else
            name='Norm-polynomial values';xl='Re N';yl='Im N';
        end
        xlabel(ax,xl); ylabel(ax,yl);
        tag=''; k=size(D.axisUnit,1);
        if ~strcmp(D.projection,'normpoly')
            if k==1
                u=D.axisUnit; names={'i','j','k'};
                [~,j]=max(abs(u));
                if abs(u(j))==1 && nnz(u)==1
                    tag=['; u=' names{j}]; if u(j)<0,tag=['; u=-' names{j}];end
                else
                    tag=sprintf('; u=(%.3g,%.3g,%.3g)',u);
                end
            else
                tag=sprintf('; %d directions',k);
            end
        end
        title(ax,sprintf('%s  |  %d members; omega=%.5g%s',name,D.sample.count,D.frequency,tag));
    end
end
if fresh
    % A padded square keeps axes comparable and shows a collinear cloud
    % without running the end markers into the frame. Always include zero.
    xa=[min([0;xp]),max([0;xp])]; ya=[min([0;yp]),max([0;yp])];
    span=max(diff(xa),diff(ya)); if span==0,span=1;end
    if isfinite(span)
        xl=xa(1)+(span*.5)*[-1.14 1.14]+diff(xa)/2;
        yl=ya(1)+(span*.5)*[-1.14 1.14]+diff(ya)/2;
        if (strcmp(mode,'zeros') || strcmp(D.projection,'class')) && yl(1) < -.06*span
            yl=yl-yl(1)-.06*span;
        end
        if all(isfinite([xl yl])) && xl(1)<xl(2) && yl(1)<yl(2)
            xlim(ax,xl); ylim(ax,yl);
        end
    end
    daspect(ax,[1 1 1]);grid(ax,'on');box(ax,'on');
end
end
function restoreHold(ax,prev)
if isgraphics(ax,'axes'),set(ax,'NextPlot',prev);end
end
