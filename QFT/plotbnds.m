function plotbnds(bdb,ptype,phase,pos)
% PLOTBNDS Plot QFT bounds.
%          PLOTBNDS(BNDS,PTYPE,PHS) plots the QFT bounds in BNDS that are
%          associated with PTYPE (closed-loop configurations) and PHS.
%          PHS is the  phase vector used for computing BNDS
%          with SISOBNDS(PTYPE,...) and GENBNDS(PTYPE,...).
%
%          PLOTBNDS(BNDS) plots all the QFT bounds in BNDS that were computed
%          using the default PHS.
%
%          PLOTBNDS(BNDS,[],PHS) plots the QFT bounds in BNDS that were computed
%          using PHS and uses default value for PTYPE (all problems).
%
%          See also SISOBNDS, GENBNDS, SECTBNDS, PLOTTMPL.

% Author: Craig Borghesani
% Date: 9/6/93
% Revised: 2/17/96 10:02 PM V1.1 updates, 7/3/03 10:05AM v2.5 updates.
% Copyright (c) 2003, Terasoft, Inc.


if ~isstr(bdb)
    if nargin==1
        [bdb,ptype,phase,axs,pos,wbs,wbs2,coora,coorb]=qplotdef(bdb,[],[],[]);
    elseif nargin==2
        [bdb,ptype,phase,axs,pos,wbs,wbs2,coora,coorb]=qplotdef(bdb,ptype,[],[]);
    elseif nargin==3
        [bdb,ptype,phase,axs,pos,wbs,wbs2,coora,coorb]=qplotdef(bdb,ptype,phase,[]);
    elseif nargin==4
        [bdb,ptype,phase,axs,pos,wbs,wbs2,coora,coorb]=qplotdef(bdb,ptype,phase,pos);
    else
        error('Too many inputs');
    end

    lf1=0.0781;
    ht=0.0417;
    tp=1-0.0625;

    [rbdb,cbdb]=size(bdb);
    lwbs2=length(wbs2);

    if lwbs2 ~= cbdb,
        if length(ptype)==1 & ptype(1)~=13,
            if ptype < 10,
                sctlt=['SISOBND',int2str(ptype),' Bound(s)'];
            else
                sctlt=['GENBND',int2str(ptype),' Bound(s)'];
            end
        elseif length(ptype)>1 & any(ptype~=13),
            sctlt='Bounds';
        else
            sctlt='Intersection of Bounds';
        end

        f1 = figure('name',sctlt,'numbertitle','off','units','norm',...
            'position',pos,'vis','off','tag','qft');

        a=gca;
        apos=get(a,'pos');
        set(a,'box','on','xgrid','on','ygrid','on',...
            'gridlinestyle',':',...
            'nextplot','add','xlim',axs(1:2),'ylim',axs(3:4));

        UICont = uicontextmenu;

        if diff(axs(1:2)) >= 360,
            XTickValues = qaxesadjust(axs(1:2));
            set(a,'xtick',XTickValues);
        end

        bnd=qplotbd(phase,bdb,coora,coorb,axs);

        % setup visibility buttons for bounds
        ct=1;
        cvec=['r';'g';'b';'c';'m'];
        clr=[cvec;cvec;cvec;cvec];
        clr=[clr;clr;clr;clr;clr];
        clr=[clr;clr;clr;clr;clr];
        [rb,cb]=size(bnd);
        bnddata2 = {};
        if rb,
            ow=bnd(rb,:); ow=sort(ow); ow(find(diff(ow)==0))=[];
            for t=1:length(ow),
                wstr=sprintf('%4.4g',ow(t));
                loc=find(wstr=='e');
                if length(loc),
                    wstr(find(wstr(loc:length(wstr))=='0')+(loc-1))=[];
                else
                    wstr(find(wstr==' '))=[];
                end
                st2=[wstr,', ',clr(t)];
                st=[wstr];
                uimenu(UICont,'label',st2,'callback','plotbnds(''Toggle'')',...
                    'userdata', t, 'checked', 'on');
                bnddata=bnd(1:(rb-2),find(ow(t)==bnd(rb,:)));
                bnddata2 = [bnddata2, {bnddata}];
                if length(get(bnddata(1),'xdata')) > 1,
                    bndhandle(ct) = bnddata(1);

                else
                    bndhandle(ct) = bnddata(2);

                end;%if

                bndlabel{ct} = st;
                %      bnd_bt(ct,1:2)={st, bnddata};
                tp=tp-ht;
                ct = ct + 1;
            end
        end
        uimenu(UICont,'label','Toggle Off','callback','plotbnds(''Toggle All'')',...
            'separator','on');
        set(a,'uicontextmenu',UICont,'userdata',bnddata2);

        hToolbar = findall(f1, 'type', 'uitoolbar');
        hChil = allchild(hToolbar);
        tooltips = get(hChil, 'tooltip');
        zoomIn = strmatch('Zoom In', tooltips, 'exact');
        zoomOut = strmatch('Zoom Out', tooltips, 'exact');
        delete(hChil(zoomOut));
        set(hChil(zoomIn),'clickedcallback','plotbnds(''Zoom Set'')');

        drawnow;
        xlabel('Open-Loop Phase (deg)');
        ylabel('Open-Loop Gain (dB)');
        legend(bndhandle, bndlabel);

        set(f1,'vis','on');
    else
        disp('Unplottable bounds.');
    end

elseif strcmp(bdb, 'Toggle'),

    tmpdata2 = get(gca,'userdata');
    CurrentItem = get(gcbo,'userdata');

    if strcmp(get(tmpdata2{CurrentItem}(1),'vis'),'on'),
        set(tmpdata2{CurrentItem},'vis','off');
        set(gcbo,'checked','off');
    else
        set(tmpdata2{CurrentItem},'vis','on');
        set(gcbo,'checked','on');
    end;%if

elseif strcmp(bdb, 'Toggle All'),

    CurrentMenu = gcbo;
    CurrentParent = get(CurrentMenu,'parent');
    CurrentChil = flipud(get(CurrentParent,'children'));
    tmpdata2 = get(gca,'userdata');

    if strcmp(get(CurrentMenu,'label'), 'Toggle Off'),
        set(CurrentMenu,'label','Toggle On');
        for k = 1:length(CurrentChil)-1,
            CurrentItem = get(CurrentChil(k),'userdata');
            set(CurrentChil(k),'checked','off');
            set(tmpdata2{k},'vis','off');
        end; %for


    else
        set(CurrentMenu,'label','Toggle Off');
        for k = 1:length(CurrentChil)-1,
            CurrentItem = get(CurrentChil(k),'userdata');
            set(CurrentChil(k),'checked','on');
            set(tmpdata2{k},'vis','on');
        end; %for

    end

elseif strcmp(bdb, 'Zoom Set'),

    f = gcf;
    if strcmp(get(f,'windowbuttondownfcn'),'plotbnds(''Zoom'')'),
        set(f,'windowbuttondownfcn','');
        set(f,'pointer','arrow');
        txt = findobj(f,'tag','uicontextholder');
        UICont = get(txt,'uicontextmenu');
        set(gca,'uicontextmenu',UICont);
        delete(txt);

    else
        txt = text(0,0,' ','vis','off','tag','uicontextholder');
        UICont = get(gca,'uicontextmenu');
        set(txt,'uicontextmenu',UICont);
        set(gca,'uicontextmenu',[]);
        set(f,'windowbuttondownfcn','plotbnds(''Zoom'')');
        setptr(f,'glass');

    end;%if

elseif strcmp(bdb,'Zoom'),

    f = gcf;
    am = gca;

    FigureSelectionType = get(f,'selectiontype');

    CurrentLimits1 = [get(am,'xlim'), get(am,'ylim')];

    switch FigureSelectionType,

        case('normal'),

            CurrentPoint1Axis1 = get(am,'currentpoint');
            fUnits = get(f,'units');
            set(f,'units','pixel');
            rbbox([get(f,'currentpoint') 0 0],...
                get(f,'currentpoint'));
            set(f,'units',fUnits);
            CurrentPoint2Axis1 = get(am,'currentpoint');

            CurrentDiff = sum(abs(CurrentPoint1Axis1(1,:) - CurrentPoint2Axis1(1,:)));

            % store the current limits
            ZoomMemoryData = getappdata(f, 'ZoomMemory');
            if ~isempty(ZoomMemoryData),
                ZoomMemoryData = [ZoomMemoryData; CurrentLimits1];

            else
                ZoomMemoryData = CurrentLimits1;

            end
            setappdata(f,'ZoomMemory',ZoomMemoryData);

            if CurrentDiff > 0.1,
                XLimits1 = [min(CurrentPoint1Axis1(1,1),CurrentPoint2Axis1(1,1)),...
                    max(CurrentPoint1Axis1(1,1),CurrentPoint2Axis1(1,1))];

                YLimits1 = [min(CurrentPoint1Axis1(1,2),CurrentPoint2Axis1(1,2)),...
                    max(CurrentPoint1Axis1(1,2),CurrentPoint2Axis1(1,2))];

            else

                % lets zoom in by 50% around the selected point
                DiffXLimits1 = diff(CurrentLimits1(1:2))/2;
                DiffYLimits1 = diff(CurrentLimits1(3:4))/2;

                XLimits1 = CurrentPoint1Axis1(1,1) + [-DiffXLimits1/2, DiffXLimits1/2];
                YLimits1 = CurrentPoint1Axis1(1,2) + [-DiffYLimits1/2, DiffYLimits1/2];

            end

            if diff(XLimits1) >= 360,
                XTickValues = qaxesadjust(XLimits1);
                set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                    'xtick',XTickValues);

            else
                set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                    'xtickmode','auto');

            end

        case('alt'),

            ZoomMemoryData = getappdata(f,'ZoomMemory');
            if length(ZoomMemoryData),
                XLimits1 = ZoomMemoryData(end, 1:2);
                YLimits1 = ZoomMemoryData(end, 3:4);
                ZoomMemoryData(end, :) = [];
                setappdata(f, 'ZoomMemory', ZoomMemoryData);
                if diff(XLimits1) >= 360,
                    XTickValues = qaxesadjust(XLimits1);
                    set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                        'xtick',XTickValues);

                else
                    set(am, 'xlim', XLimits1, 'ylim', YLimits1,...
                        'xtickmode','auto');

                end
            end

    end % FigureSelectionType

end

