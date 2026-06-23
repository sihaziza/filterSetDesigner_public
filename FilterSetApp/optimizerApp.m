classdef optimizerApp < matlab.apps.AppBase
%OPTIMIZERAPP  Archived standalone optimizer layout from OptimalFilterApp.
%   This keeps the optimizer UI available as a separate app shell while the
%   main OptimalFilterApp no longer shows an Optimizer tab.

    properties (Access = public)
        Fig
        OptMode
        OptChanDD
        OptCandList
        OptDichroicList
        OptWeight
        OptBleedWeight
        OptExCheck
        OptScoreMode
        OptSNRFields struct = struct()
        OptAFcheck
        OptResTable
        StatusLbl
    end

    methods (Access = public)
        function app = optimizerApp
            buildUI(app);
            if nargout == 0; clear app; end
        end

        function delete(app)
            if ~isempty(app.Fig) && isvalid(app.Fig)
                delete(app.Fig);
            end
        end
    end

    methods (Access = private)
        function buildUI(app)
            app.Fig = uifigure('Name','optimizerApp - archived optimizer layout', ...
                'Position',[100 100 1050 720]);
            root = uigridlayout(app.Fig,[2 1]);
            root.RowHeight = {'1x',28};
            root.Padding = [10 10 10 8];

            g = uigridlayout(root,[6 2]);
            g.RowHeight = {66,'1x',120,40,40,170};
            g.ColumnWidth = {400,'1x'};
            g.RowSpacing = 8;
            g.ColumnSpacing = 10;

            top = uigridlayout(g,[2 2]);
            top.Layout.Row = 1; top.Layout.Column = 1;
            top.ColumnWidth = {80,'1x'};
            top.RowHeight = {28,28};
            top.Padding = [0 0 0 0];
            uilabel(top,'Text','Mode:');
            app.OptMode = uidropdown(top,'Items', ...
                {'Joint (splitter + all filters)','Single channel'});
            uilabel(top,'Text','Channel:');
            app.OptChanDD = uidropdown(top,'Items',{'(compute first)'});

            p = uipanel(g,'Title','Candidate emission filters for selected channel');
            p.Layout.Row = 2; p.Layout.Column = 1;
            pg = uigridlayout(p,[2 1]);
            pg.RowHeight = {'1x',30};
            app.OptCandList = uilistbox(pg,'Multiselect','on','Items',{});
            pullBtn = uibutton(pg,'Text','+ Pull bandpass filters from FPbase near owner emission');
            pullBtn.Enable = 'off';

            pd = uipanel(g,'Title','Candidate splitter dichroics (joint mode; current always kept)');
            pd.Layout.Row = 3; pd.Layout.Column = 1;
            pdg = uigridlayout(pd,[1 1]);
            app.OptDichroicList = uilistbox(pdg,'Multiselect','on','Items',{});

            wp = uigridlayout(g,[1 5]);
            wp.Layout.Row = 4; wp.Layout.Column = 1;
            wp.ColumnWidth = {135,50,105,50,'1x'};
            wp.Padding = [0 0 0 0];
            wp.ColumnSpacing = 6;
            uilabel(wp,'Text','Crosstalk weight:');
            app.OptWeight = uieditfield(wp,'numeric','Value',5);
            uilabel(wp,'Text','Bleed weight:');
            app.OptBleedWeight = uieditfield(wp,'numeric','Value',5);
            app.OptExCheck = uicheckbox(wp,'Text','Opt. exc. filters');

            rb = uibutton(g,'Text','Run optimizer','FontWeight','bold');
            rb.Layout.Row = 5; rb.Layout.Column = 1;
            rb.Enable = 'off';

            p2 = uipanel(g,'Title','Ranked sets (best first)');
            p2.Layout.Row = [1 6]; p2.Layout.Column = 2;
            p2g = uigridlayout(p2,[2 1]);
            p2g.RowHeight = {'1x',30};
            app.OptResTable = uitable(p2g);
            app.OptResTable.ColumnName = {'Status'};
            app.OptResTable.Data = {'Archived layout only'};
            xb = uibutton(p2g,'Text','Export ranking to CSV');
            xb.Enable = 'off';

            sp = uipanel(g,'Title','Ranking score');
            sp.Layout.Row = 6; sp.Layout.Column = 1;
            spg = uigridlayout(sp,[4 4]);
            spg.RowHeight = {30,30,30,30};
            spg.ColumnWidth = {95,'1x',95,'1x'};
            spg.Padding = [6 6 6 6];
            spg.ColumnSpacing = 6;
            uilabel(spg,'Text','Score by:');
            app.OptScoreMode = uidropdown(spg,'Items', ...
                {'Optical figure of merit','Electron-domain SNR (min channel)'});
            app.OptScoreMode.Layout.Column = [2 4];
            app.OptSNRFields.power = labeledNum(spg,'Power mW',1);
            app.OptSNRFields.integ = labeledNum(spg,'Integ ms',10);
            app.OptSNRFields.conc  = labeledNum(spg,'Conc uM',10);
            app.OptSNRFields.read  = labeledNum(spg,'Read e-',1.5);
            app.OptSNRFields.NA    = labeledNum(spg,'NA',0.5);
            app.OptAFcheck = uicheckbox(spg,'Text','+tissue/fiber AF');
            app.OptAFcheck.Layout.Column = [3 4];

            app.StatusLbl = uilabel(root, ...
                'Text','Optimizer layout archived here. Compute wiring remains in OptimalFilterApp history/code if you want to restore it later.', ...
                'FontColor',[0.3 0.3 0.3]);
        end
    end
end

function h = labeledNum(parent, label, value)
    uilabel(parent,'Text',label);
    h = uieditfield(parent,'numeric','Value',value);
end
