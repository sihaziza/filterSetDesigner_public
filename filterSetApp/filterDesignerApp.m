classdef filterDesignerApp < matlab.apps.AppBase
%FILTERDESIGNERAPP  GUI for selecting optimal optical filter sets for
%   multi-colour fluorescence imaging.
%
%   Run from the filterSetApp folder (or with it on the path):
%       >> filterDesignerApp
%
%   Live tabs
%     1. Design         : define fluorophores, sources, dichroics, filters and
%                         detectors, with schematic spectra.
%     2. Metrics        : photon budget, signal efficiency and physics plots.
%     3. Compare        : side-by-side comparison of captured filter sets.
%     4. Library + Web  : browse local spectra and import/download new spectra.
%
%   The old optimizer layout is archived as standalone optimizerApp.m.
%
%   Depends on: loadSpectrum.m, fetchFPbase.m, FilterEngine.m (same folder).

    properties (Access = public)
        Fig             matlab.ui.Figure
        Lambda          double = (350:1:850)'
        LibRoot         char
        Lib             struct = struct('name',{},'category',{},'kind',{}, ...
                                        'ex',{},'em',{},'brightness',{},'file',{},'floor',{})
        % UI handles
        Tabs
        LibTree
        LibAxes
        LibYMode
        LibPlotIdx      double = []
        WebField
        WebCat
        WebResults
        WebHits         struct = struct('id',{},'name',{},'category',{},'subtype',{})
        FluorTable
        LaserTable
        ChanTable
        PrimaryDD
        DetectorDD
        RbackCheck
        BlockODField
        % schematic tab (dynamic). Box 1 excitation:
        DedupFig
        ExSrcArea
        ExSrcVals       cell = {'488 nm laser'}
        ExCleanVals     cell = {'(none)'}
        ExCombinerVals  cell = {}            % length = #sources-1 (chained)
        ExAxesTop
        ExAxesBot
        ExOpticsMode                          % Individual / Combine spectra toggle
        ExYMode
        ExXMinField
        ExXMaxField
        ExShowXtalk
        ExShowBleed
        SchXtalkPanel
        SchXtalkTable
        SchSourceLeakTable
        SchBleedPanel
        SchBleedTable
        SchDemodTable
        SchDetectorWeighted
        SchDetectorBudgetTable
        SchCrosstalkSummary
        SchXtalkSummary
        SchBleedSummary
        SchEffSummary
        SchMetricFields struct = struct()
        SchMetricTable
        SchEfficiencyTable
        SchMetricStatus
        SchStoichAxes
        SchStoichYMode                 % Y-axis display mode dropdown
        SchStoichReadout
        SchStoichSliderArea            % scrollable grid holding the dynamic knobs
        StoichTotalLines               % per-channel horizontal total-photon lines
        StoichIllumSliders = {}        % one per light source
        StoichConcSliders  = {}        % one per fluorophore
        StoichIllumMul     = []        % per-source illumination, log10 (0 = equi)
        StoichConcMul      = []        % per-fluor concentration, log10 (0 = equi)
        StoichMarkerLines              % per-fluor vertical markers
        StoichMarkerDots               % scatter of marker/curve crossings
        StoichData      struct = struct()
        SchNoiseAxes
        % Box 2 specimen:
        SchPrimaryDD
        SchFluorArea
        SchFluorVals    cell = {'mNeonGreen'}
        % Box 3 emission:
        SchChanArea
        SchEmVals       cell = {'(none)'}
        SchDetVals      cell = {'(ideal)'}
        SchOwnerVals    cell = {'mNeonGreen'}
        SchSplitterVals cell = {}            % length = #channels-1 (chained)
        ResAxes
        ResPlotMode
        ResYMode
        ResXMinField
        ResXMaxField
        ResXSlider
        SigTable
        CTTable
        BleedTable
        LastResults     struct = struct()
        LastPlotConfig  struct = struct()
        OptMode
        OptChanDD
        OptCandList
        OptDichroicList
        OptCandMap
        OptResTable
        OptWeight
        OptBleedWeight
        OptExCheck
        OptScoreMode
        OptSNRFields    struct = struct()
        OptAFcheck
        LastOptResults
        LastOptKind
        StatusLbl
        % Compare tab
        CompareCompTable
        CompareMetricTable
        CompareEntries  struct = struct('name',{},'comp',{},'metrics',{})
        Building        logical = false   % suppress recompute during startup
        MetricsDirty    logical = false   % stoich/noise sweep needs a redraw
        ResultsBuilt    logical = false   % hidden Results panel built on demand
    end

    methods (Access = public)
        function app = filterDesignerApp(libRoot)
            if nargin < 1 || isempty(libRoot)
                here = fileparts(mfilename('fullpath'));
                libRoot = fileparts(here);   % project root
            end
            [~, base] = fileparts(libRoot);
            if strcmpi(base, 'Spectra')
                libRoot = fileparts(libRoot);
            end
            app.LibRoot = libRoot;
            buildUI(app);
            % Suppress the (expensive) schematic recompute while the library
            % loads and defaults seed — the row rebuilds would otherwise each
            % trigger a full metrics + stoichiometry sweep. Recompute once.
            app.Building = true;
            scanLibrary(app);
            seedDefaults(app);
            loadStartupFilterSet(app);
            app.Building = false;
            onSchematicUpdate(app);
            if nargout == 0; clear app; end
        end

        function out = compute(app)
            %COMPUTE  Headless API: returns a struct with S, CT, eff, score,
            %   and the resolved names. Mirrors the Results tab.
            [fluors, lasers, channels, assign, detector, meta] = assembleSystem(app);
            [out.S, out.eff] = FilterEngine.signalMatrix(fluors, lasers, channels, app.Lambda, detector);
            out.CT = FilterEngine.crosstalkMatrix(out.S, assign);
            out.score = FilterEngine.systemScore(out.S, assign, optWeight(app), [], optBleedWeight(app));
            out.fluors = {fluors.name}; out.channels = {channels.name};
            out.meta = meta; out.detector = detector;
        end

        function S = spec(app, name)
            %SPEC  Public lookup of a loaded library spectrum by name.
            S = getSpec(app, name);
        end

        function plotLibrary(app, indices, append)
            %PLOTLIBRARY  Public: superimpose library entries (by index) on the
            %   Library spectral plot. append=true adds to current selection.
            if nargin < 3; append = false; end
            if append; app.LibPlotIdx = unique([app.LibPlotIdx indices(:)'], 'stable');
            else;       app.LibPlotIdx = indices(:)'; end
            redrawLibPlot(app);
        end

        function groups = duplicateGroups(app, thr)
            %DUPLICATEGROUPS  Public: cell array of index groups whose pairwise
            %   zero-lag cross-correlation is >= thr (default 0.999). Scripting.
            if nargin < 2; thr = 0.999; end
            C = correlationMatrix(app); n = size(C,1); grp = 1:n;
            for i = 1:n
                for j = i+1:n
                    if isfinite(C(i,j)) && C(i,j) >= thr; grp(grp==grp(j)) = grp(i); end
                end
            end
            [~,~,gid] = unique(grp,'stable'); groups = {};
            for g = 1:max(gid)
                m = find(gid==g); if numel(m) > 1; groups{end+1} = m; end %#ok<AGROW>
            end
        end

        function refresh(app)
            %REFRESH  Public trigger for the Results computation (same as the
            %   System Builder "Compute results" button).
            onCompute(app);
        end

        function runExcitation(app)
            %RUNEXCITATION  Public trigger for the excitation-test plot.
            onSchematicUpdate(app);
        end

        function excAddSource(app);    addExcitationSource(app);    end
        function excRemoveSource(app); removeExcitationSource(app); end
        function applySchematic(app);  applySchematicToSystem(app); end
        function findDuplicates(app);  onDedup(app); end
        function s = schematicSnapshot(app); s = schematicState(app,'snapshot'); end
        function restoreSchematic(app, s); applySchematicState(app, s); end
        function v = probeCombinerGain(app, wl)
            [~,pk] = min(abs(app.Lambda-wl)); g = combinerGain(app, pk); v = g(pk);
        end

        function runOptimizer(app)
            %RUNOPTIMIZER  Public trigger for the optimizer (= "Run optimizer").
            setStatus(app, 'Optimizer moved to the standalone optimizerApp.m archive.', true);
        end

        function addToCompare(app)
            %ADDTOCOMPARE  Public: snapshot the current schematic + its physics
            %   metrics as a new column in the Compare tab.
            onCompareCapture(app);
        end
        function clearCompare(app)
            %CLEARCOMPARE  Public: empty the Compare tab.
            onCompareClear(app);
        end

        function cfg = config(app)
            %CONFIG  Resolved system struct for the SNR app / scripting.
            [fluors, lasers, channels, assign, detector, meta] = assembleSystem(app);
            % attach FPbase photophysics (ec/qy) when available
            for i = 1:numel(fluors)
                S = getSpec(app, fluors(i).name);
                fluors(i).ec = NaN; fluors(i).qy = NaN;
                if ~isempty(S)
                    if isfield(S,'ec'); fluors(i).ec = S.ec; end
                    if isfield(S,'qy'); fluors(i).qy = S.qy; end
                end
            end
            cfg = struct('lambda',app.Lambda,'fluors',fluors,'lasers',lasers, ...
                'channels',channels,'assign',assign,'detector',detector, ...
                'Rback',meta.Rback,'blockOD',meta.blockOD);
        end
    end

    %% ---------------- UI construction ----------------
    methods (Access = private)
        function buildUI(app)
            ss = get(0,'ScreenSize');                 % use most of the screen
            w = min(1700, ss(3)*0.92); h = min(1040, ss(4)*0.90);
            x = ss(1) + (ss(3)-w)/2; y = ss(2) + (ss(4)-h)/2;
            app.Fig = uifigure('Name','Optimal Filter Set Selector', ...
                'Position',[x y w h]);
            g = uigridlayout(app.Fig,[2 1]);
            g.RowHeight = {'1x',28};
            app.Tabs = uitabgroup(g, 'SelectionChangedFcn', @(s,e)onTabChanged(app));
            app.StatusLbl = uilabel(g,'Text','Ready.','FontColor',[0 0.4 0]);

            buildSchematicTab(app);
            buildPhysicsMetricsTab(app);
            buildCompareTab(app);
            buildLibraryTab(app);
            buildWebTab(app);
            buildSystemTab(app);
            % Results panel (hidden, extra uiaxes + 3 tables) is built lazily on
            % the first Compute via ensureResultsBuilt — keeps startup lighter.
            for tt = app.Tabs.Children'
                if strcmp(tt.Title,'Design'); app.Tabs.SelectedTab = tt; break; end
            end
            bumpFonts(app, app.Fig, 2);   % +2 pt across all static components
        end

        function bumpFonts(~, root, delta)
            %BUMPFONTS  Add delta pt to every FontSize under root (idempotent
            %   per freshly-created subtree; dynamic rebuilds re-apply to their
            %   own newly-created children only).
            if isempty(root) || ~isvalid(root); return; end
            hs = findall(root, '-property', 'FontSize');
            for i = 1:numel(hs)
                try; hs(i).FontSize = hs(i).FontSize + delta; catch; end
            end
        end

        function buildLibraryTab(app)
            % Retired tab kept off-screen: only the legacy primary-dichroic /
            % detector dropdowns live here (still read by the table-based
            % assembleSystem path). The library tree, spectral plot and overlay
            % controls now live in the visible "Library & Web" tab (buildWebTab).
            t = uipanel(app.Fig,'Visible','off');
            g3 = uigridlayout(t,[2 2]); g3.ColumnWidth = {150,'1x'}; g3.RowHeight = {26,26};
            uilabel(g3,'Text','Primary dichroic:');
            app.PrimaryDD = uidropdown(g3,'Items',{'(none)'});
            uilabel(g3,'Text','Detector:');
            app.DetectorDD = uidropdown(g3,'Items',{'(ideal)'});
        end

        function buildWebTab(app)
            % Library & Web: TOP ROW = local spectra database (tree) + FPbase web
            % search, side by side. BOTTOM ROW (full width) = the spectral plot
            % showing the selected / overlaid spectra. Colour follows each curve's
            % peak wavelength; dichroics are dashed.
            t = uitab(app.Tabs,'Title','Library & Web');
            g = uigridlayout(t,[2 1]); g.RowHeight = {250,'1x'};
            g.Padding = [8 8 8 8]; g.RowSpacing = 8;

            % ===== TOP ROW: database (tree) + FPbase search =====
            topRow = uigridlayout(g,[1 2]); topRow.ColumnWidth = {'1x','1x'};
            topRow.Padding = [0 0 0 0]; topRow.ColumnSpacing = 8;

            % left: local spectra tree + overlay buttons
            leftc = uigridlayout(topRow,[2 1]); leftc.RowHeight = {'1x',30};
            leftc.Padding = [0 0 0 0]; leftc.RowSpacing = 6;
            treeP = uipanel(leftc,'Title','Local spectra database (select, then Plot / Add)');
            tpg = uigridlayout(treeP,[1 1]); tpg.Padding = [4 4 4 4];
            app.LibTree = uitree(tpg,'SelectionChangedFcn',@(s,e)onLibSelect(app));
            btns = uigridlayout(leftc,[1 4]); btns.Padding = [0 0 0 0]; btns.ColumnSpacing = 4;
            uibutton(btns,'Text','Plot','ButtonPushedFcn',@(s,e)plotSelectedLibrary(app,false));
            uibutton(btns,'Text','Add','ButtonPushedFcn',@(s,e)plotSelectedLibrary(app,true));
            uibutton(btns,'Text','Remove','ButtonPushedFcn',@(s,e)removeSelectedLibrary(app));
            uibutton(btns,'Text','Clear','ButtonPushedFcn',@(s,e)clearLibraryPlot(app));

            % right: FPbase web search
            p = uipanel(topRow,'Title','FPbase web search  (Semrock / Chroma / Omega / Zeiss + fluorophores)');
            pg = uigridlayout(p,[5 2]);
            pg.RowHeight = {28,28,28,'1x',28}; pg.ColumnWidth = {'1x','1x'};
            pg.Padding = [6 6 6 6]; pg.RowSpacing = 4; pg.ColumnSpacing = 6;
            app.WebField = uieditfield(pg,'text','Placeholder','e.g. FF01-520, Di01-R488, mNeonGreen');
            app.WebField.Layout.Row=1; app.WebField.Layout.Column=1;
            sb = uibutton(pg,'Text','Search','ButtonPushedFcn',@(s,e)onWebSearch(app));
            sb.Layout.Row=1; sb.Layout.Column=2;
            l = uilabel(pg,'Text','Spectra type:','FontWeight','bold');
            l.Layout.Row=2; l.Layout.Column=1;
            app.WebCat = uidropdown(pg,'Items', ...
                {'Fluorophore','Bandpass filter','Long/Short-pass','Dichroic (BS)', ...
                 'Light source','Detector (camera)','Any filter'});
            app.WebCat.Layout.Row=2; app.WebCat.Layout.Column=2;
            db = uibutton(pg,'Text','Download Spectrum', ...
                'ButtonPushedFcn',@(s,e)onWebDownload(app));
            db.Layout.Row=3; db.Layout.Column=1;
            ub = uibutton(pg,'Text','Upload Spectrum', ...
                'ButtonPushedFcn',@(s,e)onUserSpectrumDownload(app));
            ub.Layout.Row=3; ub.Layout.Column=2;
            app.WebResults = uilistbox(pg,'Items',{},'Multiselect','on');
            app.WebResults.Layout.Row=4; app.WebResults.Layout.Column=[1 2];
            rb = uibutton(pg,'Text','Reload local','ButtonPushedFcn',@(s,e)scanLibrary(app));
            rb.Layout.Row=5; rb.Layout.Column=1;
            ddb = uibutton(pg,'Text','Find / remove duplicates', ...
                'ButtonPushedFcn',@(s,e)onDedup(app));
            ddb.Layout.Row=5; ddb.Layout.Column=2;

            % ===== BOTTOM ROW (full width): spectral plot =====
            bottomRow = uigridlayout(g,[2 1]); bottomRow.RowHeight = {30,'1x'};
            bottomRow.Padding = [0 0 0 0]; bottomRow.RowSpacing = 4;
            ctl = uigridlayout(bottomRow,[1 3]); ctl.ColumnWidth = {28,90,'1x'};
            ctl.Padding = [0 0 0 0]; ctl.ColumnSpacing = 6;
            uilabel(ctl,'Text','Y:','FontWeight','bold');
            app.LibYMode = uidropdown(ctl,'Items',{'%T','OD'},'Value','%T', ...
                'ValueChangedFcn',@(s,e)redrawLibPlot(app));
            uilabel(ctl,'Text','Colour = peak wavelength · fluor: ex dashed / em solid · dichroics dashed · "Add" overlays multiple.', ...
                'FontColor',[0.35 0.35 0.35]);
            app.LibAxes = uiaxes(bottomRow);
            xlabel(app.LibAxes,'Wavelength (nm)'); ylabel(app.LibAxes,'%T / normalised intensity (%)');
            title(app.LibAxes,'Select a spectrum');
        end

        function buildSchematicTab(app)
            % Full optical schematic with a dropdown per element, organised as
            % four boxes: (1) fluorophores, (2) excitation, (3) specimen/primary
            % dichroic, (4) emission. Sources, fluorophores and channels are all
            % add/remove; combiner/splitter dichroics appear once there are >=2.
            t = uitab(app.Tabs,'Title','Design');
            g = uigridlayout(t,[1 2]); g.ColumnWidth = {440,'1x'};
            left = uigridlayout(g,[5 1]); left.RowHeight = {'1.35x','1.9x','0.9x','1.85x',36};
            left.Padding = [8 8 8 8]; left.RowSpacing = 8;
            left.Scrollable = 'on';

            % --- Box 1: fluorophores ---
            b1 = uipanel(left,'Title','1 · Fluorophores','FontWeight','bold','FontSize',13, ...
                'BorderType','line','BorderColor',[1 1 1],'BorderWidth',2);
            b1g = uigridlayout(b1,[2 1]); b1g.RowHeight = {'1x',30}; b1g.Padding=[8 8 8 8];
            app.SchFluorArea = uigridlayout(b1g,[1 1]); app.SchFluorArea.Padding=[0 0 0 0];
            app.SchFluorArea.Scrollable = 'on';
            c1 = uigridlayout(b1g,[1 2]); c1.Padding=[0 0 0 0];
            uibutton(c1,'Text','+ Add fluorophore','ButtonPushedFcn',@(s,e)addFluorRow(app));
            uibutton(c1,'Text','− Remove last','ButtonPushedFcn',@(s,e)removeFluorRow(app));

            % --- Box 2: excitation ---
            b2 = uipanel(left,'Title','2 · Excitation  (sources → cleanup → combiner)','FontWeight','bold','FontSize',13, ...
                'BorderType','line','BorderColor',[1 1 1],'BorderWidth',2);
            b2g = uigridlayout(b2,[2 1]); b2g.RowHeight = {'1x',30}; b2g.Padding=[8 8 8 8];
            app.ExSrcArea = uigridlayout(b2g,[1 1]); app.ExSrcArea.Padding = [0 0 0 0];
            app.ExSrcArea.Scrollable = 'on';
            c2 = uigridlayout(b2g,[1 2]); c2.Padding=[0 0 0 0];
            uibutton(c2,'Text','+ Add light source','ButtonPushedFcn',@(s,e)addExcitationSource(app));
            uibutton(c2,'Text','− Remove last','ButtonPushedFcn',@(s,e)removeExcitationSource(app));

            % --- Box 3: specimen / primary dichroic ---
            b3 = uipanel(left,'Title','3 · Specimen  (primary dichroic)','FontWeight','bold','FontSize',13, ...
                'BorderType','line','BorderColor',[1 1 1],'BorderWidth',2);
            pc = uigridlayout(b3,[1 2]); pc.RowHeight={28}; pc.ColumnWidth={'1x',150};
            pc.Padding=[8 6 8 6]; pc.ColumnSpacing=6;
            pcl = uilabel(pc,'Text','Primary dichroic:','FontSize',11);
            pcl.Layout.Row=1; pcl.Layout.Column=1;
            app.SchPrimaryDD = uidropdown(pc,'Items',{'(none)'},'FontSize',12, ...
                'ValueChangedFcn',@(s,e)onSchematicUpdate(app));
            app.SchPrimaryDD.Layout.Row=1; app.SchPrimaryDD.Layout.Column=2;

            % --- Box 4: emission (splitter + per-channel filter + detector) ---
            b4 = uipanel(left,'Title','4 · Emission  (splitter → emission filter → detector)','FontWeight','bold','FontSize',13, ...
                'BorderType','line','BorderColor',[1 1 1],'BorderWidth',2);
            b4g = uigridlayout(b4,[2 1]); b4g.RowHeight = {'1x',30}; b4g.Padding=[8 8 8 8];
            app.SchChanArea = uigridlayout(b4g,[1 1]); app.SchChanArea.Padding=[0 0 0 0];
            app.SchChanArea.Scrollable = 'on';
            c4 = uigridlayout(b4g,[1 2]); c4.Padding=[0 0 0 0];
            uibutton(c4,'Text','+ Add channel','ButtonPushedFcn',@(s,e)addChannelRow(app));
            uibutton(c4,'Text','− Remove last','ButtonPushedFcn',@(s,e)removeChannelRow(app));

            bottom = uigridlayout(left,[1 4]); bottom.Layout.Row = 5; bottom.Padding=[0 0 0 0];
            uibutton(bottom,'Text','Save filter set','ButtonPushedFcn',@(s,e)onSaveSchematic(app));
            uibutton(bottom,'Text','Load filter set','ButtonPushedFcn',@(s,e)onLoadSchematic(app));
            uibutton(bottom,'Text','Compute','FontWeight','bold', ...
                'ButtonPushedFcn',@(s,e)applySchematicToSystem(app));
            uibutton(bottom,'Text','Compare','FontWeight','bold', ...
                'Tooltip','Send this filter set + its physics metrics to the Compare tab', ...
                'ButtonPushedFcn',@(s,e)onCompareCapture(app));

            % --- right: optical schematic plots ---
            right = uigridlayout(g,[1 1]);
            right.Padding = [6 6 6 6];
            rp = uigridlayout(right,[2 2]);
            rp.RowHeight = {34,'1x'}; rp.ColumnWidth = {'1x',300};
            rp.Padding=[0 0 0 0]; rp.ColumnSpacing = 8;
            % Toolbar spans the FULL plot width (columns 1-2, same as the plots
            % below) with left/right padding matching the axes margins, so the
            % Display/tick-boxes start at the plot's left edge and the Y/X-bound
            % group ends at the plot's right edge. The two tick-boxes use 'fit'
            % columns with equal spacer columns on each side, so the gaps
            % Display->Cross-excitation and Cross-excitation->Bleedthrough match.
            % Label columns are 'fit' so each name hugs its value box (tight
            % name/value pairs); explicit 18px spacer columns separate the pairs
            % so the BETWEEN-pair gap is larger than the within-pair gap.
            tg = uigridlayout(rp,[1 15]);
            tg.Layout.Row = 1; tg.Layout.Column = [1 2];
            tg.ColumnWidth = {'fit',170,12,'fit',12,'fit','1x','fit',56,18,'fit',58,18,'fit',58};
            tg.Padding=[44 0 12 0];
            tg.ColumnSpacing = 4;
            uilabel(tg,'Text','Display:','FontWeight','bold');
            app.ExOpticsMode = uidropdown(tg,'Items', ...
                {'Individual spectra','Combine spectra'}, 'Value','Individual spectra', ...
                'ValueChangedFcn',@(s,e)redrawSchematicPlots(app));
            uilabel(tg,'Text','');   % equal spacer (Display -> Cross-excitation)
            app.ExShowXtalk = uicheckbox(tg,'Text','Cross-excitation','Value',true, ...
                'ValueChangedFcn',@(s,e)redrawSchematicPlots(app));
            uilabel(tg,'Text','');   % equal spacer (Cross-excitation -> Bleedthrough)
            app.ExShowBleed = uicheckbox(tg,'Text','Bleedthrough','Value',true, ...
                'ValueChangedFcn',@(s,e)redrawSchematicPlots(app));
            uilabel(tg,'Text','');   % flexible spacer (pushes Y/X to the right)
            uilabel(tg,'Text','Y:','FontWeight','bold');
            app.ExYMode = uidropdown(tg,'Items',{'%T','OD'},'Value','%T', ...
                'ValueChangedFcn',@(s,e)redrawSchematicPlots(app));
            uilabel(tg,'Text','');   % between-pair spacer (Y -> X min)
            uilabel(tg,'Text','X min:','FontWeight','bold');
            app.ExXMinField = uieditfield(tg,'numeric','Value',400, ...
                'Limits',[app.Lambda(1) app.Lambda(end)], ...
                'ValueChangedFcn',@(s,e)onExXRange(app));
            uilabel(tg,'Text','');   % between-pair spacer (X min -> X max)
            uilabel(tg,'Text','X max:','FontWeight','bold');
            app.ExXMaxField = uieditfield(tg,'numeric','Value',800, ...
                'Limits',[app.Lambda(1) app.Lambda(end)], ...
                'ValueChangedFcn',@(s,e)onExXRange(app));
            pax = uigridlayout(rp,[2 1]);
            pax.Layout.Row = 2; pax.Layout.Column = [1 2];
            pax.RowHeight = {'1x','1x'};
            pax.ColumnWidth = {'1x'}; pax.Padding=[0 0 0 0]; pax.ColumnSpacing=8;
            app.ExAxesTop = uiaxes(pax);
            app.ExAxesTop.Layout.Row=1; app.ExAxesTop.Layout.Column=1;
            title(app.ExAxesTop,'Excitation optics + fluorophore absorption');
            ylabel(app.ExAxesTop,'%T / relative');
            app.ExAxesBot = uiaxes(pax);
            app.ExAxesBot.Layout.Row=2; app.ExAxesBot.Layout.Column=1;
            title(app.ExAxesBot,'Fluorophore emission + collection optics');
            xlabel(app.ExAxesBot,'Wavelength (nm)'); ylabel(app.ExAxesBot,'%T / relative');
            % dynamic rows are created by populateChoosers (after the library
            % loads), so the +2pt global font bump in buildUI does not double-
            % apply to them; each rebuild bumps only its own fresh children.
        end

        function buildPhysicsMetricsTab(app)
            t = uitab(app.Tabs,'Title','Metrics');
            g = uigridlayout(t,[2 1]);
            g.RowHeight = {'1x',280};
            g.Padding = [10 10 10 10]; g.RowSpacing = 10;

            top = uigridlayout(g,[1 2]);
            top.Layout.Row = 1;
            top.RowHeight = {'1x'}; top.ColumnWidth = {540,'1x'};
            top.Padding = [0 0 0 0]; top.ColumnSpacing = 10;

            tableStack = uigridlayout(top,[4 1]);
            tableStack.Layout.Row = 1; tableStack.Layout.Column = 1;
            tableStack.RowHeight = {30,'1.05x','0.85x','0.85x'};
            tableStack.Padding = [0 0 0 0]; tableStack.RowSpacing = 8;

            ctl = uigridlayout(tableStack,[1 2]);
            ctl.Layout.Row = 1;
            ctl.ColumnWidth = {190,'1x'}; ctl.Padding = [0 0 0 0]; ctl.ColumnSpacing = 8;
            app.SchDetectorWeighted = uicheckbox(ctl, ...
                'Text','Detector-weighted tables', 'Value',true, ...
                'Tooltip','Apply the selected detector QE curves to excitation contribution and emission bleedthrough tables.', ...
                'ValueChangedFcn',@(s,e)onSchematicUpdate(app));
            app.SchDetectorWeighted.Layout.Row = 1; app.SchDetectorWeighted.Layout.Column = 1;
            desc = uilabel(ctl,'Text','Top tables use equal source power and current schematic optics.', ...
                'FontColor',[0.35 0.35 0.35], 'FontSize',11);
            desc.Layout.Row = 1; desc.Layout.Column = 2;

            app.SchXtalkPanel = uipanel(tableStack,'Title','Excitation crosstalk');
            app.SchXtalkPanel.Layout.Row = 2; app.SchXtalkPanel.Layout.Column = 1;
            xtg = uigridlayout(app.SchXtalkPanel,[3 1]);
            xtg.RowHeight = {18,'1x',32}; xtg.Padding=[4 4 4 4]; xtg.RowSpacing=2;
            uilabel(xtg,'Text','Source -> detected fluor signal (%; equal power)','FontSize',10,'FontWeight','bold');
            app.SchXtalkTable = uitable(xtg); app.SchXtalkTable.Layout.Row = 2;
            app.SchXtalkSummary = uitextarea(xtg, 'Editable','off', 'FontSize',10, ...
                'Value', {'Cross-excitation summary updates with the schematic.'});
            app.SchXtalkSummary.Layout.Row = 3;

            app.SchBleedPanel = uipanel(tableStack,'Title','Emission bleedthrough');
            app.SchBleedPanel.Layout.Row = 3; app.SchBleedPanel.Layout.Column = 1;
            blg = uigridlayout(app.SchBleedPanel,[3 1]);
            blg.RowHeight = {18,'1x',32}; blg.Padding=[4 4 4 4]; blg.RowSpacing=2;
            uilabel(blg,'Text','Detection channel -> fluor emission (% of fluor home ch)','FontSize',10,'FontWeight','bold');
            app.SchBleedTable = uitable(blg); app.SchBleedTable.Layout.Row = 2;
            app.SchBleedSummary = uitextarea(blg, 'Editable','off', 'FontSize',10, ...
                'Value', {'Bleedthrough summary updates with the schematic.'});
            app.SchBleedSummary.Layout.Row = 3;

            ep = uipanel(tableStack,'Title','Fluor signal efficiency vs ideal');
            ep.Layout.Row = 4; ep.Layout.Column = 1;
            eg = uigridlayout(ep,[2 1]); eg.RowHeight = {'1x',32};
            eg.Padding = [4 4 4 4]; eg.RowSpacing = 2;
            app.SchEfficiencyTable = uitable(eg, 'ColumnName', ...
                {'Fluor','Best ch','Current e-/frame','Ideal e-/frame','Efficiency %'}, ...
                'ColumnEditable', false(1,5), 'Data', cell(0,5));
            app.SchEfficiencyTable.Layout.Row = 1;
            app.SchEffSummary = uitextarea(eg, 'Editable','off', 'FontSize',10, ...
                'Value', {'Efficiency summary updates with the schematic.'});
            app.SchEffSummary.Layout.Row = 2;

            % stoichiometry plot + right control panel (imbalance cursors + readout)
            stoichWrap = uigridlayout(top,[1 2]);
            stoichWrap.Layout.Row = 1; stoichWrap.Layout.Column = 2;
            stoichWrap.ColumnWidth = {'1x',190};
            stoichWrap.Padding = [0 0 0 0]; stoichWrap.ColumnSpacing = 6;

            app.SchStoichAxes = uiaxes(stoichWrap);
            app.SchStoichAxes.Layout.Row = 1; app.SchStoichAxes.Layout.Column = 1;
            title(app.SchStoichAxes,'Signal vs stoichiometry');
            xlabel(app.SchStoichAxes,'Fluorophore photon flux (x equi)');
            ylabel(app.SchStoichAxes,'Detected photons (e^- / frame)');
            app.SchStoichAxes.XScale = 'log';

            sp = uigridlayout(stoichWrap,[5 1]);
            sp.Layout.Row = 1; sp.Layout.Column = 2;
            sp.RowHeight = {18,'1x',26,26,58}; sp.Padding = [2 6 2 6]; sp.RowSpacing = 4;
            uilabel(sp,'Text','Operating point (x equi):','FontSize',10,'FontWeight','bold');
            % dynamic knobs: one illumination slider per source + one
            % concentration slider per fluorophore (rebuilt in updateStoichPlot).
            app.SchStoichSliderArea = uigridlayout(sp,[1 1]);
            app.SchStoichSliderArea.Padding = [0 0 0 0]; app.SchStoichSliderArea.Scrollable = 'on';
            app.SchStoichYMode = uidropdown(sp, ...
                'Items',{'Detected photons (e-/frame)','Detected power (nW)','Fraction of channel (%)'}, ...
                'Value','Detected photons (e-/frame)', ...
                'Tooltip','What the y-axis shows for each curve.', ...
                'ValueChangedFcn',@(s,e)redrawStoichFromData(app));
            uibutton(sp,'Text','Reset cursors', ...
                'ButtonPushedFcn',@(s,e)resetStoichCursors(app));
            app.SchStoichReadout = uilabel(sp,'Text','Knobs scale each source / fluorophore from equi (1x).', ...
                'WordWrap','on','FontSize',10,'FontColor',[0.2 0.2 0.2],'VerticalAlignment','top');

            bottom = uigridlayout(g,[1 3]);
            bottom.Layout.Row = 2;
            bottom.ColumnWidth = {'1.25x','1x','1.05x'};
            bottom.Padding = [0 0 0 0]; bottom.ColumnSpacing = 10;

            mp = uipanel(bottom,'Title','Photon budget inputs');
            mp.Layout.Row = 1; mp.Layout.Column = 1;
            mg = uigridlayout(mp,[2 1]);
            mg.RowHeight = {'1x',36}; mg.Padding = [8 8 8 8]; mg.RowSpacing = 6;
            inputs = uigridlayout(mg,[6 4]);
            inputs.RowHeight = repmat({24},1,6);
            inputs.ColumnWidth = {'1x',76,'1x',118};
            inputs.RowSpacing = 4; inputs.ColumnSpacing = 6; inputs.Padding = [0 0 0 0];
            app.SchMetricFields.powerUW = metricNum(app, inputs, 'Output power (mW)', 0.25, ...
                @(s,e)updateMetricPowerFields(app,true));
            app.SchMetricFields.diamMM = metricNum(app, inputs, 'FoV/fiber diam. (mm)', 0.4, ...
                @(s,e)updateMetricPowerFields(app,true));
            app.SchMetricFields.area = metricNum(app, inputs, 'Illum. area (mm2)', 0, [], false);
            app.SchMetricFields.powerDensity = metricNum(app, inputs, 'Power density (mW/mm2)', 0, [], false);
            app.SchMetricFields.conc = metricNum(app, inputs, 'Total fluor conc. (nM)', 1);
            app.SchMetricFields.intMs = metricNum(app, inputs, 'Integration (ms)', 10);
            app.SchMetricFields.NA = metricNum(app, inputs, 'Objective NA', 0.5);
            app.SchMetricFields.read = metricNum(app, inputs, 'Read noise (e-)', 1.5);
            uilabel(inputs,'Text','Tissue AF','FontSize',11);
            app.SchMetricFields.tissueAF = uicheckbox(inputs,'Text','include','Value',true, ...
                'ValueChangedFcn',@(s,e)onSchematicUpdate(app));
            uilabel(inputs,'Text','Fiber AF','FontSize',11);
            app.SchMetricFields.fiberAF = uicheckbox(inputs,'Text','include','Value',true, ...
                'ValueChangedFcn',@(s,e)onSchematicUpdate(app));
            uilabel(inputs,'Text','Laser back-refl.','FontSize',11);
            app.RbackCheck = uicheckbox(inputs,'Text','include (0.5%)','Value',true, ...
                'Tooltip','Include laser back-reflection (R = 0.5%%) in the photon budget', ...
                'ValueChangedFcn',@(s,e)onSchematicUpdate(app));
            uilabel(inputs,'Text','Element blocking (OD)','FontSize',11);
            app.BlockODField = uieditfield(inputs,'numeric','Value',6,'Limits',[0 12], ...
                'FontSize',11,'ValueChangedFcn',@(s,e)onSchematicUpdate(app), ...
                'Tooltip',['Cap on each element''s out-of-band blocking floor for the ' ...
                'laser back-reflection estimate: effective floor = max(measured, 10^-OD).']);
            app.SchMetricFields.conc.Limits = [1 Inf];
            updateMetricPowerFields(app,false);
            app.SchMetricStatus = uilabel(mg,'Text','Metrics update with schematic selections.', ...
                'WordWrap','on','FontColor',[0.25 0.25 0.25]);
            app.SchMetricStatus.Layout.Row = 2;

            tp = uipanel(bottom,'Title','Per-channel photon budget');
            tp.Layout.Row = 1; tp.Layout.Column = 2;
            tg = uigridlayout(tp,[1 1]); tg.Padding = [4 4 4 4];
            app.SchMetricTable = uitable(tg, 'ColumnName', ...
                {'Metric'}, 'ColumnEditable', false, 'Data', cell(0,1));

            detPanel = uipanel(bottom,'Title','Detector photon origin');
            detPanel.Layout.Row = 1; detPanel.Layout.Column = 3;
            detg = uigridlayout(detPanel,[1 1]); detg.Padding=[4 4 4 4];
            app.SchDetectorBudgetTable = uitable(detg);
        end

        function rebuildFluorRows(app)
            if isempty(app.SchFluorArea) || ~isvalid(app.SchFluorArea); return; end
            delete(app.SchFluorArea.Children);
            n = numel(app.SchFluorVals);
            app.SchFluorArea.RowHeight = repmat({44}, 1, max(n,1));
            app.SchFluorArea.ColumnWidth = {'1x'}; app.SchFluorArea.RowSpacing = 4;
            fItems = fluorNames(app);
            app.SchFluorVals = sanitizeChoices(app.SchFluorVals, fItems, fItems{1});
            for i = 1:n
                gridDDv(app, app.SchFluorArea, i,1, sprintf('Fluorophore %d',i), ...
                    fItems, app.SchFluorVals{i}, @(v)setFluor(app,i,v));
            end
            bumpFonts(app, app.SchFluorArea, 2);
            onSchematicUpdate(app);
        end

        function rebuildChannelRows(app)
            if isempty(app.SchChanArea) || ~isvalid(app.SchChanArea); return; end
            delete(app.SchChanArea.Children);
            n = numel(app.SchEmVals); nSplit = max(0, n-1);
            app.SchSplitterVals = padCell(app.SchSplitterVals, nSplit, '(none)');
            nrows = n + nSplit;
            app.SchChanArea.RowHeight = repmat({54}, 1, max(nrows,1));
            app.SchChanArea.ColumnWidth = {'1x','1x'};
            app.SchChanArea.RowSpacing = 6; app.SchChanArea.ColumnSpacing = 6;
            filtItems = ['(none)' roleNames(app,'filter')];
            detItems  = ['(ideal)' detectorPreset('list',app.Lambda) roleNames(app,'detector')];
            dicItems  = ['(none)' roleNames(app,'dichroic')];
            app.SchEmVals = sanitizeChoices(app.SchEmVals, filtItems, '(none)');
            app.SchDetVals = sanitizeChoices(app.SchDetVals, detItems, '(ideal)');
            app.SchSplitterVals = sanitizeChoices(app.SchSplitterVals, dicItems, '(none)');
            for k = 1:n
                lblmode = '';
                if nSplit >= 1; lblmode = ternStr(k<n,' (reflect)',' (transmit)'); end
                gridDDv(app, app.SchChanArea, k,1, sprintf('Ch %d emission filter%s',k,lblmode), ...
                    filtItems, app.SchEmVals{k}, @(v)setChan(app,'em',k,v));
                gridDDv(app, app.SchChanArea, k,2, 'Detector', ...
                    detItems, app.SchDetVals{k}, @(v)setChan(app,'det',k,v));
            end
            for c = 1:nSplit
                gridDDv(app, app.SchChanArea, n+c,1, sprintf('Splitter dichroic %d',c), ...
                    dicItems, app.SchSplitterVals{c}, @(v)setChan(app,'split',c,v));
            end
            bumpFonts(app, app.SchChanArea, 2);
            onSchematicUpdate(app);
        end

        function setFluor(app, i, v)
            app.SchFluorVals{i} = v;
            % keep any channel owners that pointed at the old name in sync-ish
            rebuildChannelRows(app);   % owner dropdowns list current fluors
            onSchematicUpdate(app);
        end

        function setChan(app, which, k, v)
            switch which
                case 'em';    app.SchEmVals{k} = v;
                case 'det';   app.SchDetVals{k} = v;
                case 'owner'; app.SchOwnerVals{k} = v;
                case 'split'; app.SchSplitterVals{k} = v;
            end
            onSchematicUpdate(app);
        end

        function addFluorRow(app)
            cand = {'mRuby3','cyOFP','tdTomato','EGFP','iRFP670'};
            nw = cand{find(~ismember(cand,app.SchFluorVals),1)};
            if isempty(nw); nw = fluorNames(app); nw = nw{1}; end
            app.SchFluorVals{end+1} = nw;
            rebuildFluorRows(app); rebuildChannelRows(app);
        end
        function removeFluorRow(app)
            if numel(app.SchFluorVals) <= 1; setStatus(app,'Keep at least one fluorophore.',true); return; end
            app.SchFluorVals(end) = []; rebuildFluorRows(app); rebuildChannelRows(app);
        end

        function addChannelRow(app)
            app.SchEmVals{end+1} = '(none)';
            app.SchDetVals{end+1} = app.SchDetVals{end};
            fn = fluorNames(app);
            own = fn{min(numel(app.SchEmVals), numel(fn))};
            app.SchOwnerVals{end+1} = own;
            def = '(none)';
            d = roleNames(app,'dichroic');
            if any(strcmp(d,'FF564-Di01')); def = 'FF564-Di01'; elseif ~isempty(d); def = d{1}; end
            app.SchSplitterVals{end+1} = def;
            rebuildChannelRows(app);
        end
        function removeChannelRow(app)
            if numel(app.SchEmVals) <= 1; setStatus(app,'Keep at least one channel.',true); return; end
            app.SchEmVals(end) = []; app.SchDetVals(end) = []; app.SchOwnerVals(end) = [];
            if ~isempty(app.SchSplitterVals); app.SchSplitterVals(end) = []; end
            rebuildChannelRows(app);
        end

        function rebuildExcitationRows(app)
            if isempty(app.ExSrcArea) || ~isvalid(app.ExSrcArea); return; end
            delete(app.ExSrcArea.Children);
            n = numel(app.ExSrcVals); nComb = max(0, n-1);
            app.ExCombinerVals = padCell(app.ExCombinerVals, nComb, '(none)');
            nrows = n + nComb;
            app.ExSrcArea.RowHeight = repmat({54}, 1, max(nrows,1));
            app.ExSrcArea.ColumnWidth = {'1x',26,'1x'};
            app.ExSrcArea.RowSpacing = 6; app.ExSrcArea.ColumnSpacing = 6;
            filtItems = ['(none)' roleNames(app,'filter')];
            dicItems  = ['(none)' roleNames(app,'dichroic')];
            lasers = {'405 nm laser','445 nm laser','488 nm laser','505 nm laser','514 nm laser', ...
                      '561 nm laser','594 nm laser','633 nm laser'};
            srcItems = [lasers, roleNames(app,'source')];
            app.ExSrcVals = sanitizeChoices(app.ExSrcVals, srcItems, lasers{1});
            app.ExCleanVals = sanitizeChoices(app.ExCleanVals, filtItems, '(none)');
            app.ExCombinerVals = sanitizeChoices(app.ExCombinerVals, dicItems, '(none)');
            for i = 1:n
                gridDDv(app, app.ExSrcArea, i,1, sprintf('Light source %d',i), ...
                    srcItems, app.ExSrcVals{i}, @(v)setExc(app,'src',i,v));
                gridArrow(app, app.ExSrcArea, i,2, '→');
                gridDDv(app, app.ExSrcArea, i,3, 'Cleanup filter', ...
                    filtItems, app.ExCleanVals{i}, @(v)setExc(app,'clean',i,v));
            end
            for c = 1:nComb
                gridDDv(app, app.ExSrcArea, n+c,1, sprintf('Combiner dichroic %d',c), ...
                    dicItems, app.ExCombinerVals{c}, @(v)setExc(app,'comb',c,v));
            end
            bumpFonts(app, app.ExSrcArea, 2);
            onSchematicUpdate(app);
        end

        function setExc(app, which, i, v)
            switch which
                case 'src';   app.ExSrcVals{i} = v;
                case 'clean'; app.ExCleanVals{i} = v;
                case 'comb';  app.ExCombinerVals{i} = v;
            end
            onSchematicUpdate(app);
        end

        function addExcitationSource(app)
            used = app.ExSrcVals;
            cand = {'561 nm laser','488 nm laser','633 nm laser','445 nm laser','594 nm laser'};
            newSrc = cand{find(~ismember(cand,used),1)};
            if isempty(newSrc); newSrc = '488 nm laser'; end
            app.ExSrcVals{end+1} = newSrc; app.ExCleanVals{end+1} = '(none)';
            % new combiner for the added source (default FF_511_Di01 if present)
            def = '(none)';
            if any(strcmp(roleNames(app,'dichroic'),'FF_511_Di01')); def = 'FF_511_Di01'; end
            app.ExCombinerVals{end+1} = def;
            rebuildExcitationRows(app);
        end

        function removeExcitationSource(app)
            if numel(app.ExSrcVals) <= 1
                setStatus(app,'Keep at least one light source.',true); return;
            end
            app.ExSrcVals(end) = []; app.ExCleanVals(end) = [];
            if ~isempty(app.ExCombinerVals); app.ExCombinerVals(end) = []; end
            rebuildExcitationRows(app);
        end

        function gridDDv(app, parent, row, col, labelText, items, val, cb)
            c = uigridlayout(parent,[2 1]); c.Layout.Row = row; c.Layout.Column = col;
            c.RowHeight = {16,26}; c.Padding = [2 2 2 2]; c.RowSpacing = 2;
            uilabel(c,'Text',labelText,'FontSize',11);
            if ~any(strcmp(items,val)); items = [items, {val}]; end
            uidropdown(c,'Items',items,'Value',val,'FontSize',12, ...
                'ValueChangedFcn',@(s,e)cb(s.Value));
        end

        function gridArrow(~, parent, row, col, txt)
            l = uilabel(parent,'Text',txt,'FontSize',18,'HorizontalAlignment','center', ...
                'VerticalAlignment','center');
            l.Layout.Row = row; l.Layout.Column = col;
        end

        function h = metricNum(app, parent, label, val, cb, editable)
            if nargin < 5 || isempty(cb); cb = @(s,e)onSchematicUpdate(app); end
            if nargin < 6; editable = true; end
            uilabel(parent,'Text',label,'FontSize',11);
            h = uieditfield(parent,'numeric','Value',val,'FontSize',11, ...
                'Editable',editable,'ValueChangedFcn',cb);
        end

        function v = specVec(app, name, field, default)
            v = default; if strcmp(name,'(none)')||strcmp(name,'(ideal)'); return; end
            S = getSpec(app, name);
            if ~isempty(S) && isfield(S,field) && ~isempty(S.(field)); v = S.(field)(:); end
        end

        function qe = detVec(app, name)
            % Resolve a detector name to a QE vector: (ideal)->1, a built-in
            % preset, or a library Detectors spectrum.
            one = ones(numel(app.Lambda),1);
            if isempty(name) || strcmp(name,'(ideal)'); qe = one; return; end
            if any(strcmp(detectorPreset('list',app.Lambda), name))
                qe = detectorPreset(name, app.Lambda); return;
            end
            qe = specVec(app, name, 'ex', one);
        end

        function g = combinerGain(app, peakIdx)
            % Fraction of a source (peaking at index peakIdx) delivered through
            % the combiner chain. At each combiner the source takes whichever
            % port passes it best at its peak: transmit (T) if T(peak)>=0.5,
            % otherwise reflect (1-T). A well-chosen combiner => gain ~ 1.
            % With multiple sources, a missing combiner is a broken optical
            % path, not a free pass-through.
            g = ones(numel(app.Lambda),1);
            for c = 1:numel(app.ExCombinerVals)
                T = specVec(app, app.ExCombinerVals{c}, 'ex', []);
                if isempty(T)
                    if numel(app.ExSrcVals) > 1; g(:) = 0; end
                    continue;
                end
                if T(peakIdx) >= 0.5; g = g .* T; else; g = g .* (1 - T); end
            end
        end

        function r = primaryExcitationGain(app)
            if strcmp(app.SchPrimaryDD.Value,'(none)')
                r = zeros(numel(app.Lambda),1);
            else
                r = 1 - specVec(app, app.SchPrimaryDD.Value, 'ex', ones(numel(app.Lambda),1));
            end
        end

        function tf = isCombineSpectraMode(app)
            tf = strcmp(normalizeOpticsMode(app, app.ExOpticsMode.Value), 'Combine spectra');
        end

        function mode = normalizeOpticsMode(~, mode)
            if isstring(mode); mode = char(mode); end
            if strcmp(mode, 'Combined result')
                mode = 'Combine spectra';
            end
        end

        function setCombinedOverlayControls(app, enabled)
            state = ternStr(enabled, 'on', 'off');
            if ~isempty(app.ExShowXtalk) && isvalid(app.ExShowXtalk)
                app.ExShowXtalk.Enable = state;
            end
            if ~isempty(app.ExShowBleed) && isvalid(app.ExShowBleed)
                app.ExShowBleed.Enable = state;
            end
        end

        function onSchematicUpdate(app)
            % Full schematic refresh: redraw plots + recompute physics tables.
            % Selection/value changes that affect the OPTICS call this.
            if app.Building; return; end   % skip recompute during startup
            redrawSchematicPlots(app);
            updateSchematicMatrices(app);
            updateSchematicMetrics(app);
        end

        function redrawSchematicPlots(app)
            % Draw the two schematic plots ONLY (cheap: spectra lookups + trapz,
            % no snrModel). Display-only controls (%T/OD, X bounds, overlay
            % toggles, optics view) call this directly so they never trigger the
            % electron-domain physics / stoichiometry sweep.
            if app.Building; return; end
            if isempty(app.ExAxesTop) || ~isvalid(app.ExAxesTop); return; end
            lam = app.Lambda; one = ones(numel(lam),1);
            combined = isCombineSpectraMode(app);
            setCombinedOverlayControls(app, combined);
            % primary dual-band dichroic: reflects excitation, transmits emission.
            % No primary selected means excitation does not reach the specimen.
            if strcmp(app.SchPrimaryDD.Value,'(none)')
                DpRaw = zeros(numel(lam),1); DpEx = zeros(numel(lam),1); DpEm = one;
            else
                Dp = specVec(app, app.SchPrimaryDD.Value, 'ex', one);
                DpRaw = Dp;
                DpEx = 1 - Dp;
                DpEm = Dp;
            end
            % per-source combiner gain (each source uses the combiner band that
            % passes it: transmit where T>=0.5 at its peak, else reflect)
            ns = numel(app.ExSrcVals);
            cf = cell(1,ns);
            for i = 1:ns
                [~,raw] = excSourceCurve(app, app.ExSrcVals{i}, app.ExCleanVals{i});
                [~,pk] = max(raw); cf{i} = combinerGain(app, pk);
            end

            % ===================== TOP: excitation =====================
            ax = app.ExAxesTop; clearSchematicPlotAxes(ax); hold(ax,'on'); leg = {};
            exAtSample = zeros(numel(lam),1);
            exBeforePrimary = zeros(numel(lam),1);
            exBySource = cell(1,ns);
            srcFluorEff = [];
            srcFluorCurve = {};
            mainSrcForFluor = [];
            for i = 1:ns
                [y,~] = excSourceCurve(app, app.ExSrcVals{i}, app.ExCleanVals{i});
                exBeforePrimary = exBeforePrimary + y;
                deliv = y .* cf{i} .* DpEx;          % through combiner chain + primary reflect
                exBySource{i} = deliv;
                exAtSample = exAtSample + deliv;
            end
            if ~combined
                % specimen/primary dichroic on the excitation plot: show what it
                % does to the excitation light — REFLECTION (1-T) toward the
                % specimen — not its transmission curve.
                plot(ax, lam, schemY(app,DpEx), '--', 'Color',[0.86 0.86 0.86], 'LineWidth',1.2);
                leg = [leg, {'primary dichroic (reflection)'}];
                for c = 1:numel(app.ExCombinerVals)
                    Tc = specVec(app, app.ExCombinerVals{c}, 'ex', []);
                    if isempty(Tc); continue; end
                    plot(ax, lam, schemY(app,Tc), '--', 'Color',dichroicColor(lam,Tc), 'LineWidth',1.8);
                    leg = [leg, {sprintf('combiner dichroic %d (T)',c)}]; %#ok<AGROW>
                end
                for i = 1:ns
                    [y,~] = excSourceCurve(app, app.ExSrcVals{i}, app.ExCleanVals{i});
                    [~,ix] = max(y); col = wl2rgb(lam(ix));
                    % light sources: dash-dot (distinct from dichroics' dashes;
                    % only fluorophore spectra are drawn solid)
                    plot(ax, lam, schemY(app,y.*cf{i}.*DpEx), '-.', 'Color',col, 'LineWidth',1.8);
                    leg = [leg, {[app.ExSrcVals{i} ' at sample']}]; %#ok<AGROW>
                end
                for i = 1:numel(app.SchFluorVals)
                    ex = specVec(app, app.SchFluorVals{i}, 'ex', []); if isempty(ex); continue; end
                    ex = ex/max([max(ex) eps]); [~,ix]=max(ex); col = wl2rgb(lam(ix));
                    plot(ax, lam, schemY(app,ex), '-', 'Color',col, 'LineWidth',1.2);
                    leg = [leg, {[app.SchFluorVals{i} ' abs']}]; %#ok<AGROW>
                end
            else
                % one trace per fluorophore: its absorption actually excited by
                % the delivered light = abs(λ) x delivered-excitation shape
                nfTop = numel(app.SchFluorVals);
                srcFluorEff = zeros(ns,nfTop);
                srcFluorCurve = cell(ns,nfTop);
                normDen = max([max(exBeforePrimary) eps]);
                for i = 1:nfTop
                    ex = specVec(app, app.SchFluorVals{i}, 'ex', []); if isempty(ex); continue; end
                    ex = ex/max([max(ex) eps]);
                    emForColor = specVec(app, app.SchFluorVals{i}, 'em', []);
                    if isempty(emForColor)
                        [~,ix] = max(ex);
                    else
                        [~,ix] = max(emForColor);
                    end
                    col = wl2rgb(lam(ix));
                    for j = 1:ns
                        srcFluorCurve{j,i} = ex .* (exBySource{j}/normDen);
                        srcFluorEff(j,i) = trapz(lam, ex .* exBySource{j});
                    end
                    [mainEff, mainSrc] = max(srcFluorEff(:,i));
                    mainSrcForFluor(i) = mainSrc; %#ok<AGROW>
                    for j = 1:ns
                        curve = srcFluorCurve{j,i};
                        if isempty(curve) || max(curve) <= eps; continue; end
                        isMain = j == mainSrc;
                        if ~isMain && ~app.ExShowXtalk.Value; continue; end  % hide cross-excitation
                        style = ternStr(isMain, '-', '--');
                        lw = ternStr(isMain, 2.1, 1.5);
                        plot(ax, lam, schemY(app,curve), style, 'Color',col, 'LineWidth',lw);
                        leg = [leg, {sprintf('%s / %s', app.SchFluorVals{i}, ...
                            sourcePairLabel(app, app.ExSrcVals{j}))}]; %#ok<AGROW>
                    end
                end
            end
            ylabel(ax,schemAxisLabel(app,'specimen')); finishAx(app, ax, leg);

            % ===================== BOTTOM: emission =====================
            ax = app.ExAxesBot; clearSchematicPlotAxes(ax); hold(ax,'on'); leg = {};
            nc = numel(app.SchEmVals);
            splT = cell(1,numel(app.SchSplitterVals));
            for c = 1:numel(splT); splT{c} = specVec(app, app.SchSplitterVals{c}, 'ex', []); end
            sf = FilterEngine.splitterFactors(splT, numel(lam));
            % excitation efficiency of each fluorophore (0-1): absorption-weighted
            % fraction of the delivered excitation it captures
            denomE = trapz(lam, exBeforePrimary) + eps;
            nf = numel(app.SchFluorVals); excEff = zeros(1,nf); emN = cell(1,nf);
            for i = 1:nf
                ex = specVec(app, app.SchFluorVals{i}, 'ex', zeros(numel(lam),1));
                ex = ex/max([max(ex) eps]);
                excEff(i) = trapz(lam, ex .* exAtSample) / denomE;
                em = specVec(app, app.SchFluorVals{i}, 'em', zeros(numel(lam),1));
                emN{i} = em/max([max(em) eps]);
            end
            if ~combined
                plot(ax, lam, schemY(app,DpRaw), '--', 'Color',[0.86 0.86 0.86], 'LineWidth',1.2);
                leg = [leg, {'primary dichroic'}];
                for c = 1:numel(splT)
                    if isempty(splT{c}); continue; end
                    plot(ax, lam, schemY(app,splT{c}), '--', 'Color',dichroicColor(lam,splT{c}), 'LineWidth',1.8);
                    leg = [leg, {sprintf('splitter dichroic %d (T)',c)}]; %#ok<AGROW>
                end
                for i = 1:nf
                    [~,ix] = max(emN{i}); col = wl2rgb(lam(ix));
                    plot(ax, lam, schemY(app,emN{i}), '-', 'Color',col, 'LineWidth',1.4);
                    leg = [leg, {[app.SchFluorVals{i} ' em']}]; %#ok<AGROW>
                end
                for k = 1:nc
                    emF = specVec(app, app.SchEmVals{k}, 'ex', one);
                    thru = DpEm .* sf{k} .* emF .* detVec(app, app.SchDetVals{k});
                    [~,ix] = max(emF); col = wl2rgb(lam(ix));
                    % detection-path throughput (filters/dichroics/detector):
                    % dash-dot, distinct from solid fluorophore emission
                    plot(ax, lam, schemY(app,thru), '-.', 'Color',col, 'LineWidth',2.2);
                    leg = [leg, {sprintf('Ch%d collection x detector',k)}]; %#ok<AGROW>
                end
                ylabel(ax,schemAxisLabel(app,'detector'));
            else
                % Detected fluorescence flux per channel (% of ideal):
                %   solid  = the channel's own fluorophore (signal)
                %   dashed = optional diagnostic overlays:
                %            - cross-excited emission from a non-primary source
                %            - other fluorophores bleeding through into this channel
                thruK = cell(1,nc); emColl = zeros(nf,nc);
                for k = 1:nc
                    thruK{k} = DpEm .* sf{k} .* specVec(app, app.SchEmVals{k}, 'ex', one) ...
                               .* detVec(app, app.SchDetVals{k});
                    for i = 1:nf; emColl(i,k) = trapz(lam, emN{i} .* thruK{k}); end
                end
                homeI = ones(1,nf);                          % each fluor's own (best) channel
                for i = 1:nf; [~,homeI(i)] = max(emColl(i,:)); end
                ownerK = ones(1,nc);
                for k = 1:nc; [~,ownerK(k)] = max(excEff(:).*emColl(:,k)); end
                for k = 1:nc
                    o = ownerK(k);
                    detected = excEff(o) * emN{o} .* thruK{k};
                    [~,ix] = max(emN{o}); col = wl2rgb(lam(ix));
                    plot(ax, lam, schemY(app,detected), '-', 'Color',col, 'LineWidth',2.4);
                    leg = [leg, {sprintf('Ch%d %s', k, app.SchFluorVals{o})}]; %#ok<AGROW>  % numbers in table
                    if ~isempty(srcFluorEff) && app.ExShowXtalk.Value && o <= numel(mainSrcForFluor)
                        for j = 1:ns
                            if j == mainSrcForFluor(o); continue; end
                            xtalkEmission = (srcFluorEff(j,o) / denomE) * emN{o} .* thruK{k};
                            if max(xtalkEmission) <= 1e-6; continue; end
                            plot(ax, lam, schemY(app,xtalkEmission), '--', 'Color',col, 'LineWidth',1.5);
                            leg = [leg, {sprintf('Ch%d_%s/%s_xtalk', k, app.SchFluorVals{o}, ...
                                sourcePairLabel(app, app.ExSrcVals{j}))}]; %#ok<AGROW>
                        end
                    end
                    if app.ExShowBleed.Value
                        for i = 1:nf
                            if i == o; continue; end
                            contam = excEff(i) * emN{i} .* thruK{k};
                            if max(contam) <= 1e-6; continue; end
                            [~,ix] = max(emN{i}); col = wl2rgb(lam(ix));
                            plot(ax, lam, schemY(app,contam), '--', 'Color',col, 'LineWidth',1.5);
                            leg = [leg, {sprintf('Ch%d_%s_bleedthrough', k, app.SchFluorVals{i})}]; %#ok<AGROW>
                        end
                    end
                end
                ylabel(ax,schemAxisLabel(app,'detector'));
            end
            finishAx(app, ax, leg);
            xlabel(ax,'Wavelength (nm)');
            function finishAx(~, axx, lg)
                hold(axx,'off'); grid(axx,'on');
                xlo = 400; xhi = 800;
                if ~isempty(app.ExXMinField) && isvalid(app.ExXMinField)
                    xlo = app.ExXMinField.Value; xhi = app.ExXMaxField.Value;
                end
                if ~(xhi > xlo); xlo = 400; xhi = 800; end
                xlim(axx,[xlo xhi]);
                if strcmp(app.ExYMode.Value,'OD')
                    ylim(axx,[0 8]); axx.YDir = 'reverse';
                else
                    ylim(axx,[0 105]); axx.YDir = 'normal';
                end
                if ~isempty(lg)
                    L = legend(axx, lg, 'Location','northoutside','Interpreter','none', ...
                        'Orientation','horizontal','NumColumns',3);
                    L.FontSize = 9;
                end
            end
        end

        function onExXRange(app)
            % Validate the schematic plot x-axis bounds, then redraw plots only
            % (axis rescale doesn't change the physics tables).
            lo = app.ExXMinField.Value; hi = app.ExXMaxField.Value;
            if hi <= lo
                hi = min(app.Lambda(end), lo + 10);
                app.ExXMaxField.Value = hi;
            end
            redrawSchematicPlots(app);
        end

        function tf = isMetricsVisible(app)
            tf = ~isempty(app.Tabs) && isvalid(app.Tabs) && ...
                ~isempty(app.Tabs.SelectedTab) && ...
                strcmp(app.Tabs.SelectedTab.Title,'Metrics');
        end

        function onTabChanged(app)
            % Recompute the deferred Metrics stoichiometry/noise sweep the first
            % time the Metrics tab is shown after a schematic change.
            if app.Building; return; end
            if isMetricsVisible(app) && app.MetricsDirty
                updateSchematicMetrics(app);
            end
        end

        function updateSchematicMatrices(app)
            % Fill the two side tables:
            %   cross-excitation %  : rows = light sources, cols = fluorophores
            %                         (absorption-weighted % of each source absorbed)
            %   bleedthrough %      : rows = fluorophores,  cols = channels
            %                         (% of a fluor's detected emission landing in
            %                          each channel; its own channel = 100%)
            if isempty(app.SchXtalkTable) || ~isvalid(app.SchXtalkTable); return; end
            if ~isempty(app.SchXtalkPanel) && isvalid(app.SchXtalkPanel); app.SchXtalkPanel.Visible = 'on'; end
            if ~isempty(app.SchBleedPanel) && isvalid(app.SchBleedPanel); app.SchBleedPanel.Visible = 'on'; end
            M = schematicCrosstalkMetrics(app);
            app.SchXtalkTable.RowName = app.ExSrcVals;
            app.SchXtalkTable.ColumnName = app.SchFluorVals;
            app.SchXtalkTable.Data = formatPctTable(app, M.sourceFreqPct');
            if ~isempty(app.SchSourceLeakTable) && isvalid(app.SchSourceLeakTable)
                app.SchSourceLeakTable.RowName = app.ExSrcVals;
                app.SchSourceLeakTable.ColumnName = app.SchFluorVals;
                app.SchSourceLeakTable.Data = formatPctTable(app, M.sourceFreqPct');
            end
            app.SchBleedTable.RowName = arrayfun(@(k)sprintf('Ch%d',k),1:numel(app.SchEmVals),'UniformOutput',false);
            app.SchBleedTable.ColumnName = app.SchFluorVals;
            app.SchBleedTable.Data = formatPctTable(app, M.bleedPct');
            if ~isempty(app.SchDemodTable) && isvalid(app.SchDemodTable)
                app.SchDemodTable.RowName = strcat('target ', app.SchFluorVals);
                app.SchDemodTable.ColumnName = strcat('from ', app.SchFluorVals);
                app.SchDemodTable.Data = formatDemodTable(app, M.demodPct);
            end
            updateCrosstalkSummary(app, M);
        end

        function M = schematicCrosstalkMetrics(app)
            % Optical-domain tables analogous to uSMAART bleedthrough /
            % frequency-demodulation reporting:
            %   crossAbsPct   : source spectrum overlap with each fluorophore.
            %   bleedPct      : where each fluor's emission is collected.
            %   sourceFreqPct : which source-frequency generated the collected
            %                   fluorescence in that fluor's home channel.
            %   demodPct      : contaminant fluor surviving demodulation at the
            %                   target fluor's source frequency, in the target
            %                   detection channel, relative to target signal.
            lam = app.Lambda; one = ones(numel(lam),1);
            ns = numel(app.ExSrcVals); nf = numel(app.SchFluorVals); nc = numel(app.SchEmVals);
            if strcmp(app.SchPrimaryDD.Value,'(none)')
                Dp = zeros(numel(lam),1);
            else
                Dp = specVec(app, app.SchPrimaryDD.Value, 'ex', one);   % T transmits emission
            end
            DpEx = 1 - Dp;
            exN = zeros(numel(lam),nf);
            emPeak = zeros(numel(lam),nf);
            bright = ones(1,nf);
            for i = 1:nf
                ex = specVec(app, app.SchFluorVals{i}, 'ex', zeros(numel(lam),1));
                exN(:,i) = normPeak(ex);
                em = specVec(app, app.SchFluorVals{i}, 'em', zeros(numel(lam),1));
                emPeak(:,i) = normPeak(em);
                S = getSpec(app, app.SchFluorVals{i});
                if ~isempty(S) && isfield(S,'brightness') && ~isempty(S.brightness) && isfinite(S.brightness)
                    bright(i) = S.brightness;
                end
            end
            srcShape = zeros(numel(lam),ns);
            for j = 1:ns
                [y,~] = excSourceCurve(app, app.ExSrcVals{j}, app.ExCleanVals{j});
                [~,pk] = max(y);
                deliv = y .* combinerGain(app,pk) .* DpEx;
                den = trapz(lam, deliv);
                if den > 0
                    srcShape(:,j) = deliv / den; % equal power at specimen
                end
            end
            sourceAbs = zeros(nf,ns);
            for i = 1:nf
                for j = 1:ns
                    sourceAbs(i,j) = trapz(lam, exN(:,i) .* srcShape(:,j));
                end
            end
            M.crossAbsPct = 100 * sourceAbs';

            splT = cell(1,numel(app.SchSplitterVals));
            for c = 1:numel(splT); splT{c} = specVec(app, app.SchSplitterVals{c}, 'ex', []); end
            sf = FilterEngine.splitterFactors(splT, numel(lam));
            emColl = zeros(nf,nc);
            useDetector = detectorWeightedCrosstalk(app);
            for k = 1:nc
                det = one;
                if useDetector
                    det = detVec(app, app.SchDetVals{k});
                end
                thru = Dp .* sf{k} .* specVec(app, app.SchEmVals{k}, 'ex', one) .* det;
                for i = 1:nf
                    emColl(i,k) = trapz(lam, emPeak(:,i) .* thru);
                end
            end
            M.bleedPct = zeros(nf,nc);
            homeCh = ones(1,nf);
            for i = 1:nf
                [homeVal, homeCh(i)] = max(emColl(i,:));
                M.bleedPct(i,:) = 100 * emColl(i,:) / max(homeVal, eps);
            end

            contrib = zeros(nf,ns,nc);
            for i = 1:nf
                for j = 1:ns
                    for k = 1:nc
                        contrib(i,j,k) = bright(i) * sourceAbs(i,j) * emColl(i,k);
                    end
                end
            end
            M.sourceFreqPct = zeros(nf,ns);
            mainSrc = ones(1,nf);
            for i = 1:nf
                [~, mainSrc(i)] = max(sourceAbs(i,:));
                den = sum(contrib(i,:,homeCh(i)));
                M.sourceFreqPct(i,:) = 100 * contrib(i,:,homeCh(i)) / max(den, eps);
            end

            M.demodPct = nan(nf,nf);
            signalScale = max(contrib(:));
            for target = 1:nf
                src = mainSrc(target);
                ch = homeCh(target);
                targetSignal = contrib(target,src,ch);
                if targetSignal <= max(signalScale*1e-9, eps)
                    continue;
                end
                for contam = 1:nf
                    M.demodPct(target,contam) = 100 * contrib(contam,src,ch) / targetSignal;
                end
            end
            M.homeCh = homeCh;
            M.mainSrc = mainSrc;
        end

        function tf = detectorWeightedCrosstalk(app)
            tf = true;
            if ~isempty(app.SchDetectorWeighted) && isvalid(app.SchDetectorWeighted)
                tf = app.SchDetectorWeighted.Value;
            end
        end

        function data = formatPctTable(app, X)
            data = cell(size(X));
            for ii = 1:numel(X)
                data{ii} = formatPercent(app, X(ii));
            end
        end

        function data = formatDemodTable(app, X)
            data = cell(size(X));
            for r = 1:size(X,1)
                for c = 1:size(X,2)
                    if r == c
                        data{r,c} = 'signal';
                    else
                        data{r,c} = formatPercent(app, X(r,c));
                    end
                end
            end
        end

        function s = formatPercent(~, x)
            if ~isfinite(x)
                s = 'n/a';
            elseif x > 0 && x < 0.01
                s = '<0.01%';
            elseif abs(x) < 0.1
                s = sprintf('%.3f%%', x);
            elseif abs(x) < 10
                s = sprintf('%.2f%%', x);
            else
                s = sprintf('%.1f%%', x);
            end
        end

        function updateCrosstalkSummary(app, M)
            % Populates the per-table blurbs (cross-excitation + bleedthrough).
            % The old combined blurb under "Detector photon origin" was removed,
            % so app.SchCrosstalkSummary may be empty — setTextArea guards that.
            lines = {};
            nf = numel(app.SchFluorVals);
            if nf == 0 || isempty(M.bleedPct)
                setTextArea(app.SchCrosstalkSummary, 'Add fluorophores/channels to summarize crosstalk.');
                setTextArea(app.SchXtalkSummary, 'Add fluorophores/sources to summarize cross-excitation.');
                setTextArea(app.SchBleedSummary, 'Add fluorophores/channels to summarize bleedthrough.');
                return;
            end
            [bleedVal, bleedIdx] = maxOffDiagLike(app, M.bleedPct, M.homeCh);
            bleedLine = 'No significant emission bleedthrough between channels.';
            if isfinite(bleedVal)
                [fi, ch] = ind2sub(size(M.bleedPct), bleedIdx);
                bleedLine = sprintf('Optical bleedthrough: %s of %s enters Ch%d (home Ch%d).', ...
                    formatPercent(app, bleedVal), app.SchFluorVals{fi}, ch, M.homeCh(fi));
                lines{end+1} = bleedLine; %#ok<AGROW>
            end
            setTextArea(app.SchBleedSummary, bleedLine);
            sourceLeak = M.sourceFreqPct;
            for i = 1:nf
                if M.mainSrc(i) >= 1 && M.mainSrc(i) <= size(sourceLeak,2)
                    sourceLeak(i,M.mainSrc(i)) = -Inf;
                end
            end
            [srcVal, srcIdx] = max(sourceLeak(:));
            srcLine = 'Each fluorophore is driven mainly by its own source.';
            if isfinite(srcVal) && srcVal > 0
                [fi, sj] = ind2sub(size(sourceLeak), srcIdx);
                srcLine = sprintf('Source-frequency leak: %s of %s home-channel fluorescence is driven by %s.', ...
                    formatPercent(app, srcVal), app.SchFluorVals{fi}, shortSourceLabel(app, app.ExSrcVals{sj}));
                lines{end+1} = srcLine; %#ok<AGROW>
            end
            setTextArea(app.SchXtalkSummary, srcLine);
            dem = M.demodPct;
            dem(eye(size(dem))>0) = -Inf;
            dem(~isfinite(dem)) = -Inf;
            [demVal, demIdx] = max(dem(:));
            if isfinite(demVal) && demVal > 0
                [target, contam] = ind2sub(size(dem), demIdx);
                lines{end+1} = sprintf('After demod: %s %s contamination remains in the %s-frequency target signal.', ...
                    formatPercent(app, demVal), app.SchFluorVals{contam}, app.SchFluorVals{target}); %#ok<AGROW>
            end
            if isempty(lines)
                lines = {'No nonzero crosstalk terms detected with the current schematic.'};
            end
            setTextArea(app.SchCrosstalkSummary, lines);
        end

        function [mx, idx] = maxOffDiagLike(~, X, homeCh)
            Y = X;
            for i = 1:numel(homeCh)
                if homeCh(i) >= 1 && homeCh(i) <= size(Y,2)
                    Y(i,homeCh(i)) = -Inf;
                end
            end
            [mx, idx] = max(Y(:));
            if isempty(mx) || ~isfinite(mx); mx = NaN; idx = 1; end
        end

        function names = fluorNames(app)
            mask = strcmp({app.Lib.category},'Proteins') & strcmp({app.Lib.kind},'fluorophore');
            names = unique({app.Lib(mask).name}, 'stable');
            if isempty(names); names = {'mNeonGreen'}; end
        end

        function updateSchematicMetrics(app)
            if isempty(app.SchMetricTable) || ~isvalid(app.SchMetricTable); return; end
            try
                [fluors, lasers, channels, assign, Rback, blockOD] = schematicSystem(app);
                nS = numel(lasers); nC = numel(channels);
                if isempty(fluors) || isempty(lasers) || isempty(channels)
                    app.SchMetricTable.Data = cell(0,7);
                    clearFluorEfficiencyTable(app);
                    app.SchMetricStatus.Text = 'Add at least one fluorophore, source, and channel.';
                    return;
                end
                pSrc_mW = metricPowerPerSource(app);
                phys = metricPhys(app, nS, pSrc_mW);
                conc = max(1, app.SchMetricFields.conc.Value) * 1e-9;
                fp = metricPhotophysics(app, fluors, conc * ones(1,numel(fluors)));
                cfg = struct('lambda',app.Lambda,'fluors',fluors,'lasers',lasers, ...
                    'channels',channels,'assign',assign,'detector',ones(numel(app.Lambda),1), ...
                    'Rback',Rback,'blockOD',blockOD);
                af = metricAutofluor(app);
                out = snrModel(cfg, phys, fp, af);
                metricNames = {'Owner'; 'Signal e-/frame'; 'SNR'; 'Xtalk %'; 'Bleed %'; 'Limit'};
                vals = cell(numel(metricNames),nC);
                for k = 1:nC
                    sig = out.signal_e(k);
                    xt = 100*out.xtalk_e(k) / max(sig, eps);
                    bl = 100*out.exBleed_e(k) / max(sig, eps);
                    vals(:,k) = {fluors(assign(k)).name; ...
                        round(sig,3,'significant'); round(out.SNR(k),3,'significant'); ...
                        round(xt,3,'significant'); round(bl,3,'significant'); ...
                        limitingSource(out, k)};
                end
                if isfield(out,'afBySrc') && ~isempty(out.afBySrc)
                    for a = 1:numel(out.afNames)
                        metricNames{end+1,1} = [out.afNames{a} ' AF e-/frame']; %#ok<AGROW>
                        vals(end+1,:) = num2cell(round(out.afBySrc(a,:),3,'significant')); %#ok<AGROW>
                    end
                end
                app.SchMetricTable.ColumnName = [{'Metric'}, {channels.name}];
                app.SchMetricTable.ColumnEditable = false(1,nC+1);
                app.SchMetricTable.Data = [metricNames, vals];
                app.SchMetricStatus.Text = sprintf('%.3g mW/source, %.3g ms, %.3g nM', ...
                    pSrc_mW, app.SchMetricFields.intMs.Value, app.SchMetricFields.conc.Value);
                app.SchMetricStatus.FontColor = [0 0.4 0];
                updateDetectorBudgetTable(app, cfg, phys, fp, af);
                updateFluorEfficiencyTable(app, cfg, phys, fp);
                % The stoichiometry sweep runs snrModel ~81x — only draw it when
                % the Metrics tab is actually visible; otherwise mark it dirty so
                % onTabChanged recomputes it the next time the tab is opened.
                if isMetricsVisible(app)
                    updatePhysicsPlots(app, cfg, phys, fp, out);
                    app.MetricsDirty = false;
                else
                    app.MetricsDirty = true;
                end
            catch ME
                app.SchMetricTable.Data = cell(0,7);
                app.SchMetricStatus.Text = ['Metrics unavailable: ' ME.message];
                app.SchMetricStatus.FontColor = [0.7 0 0];
                clearDetectorBudgetTable(app);
                clearFluorEfficiencyTable(app);
                clearPhysicsPlots(app, ME.message);
            end
        end

        function updateFluorEfficiencyTable(app, cfg, phys, fp)
            if isempty(app.SchEfficiencyTable) || ~isvalid(app.SchEfficiencyTable); return; end
            [current, bestCh] = fluorSignalMatrixElectrons(app, cfg, phys, fp);
            ideal = idealFluorSignals(app, cfg, phys, fp);
            nF = numel(cfg.fluors);
            data = cell(nF,5);
            for i = 1:nF
                if isempty(current)
                    curBest = 0; chName = '';
                else
                    curBest = current(i,bestCh(i));
                    chName = cfg.channels(bestCh(i)).name;
                end
                effPct = 100 * curBest / max(ideal(i), eps);
                data(i,:) = {cfg.fluors(i).name, chName, ...
                    sprintf('%.3g', curBest), sprintf('%.3g', ideal(i)), ...
                    sprintf('%.2f%%', effPct)};
            end
            app.SchEfficiencyTable.ColumnName = ...
                {'Fluor','Best ch','Current e-/frame','Ideal e-/frame','Efficiency %'};
            app.SchEfficiencyTable.ColumnEditable = false(1,5);
            app.SchEfficiencyTable.Data = data;
            % blurb: weakest fluorophore vs the best-possible (peak-excited,
            % ideal long-pass) signal.
            effLine = 'Efficiency vs an ideal peak-excited / long-pass setup.';
            if nF >= 1 && ~isempty(current)
                effPctAll = zeros(1,nF);
                for i = 1:nF; effPctAll(i) = 100*current(i,bestCh(i))/max(ideal(i),eps); end
                [lo,li] = min(effPctAll);
                effLine = sprintf('%s collects %.1f%% of its best-possible signal (peak excitation + ideal long-pass).', ...
                    cfg.fluors(li).name, lo);
            end
            setTextArea(app.SchEffSummary, effLine);
        end

        function clearFluorEfficiencyTable(app)
            if isempty(app.SchEfficiencyTable) || ~isvalid(app.SchEfficiencyTable); return; end
            app.SchEfficiencyTable.ColumnName = ...
                {'Fluor','Best ch','Current e-/frame','Ideal e-/frame','Efficiency %'};
            app.SchEfficiencyTable.Data = cell(0,5);
            setTextArea(app.SchEffSummary, 'Add a fluorophore to summarize collection efficiency.');
        end

        function [N, bestCh] = fluorSignalMatrixElectrons(~, cfg, phys, fp)
            lambda = cfg.lambda(:);
            hc = 1.98644586e-25;
            photonE = hc ./ (lambda*1e-9);
            nF = numel(cfg.fluors); nS = numel(cfg.lasers); nC = numel(cfg.channels);
            N = zeros(nF,nC);
            bestCh = ones(1,nF);
            if nF == 0 || nC == 0; return; end
            th = asin(min(phys.NA/phys.n, 1));
            colEff = (1 - cos(th))/2;
            PhiTot = zeros(numel(lambda),1);
            for j = 1:nS
                src = cfg.lasers(j);
                if isfield(phys,'powers_mW') && numel(phys.powers_mW) >= j
                    src.power = phys.powers_mW(j)*1e-3;
                end
                PhiTot = PhiTot + FilterEngine.sourceExcitation(src, lambda) ./ photonE;
            end
            Tk = zeros(numel(lambda),nC);
            for k = 1:nC
                Tk(:,k) = FilterEngine.pathTransmission(cfg.channels(k), lambda);
            end
            for i = 1:nF
                absFrac = 1 - 10.^(-(fp(i).ec * cfg.fluors(i).ex(:)) * fp(i).conc_M * phys.pathLength_cm);
                absorbRate = trapz(lambda, PhiTot .* absFrac);
                emN = cfg.fluors(i).em(:);
                emArea = trapz(lambda, emN);
                if emArea > 0; emN = emN / emArea; end
                for k = 1:nC
                    detFrac = trapz(lambda, emN .* Tk(:,k));
                    N(i,k) = absorbRate * fp(i).qy * colEff * detFrac * phys.tInt_s;
                end
                [~, bestCh(i)] = max(N(i,:));
            end
        end

        function ideal = idealFluorSignals(~, cfg, phys, fp)
            % Best-achievable signal for each fluorophore: excited at its own
            % excitation peak, and collected through a PERFECT dichroic +
            % emission long-pass whose edge sits IDEAL_LP_OFFSET_NM above that
            % excitation peak (T=1 above the edge, 0 below; ideal QE=1). This is
            % the most signal a real filter set could capture, since emission
            % within a few nm of the excitation can never be separated from the
            % laser line. NA-limited collection (colEff) is kept so the
            % efficiency isolates spectral/filter quality, not solid angle.
            IDEAL_LP_OFFSET_NM = 5;
            lambda = cfg.lambda(:);
            hc = 1.98644586e-25;
            photonE = hc ./ (lambda*1e-9);
            nF = numel(cfg.fluors);
            ideal = zeros(1,nF);
            th = asin(min(phys.NA/phys.n, 1));
            colEff = (1 - cos(th))/2;
            p_mW = 0;
            if isfield(phys,'powers_mW') && ~isempty(phys.powers_mW)
                p_mW = sum(max(0, phys.powers_mW));
            end
            for i = 1:nF
                ex = cfg.fluors(i).ex(:);
                if isempty(ex) || max(ex) <= 0 || p_mW <= 0; continue; end
                ex = ex / max(ex);
                [~, pk] = max(ex);
                absFrac = 1 - 10^(-(fp(i).ec * ex(pk)) * fp(i).conc_M * phys.pathLength_cm);
                photonRate = (p_mW*1e-3) / photonE(pk);
                % ideal long-pass collection fraction (emission past the edge)
                emN = cfg.fluors(i).em(:); aEm = trapz(lambda, emN);
                if aEm > 0; emN = emN / aEm; end
                edge = lambda(pk) + IDEAL_LP_OFFSET_NM;
                idealLP = double(lambda >= edge);
                detFracIdeal = trapz(lambda, emN .* idealLP);
                ideal(i) = photonRate * absFrac * fp(i).qy * colEff * phys.tInt_s * detFracIdeal;
            end
        end

        function updateDetectorBudgetTable(app, cfg, phys, fp, af)
            if isempty(app.SchDetectorBudgetTable) || ~isvalid(app.SchDetectorBudgetTable); return; end
            B = detectorPhotonBudget(app, cfg, phys, fp, af);
            labels = {'Total e-/frame'; 'Main signal %'; 'Source xtalk %'; ...
                      'Fluor bleedthrough %'; 'Brain tissue AF %'; 'Silica fiber AF %'};
            values = [B.total; B.main; B.sourceXtalk; B.fluorBleed; B.tissueAF; B.fiberAF];
            for j = 1:numel(cfg.lasers)
                labels{end+1,1} = sprintf('Back-ref %s %%', shortSourceLabel(app, cfg.lasers(j).name)); %#ok<AGROW>
                values(end+1,:) = B.backBySource(j,:); %#ok<AGROW>
            end
            data = cell(numel(labels), numel(cfg.channels)+1);
            data(:,1) = labels;
            for r = 1:numel(labels)
                for k = 1:numel(cfg.channels)
                    if r == 1
                        data{r,k+1} = sprintf('%.3g', values(r,k));
                    else
                        data{r,k+1} = formatPercent(app, values(r,k));
                    end
                end
            end
            app.SchDetectorBudgetTable.ColumnName = [{'Component'}, {cfg.channels.name}];
            app.SchDetectorBudgetTable.ColumnEditable = false(1,numel(cfg.channels)+1);
            app.SchDetectorBudgetTable.Data = data;
        end

        function label = shortSourceLabel(~, name)
            label = char(string(name));
            tok = regexp(label,'(\d{3})','tokens','once');
            if ~isempty(tok)
                label = [tok{1} ' nm'];
            elseif numel(label) > 18
                label = [label(1:15) '...'];
            end
        end

        function label = sourcePairLabel(~, name)
            label = char(string(name));
            parts = regexp(label, '[_\s]+', 'split');
            if ~isempty(parts) && ~isempty(parts{1})
                label = parts{1};
            end
            if numel(label) > 18
                label = [label(1:15) '...'];
            end
        end

        function clearDetectorBudgetTable(app)
            if isempty(app.SchDetectorBudgetTable) || ~isvalid(app.SchDetectorBudgetTable); return; end
            app.SchDetectorBudgetTable.ColumnName = {'Component'};
            app.SchDetectorBudgetTable.Data = cell(0,1);
        end

        function B = detectorPhotonBudget(app, cfg, phys, fp, af)
            lambda = cfg.lambda(:);
            hc = 1.98644586e-25;
            photonE = hc ./ (lambda*1e-9);
            QE = ones(numel(lambda),1);
            if isfield(cfg,'detector') && ~isempty(cfg.detector); QE = cfg.detector(:); end
            nF = numel(cfg.fluors); nS = numel(cfg.lasers); nC = numel(cfg.channels);
            th = asin(min(phys.NA/phys.n, 1));
            colEff = (1 - cos(th))/2;
            tInt = phys.tInt_s;
            Phi = cell(1,nS); PhiTot = zeros(numel(lambda),1);
            for j = 1:nS
                src = cfg.lasers(j);
                if isfield(phys,'powers_mW') && numel(phys.powers_mW) >= j
                    src.power = phys.powers_mW(j)*1e-3;
                end
                E = FilterEngine.sourceExcitation(src, lambda);
                Phi{j} = E ./ photonE;
                PhiTot = PhiTot + Phi{j};
            end
            Tk = zeros(numel(lambda),nC);
            TkFloor = zeros(numel(lambda),nC);
            for k = 1:nC
                Tk(:,k) = FilterEngine.pathTransmission(cfg.channels(k), lambda);
                TkFloor(:,k) = FilterEngine.pathTransmissionFloored(cfg.channels(k), lambda, cfg.blockOD);
            end
            fluorSrc = zeros(nF,nS,nC);
            absorbRate = zeros(nF,nS);
            for i = 1:nF
                ec = fp(i).ec; qy = fp(i).qy; c = fp(i).conc_M; d = phys.pathLength_cm;
                absFrac = 1 - 10.^(-(ec * cfg.fluors(i).ex(:)) * c * d);
                emN = cfg.fluors(i).em(:); a = trapz(lambda, emN); if a > 0; emN = emN/a; end
                for j = 1:nS
                    absorbRate(i,j) = trapz(lambda, Phi{j} .* absFrac);
                    for k = 1:nC
                        detFrac = trapz(lambda, emN .* Tk(:,k) .* QE);
                        fluorSrc(i,j,k) = absorbRate(i,j) * qy * colEff * detFrac * tInt;
                    end
                end
            end
            mainSrc = ones(1,nF);
            for i = 1:nF
                if nS > 0; [~,mainSrc(i)] = max(absorbRate(i,:)); end
            end
            main = zeros(1,nC); srcX = zeros(1,nC); bleed = zeros(1,nC);
            for k = 1:nC
                o = cfg.assign(k);
                if o >= 1 && o <= nF && nS > 0
                    main(k) = fluorSrc(o,mainSrc(o),k);
                    srcX(k) = sum(fluorSrc(o,:,k)) - main(k);
                end
                other = setdiff(1:nF,o);
                if ~isempty(other); bleed(k) = sum(fluorSrc(other,:,k),'all'); end
            end
            back = zeros(nS,nC);
            for j = 1:nS
                for k = 1:nC
                    back(j,k) = cfg.Rback * trapz(lambda, Phi{j} .* TkFloor(:,k) .* QE) * tInt;
                end
            end
            afByName = zeros(numel(af),nC);
            for m = 1:numel(af)
                absFrac = af(m).strength * af(m).absorb(:);
                emN = af(m).em(:); a = trapz(lambda, emN); if a > 0; emN = emN/a; end
                for j = 1:nS
                    Rabs = trapz(lambda, Phi{j} .* absFrac);
                    for k = 1:nC
                        detFrac = trapz(lambda, emN .* Tk(:,k) .* QE);
                        afByName(m,k) = afByName(m,k) + Rabs * af(m).qy * colEff * detFrac * tInt;
                    end
                end
            end
            tissue = zeros(1,nC); fiber = zeros(1,nC);
            for m = 1:numel(af)
                nm = lower(af(m).name);
                if contains(nm,'brain') || contains(nm,'tissue')
                    tissue = tissue + afByName(m,:);
                elseif contains(nm,'fiber') || contains(nm,'silica')
                    fiber = fiber + afByName(m,:);
                end
            end
            total = main + srcX + bleed + tissue + fiber + sum(back,1);
            den = max(total, eps);
            B.total = total;
            B.main = 100*main ./ den;
            B.sourceXtalk = 100*srcX ./ den;
            B.fluorBleed = 100*bleed ./ den;
            B.tissueAF = 100*tissue ./ den;
            B.fiberAF = 100*fiber ./ den;
            B.backBySource = 100*back ./ den;
        end

        function updateMetricPowerFields(app, doRedraw)
            if nargin < 2; doRedraw = true; end
            if ~isfield(app.SchMetricFields,'powerUW') || isempty(app.SchMetricFields.powerUW) || ...
                    ~isvalid(app.SchMetricFields.powerUW)
                return;
            end
            p_mW = max(0, app.SchMetricFields.powerUW.Value);
            diam = max(eps, app.SchMetricFields.diamMM.Value);
            area = pi * (diam/2)^2;
            app.SchMetricFields.area.Value = area;
            app.SchMetricFields.powerDensity.Value = p_mW / area;
            if doRedraw; onSchematicUpdate(app); end
        end

        function pSrc_mW = metricPowerPerSource(app)
            updateMetricPowerFields(app,false);
            pSrc_mW = max(0, app.SchMetricFields.powerUW.Value);
        end

        function phys = metricPhys(app, nS, pSrc_mW)
            phys = struct('NA',max(eps,app.SchMetricFields.NA.Value), ...
                'n',1.33, 'pathLength_cm',0.01, ...
                'tInt_s',max(eps,app.SchMetricFields.intMs.Value)*1e-3, ...
                'readNoise_e',max(0,app.SchMetricFields.read.Value), ...
                'darkRate_eps',100, 'powers_mW',pSrc_mW*ones(1,nS));
        end

        function fp = metricPhotophysics(app, fluors, concs)
            fp = struct('ec',{},'qy',{},'conc_M',{});
            for i = 1:numel(fluors)
                src = getSpec(app, fluors(i).name);
                ec = 50000; qy = 0.6; haveEc = false; haveQy = false;
                if ~isempty(src)
                    if isfield(src,'ec') && isfinite(src.ec); ec = src.ec; haveEc = true; end
                    if isfield(src,'qy') && isfinite(src.qy); qy = src.qy; haveQy = true; end
                end
                % If ec/qy are not both known, fold in the fluorophore's
                % brightness (FPbase brightness = ec*qy/1000) so brighter
                % fluorophores produce proportionally more detected photons.
                b = fluorBrightnessValue(app, src, fluors(i));
                if (~haveEc || ~haveQy) && isfinite(b) && b > 0
                    ecqy = b * 1000;                 % = ec * qy
                    if haveQy && qy > 0;       ec = ecqy / qy;
                    elseif haveEc && ec > 0;   qy = min(1, ecqy / ec);
                    else;                      qy = 0.6; ec = ecqy / qy;
                    end
                end
                c = 0; if numel(concs) >= i; c = concs(i); end
                fp(i) = struct('ec',ec,'qy',qy,'conc_M',c); %#ok<AGROW>
            end
        end

        function b = fluorBrightnessValue(~, src, fluor)
            % Brightness (ec*qy/1000) from the library entry, else the resolved
            % fluor struct, else NaN (unknown -> caller uses defaults).
            b = NaN;
            if ~isempty(src) && isfield(src,'brightness') && ~isempty(src.brightness) ...
                    && isfinite(src.brightness); b = src.brightness; return; end
            if nargin >= 3 && isstruct(fluor) && isfield(fluor,'brightness') ...
                    && ~isempty(fluor.brightness) && isfinite(fluor.brightness); b = fluor.brightness; end
        end

        function af = metricAutofluor(app)
            af = struct('name',{},'absorb',{},'em',{},'strength',{},'qy',{});
            if isfield(app.SchMetricFields,'tissueAF') && ~isempty(app.SchMetricFields.tissueAF) && ...
                    isvalid(app.SchMetricFields.tissueAF) && app.SchMetricFields.tissueAF.Value
                af(end+1) = autofluorPreset('Brain tissue',app.Lambda); %#ok<AGROW>
            end
            if isfield(app.SchMetricFields,'fiberAF') && ~isempty(app.SchMetricFields.fiberAF) && ...
                    isvalid(app.SchMetricFields.fiberAF) && app.SchMetricFields.fiberAF.Value
                af(end+1) = autofluorPreset('Silica fiber',app.Lambda); %#ok<AGROW>
            end
        end

        function updatePhysicsPlots(app, cfg, phys, fp, outNow)
            updateStoichPlot(app, cfg, phys);
        end

        function clearPhysicsPlots(app, msg)
            if ~isempty(app.SchStoichAxes) && isvalid(app.SchStoichAxes)
                cla(app.SchStoichAxes); title(app.SchStoichAxes,'Stoichiometry sweep unavailable');
                text(app.SchStoichAxes,0.5,0.5,msg,'Units','normalized','HorizontalAlignment','center');
            end
        end

        function updateStoichPlot(app, cfg, phys)
            % Per-fluorophore detected signal vs that fluorophore's photon flux
            % relative to equi-illumination/equi-concentration. Each fluorophore
            % gets a curve set (solid = its home channel, dashed = its leak into
            % other channels), coloured by its emission (matching the Design
            % "Combine spectra" plot). A draggable marker per fluorophore sits at
            % the photon-flux multiplier set by the per-source illumination and
            % per-fluor concentration knobs.
            ax = app.SchStoichAxes;
            if isempty(ax) || ~isvalid(ax); return; end
            app.StoichData = struct();
            srcNames = {}; if ~isempty(cfg.lasers); srcNames = {cfg.lasers.name}; end
            fNames = {}; if ~isempty(cfg.fluors); fNames = {cfg.fluors.name}; end
            if isempty(cfg.fluors) || isempty(cfg.channels) || isempty(cfg.lasers)
                cla(ax); hold(ax,'on'); ax.XScale = 'log';
                title(ax,'Add fluorophores, a source and a channel');
                text(ax,0.5,0.5,'Define the schematic in the Design tab.', ...
                    'Units','normalized','HorizontalAlignment','center');
                hold(ax,'off');
                rebuildStoichSliders(app, srcNames, fNames);
                return;
            end
            nF = numel(cfg.fluors); nC = numel(cfg.channels); nS = numel(cfg.lasers);
            lam = cfg.lambda(:);
            concPer = max(1, app.SchMetricFields.conc.Value) * 1e-9;
            % equi per-fluor x per-channel detected electrons (all conc = concPer,
            % all sources ON). Scaling a fluorophore's photon flux by f scales its
            % row by f, so the curves are N(i,k)*f.
            N = stoichTempSignalByFluor(app, cfg, phys, concPer*ones(1,nF));   % nF x nC
            % source -> fluor excitation weights (for the marker illumination mix)
            W = zeros(nS,nF);
            for s = 1:nS
                E = FilterEngine.sourceExcitation(cfg.lasers(s), lam);
                for i = 1:nF; W(s,i) = trapz(lam, E .* cfg.fluors(i).ex(:)); end
            end
            homeCh = ones(1,nF); thresh = zeros(1,nF); fcolors = zeros(nF,3); lamEm = zeros(1,nF);
            for i = 1:nF
                h = find(cfg.assign==i, 1); if isempty(h); [~,h] = max(N(i,:)); end
                homeCh(i) = h; thresh(i) = max(0.01*N(i,h), eps);   % draw leaks > 1% of home
                [~,ix] = max(cfg.fluors(i).em); fcolors(i,:) = wl2rgb(lam(ix)); lamEm(i) = lam(ix);
            end
            ccolors = zeros(nC,3);
            for k = 1:nC; ccolors(k,:) = fcolors(min(cfg.assign(k),nF),:); end
            tInt = 1e-2; if isfield(phys,'tInt_s'); tInt = phys.tInt_s; end
            app.StoichData = struct('N',N,'W',W,'homeCh',homeCh,'thresh',thresh, ...
                'fcolors',fcolors,'ccolors',ccolors,'lamEm',lamEm,'tInt',tInt, ...
                'fnames',{fNames},'cnames',{{cfg.channels.name}}, 'fAxis',logspace(-2,2,81));
            rebuildStoichSliders(app, srcNames, fNames);
            redrawStoichFromData(app);
        end

        function redrawStoichFromData(app)
            % Plot the stored per-fluor curves (N(i,k)*f) and (re)draw markers.
            ax = app.SchStoichAxes;
            if isempty(ax) || ~isvalid(ax); return; end
            D = app.StoichData;
            if ~isfield(D,'N') || isempty(D.N); return; end
            cla(ax); hold(ax,'on');
            app.StoichMarkerLines = {}; app.StoichMarkerDots = []; app.StoichTotalLines = {};
            f = D.fAxis; nF = size(D.N,1); nC = size(D.N,2); leg = {};
            fracMode = strcmp(app.SchStoichYMode.Value,'Fraction of channel (%)');
            for i = 1:nF
                col = D.fcolors(i,:);
                for k = 1:nC
                    v = D.N(i,k);
                    if v < D.thresh(i); continue; end
                    isHome = (k == D.homeCh(i));
                    sty = '-'; lw = 2.1; if ~isHome; sty = '--'; lw = 1.3; end
                    plot(ax, f, stoichCurveY(app, D, i, k, f), sty, 'Color',col, 'LineWidth',lw);
                    if isHome
                        leg{end+1} = sprintf('%s -> %s', D.fnames{i}, D.cnames{k}); %#ok<AGROW>
                    else
                        leg{end+1} = sprintf('%s -> %s (xtalk)', D.fnames{i}, D.cnames{k}); %#ok<AGROW>
                    end
                end
            end
            xline(ax, 1, '-', 'equi', 'Color',[0.6 0.6 0.6], 'LabelOrientation','horizontal');
            grid(ax,'on'); ax.XScale = 'log'; xlim(ax,[1e-2 1e2]);
            if fracMode
                ax.YScale = 'linear'; ylim(ax,[0 100]);
                ylabel(ax, 'Fraction of channel total (%)');
            else
                ax.YScale = 'log';
                ylabel(ax, app.SchStoichYMode.Value);
            end
            xlabel(ax, 'Fluorophore photon flux (x equi)');
            title(ax, 'Per-fluorophore detected signal vs photon flux');
            if ~isempty(leg)
                legend(ax, leg, 'Location','best', 'Interpreter','none', 'FontSize',8);
            end
            hold(ax,'off');
            drawStoichMarkers(app);
        end

        function yv = stoichCurveY(app, D, i, k, f)
            % Transform a curve's detected electrons N(i,k)*f to the chosen Y-mode.
            e = D.N(i,k) .* f;
            switch app.SchStoichYMode.Value
                case 'Detected power (nW)'
                    yv = e .* (1.98644586e-25 / (D.lamEm(i)*1e-9)) ./ max(D.tInt,eps) * 1e9;
                case 'Fraction of channel (%)'
                    others = sum(D.N(:,k)) - D.N(i,k);   % other fluors at equi (1x)
                    yv = 100 .* e ./ (e + others + eps);
                otherwise
                    yv = e;
            end
        end

        function rebuildStoichSliders(app, srcNames, fluorNames)
            % (Re)create one illumination slider per source + one concentration
            % slider per fluorophore, preserving values when counts are unchanged.
            area = app.SchStoichSliderArea;
            if isempty(area) || ~isvalid(area); return; end
            delete(area.Children);
            nS = numel(srcNames); nF = numel(fluorNames);
            if numel(app.StoichIllumMul) ~= nS; app.StoichIllumMul = zeros(1,nS); end
            if numel(app.StoichConcMul)  ~= nF; app.StoichConcMul  = zeros(1,nF); end
            app.StoichIllumSliders = cell(1,nS); app.StoichConcSliders = cell(1,nF);
            nrows = max(nS + nF, 1);
            area.RowHeight = repmat({44}, 1, nrows); area.ColumnWidth = {'1x'};
            area.RowSpacing = 2;
            row = 0;
            for s = 1:nS
                row = row + 1;
                app.StoichIllumSliders{s} = makeStoichKnob(app, area, row, ...
                    sprintf('%s  illum', shortSourceLabel(app, srcNames{s})), ...
                    app.StoichIllumMul(s), 'illum', s);
            end
            for i = 1:nF
                row = row + 1;
                app.StoichConcSliders{i} = makeStoichKnob(app, area, row, ...
                    sprintf('%s  conc', fluorNames{i}), app.StoichConcMul(i), 'conc', i);
            end
        end

        function h = makeStoichKnob(app, parent, row, labelText, val, kind, idx)
            c = uigridlayout(parent,[2 1]); c.Layout.Row = row; c.Layout.Column = 1;
            c.RowHeight = {14,24}; c.Padding = [0 0 0 0]; c.RowSpacing = 0;
            uilabel(c,'Text',labelText,'FontSize',9);
            h = uislider(c,'Limits',[-2 2],'Value',val, 'MajorTicks',[-2 0 2], ...
                'MajorTickLabels',{'/100','1','x100'}, ...
                'ValueChangedFcn',@(s,e)onStoichKnob(app,kind,idx,e.Value), ...
                'ValueChangingFcn',@(s,e)onStoichKnob(app,kind,idx,e.Value));
        end

        function onStoichKnob(app, kind, idx, val)
            if strcmp(kind,'illum')
                if idx <= numel(app.StoichIllumMul); app.StoichIllumMul(idx) = val; end
            else
                if idx <= numel(app.StoichConcMul); app.StoichConcMul(idx) = val; end
            end
            drawStoichMarkers(app);
        end

        function resetStoichCursors(app)
            % Return every illumination / concentration knob to equi (1x = 0 in log).
            app.StoichIllumMul = zeros(size(app.StoichIllumMul));
            app.StoichConcMul  = zeros(size(app.StoichConcMul));
            for c = 1:numel(app.StoichIllumSliders)
                h = app.StoichIllumSliders{c}; if ~isempty(h) && isvalid(h); h.Value = 0; end
            end
            for c = 1:numel(app.StoichConcSliders)
                h = app.StoichConcSliders{c}; if ~isempty(h) && isvalid(h); h.Value = 0; end
            end
            drawStoichMarkers(app);
        end

        function drawStoichMarkers(app)
            % One vertical marker per fluorophore at its photon-flux multiplier
            % (conc knob x illumination mix), with dots where it crosses that
            % fluorophore's curves. Updated in place for smooth dragging.
            ax = app.SchStoichAxes;
            if isempty(ax) || ~isvalid(ax); return; end
            D = app.StoichData;
            if ~isfield(D,'N') || isempty(D.N); return; end
            nF = size(D.N,1); nC = size(D.N,2);
            if ~iscell(app.StoichMarkerLines) || numel(app.StoichMarkerLines) ~= nF
                if iscell(app.StoichMarkerLines)
                    for c = 1:numel(app.StoichMarkerLines)
                        hh = app.StoichMarkerLines{c};
                        if ~isempty(hh) && isvalid(hh); delete(hh); end
                    end
                end
                app.StoichMarkerLines = cell(1,nF);
            end
            xmul = ones(1,nF);                       % operating photon-flux per fluor
            for i = 1:nF; xmul(i) = 10^min(2, max(-2, stoichFluxLog(app, D, i))); end
            fracMode = strcmp(app.SchStoichYMode.Value,'Fraction of channel (%)');
            xs = []; ys = []; cs = zeros(0,3); parts = {};
            for i = 1:nF
                xi = xmul(i); col = D.fcolors(i,:);
                h = app.StoichMarkerLines{i};
                if isempty(h) || ~isvalid(h)
                    app.StoichMarkerLines{i} = xline(ax, xi, '-', 'Color',col, ...
                        'LineWidth',1.4, 'Alpha',0.85, 'Tag','stoichMarker');
                else
                    h.Value = xi;
                end
                for k = 1:nC
                    if D.N(i,k) < D.thresh(i); continue; end
                    xs(end+1) = xi; ys(end+1) = stoichCurveY(app, D, i, k, xi); cs(end+1,:) = col; %#ok<AGROW>
                end
                parts{end+1} = sprintf('%s flux x%.2g -> %.3g e-/frame', ...
                    D.fnames{i}, xi, D.N(i,D.homeCh(i))*xi); %#ok<AGROW>
            end
            if isempty(app.StoichMarkerDots) || ~isvalid(app.StoichMarkerDots)
                hold(ax,'on');
                app.StoichMarkerDots = scatter(ax, xs, ys, 42, cs, 'filled', ...
                    'MarkerEdgeColor','k', 'Tag','stoichMarker');
                hold(ax,'off');
            elseif isempty(xs)
                set(app.StoichMarkerDots, 'XData',[], 'YData',[]);
            else
                set(app.StoichMarkerDots, 'XData',xs, 'YData',ys, 'CData',cs);
            end
            % horizontal TOTAL line per detector channel = sum of all fluors'
            % contributions at the current operating point (hidden in % mode).
            if ~iscell(app.StoichTotalLines) || numel(app.StoichTotalLines) ~= nC
                if iscell(app.StoichTotalLines)
                    for c = 1:numel(app.StoichTotalLines)
                        hh = app.StoichTotalLines{c};
                        if ~isempty(hh) && isvalid(hh); delete(hh); end
                    end
                end
                app.StoichTotalLines = cell(1,nC);
            end
            for k = 1:nC
                totE = 0; totP = 0;
                for i = 1:nF
                    e = D.N(i,k) * xmul(i); totE = totE + e;
                    totP = totP + e * (1.98644586e-25/(D.lamEm(i)*1e-9)) / max(D.tInt,eps) * 1e9;
                end
                switch app.SchStoichYMode.Value
                    case 'Detected power (nW)'; yt = totP;
                    case 'Fraction of channel (%)'; yt = NaN;   % total is 100% -> no line
                    otherwise; yt = totE;
                end
                hk = app.StoichTotalLines{k};
                if isnan(yt) || yt <= 0
                    if ~isempty(hk) && isvalid(hk); hk.Visible = 'off'; end
                    continue;
                end
                if isempty(hk) || ~isvalid(hk)
                    app.StoichTotalLines{k} = yline(ax, yt, ':', ...
                        sprintf('%s total', D.cnames{k}), 'Color',D.ccolors(k,:), ...
                        'LineWidth',1.4, 'LabelHorizontalAlignment','left', ...
                        'Interpreter','none', 'Tag','stoichMarker');
                else
                    hk.Value = yt; hk.Visible = 'on';
                end
            end
            if ~isempty(app.SchStoichReadout) && isvalid(app.SchStoichReadout)
                app.SchStoichReadout.Text = strjoin(parts, char(10));
            end
        end

        function xl = stoichFluxLog(app, D, i)
            % log10 photon-flux multiplier of fluorophore i = conc knob x the
            % excitation-weighted mix of its sources' illumination knobs.
            cm = 0; if i <= numel(app.StoichConcMul); cm = app.StoichConcMul(i); end
            xl = cm;
            if isfield(D,'W') && ~isempty(D.W) && i <= size(D.W,2)
                w = D.W(:,i); sw = sum(w); im = app.StoichIllumMul(:);
                if sw > 0 && numel(im) == numel(w)
                    eff = sum(w .* 10.^im) / sw;
                    xl = cm + log10(max(eff, 1e-9));
                end
            end
        end

        function cols = stoichColors(app, cfg, nC, nS)
            cols = zeros(max(nC*nS,1),3);
            for k = 1:nC
                o = cfg.assign(k);
                if o >= 1 && o <= numel(cfg.fluors)
                    [~,ix] = max(cfg.fluors(o).em);
                    base = wl2rgb(app.Lambda(ix));
                else
                    base = lines(1);
                end
                for j = 1:nS
                    idx = (k-1)*nS + j;
                    if j == 1
                        cols(idx,:) = base;
                    else
                        [~,pk] = max(FilterEngine.sourceExcitation(cfg.lasers(j), app.Lambda));
                        srcCol = wl2rgb(app.Lambda(pk));
                        cols(idx,:) = 0.65*base + 0.35*srcCol;
                    end
                end
            end
        end

        function updateNoisePlot(app, cfg, ~, ~, outNow)
            ax = app.SchNoiseAxes;
            if isempty(ax) || ~isvalid(ax); return; end
            cla(ax);
            if isempty(cfg.channels)
                title(ax,'Noise budget needs at least one channel'); return;
            end
            vals = [outNow.shot_e(:), outNow.read_e(:), sqrt(outNow.dark_e(:)), ...
                    sqrt(max(outNow.xtalk_e(:),0)), sqrt(max(outNow.exBleed_e(:),0))];
            labs = {'Photon shot','Detector read','Detector dark','Fluor crosstalk','Laser bleed'};
            if isfield(outNow,'afBySrc') && ~isempty(outNow.afBySrc)
                vals = [vals, sqrt(max(outNow.afBySrc',0))]; %#ok<AGROW>
                labs = [labs, strcat(outNow.afNames, ' AF')]; %#ok<AGROW>
            end
            hb = bar(ax, vals, 'grouped');
            grid(ax,'on');
            ax.XTick = 1:numel(cfg.channels);
            ax.XTickLabel = {cfg.channels.name};
            ylabel(ax,'Noise electrons RMS');
            title(ax,'Noise contributors at current stoichiometry');
            legend(ax, hb, labs(1:numel(hb)), 'Location','eastoutside');
        end

        function p = rbackPercent(app)
            % Laser back-reflection R as a percent: 0.5% when the Metrics
            % checkbox is ticked, else 0 (excluded from the photon budget).
            p = 0;
            if ~isempty(app.RbackCheck) && isvalid(app.RbackCheck) && app.RbackCheck.Value
                p = 0.5;
            end
        end

        function v = blockODValue(app)
            v = 6;
            if isempty(app.BlockODField); return; end
            if isnumeric(app.BlockODField)
                v = app.BlockODField;
            elseif isvalid(app.BlockODField)
                v = app.BlockODField.Value;
            end
        end

        function setBlockODValue(app, v)
            if isempty(v) || ~isfinite(v); v = 6; end
            if isempty(app.BlockODField) || isnumeric(app.BlockODField)
                app.BlockODField = v;
            elseif isvalid(app.BlockODField)
                app.BlockODField.Value = v;
            end
        end

        function w = optWeight(app)
            w = 5;
            if ~isempty(app.OptWeight) && isvalid(app.OptWeight)
                w = app.OptWeight.Value;
            end
        end

        function w = optBleedWeight(app)
            w = 5;
            if ~isempty(app.OptBleedWeight) && isvalid(app.OptBleedWeight)
                w = app.OptBleedWeight.Value;
            end
        end

        function [fluors, lasers, channels, assign, Rback, blockOD] = schematicSystem(app)
            lam = app.Lambda; one = ones(numel(lam),1);
            nf = numel(app.SchFluorVals);
            fluors = struct('name',{},'ex',{},'em',{},'brightness',{});
            for i = 1:nf
                S = getSpec(app, app.SchFluorVals{i}); ex = 0*one; em = 0*one; b = 1;
                if ~isempty(S)
                    ex = S.ex(:); em = S.em(:);
                    if isfield(S,'brightness') && ~isempty(S.brightness) && isfinite(S.brightness); b = S.brightness; end
                end
                fluors(i) = struct('name',app.SchFluorVals{i},'ex',ex,'em',em,'brightness',b); %#ok<AGROW>
            end
            ns = numel(app.ExSrcVals);
            lasers = struct('name',{},'wl',{},'power',{},'spectrum',{},'exFilter',{});
            for i = 1:ns
                lbl = app.ExSrcVals{i}; tok = regexp(lbl,'(\d{3})','tokens','once'); spec = [];
                if contains(lower(lbl),'laser') && ~isempty(tok)
                    wl = str2double(tok{1}); [~,pk] = min(abs(lam-wl));
                else
                    S = getSpec(app, lbl);
                    if ~isempty(S); spec = S.ex(:); [~,pk] = max(S.ex); wl = lam(pk); else; wl = 500; pk = 1; end
                end
                cln = specVec(app, app.ExCleanVals{i}, 'ex', one);
                lasers(i) = struct('name',lbl,'wl',wl,'power',1,'spectrum',spec, ...
                    'exFilter',cln.*combinerGain(app, pk).*primaryExcitationGain(app)); %#ok<AGROW>
            end
            nc = numel(app.SchEmVals);
            Dp = []; pf = NaN;
            if ~strcmp(app.SchPrimaryDD.Value,'(none)')
                Sp = getSpec(app, app.SchPrimaryDD.Value);
                if ~isempty(Sp); Dp = Sp.ex(:); pf = Sp.floor; end
            end
            nSplit = numel(app.SchSplitterVals); splT = cell(1,nSplit); splF = cell(1,nSplit);
            for c = 1:nSplit
                Sc = getSpec(app, app.SchSplitterVals{c});
                if ~isempty(Sc); splT{c} = Sc.ex(:); splF{c} = Sc.floor; else; splT{c} = []; splF{c} = NaN; end
            end
            channels = struct('name',{},'path',{},'emFilter',{},'emFloor',{});
            for k = 1:nc
                path = struct('T',{},'mode',{},'floor',{});
                if ~isempty(Dp); path(end+1) = struct('T',Dp,'mode','T','floor',pf); end %#ok<AGROW>
                for c = 1:(nc-1)
                    if c < k; m = 'T'; elseif c == k; m = 'R'; else; continue; end
                    if isempty(splT{c}); continue; end
                    path(end+1) = struct('T',splT{c},'mode',m,'floor',splF{c}); %#ok<AGROW>
                end
                ef = specVec(app, app.SchEmVals{k}, 'ex', one) .* detVec(app, app.SchDetVals{k});
                eflo = NaN; Sf = getSpec(app, app.SchEmVals{k});
                if ~isempty(Sf) && isfield(Sf,'floor'); eflo = Sf.floor; end
                channels(k) = struct('name',sprintf('Ch%d',k),'path',path, ...
                    'emFilter',ef,'emFloor',eflo); %#ok<AGROW>
            end
            Stmp = FilterEngine.signalMatrix(fluors, lasers, channels, lam, one);
            assign = zeros(1,nc);
            for k = 1:nc; [~,assign(k)] = max(Stmp(:,k)); end
            Rback = rbackPercent(app)/100; blockOD = blockODValue(app);
        end

        function applySchematicToSystem(app)
            % Build the resolved system directly from the schematic (chained
            % combiners/splitters + per-channel detector) and compute results.
            lam = app.Lambda; one = ones(numel(lam),1);

            % ---- fluorophores ----
            nf = numel(app.SchFluorVals);
            fluors = struct('name',{},'ex',{},'em',{},'brightness',{});
            for i = 1:nf
                S = getSpec(app, app.SchFluorVals{i}); ex = 0*one; em = 0*one; b = 1;
                if ~isempty(S)
                    ex = S.ex(:); em = S.em(:);
                    if isfield(S,'brightness') && ~isempty(S.brightness) && isfinite(S.brightness); b = S.brightness; end
                end
                fluors(i) = struct('name',app.SchFluorVals{i},'ex',ex,'em',em,'brightness',b);
            end
            fn = {fluors.name};

            % ---- sources, with the combiner chain folded into each exFilter ----
            ns = numel(app.ExSrcVals);
            lasers = struct('name',{},'wl',{},'power',{},'spectrum',{},'exFilter',{});
            for i = 1:ns
                lbl = app.ExSrcVals{i}; tok = regexp(lbl,'(\d{3})','tokens','once'); spec = [];
                if contains(lower(lbl),'laser') && ~isempty(tok)
                    wl = str2double(tok{1}); [~,pk] = min(abs(lam-wl));
                else
                    S = getSpec(app, lbl);
                    if ~isempty(S); spec = S.ex(:); [~,pk] = max(S.ex); wl = lam(pk); else; wl = 500; pk = 1; end
                end
                cln = specVec(app, app.ExCleanVals{i}, 'ex', one);
                lasers(i) = struct('name',lbl,'wl',wl,'power',1,'spectrum',spec, ...
                    'exFilter',cln.*combinerGain(app, pk).*primaryExcitationGain(app));
            end

            % ---- channels: chained splitter paths + per-channel detector ----
            nc = numel(app.SchEmVals);
            Dp = []; pf = NaN;
            if ~strcmp(app.SchPrimaryDD.Value,'(none)')
                Sp = getSpec(app, app.SchPrimaryDD.Value);
                if ~isempty(Sp); Dp = Sp.ex(:); pf = Sp.floor; end
            end
            nSplit = numel(app.SchSplitterVals); splT = cell(1,nSplit); splF = cell(1,nSplit);
            for c = 1:nSplit
                Sc = getSpec(app, app.SchSplitterVals{c});
                if ~isempty(Sc); splT{c} = Sc.ex(:); splF{c} = Sc.floor; else; splT{c} = []; splF{c} = NaN; end
            end
            channels = struct('name',{},'path',{},'emFilter',{},'emFloor',{},'detector',{});
            assign = zeros(1,nc);
            for k = 1:nc
                path = struct('T',{},'mode',{},'floor',{});
                if ~isempty(Dp); path(end+1) = struct('T',Dp,'mode','T','floor',pf); end %#ok<AGROW>
                for c = 1:(nc-1)
                    if c < k; m = 'T'; elseif c == k; m = 'R'; else; continue; end
                    if isempty(splT{c}); continue; end
                    path(end+1) = struct('T',splT{c},'mode',m,'floor',splF{c}); %#ok<AGROW>
                end
                ef = specVec(app, app.SchEmVals{k}, 'ex', one); eflo = NaN;
                Sf = getSpec(app, app.SchEmVals{k});
                if ~isempty(Sf) && isfield(Sf,'floor'); eflo = Sf.floor; end
                channels(k).name = sprintf('Ch%d',k);
                channels(k).path = path;
                channels(k).emFilter = ef;
                channels(k).emFloor = eflo;
                channels(k).detector = detVec(app, app.SchDetVals{k});   % per-channel QE
            end

            % auto-assign each channel's owner = fluorophore with the most signal
            Rback = rbackPercent(app)/100; blockOD = blockODValue(app);
            Stmp = FilterEngine.signalMatrix(fluors, lasers, channels, lam, one);
            assign = zeros(1,nc);
            for k = 1:nc; [~,assign(k)] = max(Stmp(:,k)); end
            displayResults(app, fluors, lasers, channels, assign, one, Rback, blockOD);
            setStatus(app, sprintf(['Applied schematic: %d source(s), %d fluorophore(s), ' ...
                '%d channel(s) — chained combiners/splitters + per-channel detectors.'], ns, nf, nc));
        end

        function [y, raw] = excSourceCurve(app, srcLabel, cleanLabel)
            % raw  = normalised source spectrum (laser line -> narrow Gaussian,
            %        lamp -> library spectrum); y = raw after the cleanup filter.
            lam = app.Lambda; raw = zeros(numel(lam),1);
            tok = regexp(srcLabel,'(\d{3})','tokens','once');
            if contains(lower(srcLabel),'laser') && ~isempty(tok)
                wl = str2double(tok{1});
                raw = exp(-0.5*((lam - wl)/2).^2);     % 2 nm display linewidth
            else
                S = getSpec(app, srcLabel);
                if ~isempty(S) && ~isempty(S.ex)
                    raw = S.ex(:); m = max(raw); if m > 0; raw = raw/m; end
                end
            end
            y = raw;
            if nargin >= 3 && ~strcmp(cleanLabel,'(none)')
                F = getSpec(app, cleanLabel);
                if ~isempty(F); y = raw .* F.ex(:); end
            end
        end

        function buildSystemTab(app)
            t = uipanel(app.Fig,'Visible','off');   % retired tab; widgets kept off-screen
            g = uigridlayout(t,[3 2]);
            g.RowHeight = {150,'1x',40}; g.ColumnWidth = {'1x','1x'};

            % Fluorophores
            p1 = uipanel(g,'Title','Fluorophores (Name, Brightness=EC*QY/1000)');
            p1.Layout.Row = 1; p1.Layout.Column = 1;
            g1 = uigridlayout(p1,[2 1]); g1.RowHeight = {'1x',28};
            app.FluorTable = uitable(g1,'ColumnName',{'Fluorophore','Brightness'}, ...
                'ColumnEditable',[true true],'Data',cell(0,2));
            gb = uigridlayout(g1,[1 2]);
            uibutton(gb,'Text','Add row','ButtonPushedFcn',@(s,e)addRow(app,app.FluorTable,{'',1}));
            uibutton(gb,'Text','Remove last','ButtonPushedFcn',@(s,e)delRow(app,app.FluorTable));

            % Lasers
            p2 = uipanel(g,'Title','Excitation sources (laser line or lamp + cleanup filter)');
            p2.Layout.Row = 1; p2.Layout.Column = 2;
            g2 = uigridlayout(p2,[2 1]); g2.RowHeight = {'1x',28};
            app.LaserTable = uitable(g2, ...
                'ColumnName',{'Name','Center (nm)','Power','Source','Excitation filter'}, ...
                'ColumnEditable',[true true true true true],'Data',cell(0,5));
            gb2 = uigridlayout(g2,[1 2]);
            uibutton(gb2,'Text','Add row','ButtonPushedFcn',@(s,e)addRow(app,app.LaserTable,{'',488,1,'(laser line)','(none)'}));
            uibutton(gb2,'Text','Remove last','ButtonPushedFcn',@(s,e)delRow(app,app.LaserTable));

            % Channels (now spans both columns; dichroic/detector moved to tab 1)
            p4 = uipanel(g,'Title','Detection channels  (set primary dichroic / detector on the Build tab)');
            p4.Layout.Row = 2; p4.Layout.Column = [1 2];
            g4 = uigridlayout(p4,[2 1]); g4.RowHeight = {'1x',28};
            app.ChanTable = uitable(g4, ...
                'ColumnName',{'Channel','Splitter dichroic','Mode','Emission filter','Owner fluorophore'}, ...
                'ColumnEditable',[true true true true true],'Data',cell(0,5));
            gb4 = uigridlayout(g4,[1 2]);
            uibutton(gb4,'Text','Add channel','ButtonPushedFcn',@(s,e)addRow(app,app.ChanTable,{'','(none)','none','','' }));
            uibutton(gb4,'Text','Remove last','ButtonPushedFcn',@(s,e)delRow(app,app.ChanTable));

            cb = uibutton(g,'Text','Compute results  ▶','FontWeight','bold', ...
                'ButtonPushedFcn',@(s,e)onCompute(app));
            cb.Layout.Row = 3; cb.Layout.Column = [1 2];
        end

        function t = buildResultsTab(app)
            t = uipanel(app.Fig,'Visible','off');   % retired tab; widgets kept off-screen
            g = uigridlayout(t,[3 2]);
            g.RowHeight = {'1x','1x','1.2x'}; g.ColumnWidth = {'1x','1x'};
            p1 = uipanel(g,'Title','Signal matrix  S(fluor, channel)');
            p1.Layout.Row=1; p1.Layout.Column=1;
            app.SigTable = uitable(uigridlayout(p1,[1 1]));
            p2 = uipanel(g,'Title','Crosstalk matrix %  (col-normalised to owner)');
            p2.Layout.Row=1; p2.Layout.Column=2;
            app.CTTable = uitable(uigridlayout(p2,[1 1]));
            p3 = uipanel(g,'Title','Laser back-reflection background / signal (%) per channel & source');
            p3.Layout.Row=2; p3.Layout.Column=[1 2];
            bg = uigridlayout(p3,[1 3]); bg.ColumnWidth = {'1x',150,150};
            app.BleedTable = uitable(bg);
            uibutton(bg,'Text','Export results → CSV', ...
                'ButtonPushedFcn',@(s,e)onExportResults(app));
            uibutton(bg,'Text','Save config → SNR app', ...
                'ButtonPushedFcn',@(s,e)onSaveConfig(app));
            rp = uipanel(g,'Title','Spectral plot');
            rp.Layout.Row=3; rp.Layout.Column=[1 2];
            rg = uigridlayout(rp,[2 1]);
            rg.RowHeight = {46,'1x'};
            rg.Padding = [10 8 10 8];
            ctl = uigridlayout(rg,[1 10]);
            ctl.ColumnWidth = {45,185,30,90,45,80,45,80,58,'1x'};
            ctl.Padding = [6 4 6 4];
            ctl.ColumnSpacing = 10;
            uilabel(ctl,'Text','Plot:','FontSize',12,'FontWeight','bold');
            app.ResPlotMode = uidropdown(ctl,'Items', ...
                {'Component superposition','Channel overlay','Raw vs filtered'}, ...
                'Value','Component superposition','FontSize',12, ...
                'ValueChangedFcn',@(s,e)redrawResultsPlot(app));
            uilabel(ctl,'Text','Y:','FontSize',12,'FontWeight','bold');
            app.ResYMode = uidropdown(ctl,'Items',{'%T','OD'}, ...
                'Value','%T','FontSize',12,'ValueChangedFcn',@(s,e)redrawResultsPlot(app));
            uilabel(ctl,'Text','X min','FontSize',12);
            app.ResXMinField = uieditfield(ctl,'numeric','Value',450, ...
                'Limits',[350 850],'FontSize',12,'ValueChangedFcn',@(s,e)onXRangeEdit(app));
            uilabel(ctl,'Text','X max','FontSize',12);
            app.ResXMaxField = uieditfield(ctl,'numeric','Value',750, ...
                'Limits',[350 850],'FontSize',12,'ValueChangedFcn',@(s,e)onXRangeEdit(app));
            uilabel(ctl,'Text','Scroll','FontSize',12);
            app.ResXSlider = uislider(ctl,'Limits',[350 850],'Value',600, ...
                'ValueChangingFcn',@(s,e)onXSlider(app,e.Value), ...
                'ValueChangedFcn',@(s,e)onXSlider(app,s.Value));
            app.ResAxes = uiaxes(rg);
            xlabel(app.ResAxes,'Wavelength (nm)'); ylabel(app.ResAxes,'Throughput / emission');
        end

        function buildOptimizerTab(app)
            app.OptCandMap = containers.Map('KeyType','char','ValueType','any');
            t = uitab(app.Tabs,'Title','Optimizer');
            g = uigridlayout(t,[6 2]);
            g.RowHeight = {66,'1x',120,40,40,170}; g.ColumnWidth = {400,'1x'};

            % mode + channel selector (two rows so the long labels fit)
            top = uigridlayout(g,[2 2]); top.Layout.Row=1; top.Layout.Column=1;
            top.ColumnWidth = {80,'1x'}; top.RowHeight = {28,28}; top.Padding=[0 0 0 0];
            uilabel(top,'Text','Mode:');
            app.OptMode = uidropdown(top,'Items', ...
                {'Joint (splitter + all filters)','Single channel'});
            uilabel(top,'Text','Channel:');
            app.OptChanDD = uidropdown(top,'Items',{'(compute first)'}, ...
                'ValueChangedFcn',@(s,e)onOptChanChange(app));

            % candidate emission filters for the selected channel
            p = uipanel(g,'Title','Candidate emission filters for selected channel');
            p.Layout.Row=2; p.Layout.Column=1;
            pg = uigridlayout(p,[2 1]); pg.RowHeight = {'1x',30};
            app.OptCandList = uilistbox(pg,'Multiselect','on','Items',{}, ...
                'ValueChangedFcn',@(s,e)onCandSelChange(app));
            uibutton(pg,'Text','+ Pull bandpass filters from FPbase near owner emission', ...
                'ButtonPushedFcn',@(s,e)onPullCandidates(app));

            % candidate splitter dichroics (joint mode)
            pd = uipanel(g,'Title','Candidate splitter dichroics (joint mode; current always kept)');
            pd.Layout.Row=3; pd.Layout.Column=1;
            pdg = uigridlayout(pd,[1 1]);
            app.OptDichroicList = uilistbox(pdg,'Multiselect','on','Items',{});

            wp = uigridlayout(g,[1 5]); wp.Layout.Row=4; wp.Layout.Column=1;
            wp.ColumnWidth = {135,50,105,50,'1x'}; wp.Padding=[0 0 0 0]; wp.ColumnSpacing=6;
            uilabel(wp,'Text','Crosstalk weight:');
            app.OptWeight = uieditfield(wp,'numeric','Value',5);
            uilabel(wp,'Text','Bleed weight:');
            app.OptBleedWeight = uieditfield(wp,'numeric','Value',5);
            app.OptExCheck = uicheckbox(wp,'Text','Opt. exc. filters');

            rb = uibutton(g,'Text','Run optimizer  ▶','FontWeight','bold', ...
                'ButtonPushedFcn',@(s,e)onOptimize(app));
            rb.Layout.Row=5; rb.Layout.Column=1;

            p2 = uipanel(g,'Title','Ranked sets (best first)');
            p2.Layout.Row=[1 6]; p2.Layout.Column=2;
            p2g = uigridlayout(p2,[2 1]); p2g.RowHeight = {'1x',30};
            app.OptResTable = uitable(p2g);
            uibutton(p2g,'Text','Export ranking → CSV', ...
                'ButtonPushedFcn',@(s,e)onExportOpt(app));

            sp = uipanel(g,'Title','Ranking score'); sp.Layout.Row=6; sp.Layout.Column=1;
            spg = uigridlayout(sp,[4 4]); spg.RowHeight={30,30,30,30};
            spg.ColumnWidth={95,'1x',95,'1x'}; spg.Padding=[6 6 6 6]; spg.ColumnSpacing=6;
            uilabel(spg,'Text','Score by:');
            app.OptScoreMode = uidropdown(spg,'Items', ...
                {'Optical figure of merit','Electron-domain SNR (min channel)'});
            app.OptScoreMode.Layout.Column=[2 4];
            app.OptSNRFields.power = labeledNumO(spg,'Power mW',1);
            app.OptSNRFields.integ = labeledNumO(spg,'Integ ms',10);
            app.OptSNRFields.conc  = labeledNumO(spg,'Conc µM',10);
            app.OptSNRFields.read  = labeledNumO(spg,'Read e-',1.5);
            app.OptSNRFields.NA    = labeledNumO(spg,'NA',0.5);
            app.OptAFcheck = uicheckbox(spg,'Text','+tissue/fiber AF');
            app.OptAFcheck.Layout.Column=[3 4];
        end

        function buildCompareTab(app)
            % Side-by-side comparison of filter sets. Each "→ Compare" click on
            % the Schematic tab adds one column: the top table lists the chosen
            % element per category (Fluorophores, Excitation, dichroics, Emission
            % filters, Detectors); the bottom table lists the physics metrics.
            t = uitab(app.Tabs,'Title','Compare');
            g = uigridlayout(t,[3 1]); g.RowHeight = {30,'1.1x','1x'};
            g.Padding = [10 10 10 10]; g.RowSpacing = 8;

            ctl = uigridlayout(g,[1 4]); ctl.Layout.Row = 1;
            ctl.ColumnWidth = {'1x',150,150,160}; ctl.Padding=[0 0 0 0]; ctl.ColumnSpacing=8;
            uilabel(ctl,'Text', ['Click "Compare" on the Design tab to add the current ' ...
                'filter set as a column.'], 'FontColor',[0.3 0.3 0.3]);
            uibutton(ctl,'Text','Add current schematic', ...
                'ButtonPushedFcn',@(s,e)onCompareCapture(app));
            uibutton(ctl,'Text','Remove last','ButtonPushedFcn',@(s,e)onCompareRemoveLast(app));
            uibutton(ctl,'Text','Clear all','ButtonPushedFcn',@(s,e)onCompareClear(app));

            % Keep Compare headers visually stronger than the field contents.
            % The table itself stays large/bold for row/column headers; data
            % cells get normal-weight uistyle in refreshCompareTables.
            p1 = uipanel(g,'Title','Filter set composition  (rows = element, columns = filter set)');
            p1.Layout.Row = 2;
            app.CompareCompTable = uitable(uigridlayout(p1,[1 1]), ...
                'RowName',{}, 'ColumnName',{}, 'Data',cell(0,0), ...
                'FontSize',14, 'FontWeight','bold');

            p2 = uipanel(g,'Title','Physics metrics  (rows = metric, columns = filter set)');
            p2.Layout.Row = 3;
            cg = uigridlayout(p2,[2 1]); cg.RowHeight = {'1x',30};
            app.CompareMetricTable = uitable(cg, ...
                'RowName',{}, 'ColumnName',{}, 'Data',cell(0,0), ...
                'FontSize',14, 'FontWeight','bold');
            uibutton(cg,'Text','Export comparison → CSV', ...
                'ButtonPushedFcn',@(s,e)onCompareExport(app));
        end

        function onCompareCapture(app)
            try
                name = sprintf('Set %d', numel(app.CompareEntries)+1);
                comp = schematicComposition(app, name);
                metrics = gatherMetrics(app);
            catch ME
                setStatus(app, ['Compare: ' ME.message], true); return;
            end
            e = struct('name',name,'comp',{comp},'metrics',{metrics});
            if isempty(app.CompareEntries); app.CompareEntries = e;
            else; app.CompareEntries(end+1) = e; end
            refreshCompareTables(app);
            for tt = app.Tabs.Children'
                if strcmp(tt.Title,'Compare'); app.Tabs.SelectedTab = tt; break; end
            end
            setStatus(app, sprintf('Added "%s" to Compare (%d filter set(s)).', ...
                name, numel(app.CompareEntries)));
        end

        function onCompareRemoveLast(app)
            if isempty(app.CompareEntries); return; end
            app.CompareEntries(end) = [];
            refreshCompareTables(app);
            setStatus(app, sprintf('Removed last Compare entry (%d remaining).', ...
                numel(app.CompareEntries)));
        end

        function onCompareClear(app)
            app.CompareEntries = struct('name',{},'comp',{},'metrics',{});
            refreshCompareTables(app);
            setStatus(app, 'Cleared Compare.');
        end

        function comp = schematicComposition(app, name)
            % Cell {label, value(char)} describing every element of the current
            % schematic, in canonical category order.
            comp = {'Filter set name', name};
            for i = 1:numel(app.SchFluorVals)
                comp(end+1,:) = {sprintf('Fluorophore %d',i), app.SchFluorVals{i}}; %#ok<AGROW>
            end
            for i = 1:numel(app.ExSrcVals)
                comp(end+1,:) = {sprintf('Light source %d',i), app.ExSrcVals{i}}; %#ok<AGROW>
                comp(end+1,:) = {sprintf('Cleanup filter %d',i), app.ExCleanVals{i}}; %#ok<AGROW>
            end
            for c = 1:numel(app.ExCombinerVals)
                comp(end+1,:) = {sprintf('Combiner dichroic %d',c), app.ExCombinerVals{c}}; %#ok<AGROW>
            end
            comp(end+1,:) = {'Primary dichroic', app.SchPrimaryDD.Value};
            for c = 1:numel(app.SchSplitterVals)
                comp(end+1,:) = {sprintf('Splitter dichroic %d',c), app.SchSplitterVals{c}}; %#ok<AGROW>
            end
            for k = 1:numel(app.SchEmVals)
                comp(end+1,:) = {sprintf('Emission filter %d',k), app.SchEmVals{k}}; %#ok<AGROW>
                comp(end+1,:) = {sprintf('Detector %d',k), app.SchDetVals{k}}; %#ok<AGROW>
            end
            % (Back-reflection R and element blocking OD are SNR-model inputs, not
            % part of the optical filter-set composition, so they are not listed.)
        end

        function M = gatherMetrics(app)
            % Cell {label, value} of comparison metrics for the current schematic.
            % Combines the optical figure of merit / crosstalk / laser bleed with
            % the electron-domain SNR model used by the Physics metrics tab.
            lam = app.Lambda; one = ones(numel(lam),1);
            [fluors, lasers, channels, assign, Rback, blockOD] = schematicSystem(app);
            if isempty(fluors) || isempty(lasers) || isempty(channels)
                error('Add at least one fluorophore, source, and channel.');
            end
            nC = numel(channels); nS = numel(lasers);
            % --- optical domain (detector QE already folded into emFilter) ---
            [S, eff] = FilterEngine.signalMatrix(fluors, lasers, channels, lam, one);
            CT = FilterEngine.crosstalkMatrix(S, assign);
            bleed = FilterEngine.laserBleed(channels, lasers, lam, Rback, one, blockOD);
            score = FilterEngine.systemScore(S, assign, optWeight(app), ...
                bleed, optBleedWeight(app));
            xtPct = zeros(1,nC); bsr = zeros(1,nC); effCh = zeros(1,nC);
            for k = 1:nC
                col = CT(:,k); col(assign(k)) = 0; xtPct(k) = 100*sum(col);
                bsr(k) = 100*sum(bleed(k,:)) / (S(assign(k),k) + eps);
                effCh(k) = 100*eff(assign(k),k);
            end
            % --- electron domain (same operating point as Physics metrics) ---
            pSrc_mW = metricPowerPerSource(app);
            phys = metricPhys(app, nS, pSrc_mW);
            conc = max(1, app.SchMetricFields.conc.Value) * 1e-9;
            fp = metricPhotophysics(app, fluors, conc*ones(1,numel(fluors)));
            cfg = struct('lambda',lam,'fluors',fluors,'lasers',lasers,'channels',channels, ...
                'assign',assign,'detector',one,'Rback',Rback,'blockOD',blockOD);
            out = snrModel(cfg, phys, fp, metricAutofluor(app));
            afTot = 0;
            if isfield(out,'afBySrc') && ~isempty(out.afBySrc); afTot = sum(out.afBySrc(:)); end
            sig3 = @(x) round(x,3,'significant');
            M = { ...
                'Optical figure of merit (max)', sig3(score); ...
                'Mean collection eff. % (max)',  sig3(mean(effCh)); ...
                'Mean signal e-/frame (max)',    sig3(mean(out.signal_e)); ...
                'Min channel SNR (max)',         sig3(min(out.SNR)); ...
                'Mean channel SNR (max)',        sig3(mean(out.SNR)); ...
                'Max crosstalk % (min)',         sig3(max(xtPct)); ...
                'Mean crosstalk % (min)',        sig3(mean(xtPct)); ...
                'Max laser bleed % (min)',       sig3(max(bsr)); ...
                'Mean laser bleed % (min)',      sig3(mean(bsr)); ...
                'Total autofluor. e- (min)',     sig3(afTot); ...
                'Sources / fluorophores / channels', sprintf('%d / %d / %d', nS, numel(fluors), nC)};
        end

        function refreshCompareTables(app)
            E = app.CompareEntries; n = numel(E);
            if n == 0
                app.CompareCompTable.RowName = {}; app.CompareCompTable.ColumnName = {};
                app.CompareCompTable.Data = cell(0,0);
                app.CompareMetricTable.RowName = {}; app.CompareMetricTable.ColumnName = {};
                app.CompareMetricTable.Data = cell(0,0);
                styleCompareTableContents(app, app.CompareCompTable);
                styleCompareTableContents(app, app.CompareMetricTable);
                return;
            end
            names = {E.name};
            % --- composition (union of element labels, canonical order) ---
            labels = {};
            for e = 1:n; labels = [labels, E(e).comp(:,1)']; end %#ok<AGROW>
            labels = compSortLabels(uniqueStableC(labels));
            D = repmat({''}, numel(labels), n);
            for e = 1:n
                for r = 1:size(E(e).comp,1)
                    ri = find(strcmp(labels, E(e).comp{r,1}), 1);
                    D{ri,e} = E(e).comp{r,2};
                end
            end
            % labels live in a first DATA column (not RowName) so they can be
            % right-justified; the value columns are left-justified.
            app.CompareCompTable.RowName = {};
            app.CompareCompTable.ColumnName = [{'Element'}, names];
            app.CompareCompTable.Data = [labels(:), D];
            app.CompareCompTable.ColumnWidth = [{240}, repmat({130}, 1, n)];
            styleCompareTableContents(app, app.CompareCompTable);
            % --- metrics (labels stable & shared; union in entry order) ---
            mlabels = E(1).metrics(:,1)';
            for e = 2:n; mlabels = uniqueStableC([mlabels, E(e).metrics(:,1)']); end
            MD = repmat({[]}, numel(mlabels), n);
            for e = 1:n
                for r = 1:size(E(e).metrics,1)
                    ri = find(strcmp(mlabels, E(e).metrics{r,1}), 1);
                    MD{ri,e} = E(e).metrics{r,2};
                end
            end
            app.CompareMetricTable.RowName = {};
            app.CompareMetricTable.ColumnName = [{'Metric'}, names];
            app.CompareMetricTable.Data = [mlabels(:), MD];
            app.CompareMetricTable.ColumnWidth = [{240}, repmat({130}, 1, n)];
            styleCompareTableContents(app, app.CompareMetricTable);
        end

        function styleCompareTableContents(~, tbl)
            try
                removeStyle(tbl);
            catch
            end
            try
                [nr,nc] = size(tbl.Data);
                if nr == 0 || nc == 0; return; end
                % column 1 = row labels: bold, RIGHT-justified
                addStyle(tbl, uistyle('FontWeight','bold','HorizontalAlignment','right'), ...
                    'column', 1);
                % value columns: normal weight, LEFT-justified (numbers no longer
                % pushed to the far right edge)
                if nc >= 2
                    addStyle(tbl, uistyle('FontWeight','normal','HorizontalAlignment','left'), ...
                        'column', 2:nc);
                end
            catch
                % Older MATLAB releases may not support these table styles;
                % leave the table usable rather than failing refresh.
            end
        end

        function onCompareExport(app)
            if isempty(app.CompareEntries)
                setStatus(app,'Nothing to export — add a filter set first.',true); return;
            end
            [f,p] = uiputfile('*.csv','Export comparison', 'filter_set_comparison.csv');
            if isequal(f,0); return; end
            fid = fopen(fullfile(p,f),'w'); if fid<0; setStatus(app,'Cannot write file.',true); return; end
            cleaner = onCleanup(@()fclose(fid));
            names = {app.CompareEntries.name};
            % labels now live in Data column 1, values in columns 2..end
            writeCsvRow(fid, ['Composition', names]);
            ct = app.CompareCompTable;
            for r = 1:size(ct.Data,1)
                writeCsvRow(fid, ct.Data(r,:));
            end
            fprintf(fid,'\n');
            writeCsvRow(fid, ['Physics metric', names]);
            mt = app.CompareMetricTable;
            for r = 1:size(mt.Data,1)
                row = cell(1,size(mt.Data,2));
                for c = 1:numel(row); row{c} = num2str(mt.Data{r,c}); end
                writeCsvRow(fid, row);
            end
            setStatus(app, sprintf('Exported comparison to %s', fullfile(p,f)));
        end
    end

    %% ---------------- Library handling ----------------
    methods (Access = private)
        function scanLibrary(app)
            cats = {'Proteins','Filters','Dichroics','Illumations','Illuminations','BFMConfig','Detectors'};
            app.Lib = struct('name',{},'category',{},'kind',{},'ex',{},'em',{}, ...
                             'brightness',{},'file',{},'floor',{});
            root = spectraDir(app);
            for c = 1:numel(cats)
                d = fullfile(root, cats{c});
                if ~isfolder(d); continue; end
                files = [dir(fullfile(d,'*.txt')); dir(fullfile(d,'*.csv')); ...
                         dir(fullfile(d,'*.TXT'))];
                % case-insensitive filesystems (Windows) match *.txt and *.TXT
                % to the same files — de-duplicate so each is loaded once.
                [~, uidx] = unique(lower({files.name}), 'stable');
                files = files(uidx);
                for f = 1:numel(files)
                    fp = fullfile(files(f).folder, files(f).name);
                    if strcmpi(cats{c},'Illumations') && strcmpi(files(f).name,'Prizmatix.csv')
                        continue; % source bundle parsed into individual UHP-T files
                    end
                    try
                        S = loadSpectrum(fp, app.Lambda);
                        fl = NaN; if isfield(S,'floor'); fl = S.floor; end
                        app.Lib(end+1) = struct('name',S.name,'category',cats{c}, ...
                            'kind',S.kind,'ex',S.ex,'em',S.em,'brightness',[],'file',fp,'floor',fl); %#ok<AGROW>
                    catch ME
                        warning('Skip %s: %s', files(f).name, ME.message);
                    end
                end
            end
            refreshLibTree(app);
            populateChoosers(app);
            setStatus(app, sprintf('Loaded %d spectra from %s', numel(app.Lib), root));
        end

        function C = correlationMatrix(app)
            % Zero-lag normalised cross-correlation between every pair of library
            % spectra of the same kind. Filters compare transmission; fluorophores
            % compare excitation and emission and keep the lower of the two peaks.
            n = numel(app.Lib); C = nan(n);
            isFlu = strcmp({app.Lib.kind},'fluorophore');
            for i = 1:n
                C(i,i) = 1;
                for j = i+1:n
                    if isFlu(i) ~= isFlu(j); continue; end
                    if isFlu(i)
                        cx = localPeakXcorr(app.Lib(i).ex, app.Lib(j).ex);
                        ce = localPeakXcorr(app.Lib(i).em, app.Lib(j).em);
                        C(i,j) = min(cx, ce);
                    else
                        C(i,j) = localPeakXcorr(app.Lib(i).ex, app.Lib(j).ex);
                    end
                    C(j,i) = C(i,j);
                end
            end
        end

        function onDedup(app)
            thr = 0.999;
            if ~isempty(app.DedupFig) && isvalid(app.DedupFig); close(app.DedupFig); end
            C = correlationMatrix(app);      % reassessed fresh on every click
            n = size(C,1);
            % group near-identical entries (union-find style by threshold)
            grp = 1:n;
            for i = 1:n
                for j = i+1:n
                    if isfinite(C(i,j)) && C(i,j) >= thr; grp(grp==grp(j)) = grp(i); end
                end
            end
            [~,~,gid] = unique(grp,'stable');
            dupGroups = {};
            for g = 1:max(gid)
                members = find(gid==g);
                if numel(members) > 1; dupGroups{end+1} = members; end %#ok<AGROW>
            end
            if isempty(dupGroups)
                setStatus(app, sprintf('No duplicates above zero-lag xcorr %.4f.', thr)); return;
            end
            showDedupDialog(app, C, dupGroups, thr);
        end

        function showDedupDialog(app, C, dupGroups, thr)
            f = uifigure('Name','Find & Remove Duplicate Spectra (xcorr matrix)','Position',[200 150 1080 580]);
            app.DedupFig = f;
            g = uigridlayout(f,[2 2]); g.RowHeight={'1x',42}; g.ColumnWidth={420,'1x'};
            ax = uiaxes(g); ax.Layout.Row=1; ax.Layout.Column=1;
            imagesc(ax, C, 'AlphaData', ~isnan(C)); axis(ax,'tight'); ax.YDir='reverse';
            colorbar(ax); clim(ax,[max(0.9,thr-0.01) 1]); title(ax,'zero-lag normalised xcorr matrix');
            xlabel(ax,'spectrum #'); ylabel(ax,'spectrum #');
            % highlight the pairs that are actually counted as duplicates (>=thr)
            hold(ax,'on');
            [gi,gj] = find(triu(C,1) >= thr);
            if ~isempty(gi)
                plot(ax, [gj;gi], [gi;gj], 'rs', 'MarkerSize',7, 'LineWidth',1.3);
            end
            hold(ax,'off');
            % Build rows with one editable "Keep" checkbox per duplicate entry.
            % Default keep still prefers non-BFMConfig, but the user can change it.
            rows = {};
            for q = 1:numel(dupGroups)
                m = dupGroups{q}(:)';                 % row vector -> iterate per element
                cats = {app.Lib(m).category};
                kpref = find(~strcmp(cats,'BFMConfig'),1); if isempty(kpref); kpref=1; end
                keep = m(kpref);
                groupC = C(m,m); groupC(logical(eye(numel(m)))) = NaN;
                gmax = round(max(groupC(:),[],'omitnan'),5);
                for r = m
                    rows(end+1,:) = { r==keep, q, r, app.Lib(r).category, ...
                        app.Lib(r).name, gmax }; %#ok<AGROW>
                end
            end
            t = uitable(g); t.Layout.Row=1; t.Layout.Column=2;
            t.ColumnName={'Keep','Group','#','Category','Spectrum','Group max xcorr'};
            t.ColumnEditable = [true false false false false false];
            t.ColumnFormat = {'logical','numeric','numeric','char','char','numeric'};
            t.Data=rows;
            t.CellEditCallback = @(s,e) enforceOneKeep(e);
            lbl = uilabel(g,'Text',sprintf(['%d duplicate groups found at zero-lag xcorr >= %.4f. ' ...
                'Check exactly one Keep box per group; unchecked duplicates are removed from this app library view.'],numel(dupGroups),thr), ...
                'WordWrap','on');
            lbl.Layout.Row=2; lbl.Layout.Column=1;
            btn = uibutton(g,'Text','Find & Remove Duplicate Spectra','FontWeight','bold');
            btn.Layout.Row=2; btn.Layout.Column=2;
            btn.ButtonPushedFcn = @(s,e) doRemove();
            bumpFonts(app, f, 2);
            function enforceOneKeep(evt)
                if isempty(evt) || evt.Indices(2) ~= 1; return; end
                data = t.Data;
                row = evt.Indices(1);
                gnum = data{row,2};
                same = find(cellfun(@(x) isequal(x,gnum), data(:,2)));
                if isequal(data{row,1}, true)
                    for rr = same(:)'; data{rr,1} = false; end
                    data{row,1} = true;
                else
                    hasKeep = any(cellfun(@(x) isequal(x,true), data(same,1)));
                    if ~hasKeep
                        data{row,1} = true;
                    end
                end
                t.Data = data;
            end
            function doRemove()
                data = t.Data;
                keepIdx = zeros(1,numel(dupGroups));
                for qq = 1:numel(dupGroups)
                    rowsQ = find(cellfun(@(x) isequal(x,qq), data(:,2)));
                    keepRows = rowsQ(cellfun(@(x) isequal(x,true), data(rowsQ,1)));
                    if numel(keepRows) ~= 1
                        uialert(f, sprintf('Choose exactly one spectrum to keep in group %d.', qq), ...
                            'Find & Remove Duplicate Spectra');
                        return;
                    end
                    keepIdx(qq) = data{keepRows(1),3};
                end
                removeMask = false(1,numel(app.Lib));
                for qq = 1:numel(dupGroups)
                    m = dupGroups{qq}; removeMask(m) = true; removeMask(keepIdx(qq)) = false;
                end
                nrem = nnz(removeMask);
                app.Lib(removeMask) = [];
                refreshLibTree(app); populateChoosers(app);
                setStatus(app, sprintf('Cleaned up %d duplicate spectra from the library view.', nrem));
                close(f);
            end
        end

        function refreshLibTree(app)
            delete(app.LibTree.Children);
            app.LibPlotIdx = [];
            cats = unique({app.Lib.category},'stable');
            for c = 1:numel(cats)
                n = uitreenode(app.LibTree,'Text',cats{c});
                idx = find(strcmp({app.Lib.category},cats{c}));
                for i = idx
                    leaf = uitreenode(n,'Text',app.Lib(i).name);
                    leaf.NodeData = i;
                end
            end
        end

        function onLibSelect(app)
            sel = app.LibTree.SelectedNodes;
            if isempty(sel) || isempty(sel.NodeData); return; end
            if isempty(app.LibPlotIdx)
                app.LibPlotIdx = sel.NodeData;
                redrawLibPlot(app);
            end
        end

        function plotSelectedLibrary(app, appendMode)
            sel = app.LibTree.SelectedNodes;
            if isempty(sel) || isempty(sel.NodeData); return; end
            idx = sel.NodeData;
            if appendMode
                app.LibPlotIdx = unique([app.LibPlotIdx idx], 'stable');
            else
                app.LibPlotIdx = idx;
            end
            redrawLibPlot(app);
        end

        function removeSelectedLibrary(app)
            % Remove the tree-selected spectrum from the overlay (if present).
            sel = app.LibTree.SelectedNodes;
            if isempty(sel) || isempty(sel.NodeData)
                setStatus(app,'Select a spectrum in the tree to remove from the plot.',true); return;
            end
            app.LibPlotIdx = app.LibPlotIdx(app.LibPlotIdx ~= sel.NodeData);
            redrawLibPlot(app);
        end

        function clearLibraryPlot(app)
            app.LibPlotIdx = [];
            cla(app.LibAxes);
            xlabel(app.LibAxes,'Wavelength (nm)');
            ylabel(app.LibAxes,'%T / normalised intensity (%)');
            app.LibAxes.YDir = 'normal';
            title(app.LibAxes,'Select a spectrum');
        end

        function redrawLibPlot(app)
            idx = app.LibPlotIdx;
            idx = idx(idx >= 1 & idx <= numel(app.Lib));
            if isempty(idx); clearLibraryPlot(app); return; end
            cla(app.LibAxes); hold(app.LibAxes,'on');
            labels = {};
            plottedY = [];
            for q = 1:numel(idx)
                S = app.Lib(idx(q));
                col = wl2rgb(libPeakWavelength(app, S));   % colour by peak wavelength
                if strcmp(S.kind,'fluorophore')
                    yEx = yPlot(app,S.ex,app.LibYMode.Value);
                    yEm = yPlot(app,S.em,app.LibYMode.Value);
                    plot(app.LibAxes, app.Lambda, yEx, ...
                        '--','Color',col,'LineWidth',1.3);   % excitation dashed
                    plot(app.LibAxes, app.Lambda, yEm, ...
                        '-','Color',col,'LineWidth',1.8);    % emission solid
                    plottedY = [plottedY; yEx; yEm]; %#ok<AGROW>
                    labels = [labels, {[S.name ' ex'], [S.name ' em']}]; %#ok<AGROW>
                else
                    isDichroic = strcmpi(getfield_safe(S,'category',''),'Dichroics');
                    sty = ternStr(isDichroic, '--', '-');    % dichroics dashed
                    y = yPlot(app,S.ex,app.LibYMode.Value);
                    plot(app.LibAxes, app.Lambda, y, ...
                        sty,'Color',col,'LineWidth',1.7);
                    plottedY = [plottedY; y]; %#ok<AGROW>
                    labels{end+1} = S.name; %#ok<AGROW>
                end
            end
            hold(app.LibAxes,'off');
            if strcmp(app.LibYMode.Value,'OD')
                ylabel(app.LibAxes,'OD / -log10(normalised)');
                odMax = max(plottedY(isfinite(plottedY)));
                if isempty(odMax) || odMax <= 0; odMax = 2; end
                ylim(app.LibAxes,[0 min(8, max(2, ceil(odMax)))]);
                app.LibAxes.YDir = 'reverse';
            else
                ylabel(app.LibAxes,'%T / normalised intensity (%)');
                ylim(app.LibAxes,[0 105]);
                app.LibAxes.YDir = 'normal';
            end
            xlabel(app.LibAxes,'Wavelength (nm)');
            xlim(app.LibAxes,[app.Lambda(1) app.Lambda(end)]);
            title(app.LibAxes, sprintf('%d spectrum entry(s)', numel(idx)), 'Interpreter','none');
            legend(app.LibAxes, labels, 'Location','northeastoutside','Interpreter','none');
            grid(app.LibAxes,'on');
        end

        function wl = libPeakWavelength(app, S)
            % Representative wavelength for colouring a library curve: a
            % fluorophore's emission peak, otherwise the transmission centroid.
            lam = app.Lambda(:);
            if strcmp(S.kind,'fluorophore') && isfield(S,'em') && ~isempty(S.em)
                [~,ix] = max(S.em(:)); wl = lam(ix); return;
            end
            t = S.ex(:); t(~isfinite(t)) = 0; t = max(t,0);
            if sum(t) > 0
                wl = sum(lam .* t) / sum(t);          % transmission centroid
            else
                [~,ix] = max(t); wl = lam(ix);
            end
        end

        function [cat, sub] = webCatToQuery(app)
            % Map the friendly dropdown to FPbase (category, subtype).
            switch app.WebCat.Value
                case 'Fluorophore';        cat = 'P'; sub = '';
                case 'Bandpass filter';    cat = 'F'; sub = 'BP';
                case 'Long/Short-pass';    cat = 'F'; sub = '';   % LP/SP mixed
                case 'Dichroic (BS)';      cat = 'F'; sub = 'BS';
                case 'Light source';       cat = 'L'; sub = '';
                case 'Detector (camera)';  cat = 'C'; sub = '';
                otherwise;                 cat = 'F'; sub = '';    % Any filter
            end
        end

        function onWebSearch(app)
            q = strtrim(app.WebField.Value);
            if isempty(q); return; end
            [cat, sub] = webCatToQuery(app);
            setStatus(app, ['Searching FPbase for "' q '" ...']); drawnow;
            try
                hits = FPbase.search(q, cat, sub);
            catch ME
                setStatus(app, ['FPbase error: ' ME.message], true); return;
            end
            app.WebHits = hits;
            if isempty(hits)
                app.WebResults.Items = {}; setStatus(app,'No matches.',true); return;
            end
            labels = arrayfun(@(h) sprintf('%s  [%s]', h.name, h.subtype), ...
                hits, 'UniformOutput', false);
            app.WebResults.Items = labels;
            app.WebResults.ItemsData = num2cell(1:numel(hits));
            setStatus(app, sprintf('%d matches — select and download.', numel(hits)));
        end

        function onWebDownload(app)
            idx = app.WebResults.Value;
            if isempty(idx); setStatus(app,'Select item(s) first.',true); return; end
            if ~iscell(idx); idx = {idx}; end
            [cat, ~] = webCatToQuery(app);
            n = 0;
            for k = 1:numel(idx)
                h = app.WebHits(idx{k});
                setStatus(app, ['Downloading ' h.name ' ...']); drawnow;
                try
                    if strcmp(cat,'P')
                        % proteins: use REST API to also get ex+em+brightness
                        S = fetchFPbase(h.name, app.Lambda);
                        % persist a formatted .txt in the local DB so the
                        % download survives leaving the app (silent overwrite).
                        outFile = fullfile(spectraDir(app), 'Proteins', [cleanFileName(S.name) '.txt']);
                        writeUserSpectrumFile(app, outFile, S.name, S.ex, S.em, true);
                        addToLib(app, S.name, 'Proteins', 'fluorophore', S.ex, S.em, ...
                            getfield_safe(S,'brightness',[]), outFile);
                    else
                        S = FPbase.spectrum(h.id, app.Lambda);
                        sub = getfield_safe(S,'subtype','');
                        if strcmp(cat,'C');        category = 'Detectors';
                        elseif strcmp(sub,'BS');   category = 'Dichroics';
                        else;                      category = 'Filters'; end
                        outFile = fullfile(spectraDir(app), category, [cleanFileName(S.name) '.txt']);
                        writeUserSpectrumFile(app, outFile, S.name, S.ex, [], true);
                        addToLib(app, S.name, category, 'filter', S.ex, [], [], outFile, ...
                            getfield_safe(S,'floor',NaN));
                    end
                    n = n + 1;
                catch ME
                    setStatus(app, ['Download failed for "' h.name '": ' ME.message], true);
                end
            end
            refreshLibTree(app); populateChoosers(app);
            if n > 0
                setStatus(app, sprintf('Downloaded %d spectrum(s) and saved to %s', ...
                    n, spectraDir(app)));
            end
        end

        function onUserSpectrumDownload(app)
            try
            [folderName, typeLabel, isFluor, normCurve] = userSpectrumTarget(app);
            setStatus(app, sprintf(['Importing user %s. Accepted formats: single spectra use ' ...
                '2 columns (wavelength/value); fluorophores use either one 3-column file ' ...
                '(wavelength/excitation/emission) or two 2-column files.'], lower(typeLabel)));
            drawnow;
            if isFluor
                route = questdlg(sprintf(['Fluorophore format:\n\n' ...
                    'One 3-column file: wavelength / excitation / emission.\n' ...
                    'Two 2-column files: wavelength / excitation, then wavelength / emission.\n\n' ...
                    'The app will resample to 350-850 nm at 1 nm steps.']), ...
                    'Fluorophore import', ...
                    'One 3-column file', 'Two 2-column files', 'Cancel', ...
                    'One 3-column file');
                if isempty(route) || strcmp(route,'Cancel')
                    setStatus(app,'User spectrum import cancelled.',true); return;
                end
                if strcmp(route,'One 3-column file')
                    [inFile, inPath] = uigetfile({'*.txt;*.csv;*.tsv','Spectrum files (*.txt, *.csv, *.tsv)'; '*.*','All files'}, ...
                        'Select fluorophore spectrum (wavelength, excitation, emission)');
                    if isequal(inFile,0); setStatus(app,'User spectrum import cancelled.',true); return; end
                    [~, base] = fileparts(inFile);
                    S = loadSpectrum(fullfile(inPath, inFile), app.Lambda);
                    if ~strcmp(S.kind,'fluorophore') || isempty(S.ex) || isempty(S.em)
                        error('Fluorophore single-file import needs three numeric columns: wavelength, excitation, emission.');
                    end
                    ex = normPeak(S.ex(:));
                    em = normPeak(S.em(:));
                else
                    [exFile, exPath] = uigetfile({'*.txt;*.csv;*.tsv','Spectrum files (*.txt, *.csv, *.tsv)'; '*.*','All files'}, ...
                        'Select excitation spectrum (2 columns)');
                    if isequal(exFile,0); setStatus(app,'User spectrum import cancelled.',true); return; end
                    [emFile, emPath] = uigetfile({'*.txt;*.csv;*.tsv','Spectrum files (*.txt, *.csv, *.tsv)'; '*.*','All files'}, ...
                        'Select emission spectrum (2 columns)');
                    if isequal(emFile,0); setStatus(app,'User spectrum import cancelled.',true); return; end
                    [~, base] = fileparts(exFile);
                    ex = readUserTwoColumnCurve(app, fullfile(exPath, exFile), true);
                    em = readUserTwoColumnCurve(app, fullfile(emPath, emFile), true);
                end
                answer = inputdlg({'Local fluorophore name:'}, 'Save user spectra', [1 48], {base});
                if isempty(answer); return; end
                name = cleanFileName(strtrim(answer{1}));
                if isempty(name); setStatus(app,'Spectrum name cannot be empty.',true); return; end
                outFile = fullfile(spectraDir(app), folderName, [name '.txt']);
                writeUserSpectrumFile(app, outFile, name, ex, em);
            else
                [inFile, inPath] = uigetfile({'*.txt;*.csv;*.tsv','Spectrum files (*.txt, *.csv, *.tsv)'; '*.*','All files'}, ...
                    'Select user spectrum (2 columns)');
                if isequal(inFile,0); setStatus(app,'User spectrum import cancelled.',true); return; end
                [~, base] = fileparts(inFile);
                answer = inputdlg({'Local spectrum name:'}, 'Save user spectrum', [1 48], {base});
                if isempty(answer); return; end
                name = cleanFileName(strtrim(answer{1}));
                if isempty(name); setStatus(app,'Spectrum name cannot be empty.',true); return; end
                y = readUserTwoColumnCurve(app, fullfile(inPath, inFile), normCurve);
                outFile = fullfile(spectraDir(app), folderName, [name '.txt']);
                writeUserSpectrumFile(app, outFile, name, y, []);
            end
            Snew = loadSpectrum(outFile, app.Lambda);
            fl = NaN; if isfield(Snew,'floor'); fl = Snew.floor; end
            addToLib(app, name, folderName, Snew.kind, Snew.ex, Snew.em, [], outFile, fl);
            refreshLibTree(app);
            populateChoosers(app);
            setStatus(app, sprintf('Imported user %s "%s" into %s.', lower(typeLabel), name, outFile));
            catch ME
                setStatus(app, ['User spectrum import failed: ' ME.message], true);
            end
        end

        function [folderName, typeLabel, isFluor, normCurve] = userSpectrumTarget(app)
            isFluor = false; normCurve = false;
            switch app.WebCat.Value
                case 'Fluorophore'
                    folderName = 'Proteins'; typeLabel = 'Fluorophore'; isFluor = true; normCurve = true;
                case {'Bandpass filter','Long/Short-pass','Any filter'}
                    folderName = 'Filters'; typeLabel = 'Filter';
                case 'Dichroic (BS)'
                    folderName = 'Dichroics'; typeLabel = 'Dichroic';
                case 'Light source'
                    folderName = 'Illumations'; typeLabel = 'Light source'; normCurve = true;
                case 'Detector (camera)'
                    folderName = 'Detectors'; typeLabel = 'Detector';
                otherwise
                    folderName = 'Filters'; typeLabel = 'Filter';
            end
        end

        function y = readUserTwoColumnCurve(app, fp, normalizePeak)
            S = loadSpectrum(fp, app.Lambda);
            if strcmp(S.kind,'fluorophore')
                error('User spectra must be exactly two numeric columns: wavelength and value.');
            end
            y = S.ex(:);
            if normalizePeak
                y = normPeak(y);
            end
        end

        function writeUserSpectrumFile(app, outFile, name, y1, y2, silent)
            if nargin < 6; silent = false; end
            outDir = fileparts(outFile);
            if ~isfolder(outDir); mkdir(outDir); end
            if isfile(outFile) && ~silent
                choice = questdlg(sprintf('Replace existing local spectrum "%s"?', name), ...
                    'Confirm replace', 'Replace', 'Cancel', 'Cancel');
                if ~strcmp(choice,'Replace'); error('User cancelled replacing existing spectrum.'); end
            end
            if nargin >= 5 && ~isempty(y2)
                writeSpectrumTable(outFile, app.Lambda, y1, y2);
            else
                writeSpectrumTable(outFile, app.Lambda, y1);
            end
        end

        function addToLib(app, name, category, kind, ex, em, brightness, file, floor)
            if nargin < 9; floor = NaN; end
            % replace if a same-name entry already exists in that category
            old = strcmp({app.Lib.name},name) & strcmp({app.Lib.category},category);
            app.Lib(old) = [];
            app.Lib(end+1) = struct('name',name,'category',category,'kind',kind, ...
                'ex',ex,'em',em,'brightness',brightness,'file',file,'floor',floor);
        end

        function populateChoosers(app)
            filtNames = ['(none)' roleNames(app,'filter')];
            dicNames  = ['(none)' roleNames(app,'dichroic')];
            app.PrimaryDD.Items = keepValue(app.PrimaryDD, dicNames);
            app.ChanTable.ColumnFormat = {'char', dicNames, {'T','R','none'}, filtNames, 'char'};
            srcNames = ['(laser line)' roleNames(app,'source')];
            app.LaserTable.ColumnFormat = {'char','numeric','numeric', srcNames, filtNames};
            if ~isempty(app.OptCandList) && isvalid(app.OptCandList)
                app.OptCandList.Items = roleNames(app,'filter');
            end
            if ~isempty(app.OptDichroicList) && isvalid(app.OptDichroicList)
                app.OptDichroicList.Items = roleNames(app,'dichroic');
            end
            detNames = ['(ideal)' roleNames(app,'detector')];
            app.DetectorDD.Items = keepValue(app.DetectorDD, detNames);
            % schematic tab: rebuild dynamic rows with the full library
            if ~isempty(app.ExSrcArea) && isvalid(app.ExSrcArea)
                if ~isempty(app.SchPrimaryDD) && isvalid(app.SchPrimaryDD)
                    app.SchPrimaryDD.Items = dicNames;
                    if ~ismember(app.SchPrimaryDD.Value, dicNames)
                        app.SchPrimaryDD.Value = '(none)';
                    end
                end
                rebuildExcitationRows(app);
                rebuildFluorRows(app);
                rebuildChannelRows(app);
            end
            function items = keepValue(dd, items)
                if ~isempty(dd.Value) && ~ismember(dd.Value, items)
                    items = [items, {dd.Value}];   % don't drop a current selection
                end
            end
        end

        function names = libNames(app, cats)
            mask = ismember({app.Lib.category}, cats);
            names = {app.Lib(mask).name};
        end

        function names = roleNames(app, role)
            names = {};
            for ii = 1:numel(app.Lib)
                if strcmp(spectrumRole(app, app.Lib(ii)), role)
                    names{end+1} = app.Lib(ii).name; %#ok<AGROW>
                end
            end
            names = unique(names, 'stable');
        end

        function role = spectrumRole(~, S)
            role = 'filter';
            cat = ''; if isfield(S,'category'); cat = lower(S.category); end
            if contains(cat,'protein')
                role = 'fluorophore'; return;
            elseif contains(cat,'dichroic')
                role = 'dichroic'; return;
            elseif contains(cat,'detector')
                role = 'detector'; return;
            elseif contains(cat,'illumin')
                role = 'source'; return;
            elseif contains(cat,'filter')
                role = 'filter'; return;
            elseif isfield(S,'kind') && strcmp(S.kind,'fluorophore')
                role = 'fluorophore'; return;
            end
            fp = ''; if isfield(S,'file') && ~isempty(S.file); fp = S.file; end
            txt = lower(strjoin({S.name, cat, fp}, ' '));
            if contains(txt,'dichroic') || contains(txt,'splitter') || contains(txt,'combiner') || ...
                    ~isempty(regexp(txt, '(^|[_\-\s])di\d|di\d($|[_\-\s])|[_\-\s]di[_\-\s]', 'once'))
                role = 'dichroic';
            elseif contains(txt,'detector') || contains(txt,'sensor') || contains(txt,' qe') || ...
                    contains(txt,'apd') || contains(txt,'orca') || contains(txt,'kinetix') || ...
                    contains(txt,'hamamatsu') || contains(txt,'photometrics')
                role = 'detector';
            elseif contains(txt,'laser') || contains(txt,'led') || contains(txt,'uhp') || ...
                    contains(txt,'prizmatix') || contains(txt,'lumencor') || ...
                    contains(txt,'illumination') || contains(txt,'illumations') || ...
                    ~isempty(regexp(txt, '(^|[_\-\s])m\d{3}', 'once'))
                role = 'source';
            else
                role = 'filter';
            end
        end

        function S = getSpec(app, name)
            i = find(strcmp({app.Lib.name}, name), 1);
            if isempty(i); S = []; else; S = app.Lib(i); end
        end
    end

    %% ---------------- Defaults / table helpers ----------------
    methods (Access = private)
        function seedDefaults(app)
            % A working example mirroring the uSMAART config.
            app.FluorTable.Data = {'GFP',33.5; 'cyOFP',30.4; 'mRuby3',42.9};
            app.LaserTable.Data = {'488 laser',488,1,'(laser line)','(none)'; ...
                                   '561 laser',561,1,'(laser line)','(none)'};
            if any(strcmp(app.PrimaryDD.Items,'Di01-R488_561'))
                app.PrimaryDD.Value = 'Di01-R488_561';
            end
            app.ChanTable.Data = {
                'Green','FF562-Di03','R','FF02-520_28','GFP';
                'Red',  'FF562-Di03','T','FF01_630_92','mRuby3'};
            % schematic primary dichroic default
            if ~isempty(app.SchPrimaryDD) && isvalid(app.SchPrimaryDD) ...
                    && any(strcmp(app.SchPrimaryDD.Items,'Di01-R488_561'))
                app.SchPrimaryDD.Value = 'Di01-R488_561';
            end
            onSchematicUpdate(app);  % splitters/combiners appear when 2nd ch/src added
        end

        function loadStartupFilterSet(app)
            % Open the app on the bundled uSMAART ASAP3 / mRuby3 filter set when
            % present (self-contained folder with its own spectra). Falls back to
            % the seeded defaults if it is missing or fails to load, and reports
            % the outcome in the status bar so the load is visible.
            p = fullfile(filterSetDir(app), 'filterset_uSMAART_ASAP3_mRuby3');
            if ~isfolder(p)
                setStatus(app, ['Startup filter set folder not found: ' p ...
                    '  — using seeded defaults.'], true);
                return;
            end
            try
                [sch, n] = loadSchematicFilterSetFolder(app, p);
                nm = 'filterset_uSMAART_ASAP3_mRuby3';
                if isfield(sch,'name') && ~isempty(sch.name); nm = sch.name; end
                setStatus(app, sprintf('Opened on filter set "%s" (%d spectra resolved).', nm, n));
            catch ME
                setStatus(app, ['Startup filter set load failed (using defaults): ' ME.message], true);
            end
        end

        function addRow(~, tbl, template)
            tbl.Data = [tbl.Data; template];
        end
        function delRow(~, tbl)
            if ~isempty(tbl.Data); tbl.Data(end,:) = []; end
        end

        function onSaveFilterSet(app)
            d = filterSetDir(app);
            defaultName = ['filter_set_' datestr(now,'yyyymmdd_HHMMSS')];
            answer = inputdlg({'Filter set name:'}, 'Save filter set', [1 48], {defaultName});
            if isempty(answer); return; end
            name = cleanFileName(strtrim(answer{1}));
            if isempty(name); setStatus(app,'Filter set name cannot be empty.',true); return; end
            filterSet = collectFilterSetState(app, name);
            sch = filterSetToSchematic(app, filterSet);
            setDir = fullfile(d, name);
            if isfolder(setDir)
                choice = questdlg(sprintf('Replace existing filter set "%s"?', name), ...
                    'Confirm replace', 'Replace', 'Cancel', 'Cancel');
                if ~strcmp(choice,'Replace'); return; end
            end
            if ~isfolder(setDir); mkdir(setDir); end
            clearLegacyFilterSetSpectra(setDir);
            writeSchematicManifest(fullfile(setDir,'filter_set.txt'), sch);
            nRefs = writeSchematicPointerManifest(app, fullfile(setDir,'spectra_pointers.json'), sch);
            setStatus(app, sprintf('Saved filter set "%s" with %d local spectrum pointer(s) to %s', name, nRefs, setDir));
        end

        function onLoadFilterSet(app)
            d = filterSetDir(app);
            p = uigetdir(d, 'Load filter set folder');
            if isequal(p,0); return; end
            try
                [sch,n] = loadSchematicFilterSetFolder(app, p);
                setStatus(app, sprintf('Loaded "%s" and resolved %d spectrum pointer(s).', sch.name, n));
            catch ME
                setStatus(app, ['Load filter set failed: ' ME.message], true);
            end
        end

        function d = filterSetDir(app)
            candidates = { ...
                fullfile(app.LibRoot, 'filterSetApp', 'filterSets'), ...
                fullfile(app.LibRoot, 'filterSetApp', 'FilterSets'), ...
                fullfile(app.LibRoot, 'FilterSetApp', 'filterSets'), ...
                fullfile(app.LibRoot, 'FilterSetApp', 'FilterSets'), ...
                fullfile(app.LibRoot, 'filterSets'), ...
                fullfile(app.LibRoot, 'FilterSets')};
            d = candidates{1};
            for ii = 1:numel(candidates)
                if isfolder(candidates{ii})
                    d = candidates{ii};
                    break;
                end
            end
            if ~isfolder(d); mkdir(d); end
        end

        function d = spectraDir(app)
            candidates = { ...
                fullfile(app.LibRoot, 'spectra'), ...
                fullfile(app.LibRoot, 'Spectra')};
            d = app.LibRoot;
            for ii = 1:numel(candidates)
                if isfolder(candidates{ii})
                    d = candidates{ii};
                    break;
                end
            end
        end

        function s = schematicState(app, name)
            s = struct('version',2,'name',name,'savedOn',datestr(now));
            s.srcVals = app.ExSrcVals;   s.cleanVals = app.ExCleanVals;
            s.combinerVals = app.ExCombinerVals;
            s.primary = app.SchPrimaryDD.Value;
            s.fluorVals = app.SchFluorVals;
            s.emVals = app.SchEmVals; s.detVals = app.SchDetVals;
            s.ownerVals = app.SchOwnerVals; s.splitterVals = app.SchSplitterVals;
            s.Rback = rbackPercent(app); s.blockOD = blockODValue(app);
            s.opticsMode = app.ExOpticsMode.Value;
        end

        function applySchematicState(app, s)
            app.ExSrcVals     = ensureRow(s.srcVals, {'488 nm laser'});
            app.ExCleanVals   = padCell(ensureRow(s.cleanVals,{}), numel(app.ExSrcVals), '(none)');
            app.ExCombinerVals= padCell(ensureRow(getff(s,'combinerVals',{}),{}), numel(app.ExSrcVals)-1, '(none)');
            app.SchFluorVals  = ensureRow(s.fluorVals, {'mNeonGreen'});
            app.SchEmVals     = ensureRow(s.emVals, {'(none)'});
            app.SchDetVals    = padCell(ensureRow(getff(s,'detVals',{}),{}), numel(app.SchEmVals), '(ideal)');
            app.SchOwnerVals  = padCell(ensureRow(getff(s,'ownerVals',{}),{}), numel(app.SchEmVals), app.SchFluorVals{1});
            app.SchSplitterVals = padCell(ensureRow(getff(s,'splitterVals',{}),{}), numel(app.SchEmVals)-1, '(none)');
            if isfield(s,'Rback') && ~isempty(app.RbackCheck) && isvalid(app.RbackCheck)
                app.RbackCheck.Value = s.Rback > 0;
            end
            if isfield(s,'blockOD'); setBlockODValue(app, s.blockOD); end
            if isfield(s,'opticsMode')
                mode = normalizeOpticsMode(app, s.opticsMode);
                if any(strcmp(app.ExOpticsMode.Items, mode))
                    app.ExOpticsMode.Value = mode;
                end
            end
            if isfield(s,'primary') && ~isempty(app.SchPrimaryDD) && isvalid(app.SchPrimaryDD)
                setDropValue(app.SchPrimaryDD, s.primary); end
            rebuildExcitationRows(app); rebuildFluorRows(app); rebuildChannelRows(app);
            onSchematicUpdate(app);
            function v = getff(st,f,d); if isfield(st,f); v=st.(f); else; v=d; end; end
            function c = ensureRow(c,def); if ~iscell(c)||isempty(c); c=def; else; c=c(:)'; end; end
        end

        function onSaveSchematic(app)
            d = filterSetDir(app);
            answer = inputdlg({'Filter set name:'}, 'Save filter set', [1 48], ...
                {['filterset_' datestr(now,'yyyymmdd_HHMMSS')]});
            if isempty(answer); return; end
            name = cleanFileName(strtrim(answer{1}));
            if isempty(name); setStatus(app,'Name cannot be empty.',true); return; end
            sch = schematicState(app, name);
            setDir = fullfile(d, name);
            if isfolder(setDir)
                if ~strcmp(questdlg(sprintf('Replace "%s"?',name),'Confirm','Replace','Cancel','Cancel'),'Replace'); return; end
            end
            if ~isfolder(setDir); mkdir(setDir); end
            clearLegacyFilterSetSpectra(setDir);
            writeSchematicManifest(fullfile(setDir,'filter_set.txt'), sch);
            nRefs = writeSchematicPointerManifest(app, fullfile(setDir,'spectra_pointers.json'), sch);
            setStatus(app, sprintf('Saved filter set "%s" with %d local spectrum pointer(s) to %s.', name, nRefs, setDir));
        end

        function onLoadSchematic(app)
            d = filterSetDir(app);
            p = uigetdir(d, 'Load filter set folder');
            if isequal(p,0); return; end
            try
                [sch,n] = loadSchematicFilterSetFolder(app, p);
                setStatus(app, sprintf('Loaded filter set "%s" and resolved %d spectrum pointer(s).', getNameSafe(sch,p), n));
            catch ME
                setStatus(app, ['Load filter set failed: ' ME.message], true);
                return;
            end
            % bring the user to the schematic tab
            for tt = app.Tabs.Children'
                if strcmp(tt.Title,'Design'); app.Tabs.SelectedTab = tt; break; end
            end
            function nm = getNameSafe(st,fb); if isfield(st,'name')&&~isempty(st.name); nm=st.name; else; [~,nm]=fileparts(fb); end; end
        end

        function [sch,n] = loadSchematicFilterSetFolder(app, p)
            oldLib = app.Lib;
            oldState = schematicState(app, 'before_load');
            try
                mf = fullfile(p,'filter_set.txt');
                if ~isfile(mf)
                    [~,folderName] = fileparts(p);
                    error('"%s" has spectra but no filter_set.txt manifest, so it is not a displayable saved filter set.', folderName);
                end
                sch = readSchematicManifest(mf);
                ptr = fullfile(p, 'spectra_pointers.json');
                if isfile(ptr)
                    n = importFilterSetPointers(app, ptr, true);
                else
                    n = importFilterSetSpectra(app, p, true);
                end
                assertSchematicRefsAvailable(app, sch);
                applySchematicState(app, sch);
            catch ME
                app.Lib = oldLib;
                try; refreshLibTree(app); catch; end
                try; applySchematicState(app, oldState); catch; end
                rethrow(ME);
            end
        end

        function assertSchematicRefsAvailable(app, sch)
            missing = {};
            checkList(getCellField(sch,'fluorVals',{}), 'fluorophore');
            checkList(getCellField(sch,'srcVals',{}), 'source');
            checkList(getCellField(sch,'cleanVals',{}), 'filter');
            checkList(getCellField(sch,'combinerVals',{}), 'dichroic');
            checkList({getTextField(sch,'primary','(none)')}, 'dichroic');
            checkList(getCellField(sch,'emVals',{}), 'filter');
            checkList(getCellField(sch,'detVals',{}), 'detector');
            checkList(getCellField(sch,'splitterVals',{}), 'dichroic');
            if ~isempty(missing)
                error('Missing or wrong-type spectra: %s', strjoin(unique(missing,'stable'), ', '));
            end
            function checkList(vals, role)
                vals = vals(:)';
                for ii = 1:numel(vals)
                    name = vals{ii};
                    if isPresetPlaceholder(name); continue; end
                    if strcmp(role,'source') && ~isempty(regexp(name,'^\d{3}\s+nm laser$','once')); continue; end
                    if strcmp(role,'detector') && any(strcmp(detectorPreset('list',app.Lambda), name)); continue; end
                    S = getSpec(app, name);
                    if isempty(S) || ~strcmp(spectrumRole(app, S), role)
                        missing{end+1} = sprintf('%s (%s)', name, role); %#ok<AGROW>
                    end
                end
            end
        end

        function n = writeSchematicPointerManifest(app, fp, sch)
            refs = schematicSpectrumRefs(app, sch);
            emptyEntry = struct('field','','role','','name','','category','','kind','', ...
                                'path','','file','','preset',false);
            entries = repmat(emptyEntry, 0, 1);
            for ii = 1:numel(refs)
                r = refs(ii);
                e = emptyEntry;
                e.field = r.field; e.role = r.role; e.name = r.name;
                [S,~] = getSpecByRole(app, r.name, r.role);
                if isempty(S)
                    if isGeneratedPreset(app, r.name, r.role)
                        e.preset = true;
                        entries(end+1) = e; %#ok<AGROW>
                        continue;
                    end
                    error('Cannot save pointer for missing spectrum "%s" (%s).', r.name, r.role);
                end
                e.category = S.category;
                e.kind = S.kind;
                e.file = char(S.file);
                localFile = canonicalSpectrumFile(app, S, r.role);
                if isempty(localFile)
                    error('Cannot find a local file for "%s" (%s).', r.name, r.role);
                end
                e.path = relativePath(app, localFile);
                entries(end+1) = e; %#ok<AGROW>
            end
            P = struct('type','spectra_pointers','version',1,'savedOn',datestr(now), ...
                       'filterSet',getTextField(sch,'name','filter_set'), 'entries',entries);
            fid = fopen(fp, 'w');
            if fid < 0; error('Could not write %s', fp); end
            cleaner = onCleanup(@() fclose(fid));
            fprintf(fid, '%s\n', jsonencode(P));
            clear cleaner;
            n = numel(entries);
        end

        function refs = schematicSpectrumRefs(app, sch)
            defs = { ...
                'srcVals',       'source'; ...
                'cleanVals',     'filter'; ...
                'combinerVals',  'dichroic'; ...
                'primary',       'dichroic'; ...
                'fluorVals',     'fluorophore'; ...
                'emVals',        'filter'; ...
                'detVals',       'detector'; ...
                'splitterVals',  'dichroic'};
            refs = struct('field',{},'role',{},'name',{});
            seen = {};
            for dd = 1:size(defs,1)
                field = defs{dd,1}; role = defs{dd,2};
                if strcmp(field,'primary')
                    vals = {getTextField(sch, field, '(none)')};
                else
                    vals = getCellField(sch, field, {});
                end
                vals = vals(:)';
                for ii = 1:numel(vals)
                    name = vals{ii};
                    if isPresetPlaceholder(name) || isGeneratedPreset(app, name, role); continue; end
                    key = [role '|' name];
                    if any(strcmp(seen, key)); continue; end
                    seen{end+1} = key; %#ok<AGROW>
                    refs(end+1) = struct('field',field,'role',role,'name',name); %#ok<AGROW>
                end
            end
        end

        function [S, idx] = getSpecByRole(app, name, role)
            idx = [];
            S = [];
            hits = find(strcmp({app.Lib.name}, name));
            for hh = 1:numel(hits)
                cand = app.Lib(hits(hh));
                if isempty(role) || strcmp(spectrumRole(app, cand), role)
                    S = cand;
                    idx = hits(hh);
                    return;
                end
            end
        end

        function fp = canonicalSpectrumFile(app, S, role)
            fp = '';
            if isfield(S,'file') && ischar(S.file) && isfile(S.file) ...
                    && ~isInsideDir(S.file, filterSetDir(app))
                fp = S.file;
                return;
            end
            fp = findLibrarySpectrumFile(app, S.name, role);
            if isempty(fp) && isfield(S,'file') && ischar(S.file) && isfile(S.file)
                fp = S.file;
            end
        end

        function fp = findLibrarySpectrumFile(app, name, role)
            fp = '';
            cats = {'Proteins','Filters','Dichroics','Illumations','Illuminations','Detectors'};
            root = spectraDir(app);
            for cc = 1:numel(cats)
                d = fullfile(root, cats{cc});
                if ~isfolder(d); continue; end
                files = [dir(fullfile(d,'*.txt')); dir(fullfile(d,'*.csv')); dir(fullfile(d,'*.TXT'))];
                [~, uidx] = unique(lower({files.name}), 'stable'); files = files(uidx);
                for ff = 1:numel(files)
                    candidate = fullfile(files(ff).folder, files(ff).name);
                    if strcmpi(cats{cc},'Illumations') && strcmpi(files(ff).name,'Prizmatix.csv')
                        continue;
                    end
                    try
                        S2 = loadSpectrum(candidate, app.Lambda);
                        probe = struct('name',S2.name,'category',cats{cc},'kind',S2.kind,'file',candidate);
                        if strcmp(S2.name, name) && strcmp(spectrumRole(app, probe), role)
                            fp = candidate;
                            return;
                        end
                    catch
                    end
                end
            end
        end

        function rel = relativePath(app, fp)
            root = canonicalPath(app.LibRoot);
            fpc = canonicalPath(fp);
            if startsWithPath(fpc, root)
                rel = fpc(numel(root)+2:end);
            else
                rel = fpc;
            end
            rel = strrep(rel, filesep, '/');
        end

        function fp = pointerPath(app, relPath)
            relPath = char(relPath);
            if isempty(relPath); fp = ''; return; end
            relPath = strrep(relPath, '/', filesep);
            if isAbsolutePath(relPath)
                fp = relPath;
            else
                fp = fullfile(app.LibRoot, relPath);
            end
        end

        function tf = isGeneratedPreset(app, name, role)
            tf = false;
            if isPresetPlaceholder(name); tf = true; return; end
            if strcmp(role,'source') && ~isempty(regexp(name,'^\d{3}\s+nm laser$','once'))
                tf = true;
            elseif strcmp(role,'detector') && any(strcmp(detectorPreset('list',app.Lambda), name)) ...
                    && isempty(getSpecByRole(app, name, role))
                tf = true;
            end
        end

        function n = importFilterSetPointers(app, pointerFile, strict)
            if nargin < 3; strict = false; end
            P = jsondecode(fileread(pointerFile));
            if ~isfield(P,'entries') || isempty(P.entries); n = 0; return; end
            entries = P.entries;
            n = 0;
            for ii = 1:numel(entries)
                e = entries(ii);
                name = jsonText(e, 'name', '');
                role = jsonText(e, 'role', '');
                if isPresetPlaceholder(name) || isGeneratedPreset(app, name, role)
                    continue;
                end
                relPath = jsonText(e, 'path', '');
                fp = pointerPath(app, relPath);
                if isempty(fp)
                    if strict; error('No local spectrum pointer recorded for "%s".', name); end
                    continue;
                end
                if ~isfile(fp)
                    if strict; error('Pointer target missing for "%s": %s', name, relPath); end
                    continue;
                end
                try
                    S = loadSpectrum(fp, app.Lambda);
                    if ~strcmp(S.name, name)
                        error('expected "%s", found "%s"', name, S.name);
                    end
                    category = jsonText(e, 'category', '');
                    if isempty(category); category = inferFilterSetCategory(app, S, fp); end
                    probe = struct('name',S.name,'category',category,'kind',S.kind,'file',fp);
                    if ~isempty(role) && ~strcmp(spectrumRole(app, probe), role)
                        error('expected role "%s", found "%s"', role, spectrumRole(app, probe));
                    end
                    fl = NaN; if isfield(S,'floor'); fl = S.floor; end
                    addToLib(app, S.name, category, S.kind, S.ex, S.em, [], fp, fl);
                    n = n + 1;
                catch ME
                    if strict
                        error('Could not resolve pointer "%s": %s', name, ME.message);
                    end
                    warning('Skip pointer %s: %s', name, ME.message);
                end
            end
            refreshLibTree(app);
            populateChoosers(app);
        end

        function s = filterSetToSchematic(app, fs)
            % Map a legacy System-Builder preset onto schematic state.
            s = struct('version',2,'name',getf(fs,'name','legacy'));
            ld = getf(fs,'laserTable',{}); ns = size(ld,1);
            s.srcVals = {}; s.cleanVals = {};
            for i = 1:ns
                wl = ld{i,2}; if ischar(wl); wl = str2double(wl); end
                s.srcVals{i} = sprintf('%g nm laser', wl);
                if size(ld,2) >= 5; s.cleanVals{i} = ld{i,5}; else; s.cleanVals{i} = '(none)'; end
            end
            if isempty(s.srcVals); s.srcVals = {'488 nm laser'}; s.cleanVals = {'(none)'}; end
            s.combinerVals = repmat({'(none)'}, 1, max(0,numel(s.srcVals)-1));
            fd = getf(fs,'fluorTable',{}); s.fluorVals = {};
            for i = 1:size(fd,1); s.fluorVals{i} = fd{i,1}; end
            if isempty(s.fluorVals); s.fluorVals = {'mNeonGreen'}; end
            cd = getf(fs,'chanTable',{}); s.emVals={}; s.ownerVals={}; s.splitterVals={};
            for k = 1:size(cd,1)
                s.emVals{k} = cd{k,4}; s.ownerVals{k} = cd{k,5};
                if k >= 2; s.splitterVals{k-1} = cd{k,2}; end
            end
            if isempty(s.emVals); s.emVals = {'(none)'}; s.ownerVals = s.fluorVals(1); end
            s.detVals = repmat({getf(fs,'detector','(ideal)')}, 1, numel(s.emVals));
            s.primary = getf(fs,'primary','(none)');
            s.Rback = getf(fs,'RbackPercent',0.5); s.blockOD = getf(fs,'blockOD',6);
            function v = getf(st,f,dv); if isfield(st,f)&&~isempty(st.(f)); v=st.(f); else; v=dv; end; end
        end

        function filterSet = collectFilterSetState(app, name)
            filterSet = struct();
            filterSet.version = 1;
            filterSet.name = name;
            filterSet.savedOn = datestr(now);
            filterSet.fluorTable = app.FluorTable.Data;
            filterSet.laserTable = app.LaserTable.Data;
            filterSet.chanTable = app.ChanTable.Data;
            filterSet.primary = app.PrimaryDD.Value;
            filterSet.detector = app.DetectorDD.Value;
            filterSet.RbackPercent = rbackPercent(app);
            filterSet.blockOD = blockODValue(app);
        end

        function applyFilterSetState(app, filterSet)
            populateChoosers(app);
            app.FluorTable.Data = filterSet.fluorTable;
            app.LaserTable.Data = filterSet.laserTable;
            app.ChanTable.Data = filterSet.chanTable;
            setDropValue(app.PrimaryDD, filterSet.primary);
            setDropValue(app.DetectorDD, filterSet.detector);
            if ~isempty(app.RbackCheck) && isvalid(app.RbackCheck)
                app.RbackCheck.Value = filterSet.RbackPercent > 0;
            end
            setBlockODValue(app, filterSet.blockOD);
            app.LastResults = struct();
            app.LastPlotConfig = struct();
        end

        function n = importFilterSetSpectra(app, setDir, strict)
            if nargin < 3; strict = false; end
            files = [dir(fullfile(setDir,'**','*.txt')); dir(fullfile(setDir,'**','*.TXT')); ...
                     dir(fullfile(setDir,'**','*.csv')); dir(fullfile(setDir,'**','*.CSV'))];
            if isempty(files); n = 0; return; end
            fullNames = lower(string(fullfile({files.folder}, {files.name})));
            [~, uidx] = unique(fullNames, 'stable');
            files = files(uidx);
            n = 0;
            for ii = 1:numel(files)
                fp = fullfile(files(ii).folder, files(ii).name);
                if strcmpi(files(ii).name, 'filter_set.txt'); continue; end
                try
                    S = loadSpectrum(fp, app.Lambda);
                    category = inferFilterSetCategory(app, S, fp);
                    probe = struct('name',S.name,'category',category,'kind',S.kind,'file',fp);
                    role = spectrumRole(app, probe);
                    if ~isempty(getSpecByRole(app, S.name, role))
                        n = n + 1;
                        continue;
                    end
                    fl = NaN; if isfield(S,'floor'); fl = S.floor; end
                    addToLib(app, S.name, category, S.kind, S.ex, S.em, [], fp, fl);
                    n = n + 1;
                catch ME
                    if strict
                        error('Could not import "%s": %s', files(ii).name, ME.message);
                    end
                    warning('Skip %s: %s', fp, ME.message);
                end
            end
            refreshLibTree(app);
            populateChoosers(app);
        end

        function exportFilterSetSpectra(app, setDir, sch)
            names = {};
            fields = {'srcVals','cleanVals','combinerVals','primary','fluorVals', ...
                      'emVals','detVals','splitterVals'};
            for ii = 1:numel(fields)
                if ~isfield(sch, fields{ii}); continue; end
                v = sch.(fields{ii});
                if ischar(v) || isstring(v); v = cellstr(v); end
                if iscell(v); names = [names, v(:)']; end %#ok<AGROW>
            end
            names = unique(names(~cellfun(@isPresetPlaceholder, names)), 'stable');
            for ii = 1:numel(names)
                S = getSpec(app, names{ii});
                if isempty(S); continue; end
                fn = [cleanFileName(S.name) '.txt'];
                fp = fullfile(setDir, fn);
                if strcmp(S.kind,'fluorophore') && ~isempty(S.em)
                    writeSpectrumTable(fp, app.Lambda, S.ex, S.em);
                else
                    writeSpectrumTable(fp, app.Lambda, S.ex);
                end
            end
        end

        function category = inferFilterSetCategory(app, S, fp)
            hit = find(strcmp({app.Lib.name}, S.name), 1);
            if ~isempty(hit)
                category = app.Lib(hit).category;
                return;
            end
            T = struct('name',S.name,'category','','kind','filter','file',fp);
            switch spectrumRole(app, T)
                case 'detector'; category = 'Detectors';
                case 'source';   category = 'Illumations';
                case 'dichroic'; category = 'Dichroics';
                otherwise
                    if strcmp(S.kind,'fluorophore')
                        category = 'Proteins';
                    else
                        category = 'Filters';
                    end
            end
        end
    end

    %% ---------------- Compute / results ----------------
    methods (Access = private)
        function [fluors, lasers, channels, assign, detector, meta] = assembleSystem(app)
            % Fluorophores
            fd = app.FluorTable.Data; fluors = struct([]);
            for i = 1:size(fd,1)
                S = getSpec(app, fd{i,1});
                if isempty(S); error('Unknown fluorophore "%s"', fd{i,1}); end
                fluors(i).name = S.name; fluors(i).ex = S.ex; fluors(i).em = S.em;
                b = fd{i,2}; if ischar(b); b = str2double(b); end
                fluors(i).brightness = b;
            end
            % Excitation sources
            ld = app.LaserTable.Data; lasers = struct([]);
            for j = 1:size(ld,1)
                lasers(j).name = ld{j,1};
                lasers(j).wl = num(ld{j,2}); lasers(j).power = num(ld{j,3});
                lasers(j).spectrum = []; lasers(j).exFilter = [];
                if size(ld,2) >= 4 && ~isempty(ld{j,4}) && ~strcmp(ld{j,4},'(laser line)')
                    ss = getSpec(app, ld{j,4});
                    if ~isempty(ss); lasers(j).spectrum = ss.ex; end
                end
                if size(ld,2) >= 5 && ~isempty(ld{j,5}) && ~strcmp(ld{j,5},'(none)')
                    ef = getSpec(app, ld{j,5});
                    if ~isempty(ef); lasers(j).exFilter = ef.ex; end
                end
            end
            % Channels
            cd = app.ChanTable.Data; channels = struct([]); assign = [];
            primary = [];
            if ~strcmp(app.PrimaryDD.Value,'(none)')
                primary = getSpec(app, app.PrimaryDD.Value);
            end
            fnames = {fluors.name};
            modes = {}; splitNameShared = ''; splitFloorShared = NaN;
            primaryFloor = NaN; if ~isempty(primary) && isfield(primary,'floor'); primaryFloor = primary.floor; end
            for k = 1:size(cd,1)
                channels(k).name = cd{k,1};
                path = struct('T',{},'mode',{},'floor',{});
                if ~isempty(primary)
                    path(end+1) = struct('T',primary.ex,'mode','T','floor',primaryFloor); %#ok<AGROW>
                end
                splitName = cd{k,2}; splitMode = cd{k,3};
                modes{k} = splitMode; %#ok<AGROW>
                if ~strcmp(splitName,'(none)') && ~strcmpi(splitMode,'none')
                    sp = getSpec(app, splitName);
                    if ~isempty(sp)
                        spf = NaN; if isfield(sp,'floor'); spf = sp.floor; end
                        path(end+1) = struct('T',sp.ex,'mode',upper(splitMode),'floor',spf); %#ok<AGROW>
                        if isempty(splitNameShared); splitNameShared = splitName; splitFloorShared = spf; end
                    end
                end
                channels(k).path = path;
                ef = getSpec(app, cd{k,4});
                if isempty(ef)
                    channels(k).emFilter = ones(numel(app.Lambda),1); channels(k).emFloor = NaN;
                else
                    channels(k).emFilter = ef.ex;
                    channels(k).emFloor = NaN; if isfield(ef,'floor'); channels(k).emFloor = ef.floor; end
                end
                a = find(strcmp(fnames, cd{k,5}),1);
                if isempty(a); a = min(k, numel(fluors)); end
                assign(k) = a; %#ok<AGROW>
            end

            % detector QE
            detector = ones(numel(app.Lambda),1);
            if ~strcmp(app.DetectorDD.Value,'(ideal)')
                ds = getSpec(app, app.DetectorDD.Value);
                if ~isempty(ds); detector = ds.ex(:); end
            end

            % meta for the joint optimiser
            if ~isempty(primary); primaryT = primary.ex(:); else; primaryT = []; end
            meta = struct('names',{{channels.name}}, 'primaryT',primaryT, ...
                'primaryFloor',primaryFloor, 'splitterName',splitNameShared, ...
                'splitterFloor',splitFloorShared, 'modes',{modes}, ...
                'Rback',rbackPercent(app)/100, 'blockOD',blockODValue(app));
            function v = num(x); if ischar(x); v = str2double(x); else; v = x; end; end
        end

        function onCompute(app)
            try
                [fluors, lasers, channels, assign, detector, meta] = assembleSystem(app);
            catch ME
                setStatus(app, ME.message, true); return;
            end
            displayResults(app, fluors, lasers, channels, assign, detector, meta.Rback, meta.blockOD);
        end

        function ensureResultsBuilt(app)
            % Build the hidden Results panel on first use (lazy — see buildUI).
            if app.ResultsBuilt; return; end
            p = buildResultsTab(app);
            bumpFonts(app, p, 2);   % match the global +2 pt applied in buildUI
            app.ResultsBuilt = true;
        end

        function displayResults(app, fluors, lasers, channels, assign, detector, Rback, blockOD)
            ensureResultsBuilt(app);
            [S, eff] = FilterEngine.signalMatrix(fluors, lasers, channels, app.Lambda, detector);
            CT = FilterEngine.crosstalkMatrix(S, assign);
            bleed = FilterEngine.laserBleed(channels, lasers, app.Lambda, Rback, detector, blockOD);
            fn = {fluors.name}; cn = {channels.name}; ln = {lasers.name};

            app.SigTable.RowName = fn; app.SigTable.ColumnName = cn;
            app.SigTable.Data = round(S,3);
            app.CTTable.RowName = fn;  app.CTTable.ColumnName = cn;
            app.CTTable.Data = round(100*CT,2);   % percent

            % laser back-reflection background / signal (%), per channel x source
            bsr = zeros(numel(channels), numel(lasers));
            for k = 1:numel(channels)
                bsr(k,:) = 100 * bleed(k,:) / (S(assign(k),k) + eps);
            end
            app.BleedTable.RowName = cn; app.BleedTable.ColumnName = ln;
            app.BleedTable.Data = bsr;

            app.LastPlotConfig = struct('fluors',fluors,'lasers',lasers, ...
                'channels',channels,'assign',assign,'detector',detector);
            redrawResultsPlot(app);

            if ~isempty(app.OptChanDD) && isvalid(app.OptChanDD)
                app.OptChanDD.Items = cn; onOptChanChange(app);
            end
            sc = FilterEngine.systemScore(S, assign, optWeight(app), bleed, optBleedWeight(app));
            app.LastResults = struct('S',S,'CT',CT,'eff',eff,'bleed',bleed,'bsr',bsr, ...
                'fluors',{fn},'channels',{cn},'lasers',{ln},'score',sc);
            setStatus(app, sprintf('Computed. System score = %.3f (incl. crosstalk + laser bleed).', sc));
            for tt = app.Tabs.Children'
                if strcmp(tt.Title,'Results'); app.Tabs.SelectedTab = tt; break; end
            end
        end

        function onXRangeEdit(app)
            xmin = app.ResXMinField.Value; xmax = app.ResXMaxField.Value;
            if xmax <= xmin
                xmax = min(app.Lambda(end), xmin + 10);
                app.ResXMaxField.Value = xmax;
            end
            app.ResXSlider.Value = mean([xmin xmax]);
            redrawResultsPlot(app);
        end

        function onXSlider(app, ctr)
            if nargin < 2; ctr = app.ResXSlider.Value; end
            xmin = app.ResXMinField.Value; xmax = app.ResXMaxField.Value;
            span = max(10, xmax - xmin);
            lo = ctr - span/2; hi = ctr + span/2;
            if lo < app.Lambda(1); lo = app.Lambda(1); hi = lo + span; end
            if hi > app.Lambda(end); hi = app.Lambda(end); lo = hi - span; end
            app.ResXMinField.Value = max(app.Lambda(1), lo);
            app.ResXMaxField.Value = min(app.Lambda(end), hi);
            redrawResultsPlot(app);
        end

        function redrawResultsPlot(app)
            if isempty(app.ResAxes) || ~isvalid(app.ResAxes); return; end
            cla(app.ResAxes);
            if ~isfield(app.LastPlotConfig,'channels') || isempty(app.LastPlotConfig)
                xlabel(app.ResAxes,'Wavelength (nm)');
                ylabel(app.ResAxes,'Throughput / emission');
                title(app.ResAxes,'Compute results to plot spectra');
                return;
            end
            cfg = app.LastPlotConfig;
            lam = app.Lambda; channels = cfg.channels; fluors = cfg.fluors; lasers = cfg.lasers;
            co = lines(max(numel(channels),numel(fluors)));
            hold(app.ResAxes,'on'); labels = {}; plottedY = [];
            switch app.ResPlotMode.Value
                case 'Channel overlay'
                    for k = 1:numel(channels)
                        t = FilterEngine.pathTransmission(channels(k), lam) .* cfg.detector(:);
                        y = yPlot(app,t);
                        plot(app.ResAxes, lam, y, 'Color',co(k,:),'LineWidth',2);
                        plottedY = [plottedY; y]; %#ok<AGROW>
                        labels{end+1} = channels(k).name; %#ok<AGROW>
                    end
                    for i = 1:numel(fluors)
                        em = normPeak(fluors(i).em);
                        y = yPlot(app,em);
                        plot(app.ResAxes, lam, y, ':','Color',[.25 .25 .25],'LineWidth',1);
                        plottedY = [plottedY; y]; %#ok<AGROW>
                        labels{end+1} = [fluors(i).name ' em']; %#ok<AGROW>
                    end

                case 'Raw vs filtered'
                    for k = 1:numel(channels)
                        owner = cfg.assign(k);
                        raw = normPeak(fluors(owner).em);
                        t = FilterEngine.pathTransmission(channels(k), lam) .* cfg.detector(:);
                        filt = normPeak(fluors(owner).em(:) .* t);
                        yRaw = yPlot(app,raw);
                        yFilt = yPlot(app,filt);
                        plot(app.ResAxes, lam, yRaw, '--','Color',co(k,:),'LineWidth',1.2);
                        plot(app.ResAxes, lam, yFilt, '-','Color',co(k,:),'LineWidth',2);
                        plottedY = [plottedY; yRaw; yFilt]; %#ok<AGROW>
                        labels{end+1} = [channels(k).name ' raw ' fluors(owner).name]; %#ok<AGROW>
                        labels{end+1} = [channels(k).name ' filtered']; %#ok<AGROW>
                    end

                otherwise % Component superposition
                    for i = 1:numel(fluors)
                        yEx = yPlot(app,normPeak(fluors(i).ex));
                        yEm = yPlot(app,normPeak(fluors(i).em));
                        plot(app.ResAxes, lam, yEx, '--', ...
                            'Color',co(i,:),'LineWidth',1);
                        plot(app.ResAxes, lam, yEm, ':', ...
                            'Color',co(i,:),'LineWidth',1.2);
                        plottedY = [plottedY; yEx; yEm]; %#ok<AGROW>
                        labels{end+1} = [fluors(i).name ' ex']; %#ok<AGROW>
                        labels{end+1} = [fluors(i).name ' em']; %#ok<AGROW>
                    end
                    for k = 1:numel(channels)
                        t = FilterEngine.pathTransmission(channels(k), lam) .* cfg.detector(:);
                        y = yPlot(app,t);
                        plot(app.ResAxes, lam, y, 'Color',co(k,:),'LineWidth',2.2);
                        plottedY = [plottedY; y]; %#ok<AGROW>
                        labels{end+1} = [channels(k).name ' total']; %#ok<AGROW>
                        for p = 1:numel(channels(k).path)
                            el = channels(k).path(p).T(:);
                            if strcmpi(channels(k).path(p).mode,'R'); el = 1 - el; end
                            pale = 0.65*co(k,:) + 0.35*[1 1 1];
                            y = yPlot(app,el);
                            plot(app.ResAxes, lam, y, '-', ...
                                'Color',pale,'LineWidth',0.8);
                            plottedY = [plottedY; y]; %#ok<AGROW>
                            labels{end+1} = sprintf('%s path %d%s',channels(k).name,p,channels(k).path(p).mode); %#ok<AGROW>
                        end
                    end
            end
            for j = 1:numel(lasers)
                xline(app.ResAxes, lasers(j).wl, '-.', 'Color',[.4 .4 .4]);
            end
            hold(app.ResAxes,'off');
            xmin = app.ResXMinField.Value; xmax = app.ResXMaxField.Value;
            xlim(app.ResAxes,[xmin xmax]);
            if strcmp(app.ResYMode.Value,'OD')
                ylabel(app.ResAxes,'Optical density / -log10(normalised)');
                odMax = max(plottedY(isfinite(plottedY)));
                if isempty(odMax) || odMax <= 0; odMax = 2; end
                ylim(app.ResAxes,[0 min(8, max(2, ceil(odMax)))]);
                app.ResAxes.YDir = 'reverse';
            else
                ylabel(app.ResAxes,'%T / normalised emission (%)');
                ylim(app.ResAxes,[0 105]);
                app.ResAxes.YDir = 'normal';
            end
            xlabel(app.ResAxes,'Wavelength (nm)');
            title(app.ResAxes, app.ResPlotMode.Value);
            if ~isempty(labels)
                legend(app.ResAxes, labels, 'Location','northeastoutside','Interpreter','none');
            end
            grid(app.ResAxes,'on');
        end

        function onExportResults(app)
            if ~isfield(app.LastResults,'S') || isempty(app.LastResults.S)
                setStatus(app,'Compute results first.',true); return;
            end
            [f,p] = uiputfile('*.csv','Export results', 'filter_set_results.csv');
            if isequal(f,0); return; end
            R = app.LastResults; fid = fopen(fullfile(p,f),'w');
            c = onCleanup(@() fclose(fid));
            fprintf(fid,'Optimal Filter Set - results export\n');
            fprintf(fid,'System score,%.4f\n\n', R.score);
            writeMat(fid,'Signal matrix S(fluor,channel)', R.S, R.fluors, R.channels);
            writeMat(fid,'Crosstalk % (col-normalised to owner)', 100*R.CT, R.fluors, R.channels);
            writeMat(fid,'Collection efficiency %', 100*R.eff, R.fluors, R.channels);
            writeMat(fid,'Laser back-reflection background/signal %', R.bsr, R.channels, R.lasers);
            setStatus(app, sprintf('Exported results to %s', fullfile(p,f)));
            function writeMat(fid, title, M, rn, cn)
                fprintf(fid,'%s\n', title);
                fprintf(fid,',%s', strjoin(cn,',')); fprintf(fid,'\n');
                for r = 1:size(M,1)
                    fprintf(fid,'%s', rn{r});
                    fprintf(fid,',%.4g', M(r,:)); fprintf(fid,'\n');
                end
                fprintf(fid,'\n');
            end
        end

        function onSaveConfig(app)
            try; cfg = config(app); catch ME; setStatus(app,ME.message,true); return; end
            [f,p] = uiputfile('*.mat','Save filter-set configuration','filter_config.mat');
            if isequal(f,0); return; end
            save(fullfile(p,f),'cfg');
            setStatus(app, sprintf('Saved config (%d fluorophores, %d channels) to %s. Open it in FilterSetSNRApp.', ...
                numel(cfg.fluors), numel(cfg.channels), fullfile(p,f)));
        end

        function onExportOpt(app)
            if isempty(app.LastOptResults) || isempty(app.OptResTable.Data)
                setStatus(app,'Run the optimizer first.',true); return;
            end
            [f,p] = uiputfile('*.csv','Export ranking','filter_ranking.csv');
            if isequal(f,0); return; end
            cols = app.OptResTable.ColumnName; data = app.OptResTable.Data;
            fid = fopen(fullfile(p,f),'w'); c = onCleanup(@() fclose(fid));
            fprintf(fid,'%s\n', strjoin(cols(:)',','));
            for r = 1:size(data,1)
                cell2 = cell(1,size(data,2));
                for k = 1:size(data,2)
                    v = data{r,k};
                    if isnumeric(v); cell2{k} = num2str(v); else; cell2{k} = char(v); end
                end
                fprintf(fid,'%s\n', strjoin(cell2,','));
            end
            setStatus(app, sprintf('Exported ranking (%d rows) to %s', size(data,1), fullfile(p,f)));
        end

        function onOptChanChange(app)
            % Load the saved candidate selection for the now-selected channel.
            ch = app.OptChanDD.Value;
            if isKey(app.OptCandMap, ch)
                app.OptCandList.Value = intersect(app.OptCandMap(ch), ...
                    app.OptCandList.Items, 'stable');
            else
                app.OptCandList.Value = {};
            end
        end

        function onCandSelChange(app)
            % Persist the current listbox selection for this channel.
            v = app.OptCandList.Value; if ischar(v); v = {v}; end
            app.OptCandMap(app.OptChanDD.Value) = v;
        end

        function arr = candArray(app, names)
            arr = struct('name',{},'ex',{},'floor',{});
            for q = 1:numel(names)
                S = getSpec(app, names{q});
                if isempty(S); continue; end
                fl = NaN; if isfield(S,'floor'); fl = S.floor; end
                arr(end+1) = struct('name',S.name,'ex',S.ex(:),'floor',fl); %#ok<AGROW>
            end
        end

        function scoreFcn = buildSNRScorer(app, fluors, assign, detector, meta, snrMode, targetCh)
            % Returns [] for optical scoring, or a closure scoreFcn(ch,las) that
            % ranks a candidate set by electron-domain SNR (uses snrModel, so it
            % shares the SNR app's physics exactly). targetCh = channel index to
            % score (single mode); [] => minimum over all channels (joint mode).
            if nargin < 7; targetCh = []; end
            scoreFcn = [];
            if ~snrMode; return; end
            % photophysics per fluorophore (ec/qy from FPbase when available)
            conc = app.OptSNRFields.conc.Value * 1e-6;
            fp = struct('ec',{},'qy',{},'conc_M',{});
            for i = 1:numel(fluors)
                S = getSpec(app, fluors(i).name);
                ec = 50000; qy = 0.6;
                if ~isempty(S)
                    if isfield(S,'ec') && isfinite(S.ec); ec = S.ec; end
                    if isfield(S,'qy') && isfinite(S.qy); qy = S.qy; end
                end
                fp(i) = struct('ec',ec,'qy',qy,'conc_M',conc);
            end
            phys = struct('NA',app.OptSNRFields.NA.Value,'n',1.33, ...
                'pathLength_cm',0.01,'tInt_s',app.OptSNRFields.integ.Value*1e-3, ...
                'readNoise_e',app.OptSNRFields.read.Value,'darkRate_eps',100, ...
                'powers_mW',app.OptSNRFields.power.Value*ones(1,numel(meta.names)));
            af = struct('name',{},'absorb',{},'em',{},'strength',{},'qy',{});
            if app.OptAFcheck.Value
                af = [autofluorPreset('Brain tissue',app.Lambda), ...
                      autofluorPreset('Silica fiber',app.Lambda)];
            end
            lam = app.Lambda; flBase = fluors; asg = assign; det = detector;
            Rb = meta.Rback; bOD = meta.blockOD;
            scoreFcn = @(ch, las) localSNR(ch, las);
            function v = localSNR(ch, las)
                cfg = struct('lambda',lam,'fluors',flBase,'lasers',las, ...
                    'channels',ch,'assign',asg,'detector',det,'Rback',Rb,'blockOD',bOD);
                out = snrModel(cfg, phys, fp, af);
                if isempty(targetCh); v = min(out.SNR);   % joint: worst channel
                else; v = out.SNR(targetCh); end          % single: that channel
            end
        end

        function onOptimize(app)
            try
                [fluors, lasers, channels, assign, detector, meta] = assembleSystem(app);
            catch ME
                setStatus(app, ME.message, true); return;
            end
            cn = {channels.name}; nc = numel(channels);
            ck = find(strcmp(cn, app.OptChanDD.Value),1);
            if isempty(ck); setStatus(app,'Compute first, then pick a channel.',true); return; end
            joint = startsWith(app.OptMode.Value,'Joint');
            snrMode = startsWith(app.OptScoreMode.Value,'Electron');
            if joint; tgt = []; else; tgt = ck; end
            scoreFcn = buildSNRScorer(app, fluors, assign, detector, meta, snrMode, tgt);

            % per-channel candidate filter sets (map; empty => keep current filter)
            filterCands = cell(1,nc);
            for k = 1:nc
                names = {};
                if isKey(app.OptCandMap, cn{k}); names = app.OptCandMap(cn{k}); end
                if k == ck   % include live listbox selection for active channel
                    v = app.OptCandList.Value; if ischar(v); v = {v}; end
                    names = unique([names, v], 'stable');
                end
                if isempty(names)
                    filterCands{k} = struct('name',channels(k).name,'ex',channels(k).emFilter);
                else
                    filterCands{k} = candArray(app, names);
                end
            end
            if isempty(filterCands{ck}) || (~joint && isscalar(filterCands{ck}) && ...
                    strcmp(filterCands{ck}.name, channels(ck).name))
                setStatus(app,'Select candidate filters for the channel.',true); return;
            end

            if ~joint
                % ---- single channel: fix all others ----
                cand = cell(1,nc);
                for k = 1:nc
                    if k==ck; cand{k} = filterCands{k};
                    else; cand{k} = struct('name',channels(k).name,'ex',channels(k).emFilter); end
                end
                R = FilterEngine.optimize(fluors, lasers, channels, cand, ...
                    assign, app.Lambda, app.OptWeight.Value, detector, ...
                    meta.Rback, app.OptBleedWeight.Value, meta.blockOD, scoreFcn);
                owner = assign(ck); others = setdiff(1:numel(fluors), owner);
                scoreHdr = ternStr(snrMode, sprintf('SNR (%s)',cn{ck}), 'Score');
                rows = cell(numel(R),5);
                for r = 1:numel(R)
                    bsr = 100*sum(R(r).bleed(ck,:)) / (R(r).S(owner,ck)+eps);
                    rows(r,:) = {R(r).filters{ck}, round(R(r).score,3,'significant'), ...
                        round(R(r).eff(owner,ck)*100,1), round(100*sum(R(r).CT(others,ck)),2), ...
                        round(bsr,3,'significant')};
                end
                app.OptResTable.ColumnName = {'Emission filter',scoreHdr, ...
                    'Own collection %','Crosstalk into ch %','Laser bleed %'};
                app.OptResTable.Data = rows;
                app.LastOptResults = R; app.LastOptKind = 'single';
                setStatus(app, sprintf('Single-channel "%s": %d candidates ranked.', cn{ck}, numel(R)));
                return;
            end

            % ---- joint: co-optimise shared splitter + every channel filter ----
            splNames = app.OptDichroicList.Value; if ischar(splNames); splNames = {splNames}; end
            if ~isempty(meta.splitterName); splNames = unique([{meta.splitterName}, splNames],'stable'); end
            splNames = splNames(~cellfun(@isempty,splNames));
            if isempty(splNames)
                setStatus(app,'Joint mode needs a splitter dichroic (set one in System Builder).',true); return;
            end
            splitterCands = candArray(app, splNames);

            % optional excitation/cleanup filter candidates per source (auto-pull
            % laser-line filters near each source centre)
            exFilterCands = {}; nl = numel(lasers); exShown = false;
            if app.OptExCheck.Value
                exFilterCands = cell(1,nl);
                for j = 1:nl
                    exFilterCands{j} = pullExCandidates(app, lasers(j).wl);
                end
                exShown = true;
            end

            nCombos = numel(splitterCands);
            for k=1:nc; nCombos = nCombos*numel(filterCands{k}); end
            for j=1:numel(exFilterCands); nCombos = nCombos*max(1,numel(exFilterCands{j})); end
            if nCombos > 20000
                setStatus(app, sprintf('%d combinations is too many — trim candidates.',nCombos),true); return;
            end
            setStatus(app, sprintf('Searching %d combinations ...', nCombos)); drawnow;
            R = FilterEngine.optimizeJoint(fluors, lasers, meta.names, meta.primaryT, ...
                splitterCands, meta.modes, filterCands, assign, app.Lambda, ...
                app.OptWeight.Value, detector, meta.Rback, app.OptBleedWeight.Value, ...
                exFilterCands, meta.blockOD, meta.primaryFloor, scoreFcn);

            nShow = min(80,numel(R));
            extraCols = {}; if exShown; extraCols = strcat('ExF:', {lasers.name}); end
            rows = cell(nShow, nc+numel(extraCols)+3);
            for r = 1:nShow
                col = 1; rows{r,col} = R(r).dichroic; col = col+1;
                for k = 1:nc; rows{r,col} = R(r).filters{k}; col = col+1; end
                if exShown
                    for j = 1:nl; rows{r,col} = R(r).exfilters{j}; col = col+1; end
                end
                rows{r,col} = round(R(r).score,3,'significant'); col = col+1;
                mx = 0; bz = 0;
                for k = 1:nc
                    o = setdiff(1:numel(fluors), assign(k));
                    mx = max(mx, 100*sum(R(r).CT(o,k)));
                    bz = max(bz, 100*sum(R(r).bleed(k,:))/(R(r).S(assign(k),k)+eps));
                end
                rows{r,col} = round(mx,2); rows{r,col+1} = round(bz,3,'significant');
            end
            scoreHdr = ternStr(snrMode,'min-channel SNR','Score');
            app.OptResTable.ColumnName = [{'Splitter'}, cn, extraCols, ...
                {scoreHdr,'Max xtalk %','Max bleed %'}];
            app.OptResTable.Data = rows;
            app.LastOptResults = R; app.LastOptKind = 'joint';
            setStatus(app, sprintf('Joint optimise: %d combinations ranked (showing %d).', ...
                numel(R), nShow));
        end

        function arr = pullExCandidates(app, wl)
            % Download a few FPbase bandpass/laser-line filters near a source
            % wavelength to use as excitation-cleanup candidates ('(none)' kept).
            arr = struct('name',{'(none)'},'ex',{ones(numel(app.Lambda),1)},'floor',{NaN});
            try
                hits = FPbase.search('', 'F', 'BP', 5000);
            catch; return; end
            cen = nan(1,numel(hits));
            for i = 1:numel(hits)
                tk = regexp(hits(i).name,'(\d{3})\s*[\/-]\s*(\d{1,3})','tokens','once');
                if ~isempty(tk); cen(i) = str2double(tk{1}); end
            end
            near = find(abs(cen-wl) <= 12);
            [~,o] = sort(abs(cen(near)-wl)); near = near(o(1:min(6,numel(near))));
            for i = near
                try
                    S = FPbase.spectrum(hits(i).id, app.Lambda);
                    addToLib(app, S.name, 'Filters', 'filter', S.ex, [], [], sprintf('FPbase:%d',S.id), S.floor);
                    arr(end+1) = struct('name',S.name,'ex',S.ex(:),'floor',S.floor); %#ok<AGROW>
                catch; end
            end
        end

        function onPullCandidates(app)
            % Smart candidate generation: find the emission peak of the owner
            % fluorophore of the selected channel, then download only FPbase
            % bandpass filters whose name-encoded centre lies near that peak.
            try
                [fluors, ~, channels, assign] = assembleSystem(app);
            catch ME; setStatus(app, ME.message, true); return; end
            ck = find(strcmp({channels.name}, app.OptChanDD.Value),1);
            if isempty(ck); setStatus(app,'Compute first, then pick a channel.',true); return; end
            owner = assign(ck);
            [~, pk] = max(fluors(owner).em); peak = app.Lambda(pk);

            win = 40;   % +/- nm name-centre window around the emission peak
            hits = FPbase.search('', 'F', 'BP', 5000);
            cen = nan(1,numel(hits));
            for i = 1:numel(hits)
                tok = regexp(hits(i).name, '(\d{3})\s*[\/-]\s*(\d{1,3})', 'tokens', 'once');
                if ~isempty(tok); cen(i) = str2double(tok{1}); end
            end
            near = find(abs(cen - peak) <= win);
            if isempty(near)
                setStatus(app, sprintf('No FPbase bandpass within %dnm of %.0fnm.',win,peak),true);
                return;
            end
            [~, ord] = sort(abs(cen(near)-peak)); near = near(ord);
            near = near(1:min(25,numel(near)));   % cap downloads

            setStatus(app, sprintf('Downloading %d candidates near %.0fnm ...', ...
                numel(near), peak)); drawnow;
            added = {};
            for i = near
                try
                    S = FPbase.spectrum(hits(i).id, app.Lambda);
                    addToLib(app, S.name, 'Filters', 'filter', S.ex, [], [], ...
                        sprintf('FPbase:%d',S.id), S.floor);
                    added{end+1} = S.name; %#ok<AGROW>
                catch; end
            end
            refreshLibTree(app); populateChoosers(app);
            app.OptCandList.Items = unique([app.OptCandList.Items, added], 'stable');
            ch = app.OptChanDD.Value;
            prev = {}; if isKey(app.OptCandMap, ch); prev = app.OptCandMap(ch); end
            sel = unique([prev, added], 'stable');
            app.OptCandMap(ch) = sel;
            app.OptCandList.Value = sel;
            setStatus(app, sprintf(['Pulled %d bandpass filters near %.0fnm (%s emission) ' ...
                'into "%s" candidates.'], numel(added), peak, fluors(owner).name, ch));
        end

        function setStatus(app, msg, isErr)
            if nargin<3; isErr=false; end
            app.StatusLbl.Text = msg;
            if isErr; app.StatusLbl.FontColor=[0.7 0 0];
            else; app.StatusLbl.FontColor=[0 0.4 0]; end
        end
    end
end

function h = labeledNumO(parent, label, val)
    uilabel(parent,'Text',label);
    h = uieditfield(parent,'numeric','Value',val);
end

function u = uniqueStableC(c)
    % unique cellstr preserving first-seen order
    u = {}; for i = 1:numel(c); if ~any(strcmp(u,c{i})); u{end+1} = c{i}; end; end %#ok<AGROW>
end

function labels = compSortLabels(labels)
    % Order composition row labels by category, then by trailing index.
    cats = {'Filter set name','Fluorophore','Light source','Cleanup filter', ...
            'Combiner dichroic','Primary dichroic','Splitter dichroic', ...
            'Emission filter','Detector','Back-reflection','Element blocking'};
    key = zeros(numel(labels),1);
    for i = 1:numel(labels)
        r = numel(cats)+1;
        for c = 1:numel(cats); if startsWith(labels{i},cats{c}); r = c; break; end; end
        tok = regexp(labels{i},'(\d+)\s*$','tokens','once');
        idx = 0; if ~isempty(tok); idx = str2double(tok{1}); end
        key(i) = r*1000 + idx;
    end
    [~,ord] = sort(key); labels = labels(ord);
end

function writeCsvRow(fid, cells)
    % Write one CSV row, quoting and escaping each field.
    parts = cell(1,numel(cells));
    for i = 1:numel(cells)
        v = cells{i}; if ~ischar(v) && ~isstring(v); v = num2str(v); end
        v = char(v); v = strrep(v,'"','""');
        parts{i} = ['"' v '"'];
    end
    fprintf(fid, '%s\n', strjoin(parts, ','));
end

function c = localPeakXcorr(a, b)
    % Normalised xcorr peak in a narrow lag window on the app's 1 nm grid.
    % The narrow window catches duplicated files with tiny interpolation shifts
    % without merging legitimately shifted filters.
    a = a(:); b = b(:);
    n = min(numel(a), numel(b));
    a = a(1:n); b = b(1:n);
    a(~isfinite(a)) = 0; b(~isfinite(b)) = 0;
    if norm(a) == 0 || norm(b) == 0; c = NaN; return; end
    % Zero-lag only: duplicates share the same wavelength grid, so a match at a
    % nonzero shift means the spectra differ (e.g. two bandpasses at nearby
    % centres) and must NOT be treated as duplicates.
    maxLag = 0;
    c = -Inf;
    for lag = -maxLag:maxLag
        if lag < 0
            x = a(1-lag:end); y = b(1:end+lag);
        elseif lag > 0
            x = a(1:end-lag); y = b(1+lag:end);
        else
            x = a; y = b;
        end
        d = norm(x) * norm(y);
        if d > 0; c = max(c, dot(x,y) / d); end
    end
    c = min(1, max(-1, c));
end

function y = yPlot(app, y, yMode)
    if nargin < 3 || isempty(yMode); yMode = app.ResYMode.Value; end
    y = y(:);
    y(~isfinite(y)) = 0;
    y = max(0, y);
    if strcmp(yMode,'OD')
        y(y <= 0) = NaN;
        y = -log10(min(y,1));
    else
        y = 100 * y;
    end
end

function y = schemY(app, y)
    y = y(:);
    y(~isfinite(y)) = 0;
    y = max(0, y);
    if strcmp(app.ExYMode.Value,'OD')
        % Clamp to a display floor (OD 8) instead of letting tiny/zero
        % transmissions send -log10 to +Inf or NaN. That removes the
        % near-vertical spikes and gaps that made the OD plot look "sharp":
        % curves now approach OD 8 smoothly and stay bounded to the axis.
        y = min(max(y, 1e-8), 1);
        y = -log10(y);
    else
        y = 100 * y;
    end
end

function label = schemYLabel(app)
    if strcmp(app.ExYMode.Value,'OD')
        label = 'OD / -log10(relative)';
    else
        label = '%T / relative';
    end
end

function setTextArea(h, txt)
    % Safely set a uitextarea's text (cellstr), guarding validity.
    if isempty(h) || ~isvalid(h); return; end
    if ischar(txt) || isstring(txt); txt = cellstr(txt); end
    h.Value = txt;
end

function v = getfield_safe(S, f, dflt)
    if isstruct(S) && isfield(S, f) && ~isempty(S.(f)); v = S.(f); else; v = dflt; end
end

function label = schemAxisLabel(app, which)
    % Per-axis y-label: top plot = transmission delivered to the specimen,
    % bottom plot = transmission delivered to the detectors.
    if strcmp(which,'specimen')
        base = 'Transmission-to-Specimen';
    else
        base = 'Transmission-to-Detectors';
    end
    if strcmp(app.ExYMode.Value,'OD')
        label = [base ' (OD)'];
    else
        label = [base ' (%/relative)'];
    end
end

function clearSchematicPlotAxes(ax)
    delete(findall(ax,'Type','Line'));
    legend(ax,'off');
end

function y = normPeak(y)
    y = y(:);
    y(~isfinite(y)) = 0;
    m = max(y);
    if m > 0; y = y / m; end
end

function name = cleanFileName(name)
    name = regexprep(name, '[<>:"/\\|?*]', '_');
    name = regexprep(name, '\s+', '_');
    name = regexprep(name, '_+', '_');
    name = strtrim(name);
end

function writeSchematicManifest(fp, sch)
    fid = fopen(fp, 'w');
    if fid < 0; error('Could not write %s', fp); end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '# filterDesignerApp filter set\n');
    fprintf(fid, 'type=schematic\n');
    fprintf(fid, 'version=%g\n', getScalarField(sch, 'version', 2));
    fprintf(fid, 'name=%s\n', enc(getTextField(sch, 'name', 'filter_set')));
    fprintf(fid, 'savedOn=%s\n', enc(datestr(now)));
    writeCell(fid, 'srcVals', getCellField(sch, 'srcVals', {}));
    writeCell(fid, 'cleanVals', getCellField(sch, 'cleanVals', {}));
    writeCell(fid, 'combinerVals', getCellField(sch, 'combinerVals', {}));
    fprintf(fid, 'primary=%s\n', enc(getTextField(sch, 'primary', '(none)')));
    writeCell(fid, 'fluorVals', getCellField(sch, 'fluorVals', {}));
    writeCell(fid, 'emVals', getCellField(sch, 'emVals', {}));
    writeCell(fid, 'detVals', getCellField(sch, 'detVals', {}));
    writeCell(fid, 'ownerVals', getCellField(sch, 'ownerVals', {}));
    writeCell(fid, 'splitterVals', getCellField(sch, 'splitterVals', {}));
    fprintf(fid, 'Rback=%g\n', getScalarField(sch, 'Rback', 0.5));
    fprintf(fid, 'blockOD=%g\n', getScalarField(sch, 'blockOD', 6));
    fprintf(fid, 'opticsMode=%s\n', enc(getTextField(sch, 'opticsMode', 'Combine spectra')));
    clear cleaner;
end

function sch = readSchematicManifest(fp)
    raw = fileread(fp);
    lines = regexp(raw, '\r\n|\n|\r', 'split');
    vals = containers.Map('KeyType','char','ValueType','char');
    for ii = 1:numel(lines)
        line = strtrim(lines{ii});
        if isempty(line) || startsWith(line, '#'); continue; end
        eq = strfind(line, '=');
        if isempty(eq); continue; end
        key = strtrim(line(1:eq(1)-1));
        vals(key) = strtrim(line(eq(1)+1:end));
    end
    sch = struct('version', readNum(vals,'version',2), ...
                 'name', readText(vals,'name','filter_set'));
    sch.srcVals = readCell(vals, 'srcVals', {'488 nm laser'});
    sch.cleanVals = readCell(vals, 'cleanVals', {'(none)'});
    sch.combinerVals = readCell(vals, 'combinerVals', {});
    sch.primary = readText(vals, 'primary', '(none)');
    sch.fluorVals = readCell(vals, 'fluorVals', {'mNeonGreen'});
    sch.emVals = readCell(vals, 'emVals', {'(none)'});
    sch.detVals = readCell(vals, 'detVals', {'(ideal)'});
    sch.ownerVals = readCell(vals, 'ownerVals', sch.fluorVals(1));
    sch.splitterVals = readCell(vals, 'splitterVals', {});
    sch.Rback = readNum(vals, 'Rback', 0.5);
    sch.blockOD = readNum(vals, 'blockOD', 6);
    sch.opticsMode = readText(vals, 'opticsMode', 'Combine spectra');
end

function clearLegacyFilterSetSpectra(setDir)
    if ~isfolder(setDir); return; end
    files = [dir(fullfile(setDir,'*.txt')); dir(fullfile(setDir,'*.TXT')); ...
             dir(fullfile(setDir,'*.csv')); dir(fullfile(setDir,'*.CSV'))];
    for ii = 1:numel(files)
        if strcmpi(files(ii).name, 'filter_set.txt'); continue; end
        fp = fullfile(files(ii).folder, files(ii).name);
        try
            delete(fp);
        catch
        end
    end
end

function v = jsonText(s, field, default)
    if nargin < 3; default = ''; end
    if ~isfield(s, field) || isempty(s.(field))
        v = default;
    else
        v = char(s.(field));
    end
end

function fp = canonicalPath(fp)
    try
        fp = char(java.io.File(fp).getCanonicalPath());
    catch
        fp = char(fp);
    end
end

function tf = startsWithPath(fp, root)
    fp = canonicalPath(fp);
    root = canonicalPath(root);
    if ispc
        fp = lower(fp);
        root = lower(root);
    end
    tf = strcmp(fp, root) || startsWith(fp, [root filesep]);
end

function tf = isInsideDir(fp, root)
    if isempty(fp) || isempty(root); tf = false; return; end
    tf = startsWithPath(fp, root);
end

function tf = isAbsolutePath(fp)
    fp = char(fp);
    tf = startsWith(fp, filesep) || ~isempty(regexp(fp, '^[A-Za-z]:[\\/]', 'once'));
end

function writeSpectrumTable(fp, lambda, y1, y2)
    fid = fopen(fp, 'w');
    if fid < 0; error('Could not write %s', fp); end
    cleaner = onCleanup(@() fclose(fid));
    if nargin >= 4 && ~isempty(y2)
        fprintf(fid, 'Wavelength_nm\tExcitation\tEmission\n');
        data = [lambda(:), y1(:), y2(:)];
        fprintf(fid, '%.0f\t%.8g\t%.8g\n', data');
    else
        fprintf(fid, 'Wavelength_nm\tTransmission\n');
        data = [lambda(:), y1(:)];
        fprintf(fid, '%.0f\t%.8g\n', data');
    end
    clear cleaner;
end

function tf = isPresetPlaceholder(name)
    if isstring(name); name = char(name); end
    tf = isempty(name) || strcmp(name,'(none)') || strcmp(name,'(ideal)') || strcmp(name,'(laser line)');
end

function writeCell(fid, key, c)
    fprintf(fid, '%s=%s\n', key, strjoin(cellfun(@enc, rowCell(c), 'UniformOutput', false), '|'));
end

function c = readCell(vals, key, default)
    if ~isKey(vals, key) || isempty(vals(key))
        c = default;
        return;
    end
    c = regexp(vals(key), '\|', 'split');
    c = cellfun(@dec, c, 'UniformOutput', false);
end

function s = readText(vals, key, default)
    if isKey(vals, key); s = dec(vals(key)); else; s = default; end
end

function x = readNum(vals, key, default)
    if isKey(vals, key); x = str2double(vals(key)); else; x = default; end
    if ~isfinite(x); x = default; end
end

function c = getCellField(st, f, default)
    if isfield(st, f); c = rowCell(st.(f)); else; c = default; end
end

function s = getTextField(st, f, default)
    if isfield(st, f) && ~isempty(st.(f)); s = char(string(st.(f))); else; s = default; end
end

function x = getScalarField(st, f, default)
    if isfield(st, f) && isnumeric(st.(f)) && isscalar(st.(f)); x = st.(f); else; x = default; end
end

function c = rowCell(c)
    if isempty(c); c = {}; return; end
    if isstring(c); c = cellstr(c); end
    if ischar(c); c = {c}; end
    c = c(:)';
end

function s = enc(s)
    if isstring(s); s = char(s); end
    s = strrep(s, '%', '%25');
    s = strrep(s, '|', '%7C');
    s = strrep(s, newline, ' ');
end

function s = dec(s)
    s = strrep(s, '%7C', '|');
    s = strrep(s, '%25', '%');
end

function setDropValue(dd, value)
    if isempty(value); return; end
    if ~ismember(value, dd.Items)
        dd.Items = [dd.Items, {value}];
    end
    dd.Value = value;
end

function s = limitingSource(out, k)
    comps = [out.xtalk_e(k), out.exBleed_e(k), out.dark_e(k), out.read_e(k)^2];
    names = {'crosstalk','laser bleed','dark','read'};
    [mx, ix] = max(comps);
    if out.read_e(k)^2 > out.signal_e(k) + out.xtalk_e(k) + out.exBleed_e(k) + out.dark_e(k)
        s = 'read';
    elseif mx <= out.signal_e(k)
        s = 'shot';
    else
        s = names{ix};
    end
end

function s = ternStr(cond, a, b)
    if cond; s = a; else; s = b; end
end

function c = padCell(c, n, default)
    if nargin < 3; default = '(none)'; end
    if ~iscell(c); c = {}; end
    if numel(c) > n; c = c(1:n);
    else; for i = numel(c)+1:n; c{i} = default; end
    end
end

function c = sanitizeChoices(c, items, default)
    if ~iscell(c); c = {c}; end
    for i = 1:numel(c)
        if isempty(c{i}) || ~ismember(c{i}, items)
            c{i} = default;
        end
    end
end

function rgb = wl2rgb(wl)
%WL2RGB  Approximate visible-spectrum colour for a wavelength (nm) via a smooth
%   hue (HSV) mapping. The hue varies STEEPLY through the green band (490-545 nm)
%   so closely-spaced green fluorophores (e.g. GFP ~507, mNeonGreen ~517, YFP
%   ~527) get visibly different hues (teal-green -> green -> chartreuse) instead
%   of collapsing onto the same pure green.
    if wl < 380 || wl > 780; rgb = [0.6 0.6 0.6]; return; end
    % anchor wavelengths -> hue (degrees: 0 red, 60 yellow, 120 green,
    % 180 cyan, 240 blue, 280 violet). Extra anchors across 490-545 spread green.
    wlA  = [380 440 470 490 500 510 520 530 545 565 590 620 660 780];
    hueA = [282 240 212 184 165 142 118  98  76  56  32   6   0   0];
    h = interp1(wlA, hueA, wl, 'pchip');
    rgb = hsv2rgb([max(0,min(1,h/360)), 1, 1]);
    % fade only the deep violet / deep red ends so they stay on-screen
    if     wl < 420; f = 0.45 + 0.55*(wl-380)/40;
    elseif wl > 700; f = 0.45 + 0.55*(780-wl)/80;
    else;            f = 1;
    end
    rgb = min(1, max(0, rgb*f));
    % keep colours visible on a dark background
    if max(rgb) < 0.4; rgb = rgb + (0.4 - max(rgb)); end
end

function rgb = dichroicColor(lambda, T)
%DICHROICCOLOR  Colour a dichroic by its strongest spectral edge.
    lambda = lambda(:); T = T(:);
    n = min(numel(lambda), numel(T));
    lambda = lambda(1:n); T = T(1:n);
    good = isfinite(lambda) & isfinite(T);
    lambda = lambda(good); T = T(good);
    if numel(lambda) < 3
        rgb = [0.75 0.75 0.75]; return;
    end
    T = min(max(T,0),1);
    dT = abs(diff(T));
    % Ignore tiny noise and prefer real transition bands in the visible range.
    visible = lambda(1:end-1) >= 380 & lambda(1:end-1) <= 780;
    dT(~visible) = 0;
    [mx, ix] = max(dT);
    if mx <= 1e-4
        [~, ix] = min(abs(T - 0.5));
    end
    edge = lambda(min(ix, numel(lambda)));
    rgb = wl2rgb(edge);
    rgb = 0.85*rgb + 0.15*[1 1 1];  % keep dashed line readable on dark axes
end
