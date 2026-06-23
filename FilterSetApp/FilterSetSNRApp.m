classdef FilterSetSNRApp < matlab.apps.AppBase
%FILTERSETSNRAPP  Signal-to-noise calculator for a chosen fluorescence filter
%   set, following the IDEX/Semrock spectral model and extending it with
%   photon shot noise, detector read noise, dark current, and tissue / fibre
%   autofluorescence.
%
%   Companion to OptimalFilterApp: in that app design your system and click
%   "Save config → SNR app", then here click "Load config" (or launch with a
%   config struct):
%       >> FilterSetSNRApp                 % then Load config (.mat)
%       >> FilterSetSNRApp(app.config)     % pass a config directly
%
%   Depends on: snrModel.m, autofluorPreset.m, FilterEngine.m.

    properties (Access = public)
        Fig
        Cfg
        Tabs
        CfgLabel
        PhysFields      struct = struct()
        PowerTable
        FluorTable
        AFTable
        ResTable
        BarAxes
        SpecAxes
        StatusLbl
        LastOut
        SwParam
        SwMin
        SwMax
        SwSteps
        SweepAxes
        SweepAxes2
    end

    methods (Access = public)
        function app = FilterSetSNRApp(cfg)
            buildUI(app);
            if nargin >= 1 && ~isempty(cfg); loadConfig(app, cfg); end
            if nargout == 0; clear app; end
        end

        function out = run(app)
            %RUN  Public trigger for the SNR computation (= "Compute SNR").
            onCompute(app); out = app.LastOut;
        end

        function runSweep(app)
            %RUNSWEEP  Public trigger for the parameter sweep.
            onSweep(app);
        end
    end

    methods (Access = private)
        function buildUI(app)
            ss = get(0,'ScreenSize');
            w = min(1600, ss(3)*0.90); h = min(1040, ss(4)*0.90);
            x = ss(1) + (ss(3)-w)/2; y = ss(2) + (ss(4)-h)/2;
            app.Fig = uifigure('Name','Filter-Set SNR Calculator','Position',[x y w h]);
            g = uigridlayout(app.Fig,[2 1]); g.RowHeight = {'1x',28};
            app.Tabs = uitabgroup(g);
            app.StatusLbl = uilabel(g,'Text','Load a configuration saved from the filter app.', ...
                'FontColor',[0 0.4 0]);
            buildSetupTab(app); buildResultsTab(app); buildSweepTab(app);
            hs = findall(app.Fig,'-property','FontSize');   % +2 pt everywhere
            for i = 1:numel(hs); try; hs(i).FontSize = hs(i).FontSize + 2; catch; end; end
        end

        function buildSetupTab(app)
            t = uitab(app.Tabs,'Title','Setup');
            g = uigridlayout(t,[4 2]);
            g.RowHeight = {32,150,'1x',38}; g.ColumnWidth = {'1x','1x'};

            top = uigridlayout(g,[1 2]); top.Layout.Row=1; top.Layout.Column=[1 2];
            top.ColumnWidth = {140,'1x'};
            uibutton(top,'Text','Load config (.mat)','ButtonPushedFcn',@(s,e)onLoad(app));
            app.CfgLabel = uilabel(top,'Text','(no configuration loaded)','FontColor',[0.4 0.4 0.4]);

            % physics parameters
            p1 = uipanel(g,'Title','Acquisition & detector'); p1.Layout.Row=2; p1.Layout.Column=1;
            pg = uigridlayout(p1,[4 4]);
            app.PhysFields.NA          = labeledNum(pg,'NA',0.5);
            app.PhysFields.n           = labeledNum(pg,'index n',1.33);
            app.PhysFields.pathUm      = labeledNum(pg,'path (µm)',100);
            app.PhysFields.tIntMs      = labeledNum(pg,'integ. (ms)',10);
            app.PhysFields.readNoise   = labeledNum(pg,'read noise (e-)',1.5);
            app.PhysFields.darkRate    = labeledNum(pg,'dark (e-/s)',100);

            % source powers
            p2 = uipanel(g,'Title','Excitation power'); p2.Layout.Row=2; p2.Layout.Column=2;
            g2 = uigridlayout(p2,[1 1]);
            app.PowerTable = uitable(g2,'ColumnName',{'Source','Power (mW)'}, ...
                'ColumnEditable',[false true],'Data',cell(0,2));

            % fluorophore photophysics
            p3 = uipanel(g,'Title','Fluorophore photophysics'); p3.Layout.Row=3; p3.Layout.Column=1;
            g3 = uigridlayout(p3,[1 1]);
            app.FluorTable = uitable(g3, ...
                'ColumnName',{'Fluorophore','EC (M^-1cm^-1)','QY','Conc (µM)'}, ...
                'ColumnEditable',[false true true true],'Data',cell(0,4));

            % autofluorescence
            p4 = uipanel(g,'Title','Autofluorescence sources (excitation-dependent, broad emission)');
            p4.Layout.Row=3; p4.Layout.Column=2;
            g4 = uigridlayout(p4,[2 1]); g4.RowHeight={'1x',28};
            app.AFTable = uitable(g4, ...
                'ColumnName',{'Type','Strength (abs. frac.)','QY'}, ...
                'ColumnEditable',[true true true],'Data',cell(0,3));
            app.AFTable.ColumnFormat = {[autofluorPreset('list',(350:850)') {'(none)'}],'numeric','numeric'};
            gb = uigridlayout(g4,[1 2]);
            uibutton(gb,'Text','Add','ButtonPushedFcn',@(s,e)addRow(app,app.AFTable,{'Brain tissue',2e-3,0.1}));
            uibutton(gb,'Text','Remove last','ButtonPushedFcn',@(s,e)delRow(app,app.AFTable));

            cb = uibutton(g,'Text','Compute SNR  ▶','FontWeight','bold', ...
                'ButtonPushedFcn',@(s,e)onCompute(app));
            cb.Layout.Row=4; cb.Layout.Column=[1 2];
        end

        function buildResultsTab(app)
            t = uitab(app.Tabs,'Title','Results');
            g = uigridlayout(t,[2 2]); g.RowHeight={'1.1x','1x'}; g.ColumnWidth={'1x','1x'};
            p = uipanel(g,'Title','Per-channel signal, noise (electrons) and SNR');
            p.Layout.Row=1; p.Layout.Column=[1 2];
            pg = uigridlayout(p,[1 2]); pg.ColumnWidth={'1x',130};
            app.ResTable = uitable(pg);
            uibutton(pg,'Text','Export → CSV','ButtonPushedFcn',@(s,e)onExport(app));
            app.BarAxes = uiaxes(g); app.BarAxes.Layout.Row=2; app.BarAxes.Layout.Column=1;
            title(app.BarAxes,'Noise variance breakdown (e-^2)');
            app.SpecAxes = uiaxes(g); app.SpecAxes.Layout.Row=2; app.SpecAxes.Layout.Column=2;
            title(app.SpecAxes,'Autofluorescence emission vs channels');
            xlabel(app.SpecAxes,'Wavelength (nm)');
        end

        function buildSweepTab(app)
            t = uitab(app.Tabs,'Title','Sweep');
            g = uigridlayout(t,[2 1]); g.RowHeight = {44,'1x'};
            ctl = uigridlayout(g,[1 8]);
            ctl.ColumnWidth = {55,'1x',40,70,40,70,40,60};
            uilabel(ctl,'Text','Sweep:');
            app.SwParam = uidropdown(ctl,'Items',{'Integration time (ms)','Laser power (mW)', ...
                'Concentration (µM)','Read noise (e-)','NA','Autofluor strength x'});
            uilabel(ctl,'Text','min'); app.SwMin = uieditfield(ctl,'numeric','Value',0.1);
            uilabel(ctl,'Text','max'); app.SwMax = uieditfield(ctl,'numeric','Value',100);
            uilabel(ctl,'Text','pts'); app.SwSteps = uieditfield(ctl,'numeric','Value',40);
            uibutton(ctl,'Text','Run ▶','ButtonPushedFcn',@(s,e)onSweep(app));

            gg = uigridlayout(g,[1 2]); gg.Layout.Row=2;
            app.SweepAxes  = uiaxes(gg);
            app.SweepAxes2 = uiaxes(gg);
            xlabel(app.SweepAxes,'parameter'); ylabel(app.SweepAxes,'SNR');
            title(app.SweepAxes,'SNR vs parameter (per channel)');
            xlabel(app.SweepAxes2,'parameter'); ylabel(app.SweepAxes2,'noise (e-)');
            title(app.SweepAxes2,'Shot vs read noise (worst channel)');
        end

        function onSweep(app)
            if isempty(app.Cfg); setStatus(app,'Load a configuration first.',true); return; end
            cfg = app.Cfg; [phys, fp, af] = gatherParams(app);
            lo = app.SwMin.Value; hi = app.SwMax.Value; n = max(2,round(app.SwSteps.Value));
            if lo<=0 || hi<=lo; setStatus(app,'Need 0 < min < max.',true); return; end
            logx = (hi/lo) >= 50;             % auto log axis for wide ranges
            if logx; vals = logspace(log10(lo),log10(hi),n); else; vals = linspace(lo,hi,n); end
            par = app.SwParam.Value;

            nC = numel(cfg.channels); cnames = {cfg.channels.name};
            SNR = zeros(n,nC); shot = zeros(n,nC); read = zeros(n,nC);
            for q = 1:n
                [ph,fpq,afq] = applyOverride(phys, fp, af, par, vals(q));
                o = snrModel(cfg, ph, fpq, afq);
                SNR(q,:) = o.SNR; shot(q,:) = o.shot_e; read(q,:) = o.read_e;
            end

            cla(app.SweepAxes); hold(app.SweepAxes,'on'); co = lines(nC);
            for k=1:nC; plot(app.SweepAxes, vals, SNR(:,k),'-o','Color',co(k,:), ...
                    'LineWidth',1.5,'MarkerSize',3); end
            hold(app.SweepAxes,'off'); grid(app.SweepAxes,'on');
            legend(app.SweepAxes, cnames, 'Location','best','Interpreter','none');
            xlabel(app.SweepAxes, par); ylabel(app.SweepAxes,'SNR');
            if logx; app.SweepAxes.XScale='log'; else; app.SweepAxes.XScale='linear'; end

            % worst channel: shot vs read, mark crossover
            [~,wc] = min(SNR(end,:));
            cla(app.SweepAxes2); hold(app.SweepAxes2,'on');
            plot(app.SweepAxes2, vals, shot(:,wc),'-','LineWidth',1.8);
            plot(app.SweepAxes2, vals, read(:,wc),'--','LineWidth',1.8);
            d = shot(:,wc)-read(:,wc); ix = find(d(1:end-1).*d(2:end) <= 0, 1);
            xc = NaN;
            if ~isempty(ix)
                xc = interp1(d(ix:ix+1), vals(ix:ix+1), 0);
                xline(app.SweepAxes2, xc, '-.','Color',[.4 .4 .4], ...
                    'Label','read = shot','LabelOrientation','horizontal');
            end
            hold(app.SweepAxes2,'off'); grid(app.SweepAxes2,'on');
            legend(app.SweepAxes2,{'shot noise','read noise'},'Location','best');
            xlabel(app.SweepAxes2, par); ylabel(app.SweepAxes2,'noise (e-)');
            if logx; app.SweepAxes2.XScale='log'; app.SweepAxes2.YScale='log';
            else; app.SweepAxes2.XScale='linear'; app.SweepAxes2.YScale='linear'; end
            title(app.SweepAxes2, sprintf('Shot vs read (%s)', cnames{wc}));

            if isnan(xc); msg = 'no read/shot crossover in range';
            else; msg = sprintf('read=shot at %s = %.3g', par, xc); end
            setStatus(app, ['Swept ' par ' over ' num2str(n) ' points; ' msg '.']);
        end

        %% ---- config loading ----
        function onLoad(app)
            [f,p] = uigetfile('*.mat','Load filter-set configuration');
            if isequal(f,0); return; end
            S = load(fullfile(p,f));
            if ~isfield(S,'cfg'); setStatus(app,'File has no ''cfg'' variable.',true); return; end
            loadConfig(app, S.cfg);
            setStatus(app, ['Loaded ' f]);
        end

        function loadConfig(app, cfg)
            app.Cfg = cfg;
            app.CfgLabel.Text = sprintf('%d fluorophores | %d channels | %d sources | Rback=%.2g%% | OD cap=%g', ...
                numel(cfg.fluors), numel(cfg.channels), numel(cfg.lasers), 100*cfg.Rback, cfg.blockOD);
            app.CfgLabel.FontColor = [0 0 0];
            % powers
            pd = cell(numel(cfg.lasers),2);
            for j=1:numel(cfg.lasers); pd(j,:) = {cfg.lasers(j).name, 1}; end
            app.PowerTable.Data = pd;
            % fluor photophysics (defaults from FPbase ec/qy when present)
            fd = cell(numel(cfg.fluors),4);
            for i=1:numel(cfg.fluors)
                ec = 50000; qy = 0.6;
                if isfield(cfg.fluors,'ec') && isfinite(cfg.fluors(i).ec); ec = cfg.fluors(i).ec; end
                if isfield(cfg.fluors,'qy') && isfinite(cfg.fluors(i).qy); qy = cfg.fluors(i).qy; end
                fd(i,:) = {cfg.fluors(i).name, ec, qy, 10};
            end
            app.FluorTable.Data = fd;
            if isempty(app.AFTable.Data)
                app.AFTable.Data = {'Brain tissue',2e-3,0.1; 'Silica fiber',5e-4,0.05};
            end
        end

        %% ---- compute ----
        function [phys, fp, af] = gatherParams(app)
            lambda = app.Cfg.lambda;
            phys.NA = num(app,app.PhysFields.NA); phys.n = num(app,app.PhysFields.n);
            phys.pathLength_cm = num(app,app.PhysFields.pathUm)*1e-4;
            phys.tInt_s = num(app,app.PhysFields.tIntMs)*1e-3;
            phys.readNoise_e = num(app,app.PhysFields.readNoise);
            phys.darkRate_eps = num(app,app.PhysFields.darkRate);
            pd = app.PowerTable.Data; phys.powers_mW = cell2mat(pd(:,2))';
            fd = app.FluorTable.Data; fp = struct('ec',{},'qy',{},'conc_M',{});
            for i=1:size(fd,1)
                fp(i).ec = numv(fd{i,2}); fp(i).qy = numv(fd{i,3}); fp(i).conc_M = numv(fd{i,4})*1e-6;
            end
            ad = app.AFTable.Data; af = struct('name',{},'absorb',{},'em',{},'strength',{},'qy',{});
            for m=1:size(ad,1)
                nm = ad{m,1}; if strcmp(nm,'(none)'); continue; end
                s = autofluorPreset(nm, lambda);
                s.strength = numv(ad{m,2}); s.qy = numv(ad{m,3});
                af(end+1) = s; %#ok<AGROW>
            end
        end

        function onCompute(app)
            if isempty(app.Cfg); setStatus(app,'Load a configuration first.',true); return; end
            cfg = app.Cfg;
            [phys, fp, af] = gatherParams(app);
            out = snrModel(cfg, phys, fp, af);
            app.LastOut = out;
            fillResults(app, out, af);
            parts = arrayfun(@(k) sprintf('%s=%.1f', out.channels{k}, out.SNR(k)), ...
                1:numel(out.channels), 'UniformOutput', false);
            setStatus(app, ['Computed. SNR per channel: ' strjoin(parts, ', ')]);
            app.Tabs.SelectedTab = app.Tabs.Children(2);
        end

        function fillResults(app, out, af)
            cn = out.channels; nC = numel(cn);
            rowNames = {'Signal (e-)','Crosstalk (e-)','Excitation bleed (e-)', ...
                'Autofluor. (e-)','Dark (e-)','Read noise (e-)','Shot noise (e-)', ...
                'Total noise (e-)','SNR','SNR (optical only)'};
            M = [out.signal_e; out.xtalk_e; out.exBleed_e; out.af_e; out.dark_e; ...
                 out.read_e; out.shot_e; out.noise_e; out.SNR; out.SNR_optical];
            app.ResTable.RowName = rowNames; app.ResTable.ColumnName = cn;
            app.ResTable.Data = arrayfun(@(x) round(x,3,'significant'), M);

            % noise variance breakdown (e-^2)
            comps = [out.signal_e; out.xtalk_e; out.exBleed_e; out.af_e; out.dark_e; out.read_e.^2];
            cla(app.BarAxes);
            bar(app.BarAxes, categorical(cn), comps', 'stacked');
            legend(app.BarAxes, {'signal shot','crosstalk','exc. bleed','autofluor','dark','read^2'}, ...
                'Location','bestoutside'); ylabel(app.BarAxes,'variance (e-^2)');
            app.BarAxes.YScale = 'log';

            % spectral overlay: AF emission (area→peak) vs channel transmission
            cla(app.SpecAxes); hold(app.SpecAxes,'on'); lam = app.Cfg.lambda;
            co = lines(nC);
            for k=1:nC
                T = FilterEngine.pathTransmission(app.Cfg.channels(k), lam);
                plot(app.SpecAxes, lam, T, 'Color',co(k,:),'LineWidth',2);
            end
            for m=1:numel(af)
                plot(app.SpecAxes, lam, af(m).em/max(af(m).em),'--','LineWidth',1.2);
            end
            for j=1:numel(app.Cfg.lasers)
                xline(app.SpecAxes, app.Cfg.lasers(j).wl,'-.','Color',[.4 .4 .4]);
            end
            hold(app.SpecAxes,'off'); xlim(app.SpecAxes,[450 750]);
            legend(app.SpecAxes, [cn, {af.name}], 'Location','bestoutside','Interpreter','none');
        end

        function onExport(app)
            if isempty(app.LastOut); setStatus(app,'Compute first.',true); return; end
            [f,p] = uiputfile('*.csv','Export SNR results','snr_results.csv');
            if isequal(f,0); return; end
            o = app.LastOut; fid = fopen(fullfile(p,f),'w'); c = onCleanup(@()fclose(fid));
            fprintf(fid,'Filter-set SNR results\n,%s\n', strjoin(o.channels,','));
            rows = {'Signal (e-)',o.signal_e; 'Crosstalk (e-)',o.xtalk_e; ...
                'Excitation bleed (e-)',o.exBleed_e; 'Autofluorescence (e-)',o.af_e; ...
                'Dark (e-)',o.dark_e; 'Read noise (e-)',o.read_e; 'Shot noise (e-)',o.shot_e; ...
                'Total noise (e-)',o.noise_e; 'SNR',o.SNR; 'SNR optical',o.SNR_optical};
            for r=1:size(rows,1)
                fprintf(fid,'%s', rows{r,1}); fprintf(fid,',%.5g', rows{r,2}); fprintf(fid,'\n');
            end
            setStatus(app, ['Exported ' fullfile(p,f)]);
        end

        %% ---- helpers ----
        function v = num(~, h); v = h.Value; end
        function addRow(~, tbl, t); tbl.Data = [tbl.Data; t]; end
        function delRow(~, tbl); if ~isempty(tbl.Data); tbl.Data(end,:)=[]; end; end
        function setStatus(app, msg, isErr)
            if nargin<3; isErr=false; end
            app.StatusLbl.Text = msg;
            if isErr; app.StatusLbl.FontColor=[0.7 0 0]; else; app.StatusLbl.FontColor=[0 0.4 0]; end
        end
    end
end

function h = labeledNum(parent, label, val)
    uilabel(parent,'Text',label);
    h = uieditfield(parent,'numeric','Value',val);
end

function v = numv(x)
    if ischar(x) || isstring(x); v = str2double(x); else; v = x; end
end

function [phys, fp, af] = applyOverride(phys, fp, af, par, v)
    switch par
        case 'Integration time (ms)';  phys.tInt_s = v*1e-3;
        case 'Laser power (mW)';        phys.powers_mW = v*ones(size(phys.powers_mW));
        case 'Concentration (µM)';      for i=1:numel(fp); fp(i).conc_M = v*1e-6; end
        case 'Read noise (e-)';         phys.readNoise_e = v;
        case 'NA';                      phys.NA = v;
        case 'Autofluor strength x';    for m=1:numel(af); af(m).strength = af(m).strength*v; end
    end
end
