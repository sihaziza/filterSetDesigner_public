classdef FilterEngine
%FILTERENGINE  General optics math for multi-channel fluorescence filter sets.
%
%   A SYSTEM is described by:
%     fluors(i)  : struct with .name .ex .em (on common lambda) .brightness
%     lasers(j)  : struct with .name .wl (nm) .power (a.u.)  [excitation lines]
%     channels(k): struct with .name
%                            .emFilter  (transmission vector on lambda)
%                            .path      struct array of {.T (vector), .mode}
%                                       mode = 'T' (transmit) or 'R' (reflect)
%
%   All transmission/reflection vectors are on the same lambda grid.
%
%   Core outputs:
%     S(i,k)  = detected signal from fluorophore i in channel k (summed over
%               all lasers, brightness- and excitation-weighted)
%     eff(i,k)= fraction of fluorophore i's emission collected by channel k
%     CT      = crosstalk matrix (signal leaked / signal in own channel)

    methods (Static)

        function t = pathTransmission(channel, lambda)
            % Product of emission filter and each dichroic (T or 1-T).
            t = ones(numel(lambda),1);
            if isfield(channel,'emFilter') && ~isempty(channel.emFilter)
                t = t .* channel.emFilter(:);
            end
            if isfield(channel,'path')
                for n = 1:numel(channel.path)
                    p = channel.path(n);
                    if strcmpi(p.mode,'R')
                        t = t .* (1 - p.T(:));
                    else
                        t = t .* p.T(:);
                    end
                end
            end
        end

        function E = sourceExcitation(src, lambda)
            %SOURCEEXCITATION Effective excitation spectrum E(λ) hitting the
            %   sample, normalised so trapz(E)=power, after the cleanup filter.
            %   A source is either a laser LINE (.wl + .power) or a broadband
            %   spectrum (.spectrum, e.g. an LED/lamp), optionally with an
            %   excitation/cleanup filter (.exFilter).
            lambda = lambda(:); dl = mean(diff(lambda));
            if isfield(src,'E') && ~isempty(src.E); E = src.E(:); return; end
            p = 1; if isfield(src,'power') && ~isempty(src.power); p = src.power; end
            if isfield(src,'spectrum') && ~isempty(src.spectrum)
                s = src.spectrum(:); a = trapz(lambda, s); if a>0; s = s/a; end
                E = p * s;
            else
                [~, idx] = min(abs(lambda - src.wl));
                E = zeros(numel(lambda),1); E(idx) = p/dl;   % unit-area impulse
            end
            if isfield(src,'exFilter') && ~isempty(src.exFilter)
                E = E .* src.exFilter(:);
            end
        end

        function [S, eff] = signalMatrix(fluors, lasers, channels, lambda, detector)
            % detector (optional): QE vector on lambda (0-1); applied to every
            % channel. [] or omitted => ideal detector (QE = 1).
            if nargin < 5 || isempty(detector); detector = ones(numel(lambda),1); end
            detector = detector(:); lambda = lambda(:);
            nf = numel(fluors); nc = numel(channels);
            % effective excitation spectrum per source (incl. cleanup filter)
            E = cell(1,numel(lasers));
            for j = 1:numel(lasers); E{j} = FilterEngine.sourceExcitation(lasers(j), lambda); end
            S = zeros(nf, nc); eff = zeros(nf, nc);
            hasDet = isfield(channels,'detector');
            for k = 1:nc
                det = detector;                          % per-channel QE override
                if hasDet && ~isempty(channels(k).detector); det = channels(k).detector(:); end
                t = FilterEngine.pathTransmission(channels(k), lambda) .* det;
                for i = 1:nf
                    em = fluors(i).em(:);
                    emTot = trapz(lambda, em) + eps;
                    collected = trapz(lambda, em .* t);
                    eff(i,k) = collected / emTot;            % geometry-free fraction
                    exW = 0;                                 % absorbed excitation
                    for j = 1:numel(lasers)
                        exW = exW + trapz(lambda, E{j} .* fluors(i).ex(:));
                    end
                    b = 1; if isfield(fluors,'brightness') && ~isempty(fluors(i).brightness) ...
                              && isfinite(fluors(i).brightness); b = fluors(i).brightness; end
                    S(i,k) = b * exW * collected;
                end
            end
        end

        function bleed = laserBleed(channels, lasers, lambda, Rback, detector, blockOD)
            %LASERBLEED Back-reflected excitation light reaching each detector.
            %   A fraction Rback (default 0.005 = 0.5%) of every source's light
            %   back-scatters off non-AR-coated optics into the detection path,
            %   where it is attenuated by that channel's dichroics + emission
            %   filter + detector QE.  bleed(k,j) is in the same units as the
            %   signal matrix, so bleed(k)/S(owner,k) is a true
            %   background-to-signal ratio (incident flux cancels).
            %
            %   blockOD (default 6): realistic per-element blocking floor. Real
            %   hard-coated filters/dichroics are specified to a finite optical
            %   density (e.g. OD6 => 1e-6 transmission); measured spectra often
            %   round to 0 out of band, which would over-state blocking. Each
            %   element's transmission is floored at 10^-blockOD before use so
            %   the bleed estimate is not artificially optimistic.
            if nargin < 4 || isempty(Rback);    Rback = 0.005; end
            if nargin < 5 || isempty(detector); detector = ones(numel(lambda),1); end
            if nargin < 6 || isempty(blockOD);  blockOD = 6; end
            detector = detector(:); lambda = lambda(:);
            nc = numel(channels); nl = numel(lasers);
            hasDet = isfield(channels,'detector');
            bleed = zeros(nc, nl);
            for k = 1:nc
                det = detector;
                if hasDet && ~isempty(channels(k).detector); det = channels(k).detector(:); end
                t = FilterEngine.pathTransmissionFloored(channels(k), lambda, blockOD) .* det;
                for j = 1:nl
                    E = FilterEngine.sourceExcitation(lasers(j), lambda);
                    bleed(k,j) = Rback * trapz(lambda, E .* t);
                end
            end
        end

        function t = pathTransmissionFloored(channel, lambda, blockOD)
            % Like pathTransmission but each optical element's deep-blocking
            % region is floored at its OWN measured blocking depth (channel
            % .emFloor / path(n).floor, from loadSpectrum/FPbase), with the
            % global cap 10^-blockOD limiting credited optimism. Elements
            % whose floor is NaN ("no deep-blocking data") fall back to the cap.
            cap = 10^(-blockOD);
            t = ones(numel(lambda),1);
            if isfield(channel,'emFilter') && ~isempty(channel.emFilter)
                fe = FilterEngine.elementFloor(getfloor(channel,'emFloor'), cap);
                t = t .* max(channel.emFilter(:), fe);
            end
            if isfield(channel,'path')
                for n = 1:numel(channel.path)
                    p = channel.path(n);
                    if strcmpi(p.mode,'R'); el = 1 - p.T(:); else; el = p.T(:); end
                    mf = NaN; if isfield(p,'floor'); mf = p.floor; end
                    t = t .* max(el, FilterEngine.elementFloor(mf, cap));
                end
            end
            function v = getfloor(s,f); if isfield(s,f); v = s.(f); else; v = NaN; end; end
        end

        function fe = elementFloor(measured, cap)
            % measured = curve's demonstrated blocking floor (NaN if unknown);
            % cap = global optimism limit. Conservative: never deeper than data.
            if isnan(measured); fe = cap; else; fe = max(measured, cap); end
        end

        function CT = crosstalkMatrix(S, assign)
            % assign(k) = index of the fluorophore that "owns" channel k.
            % CT(i,k) = S(i,k) / S(owner_k, k) ; diagonal-of-interest = 1.
            nc = size(S,2);
            CT = zeros(size(S));
            for k = 1:nc
                denom = S(assign(k), k) + eps;
                CT(:,k) = S(:,k) / denom;
            end
        end

        function score = systemScore(S, assign, ctWeight, bleed, bleedWeight)
            % Figure of merit: wanted signal - weighted crosstalk - weighted
            % laser back-reflection background. Higher is better.
            %   bleed (optional)      : nc-vector or nc×nLaser matrix from
            %                           laserBleed (same units as S).
            %   bleedWeight (optional): penalty on bleed-to-signal ratio.
            if nargin < 3 || isempty(ctWeight); ctWeight = 5; end
            if nargin < 4; bleed = []; end
            if nargin < 5 || isempty(bleedWeight); bleedWeight = 0; end
            nc = size(S,2);
            wanted = 0; leak = 0; blk = 0;
            if ~isempty(bleed); bleedTot = sum(bleed,2); else; bleedTot = zeros(nc,1); end
            for k = 1:nc
                own = assign(k);
                colMax = max(S(:,k)) + eps;
                wanted = wanted + S(own,k)/colMax;          % normalised capture
                others = setdiff(1:size(S,1), own);
                if ~isempty(others)
                    leak = leak + sum(S(others,k)) / (S(own,k)+eps);
                end
                blk = blk + bleedTot(k) / (S(own,k)+eps);   % laser bkgnd / signal
            end
            score = wanted - ctWeight*leak - bleedWeight*blk;
        end

        function R = optimize(fluors, lasers, baseChannels, candidates, assign, lambda, ctWeight, detector, Rback, bleedWeight, blockOD, scoreFcn)
            %OPTIMIZE Brute-force search over candidate emission filters per channel.
            %  candidates : 1xnc cell, candidates{k} = struct array of filter
            %               spectra (each with .name and .ex) to try in channel k.
            %  detector   : optional QE vector (see signalMatrix).
            %  Rback,bleedWeight : laser back-reflection fraction and its score
            %               penalty (see laserBleed / systemScore).
            %  Returns R: struct array sorted by score (best first) with fields
            %    .filters (1xnc names) .score .S .CT .eff .bleed
            if nargin < 7; ctWeight = 5; end
            if nargin < 8; detector = []; end
            if nargin < 9 || isempty(Rback); Rback = 0.005; end
            if nargin < 10 || isempty(bleedWeight); bleedWeight = 0; end
            if nargin < 11 || isempty(blockOD); blockOD = 6; end
            if nargin < 12; scoreFcn = []; end   % scoreFcn(channels,lasers)->score
            nc = numel(baseChannels);
            grids = cell(1,nc);
            for k = 1:nc; grids{k} = 1:numel(candidates{k}); end
            combos = FilterEngine.cartprod(grids);
            R = struct('filters',{},'score',{},'S',{},'CT',{},'eff',{},'bleed',{});
            for c = 1:size(combos,1)
                ch = baseChannels;
                names = cell(1,nc);
                for k = 1:nc
                    f = candidates{k}(combos(c,k));
                    ch(k).emFilter = f.ex(:);
                    if isfield(f,'floor'); ch(k).emFloor = f.floor; end
                    names{k} = f.name;
                end
                [S, eff] = FilterEngine.signalMatrix(fluors, lasers, ch, lambda, detector);
                bleed = FilterEngine.laserBleed(ch, lasers, lambda, Rback, detector, blockOD);
                if isempty(scoreFcn)
                    sc = FilterEngine.systemScore(S, assign, ctWeight, bleed, bleedWeight);
                else
                    sc = scoreFcn(ch, lasers);
                end
                R(end+1) = struct('filters',{names},'score',sc,'S',S, ...
                    'CT',FilterEngine.crosstalkMatrix(S,assign),'eff',eff,'bleed',bleed); %#ok<AGROW>
            end
            [~, ord] = sort([R.score], 'descend');
            R = R(ord);
        end

        function f = combinerFactors(Ts, nLam)
            %COMBINERFACTORS Fraction of each excitation source reaching the
            %   output through a chain of combiner dichroics. Ts is a 1x(N-1)
            %   cell of transmission vectors ([] = no element / no loss). Source
            %   1 is the base (transmitted through all combiners); each later
            %   source enters by reflection at its combiner. Returns 1xN cell.
            N = numel(Ts) + 1; [Tt, Tr] = FilterEngine.stages(Ts, nLam);
            f = cell(1,N);
            f{1} = FilterEngine.prodFrom(Tt, 1, N-1, nLam);
            for i = 2:N
                f{i} = Tr{i-1} .* FilterEngine.prodFrom(Tt, i, N-1, nLam);
            end
        end

        function f = splitterFactors(Ts, nLam)
            %SPLITTERFACTORS Fraction of emission directed to each channel through
            %   a chain of splitter dichroics. Channel k is reflected at splitter
            %   k (transmitted through 1..k-1); the last channel is all-transmit.
            N = numel(Ts) + 1; [Tt, Tr] = FilterEngine.stages(Ts, nLam);
            f = cell(1,N);
            for k = 1:N-1
                f{k} = FilterEngine.prodFrom(Tt, 1, k-1, nLam) .* Tr{k};
            end
            f{N} = FilterEngine.prodFrom(Tt, 1, N-1, nLam);
        end

        function [Tt, Tr] = stages(Ts, nLam)
            m = numel(Ts); Tt = cell(1,m); Tr = cell(1,m);
            for i = 1:m
                if isempty(Ts{i}); Tt{i} = ones(nLam,1); Tr{i} = ones(nLam,1);
                else; Tt{i} = Ts{i}(:); Tr{i} = 1 - Ts{i}(:); end
            end
        end

        function p = prodFrom(Tt, a, b, nLam)
            p = ones(nLam,1);
            for i = a:b; p = p .* Tt{i}; end
        end

        function channels = makeChannels(names, primaryT, splitterT, modes, emFilters, ...
                                         emFloors, primaryFloor, splitterFloor)
            %MAKECHANNELS Assemble channel structs from shared optics.
            %  names{k}    : channel name
            %  primaryT    : primary dichroic transmission (always Transmit), or []
            %  splitterT   : shared splitter dichroic transmission, or []
            %  modes{k}    : 'T' | 'R' | 'none'  (how channel k uses the splitter)
            %  emFilters{k}: emission filter transmission vector
            %  emFloors{k}, primaryFloor, splitterFloor (optional): measured
            %    blocking floors used by laserBleed; NaN/omitted => unknown.
            nc = numel(names);
            if nargin < 6 || isempty(emFloors); emFloors = num2cell(nan(1,nc)); end
            if nargin < 7; primaryFloor = NaN; end
            if nargin < 8; splitterFloor = NaN; end
            channels = struct('name',cell(1,nc),'path',[],'emFilter',[],'emFloor',[]);
            for k = 1:nc
                path = struct('T',{},'mode',{},'floor',{});
                if ~isempty(primaryT)
                    path(end+1) = struct('T',primaryT(:),'mode','T','floor',primaryFloor); %#ok<AGROW>
                end
                if ~isempty(splitterT) && ~strcmpi(modes{k},'none')
                    path(end+1) = struct('T',splitterT(:),'mode',upper(modes{k}),'floor',splitterFloor); %#ok<AGROW>
                end
                channels(k).name = names{k};
                channels(k).path = path;
                channels(k).emFilter = emFilters{k}(:);
                channels(k).emFloor = emFloors{k};
            end
        end

        function R = optimizeJoint(fluors, lasers, names, primaryT, splitterCands, ...
                                   modes, filterCands, assign, lambda, ctWeight, ...
                                   detector, Rback, bleedWeight, exFilterCands, blockOD, primaryFloor, scoreFcn)
            %OPTIMIZEJOINT Co-optimise the shared splitter dichroic, every
            %   channel's emission filter, and (optionally) each source's
            %   excitation/cleanup filter, simultaneously.
            %  splitterCands : struct array (.name .ex) of candidate splitters.
            %  filterCands   : 1xnc cell, each a struct array (.name .ex).
            %  exFilterCands : 1xnLaser cell, each a struct array (.name .ex);
            %                  empty entry => keep that source's current filter.
            %  Returns R sorted by score with fields:
            %    .dichroic .filters (1xnc) .exfilters (1xnLaser) .score .S .CT .eff .bleed
            if nargin < 10 || isempty(ctWeight); ctWeight = 5; end
            if nargin < 11; detector = []; end
            if nargin < 12 || isempty(Rback); Rback = 0.005; end
            if nargin < 13 || isempty(bleedWeight); bleedWeight = 0; end
            if nargin < 14; exFilterCands = {}; end
            if nargin < 15 || isempty(blockOD); blockOD = 6; end
            if nargin < 16; primaryFloor = NaN; end
            if nargin < 17; scoreFcn = []; end   % scoreFcn(channels,lasers)->score
            nc = numel(names); nl = numel(lasers);
            % axes: splitter, per-channel filter, per-laser ex filter
            axes = [{1:numel(splitterCands)}, cell(1,nc), cell(1,nl)];
            for k = 1:nc; axes{1+k} = 1:numel(filterCands{k}); end
            for j = 1:nl
                if numel(exFilterCands)>=j && ~isempty(exFilterCands{j})
                    axes{1+nc+j} = 1:numel(exFilterCands{j});
                else
                    axes{1+nc+j} = 1;   % keep current
                end
            end
            combos = FilterEngine.cartprod(axes);
            R = struct('dichroic',{},'filters',{},'exfilters',{},'score',{}, ...
                       'S',{},'CT',{},'eff',{},'bleed',{});
            for c = 1:size(combos,1)
                sp = splitterCands(combos(c,1));
                spFloor = NaN; if isfield(sp,'floor'); spFloor = sp.floor; end
                emF = cell(1,nc); emFl = cell(1,nc); fnames = cell(1,nc);
                for k = 1:nc
                    f = filterCands{k}(combos(c,1+k));
                    emF{k} = f.ex(:); fnames{k} = f.name;
                    if isfield(f,'floor'); emFl{k} = f.floor; else; emFl{k} = NaN; end
                end
                las = lasers; exnames = cell(1,nl);
                for j = 1:nl
                    if numel(exFilterCands)>=j && ~isempty(exFilterCands{j})
                        ef = exFilterCands{j}(combos(c,1+nc+j));
                        las(j).exFilter = ef.ex(:); exnames{j} = ef.name;
                    else
                        exnames{j} = '(unchanged)';
                    end
                end
                ch = FilterEngine.makeChannels(names, primaryT, sp.ex, modes, emF, ...
                    emFl, primaryFloor, spFloor);
                [S, eff] = FilterEngine.signalMatrix(fluors, las, ch, lambda, detector);
                bleed = FilterEngine.laserBleed(ch, las, lambda, Rback, detector, blockOD);
                if isempty(scoreFcn)
                    sc = FilterEngine.systemScore(S, assign, ctWeight, bleed, bleedWeight);
                else
                    sc = scoreFcn(ch, las);
                end
                R(end+1) = struct('dichroic',sp.name,'filters',{fnames}, ...
                    'exfilters',{exnames},'score',sc,'S',S, ...
                    'CT',FilterEngine.crosstalkMatrix(S,assign),'eff',eff,'bleed',bleed); %#ok<AGROW>
            end
            [~, ord] = sort([R.score], 'descend');
            R = R(ord);
        end

        function out = cartprod(grids)
            % Cartesian product of index vectors -> rows of combinations.
            n = numel(grids);
            [g{1:n}] = ndgrid(grids{:});
            out = zeros(numel(g{1}), n);
            for k = 1:n; out(:,k) = g{k}(:); end
        end
    end
end
