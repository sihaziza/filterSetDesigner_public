classdef FPbase
%FPBASE  Client for the FPbase GraphQL API (filters, dichroics, light
%   sources, and fluorophore spectra). https://www.fpbase.org/graphql/
%
%   FPbase hosts the full Semrock / Chroma / Omega / Zeiss catalogues, so
%   you never need to download per-filter text files or scrape vendor sites.
%
%   Catalogue (FPbase categories):
%     F = filter      (~4100; subtype BP/LP/SP bandpass-etc, and BS = dichroic
%                      beamsplitter -- dichroics live HERE, not under D)
%     L = light source(~270)
%     C = camera      (~70)
%     P = protein  ,  D = organic dye  (both fluorophores, EX/EM subtypes)
%
%   Usage
%     list = FPbase.catalogue();                 % cached struct array
%     hits = FPbase.search('FF01-520');          % name substring (any filter)
%     hits = FPbase.search('Di01', 'F', 'BS');   % dichroics only (subtype BS)
%     S    = FPbase.spectrum(hits(1).id, lambda);% download one curve
%
%   FPbase.search/spectrum return structs compatible with loadSpectrum:
%     S.name, S.lambda, S.ex (the transmission/curve), S.em ([] for filters),
%     S.kind, S.id, S.category, S.owner.

    properties (Constant)
        ENDPOINT = 'https://www.fpbase.org/graphql/'
    end

    methods (Static)
        function C = catalogue(forceRefresh)
            %CATALOGUE  Cached list of every spectrum on FPbase (id/category/
            %   subtype/owner name). Cached in memory + a local .mat file.
            persistent CACHE
            if nargin < 1; forceRefresh = false; end
            matFile = fullfile(fileparts(mfilename('fullpath')), 'fpbase_catalogue.mat');
            if ~forceRefresh
                if ~isempty(CACHE); C = CACHE; return; end
                if isfile(matFile)
                    L = load(matFile,'C'); CACHE = L.C; C = CACHE; return;
                end
            end
            wo = weboptions('MediaType','application/json','Timeout',90);
            q  = struct('query','{ spectra { id category subtype owner { name } } }');
            R  = webwrite(FPbase.ENDPOINT, q, wo);
            raw = R.data.spectra;
            n = numel(raw);
            C = struct('id',cell(1,n),'category',[],'subtype',[],'name',[]);
            for i = 1:n
                x = raw(i); if iscell(raw); x = raw{i}; end
                C(i).id = str2double(x.id);                 % API returns char
                C(i).category = char(x.category);
                C(i).subtype  = char(string(x.subtype));    % BP/LP/SP/BS/EX/EM
                nm = '';
                if isstruct(x.owner) && isfield(x.owner,'name') && ~isempty(x.owner.name)
                    nm = char(x.owner.name);                % the product name
                end
                C(i).name  = nm;
            end
            CACHE = C;
            try; save(matFile,'C'); catch; end
        end

        function hits = search(query, category, subtype, maxHits)
            %SEARCH  Case-insensitive substring match on the product name.
            %   category: '' (any) or 'F','L','C','P','D'
            %   subtype : '' (any) or 'BP','LP','SP','BS','EX','EM' ...
            if nargin < 2; category = ''; end
            if nargin < 3; subtype  = ''; end
            if nargin < 4; maxHits  = 100; end
            C = FPbase.catalogue();
            if ~isempty(category); C = C(strcmpi({C.category}, category)); end
            if ~isempty(subtype);  C = C(strcmpi({C.subtype},  subtype));  end
            keep = contains(lower({C.name}), lower(strtrim(query)));
            hits = C(keep);
            % de-duplicate identical names (keep first), sort by name
            [~, ia] = unique({hits.name}, 'stable'); hits = hits(ia);
            [~, ord] = sort(lower({hits.name})); hits = hits(ord);
            if numel(hits) > maxHits; hits = hits(1:maxHits); end
        end

        function S = spectrum(id, lambda)
            %SPECTRUM  Download one spectrum's data and resample onto lambda.
            if nargin < 2 || isempty(lambda); lambda = (350:1:850)'; end
            lambda = lambda(:);
            wo = weboptions('MediaType','application/json','Timeout',60);
            q  = struct('query', sprintf( ...
                '{ spectrum(id: %d){ subtype category data owner { name } } }', id));
            R  = webwrite(FPbase.ENDPOINT, q, wo);
            sp = R.data.spectrum;
            d  = sp.data; if iscell(d); d = cell2mat(d); end
            w = d(:,1); v = d(:,2);
            vi = interp1(w, v, lambda, 'linear', 0);
            vi(~isfinite(vi)) = 0;

            cat = char(sp.category);
            sub = char(string(sp.subtype));
            owner = ''; if isstruct(sp.owner)&&isfield(sp.owner,'name'); owner = char(sp.owner.name); end
            S.id = id; S.category = cat; S.subtype = sub;
            S.name = owner; S.lambda = lambda;
            if any(strcmpi(sub,{'ex','ab','abs','em','2p'}))
                % fluorophore single curve (protein P or dye D)
                S.kind = 'fluorophore';
                if strcmpi(sub,'em'); S.em = normPeak(vi); S.ex = zeros(size(vi));
                else;                 S.ex = normPeak(vi); S.em = zeros(size(vi)); end
                S.floor = NaN;
            else
                S.kind = 'filter';            % BP/LP/SP/BS (BS = dichroic)
                if max(vi) > 1.5; vi = vi/100; end   % percent -> fraction
                S.ex = min(max(vi,0),1); S.em = [];
                pos = S.ex(S.ex>0); mp = min([pos; inf]);
                if mp < 1e-4; S.floor = max(mp,1e-7); else; S.floor = NaN; end
            end
            function y = normPeak(y); m=max(y); if m>0; y=y/m; end; end
        end
    end
end
