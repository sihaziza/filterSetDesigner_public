function out = snrModel(cfg, phys, fp, af)
%SNRMODEL  Absolute signal-to-noise model for a fluorescence filter set.
%
%   Implements the IDEX/Semrock "Spectral Modeling in Fluorescence Microscopy"
%   optical model (signal S, excitation-reflection noise N_E, autofluorescence
%   noise N_AF, optical SNR = S/(N_E+N_F)) and extends it into the detector
%   electron domain with PHOTON SHOT NOISE, DETECTOR READ NOISE and DARK
%   current, giving a realistic per-channel SNR.
%
%   Everything is computed in detected photo-electrons accumulated over the
%   integration time, so read noise (a fixed electron count) and shot noise
%   (sqrt of the electron count) trade off correctly.
%
%   INPUTS
%     cfg  : system from the filter-set app (assembleSystem / saveConfig):
%            .lambda .fluors(.name .ex .em) .lasers(.name .wl .power .spectrum
%            .exFilter) .channels .assign .detector .Rback .blockOD
%     phys : .NA .n .pathLength_cm .tInt_s .readNoise_e .darkRate_eps
%            .powers_mW (1 x nSource, overrides cfg laser power)
%     fp   : per-fluorophore photophysics, 1 x nFluor struct:
%            .ec (M^-1 cm^-1, peak extinction) .qy .conc_M
%     af   : autofluorescence sources, struct array:
%            .name .absorb(λ, peak 1) .em(λ, peak 1) .strength (peak absorbed
%            fraction, ~ ε·c·d) .qy
%
%   OUTPUT struct OUT (per channel k = column):
%     .channels .signal_e .xtalk_e .exBleed_e .af_e(total) .afBySrc
%     .dark_e .shot_e .read_e .noise_e .SNR (electron) .SNR_optical
%     .afNames

    lambda = cfg.lambda(:); dl = mean(diff(lambda));
    hc = 1.98644586e-25;                 % J*m
    photonE = hc ./ (lambda*1e-9);       % J per photon at each λ
    QE = ones(numel(lambda),1);
    if isfield(cfg,'detector') && ~isempty(cfg.detector); QE = cfg.detector(:); end

    nF = numel(cfg.fluors); nC = numel(cfg.channels); nS = numel(cfg.lasers);

    % ---- incident excitation photon-rate spectrum (photons/s/nm) ----
    PhiTot = zeros(numel(lambda),1);
    Phi = cell(1,nS);
    for j = 1:nS
        src = cfg.lasers(j);
        if isfield(phys,'powers_mW') && numel(phys.powers_mW)>=j
            src.power = phys.powers_mW(j)*1e-3;          % mW -> W
        end
        E = FilterEngine.sourceExcitation(src, lambda);  % W/nm, integral=power
        Phi{j} = E ./ photonE;                           % photons/s/nm
        PhiTot = PhiTot + Phi{j};
    end

    % ---- collection efficiency from NA ----
    th = asin(min(phys.NA/phys.n, 1));
    colEff = (1 - cos(th))/2;            % fraction of 4π collected
    t = phys.tInt_s;

    % ---- per-channel detection transmission (no detector) ----
    Tk = zeros(numel(lambda), nC); TkFloor = zeros(numel(lambda), nC);
    for k = 1:nC
        Tk(:,k)      = FilterEngine.pathTransmission(cfg.channels(k), lambda);
        TkFloor(:,k) = FilterEngine.pathTransmissionFloored(cfg.channels(k), lambda, cfg.blockOD);
    end

    % ---- fluorophore signal electrons N_S(i,k) ----
    N_S = zeros(nF, nC);
    for i = 1:nF
        ec = fp(i).ec; qy = fp(i).qy; c = fp(i).conc_M; d = phys.pathLength_cm;
        absFrac = 1 - 10.^(-(ec * cfg.fluors(i).ex(:)) * c * d);    % per λ
        Rabs = trapz(lambda, PhiTot .* absFrac);                    % photons/s absorbed
        emN = cfg.fluors(i).em(:); a = trapz(lambda, emN); if a>0; emN = emN/a; end
        for k = 1:nC
            detFrac = trapz(lambda, emN .* Tk(:,k) .* QE);          % collected & detected /emitted
            N_S(i,k) = Rabs * qy * colEff * detFrac * t;            % electrons
        end
    end

    % ---- excitation back-reflection electrons N_E(k) ----
    N_E = zeros(1,nC);
    for k = 1:nC
        N_E(k) = cfg.Rback * trapz(lambda, PhiTot .* TkFloor(:,k) .* QE) * t;
    end

    % ---- autofluorescence electrons N_AF(m,k) ----
    nAF = numel(af);
    AF = zeros(nAF, nC);
    AFsrc = zeros(nAF, nS, nC);
    for m = 1:nAF
        absFrac = af(m).strength * af(m).absorb(:);                 % absorbed fraction per λ
        emN = af(m).em(:); a = trapz(lambda, emN); if a>0; emN = emN/a; end
        for j = 1:nS
            Rabs = trapz(lambda, Phi{j} .* absFrac);
            for k = 1:nC
                detFrac = trapz(lambda, emN .* Tk(:,k) .* QE);
                AFsrc(m,j,k) = Rabs * af(m).qy * colEff * detFrac * t;
            end
        end
    end
    if nAF > 0
        AF = squeeze(sum(AFsrc,2));
        if nAF == 1; AF = reshape(AF,1,nC); end
    end

    % ---- assemble per-channel SNR ----
    assign = cfg.assign;
    N_dark = phys.darkRate_eps * t;
    sig = zeros(1,nC); xtk = zeros(1,nC); afTot = sum(AF,1);
    shot = zeros(1,nC); noise = zeros(1,nC); snr = zeros(1,nC); snrOpt = zeros(1,nC);
    for k = 1:nC
        o = assign(k);
        sig(k) = N_S(o,k);
        xtk(k) = sum(N_S(:,k)) - N_S(o,k);
        bg = xtk(k) + N_E(k) + afTot(k) + N_dark;
        shot(k) = sqrt(sig(k) + bg);
        noise(k) = sqrt(sig(k) + bg + phys.readNoise_e^2);
        snr(k) = sig(k) / max(noise(k), eps);
        opt = N_E(k) + afTot(k) + xtk(k);
        snrOpt(k) = sig(k) / max(opt, eps);
    end

    out.channels   = {cfg.channels.name};
    out.afNames    = {af.name};
    out.signal_e   = sig;
    out.xtalk_e    = xtk;
    out.exBleed_e  = N_E;
    out.af_e       = afTot;
    out.afBySrc    = AF;
    out.afByExcitation = AFsrc;
    out.excitationNames = {cfg.lasers.name};
    out.dark_e     = N_dark*ones(1,nC);
    out.shot_e     = shot;
    out.read_e     = phys.readNoise_e*ones(1,nC);
    out.noise_e    = noise;
    out.SNR        = snr;
    out.SNR_optical= snrOpt;
    out.colEff     = colEff;
end
