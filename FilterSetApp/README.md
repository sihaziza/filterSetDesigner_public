# Optimal Filter Set Selector

Two MATLAB GUIs for fluorescence imaging system design:

1. **`OptimalFilterApp`** — choose optimal optical filter sets (signal,
   crosstalk, laser back-reflection, optimizer) with a web connection to
   [FPbase](https://www.fpbase.org).
2. **`FilterSetSNRApp`** — take the chosen filter set and compute absolute
   signal-to-noise following the IDEX/Semrock spectral model, extended with
   photon shot noise, detector read noise, dark current, and tissue / fibre
   autofluorescence.

## Run

```matlab
cd 'Optimal Filter sets'
runApp                  % 1) design the filter set
% ... in Results, click "Save config → SNR app"
addpath('FilterSetApp')
FilterSetSNRApp         % 2) load that config and compute SNR
```

The app auto-loads spectra from the project `Spectra` folder:
`Spectra/Proteins`, `Spectra/Filters`, `Spectra/Dichroics`,
`Spectra/Illumations`, and `Spectra/Detectors`. Filter-set manifests live in
`FilterSetApp/FilterSets`. To point at a different project root:

```matlab
OptimalFilterApp('C:\path\to\project_root')
```

## Tabs

### Tab layout
- **Build Fluorescence Imaging Filter Set** — library tree, the main spectral plot, and the **primary dichroic / detector / back-reflection** controls + **Load/Save filter set**.
- **Spectra Library + Web** — FPbase web search/download + the duplicate finder.
- **System Builder** — fluorophores, excitation sources, detection channels, Compute.
- **Results**, **Optimizer** — as below.

### Spectral viewing & cleanup
- **Superposition viewer** (Build tab "Main spectral plot"): select tree entries and **Plot selected** / **Add selected** / **Remove selected** to build an overlay (fluorophores as ex-dashed / em-solid), **Clear** to reset — SearchLight-style. **Y axis** toggles **%T ↔ OD**.
- **Duplicate finder** (Spectra Library + Web tab, "Find / remove duplicate spectra (xcorr)"): computes a **zero-lag normalised cross-correlation matrix** across same-kind spectra, shows it as a heatmap with the duplicate groups, and removes extras (keeps one per group, preferring non-`BFMConfig`). Also fixes a double-load bug where case-insensitive `*.txt`/`*.TXT` globbing loaded every file twice.
- **Results spectral plot**: **Plot** dropdown switches **Component superposition / Channel overlay / Raw vs filtered**; **Y** toggles **%T ↔ OD**; **X min/max** fields + a **scroll** slider pan/zoom the wavelength axis.

1. **Library + Web** — browse local spectra; type a fluorophore name and
   **Fetch** to pull its excitation/emission spectra **and brightness
   (EC×QY/1000)** straight from FPbase.
   *Filters & dichroics too:* pick **Bandpass / Long-Short-pass / Dichroic (BS)
   / Light source** in the Type menu and search by part name (e.g. `FF01-520`,
   `Di01-R488`, `ET525`). FPbase mirrors the full **Semrock, Chroma, Omega and
   Zeiss** catalogues (~4100 filters incl. dichroics, ~270 light sources), so
   you never download a vendor text file again. The catalogue is cached to
   `fpbase_catalogue.mat` for offline reuse.
2. **System Builder** — define fluorophores (+brightness); excitation sources
   (centre λ + power, a **laser line or a lamp spectrum**, plus an optional
   **excitation/cleanup filter**); a primary dichroic; a **detector QE**
   (camera/APD, or "(ideal)"); the **back-reflection R (%)** and per-element
   **blocking OD**; and one row per detection channel (splitter dichroic +
   Transmit/Reflect + emission filter + owning fluorophore). Pre-seeded with the
   uSMAART example.
3. **Results** — signal matrix `S(fluor, channel)`, crosstalk matrix
   (%, column-normalised to each channel's owner), the **laser back-reflection
   background/signal** table (per channel × source), and a spectral overlay
   (channel throughputs, fluorophore emission, laser lines). **Export results →
   CSV** saves every matrix to a spreadsheet-friendly file.
4. **Optimizer** — two modes:
   * **Single channel** — sweep one channel's emission filter, others fixed.
   * **Joint (splitter + all filters)** — co-optimise the shared splitter
     dichroic *and* every channel's emission filter simultaneously, ranking the
     full combination by the figure of merit. Set candidate filters per channel
     via the Channel menu (each channel remembers its own list); pick candidate
     splitter dichroics in the lower list (the current one is always included).

   The **+ Pull bandpass filters from FPbase near owner emission** button
   auto-downloads only the real Semrock/Chroma/etc. bandpass filters whose
   centre is near the owning fluorophore's emission peak (parsed from the part
   name) and loads them as candidates for that channel. **Crosstalk weight** and
   **bleed weight** tune the score; tick **Opt. excitation filters** to also
   co-select each source's cleanup filter (auto-pulled near the source line).
   Detector QE and laser back-reflection from the System Builder apply
   throughout. **Export ranking → CSV** saves the ranked table.

## On "automatic" optimisation & web-scraping (realistic scope)

Direct HTML scraping of Semrock/Chroma is unnecessary and fragile (page
changes, rate limits, terms of use). FPbase already exposes their catalogues as
structured data through one GraphQL endpoint — that is what `FPbase.m` uses.

Searching *every* filter (~4100) across multiple channels is combinatorially
too large to brute-force blindly. The realistic, implemented strategy is
**prune-then-search**: parse each filter's centre wavelength from its name,
keep only those within a window of the target fluorophore's emission peak
(typically 10–30 candidates), download just those, then exhaustively score the
small set. This finds a near-optimal emission filter per channel in seconds and
is fully automatic from a light-source + fluorophore configuration.

## SNR calculator (`FilterSetSNRApp`)

Loads a configuration saved from the filter app and computes, **per channel, in
detected photo-electrons** over the integration time:

| Quantity | Meaning |
|---|---|
| Signal | owner-fluorophore electrons (Semrock S, via ε·c·d absorption, QY, NA collection, channel throughput × QE) |
| Crosstalk | electrons from the other fluorophores in that channel |
| Excitation bleed | back-reflected laser electrons (Semrock N_E; uses R and per-element OD) |
| Autofluorescence | tissue / fibre / flavin / lipofuscin background (Semrock N_AF) |
| Dark | dark-current electrons (rate × time) |
| Shot noise | √(signal + all background) — photon statistics |
| Read noise | detector read noise (e⁻ RMS, added in quadrature) |
| Total noise | √(signal + background + read²) |
| **SNR** | signal / total noise (electron-domain) |
| SNR (optical) | Semrock S/(N_E+N_F), scale-invariant |

**Autofluorescence** sources are excitation-wavelength-dependent (stronger at
shorter λ) with broad, red-tailed emission — presets `Brain tissue`,
`Silica fiber`, `Flavin (FAD)`, `Lipofuscin` (see `autofluorPreset.m`); tune
strength (peak absorbed fraction ≈ ε·c·d) and QY per source.

The **Sweep** tab plots SNR vs one parameter (integration time, laser power,
concentration, read noise, NA, or autofluorescence strength) for every channel,
and marks the shot-noise = read-noise crossover for the worst channel — useful
for finding where you become read-noise- vs background-limited.

The filter-set **Optimizer** (in `OptimalFilterApp`) can rank candidate sets by
either the optical figure of merit or the **electron-domain SNR** (Score-by
menu) — the latter calls this same `snrModel`, so the optimizer optimises the
real SNR including shot/read/dark/autofluorescence at a chosen operating point.

**Collection model.** Signal is computed for a **lumped / point detector** — the
whole excited detection volume funnels onto one detector (fiber photometry, APD,
single-channel PMT). For camera-pixel imaging, set the concentration to the
per-pixel molecule equivalent. SNR_optical is independent of absolute scale;
the electron-domain SNR (and the read-vs-shot-noise balance) depends on power,
concentration, path length and integration time, so dim/low-light regimes
correctly show read noise mattering while bright regimes are shot-limited.

## Scripting / batch use

```matlab
app = OptimalFilterApp;
r   = app.compute;      % r.S, r.CT, r.eff, r.score, r.fluors, r.channels
```

## Files

| File | Role |
|------|------|
| `OptimalFilterApp.m` | filter-set design GUI (programmatic App Designer class) |
| `FilterSetSNRApp.m` | SNR calculator GUI (loads a saved config) |
| `snrModel.m` | absolute electron-domain signal/noise/SNR model |
| `autofluorPreset.m` | tissue/fibre/flavin/lipofuscin autofluorescence presets |
| `loadSpectrum.m` | robust importer (auto delimiter/header/scale/step) |
| `fetchFPbase.m` | FPbase REST client → fluorophore spectrum + brightness |
| `FPbase.m` | FPbase GraphQL client → filter/dichroic/light-source catalogue + download (cached) |
| `FilterEngine.m` | optics math: signal, efficiency, crosstalk, optimizer |
| `selftest.m` | headless regression test against the bundled data |

## Physics summary

For channel *k* with throughput `T_k(λ) = emFilter · Π(dichroic T or 1−T)`,
detector QE, and source *j* effective excitation `E_j(λ)` (a laser line or lamp
spectrum after its cleanup filter, with `∫E_j = power_j`):

```
collected(i,k) = ∫ em_i(λ)·T_k(λ)·QE(λ) dλ
eff(i,k)       = collected / ∫ em_i(λ) dλ
S(i,k)         = brightness_i · Σ_j[∫ E_j(λ)·ex_i(λ) dλ] · collected(i,k)
CT(i,k)        = S(i,k) / S(owner_k, k)
```

**Laser back-reflection.** A fraction `R` (~0.5%) of each source back-scatters
off non-AR-coated optics into the detection path, attenuated by the same optics:

```
bleed(k,j) = R · ∫ E_j(λ)·T_k(λ)·QE(λ) dλ
bleed/signal(k) = Σ_j bleed(k,j) / S(owner_k, k)
```

The incident flux cancels in `bleed/signal`, so it is a true background-to-signal
ratio set entirely by the filter set — which is what rewards deep emission-filter
/ dichroic OD at the laser line.

**Per-element blocking floor.** Each filter/dichroic uses its *own* measured
blocking depth: `loadSpectrum`/`FPbase` record the deepest demonstrated
transmission of each curve (e.g. real Semrock files reach ~1e-7, i.e. OD7).
Curves that show no genuine deep blocking (linear-scale data that rounds to 0)
are marked "unknown" and fall back to the global **Element blocking (OD)** field.
The effective floor is `max(measured, 10^−OD)`, so the global field is a
conservative optimism cap — the estimate is never credited deeper than the data
supports, yet filters with real high-dynamic-range data use their true OD.

The **joint optimiser** evaluates
`score = Σ_k capture_k − w_ct·Σ_k leak_k − w_bleed·Σ_k bleed/signal_k`
over the Cartesian product of {candidate splitter dichroics} × {candidate
emission filter per channel} × {candidate excitation filter per source},
returning the ranked list (capped at 20 000 combinations).

`selftest.m` validates this against your GFP / cyOFP / mRuby3 uSMAART set.
