# CLAUDE.md — EPOS Au+Au 200 GeV Flow Analysis

Project documentation for Claude Code sessions. Everything below was established
empirically during the 2026-07 elliptic-flow analysis of `auau200-7`.

## Project layout

- `/media/Students/aliu/` — project root (moved here 2026-07-02 from `/media/Disk_Jin/aliu`, old path is dead)
- `auau200-N/` — datasets; `auau200-7/` is the active one: 90 files `z-auau_run_*.root`,
  1000 events each (~84,000 usable events, ~240 MB/file, ~20 GB total)
- `epos1..epos5/` (under project root) — 5 concurrent EPOS 4.0.3 generator instances
- `~/epos_scheduler.sh` → `~/cycle_epos.sh` — production automation: detached daemon loops every 7 h,
  skips instances whose `epos` process is still alive (pgrep guard), moves finished
  `z-auau_run_N.root` into the dataset with the next free index, restarts each instance
  in an isolated subshell. `RESTART=0` = harvest only.
- `scripts/` — the production scripts above, versioned (2026-08-13). They physically live
  here; `~/epos_scheduler.sh`, `~/cycle_epos.sh`, `~/run_all.sh`, `~/run.sh`, `~/startup.sh`
  and `~/epos{1..5}/startup.sh` are **symlinks into this directory**. Edit either path — same file.

## Version control

Git repo at the project root since 2026-08-13. Tracked: `*.C`, `*.h`, `*.sh`, `*.md`,
`idt.dt`/`read.txt`/`good_files.txt`/`bad_files.txt`, HEPData CSVs, and small outputs
(`cen*.root`, `*.png`) so `MakeFigure_*.C` runs from a clean checkout without a re-loop.

`.gitignore` is a **whitelist** — it ignores `*` and re-includes only the above. This is
deliberate: 150 GB of `z-*.root` sits in the same directories as the macros, so a blacklist
plus one `git add .` would commit a 93 MB event file permanently. When adding a new file
type worth tracking, add a `!pattern` line rather than removing the `*`.

Remotes: `backup` → `/home/aliu/backup/auau-analysis.git` (bare, on the other physical
disk `nvme0n1p3`; `git push backup main`).

**After a reboot**: `/media/Students` must be mounted *before* the scheduler daemon is
restarted, or every `~/…startup.sh` symlink dangles and the epos runs fail at `source`.
Check with `readlink -e ~/cycle_epos.sh` (empty output = not mounted yet). Restart with:
`cd ~ && nohup ./epos_scheduler.sh auau200-7 > scheduler.log 2>&1 &`

Known bug, unfixed: `scripts/run_all.sh` still points `JINBASE` at the dead
`/media/Disk_Jin/aliu` path. `cycle_epos.sh` (what the daemon actually calls) is correct.

## Dataset: TTree `teposevent` — CRITICAL branch semantics

These were verified empirically; several branch names are misleading:

| branch | actual meaning |
|---|---|
| `phi`  | **true reaction-plane angle Ψ_RP** — random per-event rotation, uniform in [0, 2π). USE THIS as the reaction plane. |
| `phir` | small fluctuation of the participant plane relative to the RP (peaked at 0, RMS ≈ 0.22). Participant plane = `phi + phir`. **NOT the reaction plane** — using it alone gives v₂ ≈ 0. |
| `psi2`–`psi5`, `ecci2`–`ecci5` | unfilled (~1e-4 in every event). Do not use. |
| `e`    | particle **mass**, not energy (π → 0.1396). Energy = √(m²+p²); rapidity y = ½ ln((E+pz)/(E−pz)). |
| `bim`  | impact parameter; this dataset was generated in 8.4–9.2 fm = 30–40% centrality (all events in range). |
| `ist`  | 0 = final-state particle; 1 = decayed (π⁰ are ALWAYS ist==1 — they never appear as final state). |
| `id`   | **EPOS ids, not PDG**: π± = ±120, K± = ±130, p/p̄ = ±1120, π⁰ = 110. Full table in `auau200-7/idt.dt`. |
| `np`   | tracks per event (~9,000–10,000). Charged-pion POI (0.2<pT<2, |y|<1) ≈ 110/event, max ≈ 180. |

How the `phi`/`phir` mislabel was proven: reconstruct each event's plane from its own
particles (Q-vector, Ψ_Q = ½ atan2(Σsin2φ, Σcos2φ)) and correlate against candidates:
⟨cos2(Ψ_Q−phi)⟩ = 0.44, ⟨cos2(Ψ_Q−phi−phir)⟩ = 0.47, vs 0.04 for `phir` and 0.00 for the lab axis.

## Performance — always do these

- **Branch selection**: `chain->SetBranchStatus("*",0)` then enable only used branches
  (`np,bim,phi,phir,px,py,pz,id,ist[,e]`). Skips decompressing x/y/z/t/zus etc. (~half the
  file volume) → measured **~10× wall-clock speedup**. Already in `Elliptic_KP/EP/ESS.C`;
  `Elliptic.C` predates it.
- **Line-buffered logs for background runs**: `stdbuf -oL -eL root -l -b -q 'X.C+(5,0)' > log 2>&1`
  — otherwise stdout is block-buffered and the log looks frozen for many minutes.
- A full 84k-event pass takes ~10–20 min depending on disk contention. If a job "looks hung",
  check `ps` first: state `D` + still-accumulating CPU time = disk-I/O-bound, not deadlocked.
- ACLiC artifacts (`*_C.so`, `*_C.d`, `*.pcm`) may embed pre-move paths; delete and let `.C+` rebuild.

## Analysis scripts (all in `auau200-7/`, run from inside that dir)

All follow the same pattern as the original `analyze.C`/`Centrality.C`: TChain over
`z-*.root`, TLeaf access, b-cut 8.4–9.2 fm, output `cen5.<name>_job0.root`.
Run: `root -l -b -q 'Script.C+(5,0)'`. Figures: matching `MakeFigure_*.C` reads the
output file, so plots regenerate without re-looping.

| analysis | figure macro | output file | contents |
|---|---|---|---|
| `Elliptic.C` | `MakeFigure_Elliptic.C` | `cen5.v2pT_pion_job0.root` | π±/π⁰ v₂(pT), v₂(η) vs true RP and PP; π⁰ selected with ist==1 |
| `Elliptic_KP.C` | `MakeFigure_Elliptic_KP.C` | `cen5.v2pT_kaon_proton_job0.root` | K±, p/p̄ v₂(pT), v₂(η); single pass fills both species |
| `Elliptic_EP.C` | `MakeFigure_Elliptic_EP.C` | `cen5.v2pT_EP_job0.root` | reconstructed-EP methods: full EP (autocorrelation removed by Q-vector subtraction) and η-sub EP (sub-events η<−0.05 / η>0.05, POI vs opposite hemisphere); raw profiles + `pRes` resolution terms; figure macro applies Ollitrault R(χ) correction (Bessel I₀/I₁, bisection inversion) |
| `Elliptic_ESS.C` | `MakeFigure_Elliptic_ESS.C` | `cen5.v2q2_ESS_job0.root` | event-shape selection of arXiv:2307.14997 Eqs. (9)–(11): single/pair q₂², single/pair v₂ for π± (|y|<1, 0.2<pT<2, true RP); two-pass (pass 1: cumulants v₂{2}, v₂,pair{2}; pass 2: normalized q₂² + profiles); pair angle φᵖ = angle of pair momentum sum, ALL charge combos (~5,800 pairs/event, ~5×10⁸ total) |

## Key numeric results (30–40% Au+Au 200 GeV, this dataset)

- **True-RP v₂(pT)**: pions rise from ~0.006 (pT=0.075) to plateau ≈ 0.07–0.08 above 1 GeV;
  kaons similar shape, slightly lower; protons keep climbing to ≈ 0.13 at 2 GeV.
  Mass ordering at low pT and baryon/meson splitting at high pT clearly reproduced.
  Participant-plane values ~8–15% above RP values.
- **EP resolutions** (charged π/K/p Q-vector, 0.2<pT<2, |η|<1): ⟨cos2(Ψ_A−Ψ_B)⟩ = 0.1358
  → R_sub = 0.368, R_full = 0.510 (χ_sub = 0.425). Truth cross-check ⟨cos2(Ψ_full−Ψ_RP)⟩ = 0.372:
  the sub-event method over-estimates R because A/B share PP fluctuations + non-flow;
  consequently v₂{EP} and v₂{sub-EP} sit ~30–50% above v₂{RP,true} (non-flow grows with pT).
- **ESS constants**: v₂{2} = 0.0557, v₂,pair{2} = 0.0810. Sanity: ⟨single q₂²⟩ = 0.992,
  ⟨pair q₂²⟩ = 1.068 (pair >1 is expected — static v₂,pair{2} vs event-by-event N–flow correlation).
- **ESS 4-panel** (`v2_vs_q2_ESS.png`) slopes/intercepts:
  (a) single q₂²→single v₂: 0.040 / 0.001; (b) single→pair: 0.024 / 0.008;
  (c) pair→single: 0.027 / 0.012; (d) pair→pair: 0.028 / 0.001.
  Unmixed combos (a,d) extrapolate to ≈0 at q₂²→0; mixed combos (b,c) keep positive
  intercepts — reproduces the paper's argument for mixed recipes in CME background control.
- Statistics passing cuts (84k events): 12.1M π±, 1.30M K±, 0.56M p+p̄ (|η|<1, pT>0.05).

## Reference

- ESS paper: Xu, Chan, Wang, Tang, Huang, "Event Shape Selection Method in Search of the
  Chiral Magnetic Effect in Heavy-ion Collisions", arXiv:2307.14997. Screenshot of the
  Eq. (9)–(11) page: `auau200-7/Screenshot 2026-07-02 at 3.11.20 PM.png`.
