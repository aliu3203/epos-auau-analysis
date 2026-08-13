# External reference data

`HEPData-ins2928164-v1-csv/` — STAR measurement, downloaded from HEPData.

- INSPIRE record: `ins2928164`
- HEPData record DOI: `10.17182/hepdata.159952.v1`
- Content: Au+Au, sqrt(s_NN) = 7.7–200 GeV; centrality dependence of Npart*Δγ and the
  event-shape-selected Npart*Δγ_ESS, plus the κ tables. Each CSV carries its own
  `table_doi` and description in the header comments.
- License: HEPData releases records under CC0 (public domain), so redistribution here
  is unrestricted.

## Why this is in the repo

**Nothing reads these files at runtime.** The STAR points are hardcoded as literal arrays
in the figure macros:

- `auau200-7/MakeFigure_Cent_STAR.C`
- `auau200-7/MakeFigure_Kappa_STAR.C`

(`Correlators_Cent.C` and `Correlators_mesons.C` reference the record in comments only.)

These CSVs are therefore the **provenance for those hardcoded numbers** — the thing to
diff against if a comparison plot ever looks wrong, and the thing a referee will ask for.
Do not delete them on the grounds that no code opens them.
