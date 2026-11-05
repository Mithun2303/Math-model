# Shared conventions for report-part writers (COMPLETE_RESULTS_REPORT)

Output: one Markdown part file per assignment in `results/_report_parts/`. NO title, NO front matter. Start directly with the numbered `## N. Title` heading of your first section. The assembler concatenates parts in section order, converts to HTML (python-markdown + KaTeX) and prints to PDF with a headless Chromium. There is no pandoc/LaTeX.

## Evidence rules
- Read the relevant PHASE_n_SUMMARY.md, the detailed `report/phase_*.md`, literature files, configs, experiment scripts and CSV tables (`pd.read_csv(path, comment="#")`; use `/Users/mithunkarthickvenkatesan/tumour-dormancy-relapse/.venv/bin/python`). **Open (Read tool) every PNG you describe before describing it.**
- Every number in text/tables must come from a table, report or summary you read; copy/compute it from the CSV. Do not invent or infer values from filenames. If a value is unavailable, omit it or say "not computed". Report the precision the source supports.
- Separate **proved** (analytic derivation / exact computation) from **numerical observation** from **assumption/hypothetical**. All parameter values are hypothetical; no clinical/patient-risk claims; METABRIC results are descriptive, never "validation".
- Do not rerun expensive simulations; do not modify models, data, earlier outputs.

## Format
- Headings: `## 7. Title`, `### 7.1 Title`, `#### 7.1.1 Title` with explicit numbers.
- Math: `$inline$` and display math on their own lines as `$$ ... $$` (KaTeX-compatible: use `\begin{aligned}...\end{aligned}`, `\mathrm`, `\partial`; NO `\begin{align}`, NO `\label`, NO `\eqref`). Number key equations by appending `\tag{7.1}` inside the $$ block. Escape nothing else. Notation: P, D (cell fractions of K), N=P+D, e (immune effectiveness), parameters r,q,a,δ_D (write `\delta_D`),β,ρ,η,μ,h,K as in the Phase 2 report; use the SAME symbols as `report/phase_2_mathematical_analysis.md`.
- Figures: for EACH figure assigned to you, exactly this block (path relative to results/, real filename, shown ONCE in the whole report):
  ```
  ![Figure 7.1. Descriptive title](figures/real_filename.png)

  **Figure 7.1.** Descriptive title. *Setting:* parameters/experiment. *Axes:* ... *Shows:* ... *Mathematical interpretation:* ... *Biological interpretation:* ... *Limitations:* ...
  ```
  Figure numbers are `section.counter` (e.g. 7.1, 7.2). Refer to figures in prose as "Figure 7.1". Cover all nine items (number, title, image, setting, axes, what it shows, mathematical interpretation, biological interpretation, limitations). Keep each caption paragraph 80–200 words, accurate to what you SAW in the image.
- Tables: caption line `**Table 7.1.** Caption (source: results/tables/xxx.csv).` immediately followed by a blank line and a pipe table. Numbered `section.counter`. Keep tables ≤ 9 columns (they are printed on A4 portrait); round sensibly.
- Citations: `[@bibkey]` or `[@key1; @key2]` using ONLY keys that exist in `report/references.bib` (grep to verify). Do not cite anything not in the bib. Statements of literature facts must be ones the bib entry plausibly supports; prefer to cite sparingly and accurately.
- Prose: professional research-paper style, explicit assumptions, no hype, no unsupported claims, connect to research questions. Self-contained for a faculty examiner who will not read source code. Mention config/seed/solver settings where relevant (in text or one table).
- Cross-reference other sections as "Section N".
- When done, reply with: file path(s), number of figures, number of tables, word count, any values you could not verify, and any inconsistency found between sources. Keep that reply short.
