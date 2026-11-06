from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from PIL import Image
F = 'D:/Math-model/results/figures/'
NAVY = RGBColor(0x1F, 0x3A, 0x5F); GREY = RGBColor(0x44, 0x44, 0x44)
prs = Presentation(); prs.slide_width = Inches(13.333); prs.slide_height = Inches(7.5)
blank = prs.slide_layouts[6]


def tb(s, x, y, w, h, text, size=18, bold=False, color=GREY):
    t = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h)).text_frame
    t.word_wrap = True
    items = text if isinstance(text, list) else [text]
    for i, it in enumerate(items):
        p = t.paragraphs[0] if i == 0 else t.add_paragraph()
        p.text = it; p.font.size = Pt(size); p.font.bold = bold
        p.font.color.rgb = color; p.space_after = Pt(8)


def base(title):
    s = prs.slides.add_slide(blank)
    bar = s.shapes.add_shape(1, 0, 0, prs.slide_width, Inches(1.0))
    bar.fill.solid(); bar.fill.fore_color.rgb = NAVY; bar.line.fill.background()
    tb(s, 0.5, 0.15, 12.3, 0.8, title, 28, True, RGBColor(255, 255, 255))
    return s


def text_slide(title, bul):
    s = base(title); tb(s, 0.7, 1.4, 11.9, 5.6, ['\u2022 ' + b for b in bul], 22); return s


def img_slide(title, img, bul, split=0.64):
    s = base(title)
    im = Image.open(F + img + '.png'); ar = im.width / im.height
    mw = 13.333 * split - 0.6; mh = 5.9
    w = mw; h = w / ar
    if h > mh:
        h = mh; w = h * ar
    s.shapes.add_picture(F + img + '.png', Inches(0.4), Inches(1.2 + (mh - h) / 2), Inches(w), Inches(h))
    x = 0.4 + mw + 0.3
    tb(s, x, 1.4, 13.333 - x - 0.4, 5.6, ['\u2022 ' + b for b in bul], 17); return s


s = prs.slides.add_slide(blank)
bg = s.shapes.add_shape(1, 0, 0, prs.slide_width, prs.slide_height)
bg.fill.solid(); bg.fill.fore_color.rgb = NAVY; bg.line.fill.background()
tb(s, 0.8, 2.2, 11.7, 2, 'Predicting Cancer Relapse Through Mathematical Modelling of Tumour Dormancy', 38, True, RGBColor(255, 255, 255))
tb(s, 0.8, 4.5, 11.7, 1.5, ['M.Sc. project: Phases 1-5 results', 'Proliferating-Dormant-Immune (P-D-E) ODE model, sensitivity, identifiability and METABRIC comparison'], 20, False, RGBColor(0xDD, 0xE6, 0xF0))

text_slide('Background and motivation', ['Breast cancer can relapse years or decades after apparently successful treatment', 'Dormant (growth-arrested) cells and immune control are proposed explanations for delayed relapse', 'Aim: build an interpretable mathematical model, study when relapse occurs, and check what clinical data can and cannot say', 'All model parameter values are hypothetical; results are not clinical predictions'])
text_slide('Research gap and questions', ['Dormancy-immune ODE models are mostly analysed qualitatively; few are turned into relapse-time distributions compared with censored data', 'Identifiability and endpoint mismatch (time origin, composite event) are rarely confronted explicitly', 'RQ1-RQ6: immune-control bounds, role of dormancy, relapse-time control, estimating relapse likelihood, matching group-level recurrence, identifiability of the P/D split', 'Negative or limiting results are treated as first-class outputs'])
text_slide('The P-D-E model', ['P: proliferating cells, D: dormant cells, E: effective immune activity (cell fractions of capacity K)', 'dP/dt = rP(1 - N/K) - qP + aD - kEP', 'dD/dt = qP - (a + \u03b4_D)D', 'dE/dt = s + \u03c1\u00b7P/(h+P)\u00b7E - \u03bcE - \u03b7NE,   N = P + D', 'Baseline: r=0.12, q=0.04, a=0.005, \u03b4_D=0.001, k=0.1, s=0.02, \u03bc=0.05, \u03c1=0.04, h=0.02, \u03b7=0.08, K=1', 'Relapse event: first upward crossing of N_detect = 0.1; horizon 240 months'])
img_slide('Baseline trajectories', 'baseline_trajectories', ['Deterministic solution from a small post-treatment residual burden', 'Long quiescent phase, then escape: delayed relapse emerges from slow growth near the tumour-free state', 'Shows the three compartments and the detection threshold'])
img_slide('Stability and bifurcation structure', 'fig_B_stability_bifurcation_diagram', ['Tumour-free state is stable iff r < ks/\u03bc + q\u03b4_D/(a+\u03b4_D) (proved)', 'Transcritical bifurcation at k_c = 0.28333', 'Bistability (two kinds) found numerically; no Hopf bifurcation detected', 'If \u03b4_D = 0, dormancy alone does not change tumour-free stability'])
img_slide('Equilibrium branches', 'fig_A_equilibrium_branches', ['Nonzero equilibria obtained by reduction to one scalar equation', 'Below k_c a stable nonzero equilibrium with high dormant fraction exists', 'Its burden tends to 0 as k approaches k_c'])
img_slide('Phase portraits', 'fig_E_phase_portraits', ['Trajectories converge to different attractors depending on initial state', 'Separatrix structure explains why small residual disease can be controlled or escape'])
img_slide('Relapse-time estimators', 'fig_G_relapse_time_approximation', ['Naive linear-escape estimate ln(N_detect/N0)/\u03bb+ is not supported (about 240 vs 375 months at k = 0.5 k_c)', 'Eigenvector-projection and nonlinear corrections improve the estimate', 'Relapse timing is set by residual burden and proximity to the stability boundary'])
text_slide('Phase 2 summary', ['Closed-form tumour-free stability boundary verified against eigenvalues', 'Bifurcations classified: transcritical, saddle-node, bistability; no Hopf', 'M1 normalisation (e = \u03bcE/s, \u03b2 = ks/\u03bc) removes the k-s gauge redundancy', 'Dormancy and relapse thresholds are illustrative; threshold sensitivity is tested in Phase 3'])
img_slide('Local sensitivity (Phase 3)', 'fig_p3_01_local_sensitivity_tornado', ['Local sensitivity of relapse outcome to each parameter around the baseline', 'Parameters near the immune-control balance dominate'])
img_slide('Global sensitivity (Phase 3)', 'fig_p3_02_global_sensitivity_ranking', ['Spearman, PRCC, Morris and Sobol analyses compared', 'Parameter interactions matter: effects depend on proximity to the bifurcation'])
img_slide('Basins of attraction', 'fig_p3_04_basins_tumour_free_nonzero', ['Initial conditions partition into tumour-free and nonzero outcomes', 'Event occurrence depends on both parameters and residual state'])
img_slide('Simulated cumulative relapse probability', 'fig_p3_07_cumulative_relapse_probability', ['Monte Carlo over priors P1-P7 gives model event curves', 'Prior P1 (n = 15,000): P(5 y) = 0.0005, P(10 y) = 0.091, P(20 y) = 0.436', 'Curve shape depends strongly on the prior'])
img_slide('Event probabilities at 5, 10, 20 years', 'fig_p3_08_event_probabilities_5_10_20y', ['Simulated event probabilities with Monte Carlo intervals across priors', 'Events are strongly back-loaded: very few early, many late'])
img_slide('Detection threshold sensitivity', 'fig_p3_09_detection_threshold_sensitivity', ['Moving N_detect shifts the whole event-time distribution', 'The threshold is a modelling convention, so absolute probabilities are threshold-dependent'])
img_slide('Identifiability: profile likelihoods (Phase 4)', 'fig_p4_03_profile_likelihoods', ['Some parameter combinations are well constrained; others have flat profiles', 'Total-burden data alone cannot separate the P/D split'])
img_slide('Structural vs practical identifiability', 'fig_p4_11_structural_vs_practical_identifiability', ['Structural redundancy (k-s) plus practical limits from noise and sampling', 'Pre-registered criteria: A (descriptive survival comparison) yes; B (mechanistic calibration) not supported; C (patient-level prediction) unsuitable'])
img_slide('Time-origin uncertainty', 'fig_p4_10_time_origin_uncertainty', ['METABRIC RFS origin is not documented as treatment completion', 'Offsets confound model time and clinical time, so absolute timing cannot be fitted'])
img_slide('METABRIC: overall relapse-free survival', 'fig_p5_01_metabric_km_overall', ['2,380 usable records, 964 events', 'Cumulative incidence (95% CI): 5 y 0.263 (0.245-0.281); 10 y 0.376 (0.355-0.397); 15 y 0.456 (0.433-0.480); 20 y 0.505 (0.478-0.532)', 'Kaplan-Meier with Greenwood log-log CI; no competing-risk analysis'])
img_slide('ER-positive vs ER-negative', 'fig_p5_02_er_positive_vs_negative', ['ER- has more early events (0.390 vs 0.225 at 5 y); curves converge by 15-20 y (0.497 vs 0.515)', 'ER- interval hazard falls 10.35 to 0.91 per 100 person-years; ER+ only 5.07 to 2.61', 'Log-rank p < 0.001 at 5 and 10 y'])
img_slide('Model vs METABRIC (descriptive)', 'fig_p5_05_model_vs_metabric', ['Pre-registered grid of 1,764 scenarios; nothing fitted or tuned', 'Reference scenario: simulated 0.0005 / 0.091 / 0.436 at 5 / 10 / 20 y vs METABRIC 0.263 / 0.376 / 0.505'])
img_slide('Consistency envelope', 'fig_p5_06_consistency_envelope', ['5 y: simulated incidence below the METABRIC CI in all 1,764 scenarios', '10 y below in 98.8%; late events more back-loaded in all scenarios', '20-y value inside the CI in only 11.1% of scenarios'])
img_slide('Phase 5 summary', 'fig_p5_10_clinical_comparison_summary', ['The model has no treatment era or early-recurrence mechanism, so early mismatch is expected', 'Neither validation nor refutation of the mechanism', 'Threshold and time-origin choices change magnitude but not the qualitative gaps'])
text_slide('Conclusions and limitations', ['Analytic stability boundary and bifurcation structure of the P-D-E model established', 'Relapse-time distributions are strongly back-loaded and prior- and threshold-dependent', 'Total-burden data cannot identify the P/D split; mechanistic calibration is not supported', 'METABRIC comparison is descriptive: undocumented RFS origin, composite event, curve ends at about 230 months', 'Future work: treatment-era and early-recurrence mechanisms, stochastic extinction, data with dormancy markers', '164 automated tests pass; outputs reproducible from the repository'])
prs.save('D:/Math-model/report/Project_Presentation.pptx'); print(len(prs.slides))
