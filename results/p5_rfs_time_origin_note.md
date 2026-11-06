# RFS time-origin limitation (Phase 5 descriptive METABRIC analysis)

**Verified (documentation / exact computation, Phase 4):**
- RFS_STATUS is a composite endpoint (loco-regional relapse, distant relapse or disease-specific death); the cBioPortal export does not
  provide the time of each component, nor which event came first.
- OS_MONTHS is documented as measured from initial diagnosis. **No documentation of the RFS_MONTHS origin was found** (export metadata,
  open text of the METABRIC papers, datahub history). The export has no treatment or surgery dates.

**Inference (not proof of the origin):** RFS_MONTHS equals OS_MONTHS exactly for all 745 patients not recurred and living and all 432 who
died of other causes (0 violations of RFS <= OS in 1,980 patients with both times), which suggests a common clock for the two
endpoints. This does not establish whether the clock starts at diagnosis, surgery or another date.

**Consequence for this analysis:** every Kaplan-Meier time in Phase 5 is "months since the (unknown) RFS origin". The offset between
diagnosis and the RFS origin, the diagnosis-to-treatment interval and any detection delay are nuisance quantities (Phase 4, section 5-6)
and are NOT estimated here. Results are descriptive of the exported variable; they are not comparable to model times without a
stated time-origin assumption, and they are not model validation.

Analysis sample: 2380 patients with usable RFS time and status (964 events). Cumulative incidence in this file agrees with the
Phase 1 table results/tables/metabric_km_overall.csv to within 2.6e-05.
