# lncRNA Function & Literature Report — CRC Biomarker Candidates

> Deep web search · 2026-06-20 · 18 selected lncRNAs investigated

---

## Classification Summary

| Tier | lncRNAs | Status |
|------|---------|--------|
| **Well-characterized** | MIR210HG, EIF3J-AS1, LINC00261, PCAT6, CD27-AS1, ARRDC1-AS1, TNRC6C-AS1 | Multiple publications, known mechanisms |
| **Emerging** | AC113189.5, AC114730.3, LINC00174, AC156455.1, RP11-199F11.2, RP11-277P12.20, RP11-549B18.1, RP5-1021I20.1 | Some literature, mostly signature components |
| **Novel/Uncharacterized** | CTB-25B13.12, AP006621.5, RP11-440D17.3, RP11-268J15.5, RP11-44M6.7 | No published functional studies — **pipeline discoveries** |

---

## 1. Well-Characterized lncRNAs

### 1.1 MIR210HG ⭐ Top Hit

| Property | Finding |
|----------|---------|
| **Expression in CRC** | Consistently upregulated (tissue, serum, cell lines) |
| **Our result** | HR = 1.40, p = 0.0003 (most significant) |
| **Mechanism 1** | ceRNA sponge of miR-1226-3p → proliferation ↑, apoptosis ↓ |
| **Mechanism 2** | Binds PCBP1 → inhibits ferroptosis (GPX4 ↑) |
| **Upstream** | Transcriptionally activated by CREB3 |
| **Diagnostic AUC** | 0.870 (serum), sensitivity 86.5% |
| **Prognosis** | High expression → poor OS (HR = 3.93 in multivariate) |
| **Clinical correlates** | Lymph node metastasis, advanced TNM stage, tumor size |

> **Key references**: Jiang & Zhao (2024) *Turk J Gastroenterol*; *Sci Rep* (2025) 15:871; CREB3-MIR210HG axis in *Transl Cancer Res* (2024)

---

### 1.2 EIF3J-AS1 (EIF3J-DT)

| Property | Finding |
|----------|---------|
| **Our result** | HR = 1.70, p = 0.0036 |
| **Gastric cancer** | Activates autophagy → 5-FU/OXA chemoresistance via ATG14 / miR-188-3p |
| **HCC** | Hypoxia-induced, targets miR-122-5p/CTNND2 axis |
| **Esophageal cancer** | Upregulated, correlates with TNM stage and invasion |

> **Key reference**: EIF3J-DT/ATG14 autophagy axis in drug-resistant gastric cancer

---

### 1.3 LINC00261 🛡️ Protective

| Property | Finding |
|----------|---------|
| **Our result** | HR = 0.90, p = 0.022 (protective) |
| **Expression** | Downregulated in gastric cancer |
| **Mechanism** | Binds Slug protein → enhances GSK3β-Slug interaction → Slug degradation |
| **Effect** | Suppresses EMT, invasion, and lung metastasis in vivo |
| **Diagnostic AUC** | 0.724 (tissue) |

> **Key reference**: Yu Y et al. (2017) *J Cell Mol Med* 21:955–967

---

### 1.4 PCAT6

| Property | Finding |
|----------|---------|
| **Our result** | HR = 1.40, p = 0.0025 |
| **Expression** | Significantly upregulated in CRC |
| **Mechanism** | ceRNA: sponges miR-204 → HMGA2/PI3K ↑; interacts with EZH2 |
| **CRC phenotype** | Promotes EMT, stemness, 5-FU chemoresistance |
| **Pan-cancer** | Oncogenic in prostate, gastric, breast, bladder, ovarian cancers |

> **Key reference**: Ghafouri-Fard S et al. (2021) *Biomed Pharmacother*; Wang S et al. (2021) *Front Oncol*

---

### 1.5 CD27-AS1

| Property | Finding |
|----------|---------|
| **Our result** | HR = 1.80, p = 0.0010 (highest HR among top hits) |
| **Immune link** | Antisense to CD27 (TNF receptor family), enriched in NK cell pathways |
| **AML** | ceRNA: sponges miR-224-5p → PBX3/MAPK ↑ |
| **Melanoma** | CD27-AS1-208 isoform binds STAT3 → proliferation |
| **COAD** | Component of autophagy-related lncRNA prognostic signature |

> **Key reference**: Tao Y et al. (2021) *Cell Death Dis*; Ma J et al. (2022) *Front Oncol*

---

### 1.6 ARRDC1-AS1

| Property | Finding |
|----------|---------|
| **Our result** | HR = 1.80, p = 0.0031 |
| **Pan-cancer** | Oncogenic ceRNA in glioma (miR-432-5p/PRMT5), DLBCL (miR-2355-5p/ATG5), breast (miR-4731-5p/AKT1) |
| **COAD** | Part of 7-lncRNA ferroptosis-related signature; high expression → advanced stage |
| **Pathways** | Wnt/β-catenin, mTORC1, MAPK |

> **Key reference**: Zou et al. (2021) glioma; Xu et al. (2021) DLBCL; COAD ferroptosis signature (2021)

---

### 1.7 TNRC6C-AS1

| Property | Finding |
|----------|---------|
| **Our result** | HR = 1.40, p = 0.0057 |
| **CRC** | Component of pyroptosis-associated 4-lncRNA prognostic signature |
| **CRC liver metastasis** | Key lncRNA in ceRNA network, Notch pathway |
| **Thyroid cancer** | Suppresses TNRC6C expression → proliferation, migration ↑ |

> **Key reference**: Pyroptosis signature in CRC; ceRNA network in CRC liver metastases (*Aging*, 2021)

---

## 2. Emerging lncRNAs (Limited Literature)

### 2.1 AC113189.5
- **Our result**: HR = 1.50, p = 0.0005
- Reported in colon cancer ferroptosis-related lncRNA signatures
- No dedicated functional studies — **pipeline discovery potential**

### 2.2 AC114730.3
- **Our result**: HR = 0.87, p = 0.012 (protective)
- Upregulated in tumor (logFC = 1.69)
- No published functional studies — **novel protective lncRNA candidate**

### 2.3 LINC00174
- **Our result**: HR = 1.40, p = 0.0032
- Reported in pan-cancer lncRNA expression analyses
- Limited functional characterization in CRC

### 2.4 AC156455.1
- **Our result**: HR = 1.20, p = 0.0061
- Occasionally appears in lncRNA prognostic signatures
- No dedicated functional studies

### 2.5 RP11-199F11.2
- **Our result**: HR = 1.40, p = 0.0054
- Appears in some cancer lncRNA profiling studies
- Function unknown

### 2.6 RP11-277P12.20
- **Our result**: HR = 0.87, p = 0.0040 (protective)
- Rarely reported — **novel protective lncRNA candidate**

### 2.7 RP11-549B18.1
- **Our result**: HR = 0.71, p = 0.0052 (strongest protective effect)
- No published studies — **highly novel protective lncRNA**

### 2.8 RP5-1021I20.1
- **Our result**: HR = 1.30, p = 0.0070
- Appears in gene annotation databases
- No dedicated cancer studies

---

## 3. Novel/Uncharacterized lncRNAs ⭐ Pipeline Discoveries

These lncRNAs have **no published functional studies** and represent **novel findings** from the pipeline:

| lncRNA | Our Finding | HR | p-value | Potential Significance |
|--------|------------|-----|---------|----------------------|
| **CTB-25B13.12** | 3rd most significant | 1.70 | 0.0004 | High-priority novel oncogenic candidate |
| **AP006621.5** | 2nd most significant | 1.60 | 0.0003 | High-priority novel oncogenic candidate |
| **RP11-440D17.3** | Stepwise Cox selected | 1.10 | 0.064 | Borderline significant, Cox-selected feature |
| **RP11-268J15.5** | Cox significant | 1.20 | 0.0045 | Novel risk factor |
| **RP11-44M6.7** | Cox significant | 1.50 | 0.0012 | Novel risk factor, substantial HR |

**This is a strength**: The pipeline identified lncRNAs with strong statistical association with CRC survival that have never been functionally characterized — representing genuine discovery potential.

---

## 4. Pathway & Functional Enrichment Summary

| Pathway/Process | Associated lncRNAs |
|----------------|-------------------|
| **Ferroptosis** | MIR210HG (inhibits via PCBP1), ARRDC1-AS1 |
| **EMT / Metastasis** | LINC00261 (suppresses via Slug), PCAT6 (promotes) |
| **Autophagy** | EIF3J-AS1 (activates via ATG14), ARRDC1-AS1 |
| **Immune modulation** | CD27-AS1 (NK cells, TNF pathway), PCAT6 (pDC/Treg) |
| **ceRNA sponging** | MIR210HG, PCAT6, CD27-AS1, ARRDC1-AS1, TNRC6C-AS1 |
| **Wnt/β-catenin** | ARRDC1-AS1, PCAT6 |
| **PI3K/AKT/mTOR** | PCAT6, ARRDC1-AS1 |
| **MAPK signaling** | CD27-AS1, ARRDC1-AS1 |
| **Notch signaling** | TNRC6C-AS1 |
| **5-FU resistance** | EIF3J-AS1, PCAT6 |

---

## 5. Concordance with Our Pipeline Results

| Comparison | Finding |
|------------|---------|
| **MIR210HG as top hit** | ✅ Consistent — most studied CRC lncRNA, most significant in our Cox analysis |
| **LINC00261 as protective** | ✅ Consistent — known tumor suppressor, our HR = 0.90 |
| **PCAT6 oncogenic** | ✅ Consistent — known pan-cancer oncogene, our HR = 1.40 |
| **EIF3J-AS1 drug resistance** | ✅ Consistent — known autophagy/drug resistance mediator |
| **Novel lncRNAs** | CTB-25B13.12, AP006621.5, RP11-549B18.1 are **genuinely novel findings** |

---

## 6. Research Recommendations

1. **High priority for experimental validation**:
   - CTB-25B13.12 (HR = 1.70, p = 0.0004) — novel, highly significant
   - AP006621.5 (HR = 1.60, p = 0.0003) — novel, highly significant
   - RP11-549B18.1 (HR = 0.71, p = 0.0052) — strongest protective effect, novel

2. **Mechanistic studies** for MIR210HG/PCBP1/ferroptosis axis in our CRC cohort

3. **ceRNA network construction** for TNRC6C-AS1 and CD27-AS1 in CRC

4. **Drug resistance validation** for EIF3J-AS1 and PCAT6 in CRC cell lines

5. **Multi-omics integration** — combine lncRNA signature with methylation/immune profiling

---

*Sources: PubMed, PMC, Google Scholar, Semantic Scholar, GeneCards, Lnc2Cancer, DOAJ*
