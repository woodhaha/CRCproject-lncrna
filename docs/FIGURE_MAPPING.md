# Figure Mapping — Pipeline Output → Manuscript

Copy pipeline output PDFs to manuscript figure names:

```bash
cd D:\Researching\CRCproject

# Main Figures
cp outputs/Fig_2_feature_lncRNAs.pdf          paper/Figure-1.pdf   # Study flowchart (adapt)
cp outputs/Fig_3a_crossvalidation_deviance.pdf paper/Figure-2a.pdf  # Elastic Net surface
cp outputs/Fig_3b_elastic_net_deviance.pdf     paper/Figure-2b.pdf  # CV deviance curve
cp outputs/Fig_3c_coefficient_paths.pdf        paper/Figure-2c.pdf  # Coefficient paths (lambda)
cp outputs/Fig_3d_fraction_deviance.pdf        paper/Figure-2d.pdf  # Coefficient paths (deviance)
cp outputs/Fig_3_forestplot_cox.pdf            paper/Figure-3.pdf   # Forest plot
cp outputs/Fig_6_survival_train.pdf            paper/Figure-4a.pdf  # KM curve (train)
cp outputs/Fig_6_survival_test.pdf             paper/Figure-4b.pdf  # KM curve (test)
cp outputs/Fig_9_rf_survival_training.pdf      paper/Figure-4c.pdf  # RF survival (train)
cp outputs/Fig_9_rf_survival_testing.pdf       paper/Figure-4d.pdf  # RF survival (test)
cp outputs/Fig_4_roc_1year.pdf                 paper/Figure-5a.pdf  # ROC 1yr
cp outputs/Fig_4_roc_3year.pdf                 paper/Figure-5b.pdf  # ROC 3yr
cp outputs/Fig_4_roc_5year.pdf                 paper/Figure-5c.pdf  # ROC 5yr
cp outputs/Fig_5_auc_over_time.pdf             paper/Figure-5d.pdf  # AUC over time
cp outputs/Fig_12_brier_train.pdf              paper/Figure-6a.pdf  # Brier (train)
cp outputs/Fig_12_brier_test.pdf               paper/Figure-6b.pdf  # Brier (test)
cp outputs/Fig_13_cindex_train.pdf             paper/Figure-6c.pdf  # C-index (train)
cp outputs/Fig_13_cindex_test.pdf              paper/Figure-6d.pdf  # C-index (test)
cp outputs/Fig_14_roc_train.pdf                paper/Figure-7a.pdf  # ROC classifiers (train)
cp outputs/Fig_14_roc_test.pdf                 paper/Figure-7b.pdf  # ROC classifiers (test)
cp outputs/Fig_15_lift_train.pdf               paper/Figure-7c.pdf  # Lift (train)
cp outputs/Fig_15_lift_test.pdf                paper/Figure-7d.pdf  # Lift (test)

# Supplementary
cp outputs/Fig_S1_survival_by_stage.pdf        paper/Figure-S1.pdf
cp outputs/Fig_S2_rf_vimp.pdf                  paper/Figure-S2a.pdf
cp outputs/Fig_S3_rf_vimp_bar.pdf              paper/Figure-S2b.pdf
cp outputs/Fig_S4_feature_density.pdf          paper/Figure-S4.pdf
cp outputs/Fig_S5_feature_boxplot.pdf          paper/Figure-S5.pdf
cp outputs/Fig_S6_model_tuning.pdf             paper/Figure-S6.pdf
```

Or use the `copy_figures.ps1` / `copy_figures.sh` script generated from this mapping.
