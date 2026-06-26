library(GEOquery); library(SummarizedExperiment)

cat('═══════════ CROSS-COHORT P-VALUES ═══════════\n\n')

# ── Load data ──
clin <- read.csv('../02_Data/raw/tcga_coad_clinical_followup_data.csv', stringsAsFactors=FALSE, row.names=1)
c_age <- clin$age_at_initial_pathologic_diagnosis
c_gender <- clin$gender
c_stage_raw <- tolower(gsub('stage ','',clin$pathologic_stage))
c_stage <- ifelse(grepl('^i[abc]?$',c_stage_raw)&!grepl('ii',c_stage_raw),'I',
           ifelse(grepl('^iii',c_stage_raw),'III',
           ifelse(grepl('^ii',c_stage_raw),'II',
           ifelse(grepl('^iv',c_stage_raw),'IV',NA))))
c_dead <- ifelse(clin$vital_status=='dead',1,0)

se <- readRDS('external_validation/READ/READ_SE.rds'); cd <- as.data.frame(colData(se))
r_age <- cd$age_at_index; r_gender <- cd$gender
sm <- c('Stage I'='I','Stage IA'='I','Stage II'='II','Stage IIA'='II','Stage IIB'='II',
        'Stage III'='III','Stage IIIA'='III','Stage IIIB'='III','Stage IIIC'='III',
        'Stage IV'='IV','Stage IVA'='IV')
r_stage <- unname(sm[cd$ajcc_pathologic_stage]); r_stage[is.na(r_stage)] <- NA
r_dead <- ifelse(cd$vital_status=='Dead',1,0)

gse3 <- getGEO('GSE39582', GSEMatrix=TRUE); pd3 <- pData(gse3[[1]])
g3_age <- as.numeric(pd3[['age.at.diagnosis (year):ch1']])
g3_gender <- pd3[['Sex:ch1']]
g3_stage <- pd3[['ajcc.stage:ch1']]; g3_stage[g3_stage=='0'] <- NA
g3_dead <- as.numeric(pd3[['os.event:ch1']])

gse1 <- getGEO('GSE17536', GSEMatrix=TRUE); pd1 <- pData(gse1[[1]])
g1_age <- as.numeric(pd1[['age (year):ch1']])
g1_gender <- pd1[['gender:ch1']]
g1_stage <- pd1[['ajcc_stage:ch1']]; g1_stage[g1_stage=='0'] <- NA
g1_evt <- pd1[['overall_event (death from any cause):ch1']]
g1_dead <- ifelse(grepl('death',g1_evt,ignore.case=TRUE)&!grepl('no death',g1_evt,ignore.case=TRUE),1,0)

# ── 1. AGE ──
cat('─── Age (ANOVA) ───\n')
age_df <- data.frame(
  age = c(c_age, r_age, g3_age, g1_age),
  cohort = factor(rep(c('COAD','READ','GSE39582','GSE17536'),
                c(length(c_age),length(r_age),length(g3_age),length(g1_age)))))
age_df <- age_df[!is.na(age_df$age),]
aov_fit <- summary(aov(age ~ cohort, data=age_df))
cat(sprintf('F=%.3f, p=%.4f\n', aov_fit[[1]][['F value']][1], aov_fit[[1]][['Pr(>F)']][1]))

# ── 2. GENDER ──
cat('\n─── Gender (Chi-sq) ───\n')
gender_tbl <- rbind(
  COAD=c(sum(c_gender=='male',na.rm=TRUE),sum(c_gender=='female',na.rm=TRUE)),
  READ=c(sum(r_gender=='male',na.rm=TRUE),sum(r_gender=='female',na.rm=TRUE)),
  GSE39582=c(sum(g3_gender=='Male',na.rm=TRUE),sum(g3_gender=='Female',na.rm=TRUE)),
  GSE17536=c(sum(g1_gender=='male',na.rm=TRUE),sum(g1_gender=='female',na.rm=TRUE)))
colnames(gender_tbl) <- c('Male','Female')
gt <- chisq.test(gender_tbl)
cat(sprintf('Chi-sq=%.2f, df=%d, p=%.4f\n', gt$statistic, gt$parameter, gt$p.value))

# ── 3. STAGE (COAD vs READ vs GSE17536) ──
cat('\n─── Stage (Chi-sq, 3 cohorts) ───\n')
stage_tbl <- rbind(
  COAD=c(sum(c_stage=='I',na.rm=TRUE),sum(c_stage=='II',na.rm=TRUE),sum(c_stage=='III',na.rm=TRUE),sum(c_stage=='IV',na.rm=TRUE)),
  READ=c(sum(r_stage=='I',na.rm=TRUE),sum(r_stage=='II',na.rm=TRUE),sum(r_stage=='III',na.rm=TRUE),sum(r_stage=='IV',na.rm=TRUE)),
  GSE17536=c(sum(g1_stage=='1',na.rm=TRUE),sum(g1_stage=='2',na.rm=TRUE),sum(g1_stage=='3',na.rm=TRUE),sum(g1_stage=='4',na.rm=TRUE)))
colnames(stage_tbl) <- c('I','II','III','IV')
st <- chisq.test(stage_tbl)
cat(sprintf('Chi-sq=%.2f, df=%d, p=%.4f\n', st$statistic, st$parameter, st$p.value))

# ── 4. MORTALITY ──
cat('\n─── Mortality (Chi-sq, 4 cohorts) ───\n')
death_tbl <- rbind(
  COAD=c(sum(c_dead==1,na.rm=TRUE),sum(c_dead==0,na.rm=TRUE)),
  READ=c(sum(r_dead==1,na.rm=TRUE),sum(r_dead==0,na.rm=TRUE)),
  GSE39582=c(sum(g3_dead==1,na.rm=TRUE),sum(g3_dead==0,na.rm=TRUE)),
  GSE17536=c(sum(g1_dead==1,na.rm=TRUE),sum(g1_dead==0,na.rm=TRUE)))
colnames(death_tbl) <- c('Dead','Alive')
dt <- chisq.test(death_tbl)
cat(sprintf('Chi-sq=%.2f, df=%d, p=%.4f\n', dt$statistic, dt$parameter, dt$p.value))

# ── 5. PAIRWISE: COAD vs each ──
cat('\n─── Pairwise: COAD vs each validation ───\n')
for (nm in c('READ','GSE39582','GSE17536')) {
  cat(sprintf('\nCOAD vs %s:\n', nm))
  if (nm=='READ') { va<-r_age; vg<-r_gender; vd<-r_dead; vs<-r_stage }
  else if (nm=='GSE39582') { va<-g3_age; vg<-g3_gender; vd<-g3_dead; vs<-g3_stage }
  else { va<-g1_age; vg<-g1_gender; vd<-g1_dead; vs<-g1_stage }

  tt <- t.test(c_age, va); cat(sprintf('  Age: t=%.3f, p=%.4f\n', tt$statistic, tt$p.value))

  if (nm=='GSE39582') {
    gt2 <- chisq.test(rbind(c(sum(c_gender=='male',na.rm=TRUE),sum(c_gender=='female',na.rm=TRUE)),
                            c(sum(vg=='Male',na.rm=TRUE),sum(vg=='Female',na.rm=TRUE))))
  } else {
    gt2 <- chisq.test(rbind(c(sum(c_gender=='male',na.rm=TRUE),sum(c_gender=='female',na.rm=TRUE)),
                            c(sum(vg=='male',na.rm=TRUE),sum(vg=='female',na.rm=TRUE))))
  }
  cat(sprintf('  Gender: Chi-sq=%.2f, p=%.4f\n', gt2$statistic, gt2$p.value))

  dt2 <- chisq.test(rbind(c(sum(c_dead==1,na.rm=TRUE),sum(c_dead==0,na.rm=TRUE)),
                           c(sum(vd==1,na.rm=TRUE),sum(vd==0,na.rm=TRUE))))
  cat(sprintf('  Death: Chi-sq=%.2f, p=%.4f\n', dt2$statistic, dt2$p.value))

  if (nm != 'GSE39582') {
    if (nm=='READ') st2_tbl <- rbind(
      c(sum(c_stage=='I',na.rm=TRUE),sum(c_stage=='II',na.rm=TRUE),sum(c_stage=='III',na.rm=TRUE),sum(c_stage=='IV',na.rm=TRUE)),
      c(sum(vs=='I',na.rm=TRUE),sum(vs=='II',na.rm=TRUE),sum(vs=='III',na.rm=TRUE),sum(vs=='IV',na.rm=TRUE)))
    else st2_tbl <- rbind(
      c(sum(c_stage=='I',na.rm=TRUE),sum(c_stage=='II',na.rm=TRUE),sum(c_stage=='III',na.rm=TRUE),sum(c_stage=='IV',na.rm=TRUE)),
      c(sum(vs=='1',na.rm=TRUE),sum(vs=='2',na.rm=TRUE),sum(vs=='3',na.rm=TRUE),sum(vs=='4',na.rm=TRUE)))
    st2 <- chisq.test(st2_tbl)
    cat(sprintf('  Stage: Chi-sq=%.2f, df=%d, p=%.4f\n', st2$statistic, st2$parameter, st2$p.value))
  }
}
cat('\n═══ DONE ═══\n')
