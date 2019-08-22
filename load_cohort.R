load_cohort = function(studies, k=6, readlen=150) {
  load_kmer(studies, k=k, readlen=readlen)
  
  if ('amerigut' %in% studies) {
    amerigut.kmer.humanstool = amerigut.kmer[replace_na(amerigut.meta$host == 'human' & amerigut.meta$sample_type == 'stool', replace=FALSE),]
    amerigut.meta.humanstool = amerigut.meta[replace_na(amerigut.meta$host == 'human' & amerigut.meta$sample_type == 'stool', replace=FALSE),]
    
    amerigut.kmer.healthy = amerigut.kmer.humanstool[is.na(amerigut.meta.humanstool$ibd_diagnosis) &
                                                       replace_na(amerigut.meta.humanstool$mental_illness_anorexia == "false" &
                                                                    amerigut.meta.humanstool$cdiff == 'I do not have this condition' &
                                                                    amerigut.meta.humanstool$ibs == 'I do not have this condition' &
                                                                    amerigut.meta.humanstool$collection_date >= 20170000, replace=FALSE),]
    
    amerigut.kmer.uc = amerigut.kmer.humanstool[replace_na(amerigut.meta.humanstool$ibd_diagnosis == "Ulcerative colitis", replace=FALSE),]
    amerigut.kmer.cd = amerigut.kmer.humanstool[replace_na(amerigut.meta.humanstool$ibd_diagnosis == "Crohn's disease", replace=FALSE),]
    amerigut.kmer.an = amerigut.kmer.humanstool[replace_na(amerigut.meta.humanstool$mental_illness_anorexia == "true", replace=FALSE),]
    amerigut.kmer.cdiff = amerigut.kmer.humanstool[replace_na(amerigut.meta.humanstool$cdiff == "Diagnosed by a medical professional (doctor, physician assistant)", replace=FALSE),]
    amerigut.kmer.ibs = amerigut.kmer.humanstool[replace_na(amerigut.meta.humanstool$ibs == "Diagnosed by a medical professional (doctor, physician assistant)", replace=FALSE),]
    
    amerigut.kmer.redmeatfreq = amerigut.kmer.humanstool[replace_na(amerigut.meta.humanstool$high_fat_red_meat_freq %in% c('Daily', 'Regularly (3-5 times/week)'), replace=FALSE),]
    amerigut.kmer.redmeatnever = amerigut.kmer.humanstool[replace_na(amerigut.meta.humanstool$high_fat_red_meat_freq == 'Never', replace=FALSE),]
    
    amerigut.kmer.recent = amerigut.kmer.humanstool[replace_na(amerigut.meta.humanstool$collection_date >= 20170000, replace=FALSE),]
  }
  
  if ('anorex' %in% studies) {
    assign('anorex.kmer.healthy', anorex.kmer[anorex.meta$pheno == 'healthy',], envir=.GlobalEnv)
    assign('anorex.kmer.an', anorex.kmer[anorex.meta$pheno == 'anorexia',], envir=.GlobalEnv)
  }
  
  if ('apcire' %in% studies) {
    assign('apcire.kmer.control.stool', apcire.kmer[replace_na(apcire.meta$sample_type == 'stool' & apcire.meta$pheno == 'control', replace=FALSE),], envir=.GlobalEnv)
    assign('apcire.kmer.adenoma.stool', apcire.kmer[replace_na(apcire.meta$sample_type == 'stool' & apcire.meta$pheno == 'adenoma (polyp)', replace=FALSE),], envir=.GlobalEnv)
    assign('apcire.kmer.carcinoma.stool', apcire.kmer[replace_na(apcire.meta$sample_type == 'stool' & apcire.meta$pheno == 'carcinoma', replace=FALSE),], envir=.GlobalEnv)
    
    assign('apcire.kmer.control.swab', apcire.kmer[replace_na(apcire.meta$sample_type == 'cheek swab' & apcire.meta$pheno == 'control', replace=FALSE),], envir=.GlobalEnv)
    assign('apcire.kmer.adenoma.swab', apcire.kmer[replace_na(apcire.meta$sample_type == 'cheek swab' & apcire.meta$pheno == 'adenoma (polyp)', replace=FALSE),], envir=.GlobalEnv)
    assign('apcire.kmer.carcinoma.swab', apcire.kmer[replace_na(apcire.meta$sample_type == 'cheek swab' & apcire.meta$pheno == 'carcinoma', replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('aus' %in% studies) {
    aus.kmer.kangaroo = aus.kmer[replace_na(aus.meta$host == 'kangaroo', replace=FALSE),]
    aus.kmer.rabbit = aus.kmer[replace_na(aus.meta$host == 'rabbit', replace=FALSE),]
    aus.kmer.sambar_deer = aus.kmer[replace_na(aus.meta$host == 'sambar deer', replace=FALSE),]
    aus.kmer.wombat = aus.kmer[replace_na(aus.meta$host == 'wombat', replace=FALSE),]
  }
  
  if ('cdiff311' %in% studies) {
    assign('cdiff311.kmer.recur', cdiff311.kmer[replace_na(cdiff311.meta$treatment_outcome == 'recurrence', replace=FALSE),], envir=.GlobalEnv)
    assign('cdiff311.kmer.response', cdiff311.kmer[replace_na(cdiff311.meta$treatment_outcome == 'response', replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('cdiff379' %in% studies) {
    assign('cdiff379.kmer.healthy', cdiff379.kmer[cdiff379.meta$pheno == 'healthy control',], envir=.GlobalEnv)
    assign('cdiff379.kmer.cdiff', cdiff379.kmer[cdiff379.meta$pheno == "C. difficile infection",], envir=.GlobalEnv)
    
    assign('cdiff379.kmer.ridi', cdiff379.kmer[cdiff379.meta$pheno == 'C. difficile infection' & cdiff379.meta$treatment == 'ridinilazole',], envir=.GlobalEnv)
    assign('cdiff379.kmer.vanc', cdiff379.kmer[cdiff379.meta$pheno == "C. difficile infection" & cdiff379.meta$treatment == 'vancomycin',], envir=.GlobalEnv)
    
    assign('cdiff379.kmer.day1', cdiff379.kmer[cdiff379.meta$pheno == 'C. difficile infection' & cdiff379.meta$day == "1 (baseline, start of treatment)",], envir=.GlobalEnv)
    assign('cdiff379.kmer.day5', cdiff379.kmer[cdiff379.meta$pheno == 'C. difficile infection' & cdiff379.meta$day == "5",], envir=.GlobalEnv)
    assign('cdiff379.kmer.day10', cdiff379.kmer[cdiff379.meta$pheno == 'C. difficile infection' & cdiff379.meta$day == "10 (end of treatment)",], envir=.GlobalEnv)
    assign('cdiff379.kmer.day25', cdiff379.kmer[cdiff379.meta$pheno == 'C. difficile infection' & cdiff379.meta$day == "25",], envir=.GlobalEnv)
    assign('cdiff379.kmer.day40', cdiff379.kmer[cdiff379.meta$pheno == 'C. difficile infection' & cdiff379.meta$day == "40",], envir=.GlobalEnv)
    assign('cdiff379.kmer.recur', cdiff379.kmer[cdiff379.meta$pheno == 'C. difficile infection' & cdiff379.meta$day == "recurrence",], envir=.GlobalEnv)
  }
  
  if ('cdiff427' %in% studies) {
    assign('cdiff427.kmer.healthy', cdiff427.kmer[cdiff427.meta$pheno == 'healthy',], envir=.GlobalEnv)
    assign('cdiff427.kmer.cdiff', cdiff427.kmer[cdiff427.meta$pheno == 'C. difficile infection',], envir=.GlobalEnv)
  }
  
  if ('celiac385' %in% studies) {
    assign('celiac385.kmer.healthy.stool', celiac385.kmer[celiac385.meta$pheno == 'healthy' & celiac385.meta$sample_type == 'stool',], envir=.GlobalEnv)
    assign('celiac385.kmer.celiac.stool', celiac385.kmer[celiac385.meta$pheno == 'celiac disease' & celiac385.meta$sample_type == 'stool',], envir=.GlobalEnv)
    assign('celiac385.kmer.healthy.duod', celiac385.kmer[celiac385.meta$pheno == 'healthy' & celiac385.meta$sample_type == 'duodenal biopsy',], envir=.GlobalEnv)
    assign('celiac385.kmer.celiac.duod', celiac385.kmer[celiac385.meta$pheno == 'celiac disease' & celiac385.meta$sample_type == 'duodenal biopsy',], envir=.GlobalEnv)
  }
  
  if ('celiac401' %in% studies) {
    assign('celiac401.kmer.healthy', celiac401.kmer[celiac401.meta$pheno == 'healthy' & celiac401.meta$sample_type == 'duodenal biopsy',], envir=.GlobalEnv)
    assign('celiac401.kmer.gluten', celiac401.kmer[celiac401.meta$pheno == 'non-celiac gluten sensitivity' & celiac401.meta$sample_type == 'duodenal biopsy',], envir=.GlobalEnv)
    assign('celiac401.kmer.celiac', celiac401.kmer[celiac401.meta$pheno == 'celiac disease' & celiac401.meta$sample_type == 'duodenal biopsy',], envir=.GlobalEnv)
    assign('celiac401.kmer.baseline', celiac401.kmer[celiac401.meta$diet == 'baseline' & celiac401.meta$sample_type == 'duodenal biopsy',], envir=.GlobalEnv)
    assign('celiac401.kmer.gf', celiac401.kmer[celiac401.meta$diet == '4 weeks on gluten-free diet' & celiac401.meta$sample_type == 'duodenal biopsy',], envir=.GlobalEnv)
  }
  
  if ('chimp' %in% studies) {
    assign('chimp.kmer.healthy', chimp.kmer[chimp.meta$pheno == 'healthy',], envir=.GlobalEnv)
    assign('chimp.kmer.cd', chimp.kmer[chimp.meta$pheno == "Crohn's disease",], envir=.GlobalEnv)
    assign('chimp.kmer.uc', chimp.kmer[chimp.meta$pheno == 'ulcerative colitis',], envir=.GlobalEnv)
  }
  
  if ('costello' %in% studies) {
    assign('costello.kmer.feces', costello.kmer[costello.meta$sample_type == 'feces',], envir=.GlobalEnv)
    assign('costello.kmer.skin', costello.kmer[costello.meta$sample_type == 'sebum',], envir=.GlobalEnv)
    assign('costello.kmer.hair', costello.kmer[costello.meta$sample_type == 'hair',], envir=.GlobalEnv)
    assign('costello.kmer.urine', costello.kmer[costello.meta$sample_type == 'urine',], envir=.GlobalEnv)
    assign('costello.kmer.saliva', costello.kmer[costello.meta$sample_type == 'saliva',], envir=.GlobalEnv)
    assign('costello.kmer.cerumen', costello.kmer[costello.meta$sample_type == 'cerumen',], envir=.GlobalEnv)
    assign('costello.kmer.nostril', costello.kmer[costello.meta$sample_type == 'mucus',], envir=.GlobalEnv)
  }
  
  if ('crc280' %in% studies) {
    assign('crc280.kmer.normal', crc280.kmer[replace_na(crc280.meta$sample_type == 'normal colonic tissue', replace=FALSE),], envir=.GlobalEnv)
    assign('crc280.kmer.adenoma', crc280.kmer[replace_na(crc280.meta$sample_type == 'adenoma colonic tissue', replace=FALSE),], envir=.GlobalEnv)
    assign('crc280.kmer.adenoma_adj', crc280.kmer[replace_na(crc280.meta$sample_type == 'adenoma-adjacent colonic tissue', replace=FALSE),], envir=.GlobalEnv)
    assign('crc280.kmer.carcinoma', crc280.kmer[replace_na(crc280.meta$sample_type == 'carcinoma colonic tissue', replace=FALSE),], envir=.GlobalEnv)
    assign('crc280.kmer.carcinoma_adj', crc280.kmer[replace_na(crc280.meta$sample_type == 'carcinoma-adjacent colonic tissue', replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('crc445' %in% studies) {
    assign('crc445.kmer.mucosa_apx.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (appendix)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_apx.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (appendix)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_apx.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (appendix)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.mucosa_asc.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (ascending colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_asc.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (ascending colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_asc.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (ascending colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.mucosa_des.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (descending colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_des.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (descending colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_des.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (descending colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.mucosa_rect.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (rectum)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_rect.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (rectum)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_rect.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (rectum)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.mucosa_sigmoid.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (sigmoid)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_sigmoid.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (sigmoid)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_sigmoid.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (sigmoid)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.mucosa_trans.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (transverse colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_trans.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (transverse colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.mucosa_trans.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'mucosa (transverse colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.stool_na.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_na.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_na.t',  crc445.kmer[replace_na(crc445.meta$sample_type == 'stool' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.stool_apx.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (appendix)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_apx.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (appendix)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_apx.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (appendix)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.stool_asc.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (ascending colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_asc.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (ascending colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_asc.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (ascending colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.stool_des.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (descending colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_des.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (descending colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_des.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (descending colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.stool_rect.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (rectum)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_rect.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (rectum)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_rect.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (rectum)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.stool_sigmoid.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (sigmoid)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_sigmoid.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (sigmoid)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_sigmoid.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (sigmoid)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.stool_trans.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (transverse colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_trans.a',  crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (transverse colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.stool_trans.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'stool (transverse colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_na.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_na.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_na.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_apx.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (appendix)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_apx.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (appendix)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_apx.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (appendix)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_asc.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (ascending colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_asc.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (ascending colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_asc.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (ascending colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_des.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (descending colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_des.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (descending colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_des.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (descending colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_rect.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (rectum)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_rect.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (rectum)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_rect.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (rectum)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_sigmoid.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (sigmoid)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_sigmoid.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (sigmoid)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_sigmoid.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (sigmoid)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_trans.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (transverse colon)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_trans.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (transverse colon)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_trans.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (transverse colon)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_ileum.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (ileum)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_ileum.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (ileum)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_ileum.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (ileum)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
    
    assign('crc445.kmer.tissue_cecum.n', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (cecum)' & crc445.meta$sample_pheno == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_cecum.a', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (cecum)' & crc445.meta$sample_pheno == 'adjacent', replace=FALSE),], envir=.GlobalEnv)
    assign('crc445.kmer.tissue_cecum.t', crc445.kmer[replace_na(crc445.meta$sample_type == 'tissue (cecum)' & crc445.meta$sample_pheno == 'tumor', replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('crc607' %in% studies) {
    assign('crc607.kmer.normal', crc607.kmer[crc607.meta$sample_type == 'stool' & crc607.meta$pheno == 'normal',], envir=.GlobalEnv)
    assign('crc607.kmer.cancer', crc607.kmer[crc607.meta$sample_type == 'stool' & crc607.meta$pheno == 'cancer',], envir=.GlobalEnv)
  }
  
  if ('crc678' %in% studies) {
    assign('crc678.kmer.normal', crc678.kmer[replace_na(crc678.meta$tissue == 'normal', replace=FALSE),], envir=.GlobalEnv)
    assign('crc678.kmer.cancer', crc678.kmer[replace_na(crc678.meta$tissue == 'tumor', replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('emp148' %in% studies) {
    emp148.kmer.nostril = emp148.kmer[replace_na(emp148.meta$body_habitat == 'nostril', replace=FALSE),]
    emp148.kmer.vagina = emp148.kmer[replace_na(emp148.meta$body_habitat == 'vagina', replace=FALSE),]
    emp148.kmer.skin = emp148.kmer[replace_na(emp148.meta$sample_type == 'sebum', replace=FALSE),]
    emp148.kmer.saliva = emp148.kmer[replace_na(emp148.meta$sample_type == 'saliva', replace=FALSE),]
    emp148.kmer.feces = emp148.kmer[replace_na(emp148.meta$sample_type == 'feces', replace=FALSE),]
    emp148.kmer.anal = emp148.kmer[replace_na(emp148.meta$sample_type == 'anal mucosa', replace=FALSE),]
  }
  
  if ('emp150' %in% studies) {
    emp150.kmer.nofert = emp150.kmer[emp150.meta$fertilizer == 'none',]
    emp150.kmer.fert = emp150.kmer[emp150.meta$fertilizer == 'nitrogen fertilizer',]
  }
  
  if ('emp198' %in% studies) {
    emp198.kmer.feces = emp198.kmer[emp198.meta$sample_type == 'feces',]
    emp198.kmer.saliva = emp198.kmer[emp198.meta$sample_type == 'saliva',]
    emp198.kmer.skin = emp198.kmer[emp198.meta$sample_type == 'skin',]
  }
  
  if ('emp629' %in% studies) {
    emp629.kmer.skin = emp629.kmer[emp629.meta$sample_type == 'sebum' & emp629.meta$host == 'human',]
    emp629.kmer.nostril = emp629.kmer[emp629.meta$sample_type == 'mucus' & emp629.meta$host == 'human',]
  }
  
  if ('hmp' %in% studies) {
    assign('hmp.kmer.healthy', hmp.kmer[replace_na(hmp.meta$ibd == 'healthy', replace=FALSE),], envir=.GlobalEnv)
    assign('hmp.kmer.cd', hmp.kmer[replace_na(hmp.meta$ibd == "Crohn's disease", replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('hmpbody' %in% studies) {
    assign('hmpbody.kmer.stool', hmpbody.kmer[hmpbody.meta$body_habitat == 'feces',], envir=.GlobalEnv)
    assign('hmpbody.kmer.oral_cavity', hmpbody.kmer[hmpbody.meta$body_habitat == 'oral cavity',], envir=.GlobalEnv)
    assign('hmpbody.kmer.saliva', hmpbody.kmer[hmpbody.meta$body_habitat == 'saliva',], envir=.GlobalEnv)
    assign('hmpbody.kmer.skin', hmpbody.kmer[hmpbody.meta$body_habitat == 'skin',], envir=.GlobalEnv)
    assign('hmpbody.kmer.vagina', hmpbody.kmer[hmpbody.meta$body_habitat == 'vagina',], envir=.GlobalEnv)
  }
  
  if ('ibd' %in% studies) {
    assign('ibd.kmer.healthy', ibd.kmer[ibd.meta$pheno == 'healthy',], envir=.GlobalEnv)
    assign('ibd.kmer.cd', ibd.kmer[ibd.meta$pheno == "Crohn's disease",], envir=.GlobalEnv)
    assign('ibd.kmer.uc', ibd.kmer[ibd.meta$pheno == 'ulcerative colitis',], envir=.GlobalEnv)
    assign('ibd.kmer.cc', ibd.kmer[ibd.meta$pheno == 'collagenous colitis',], envir=.GlobalEnv)
    
    ibd.kmer.ccrohn = ibd.kmer[ibd.meta$ibd_subtype == "colonic Crohn's disease",]
    ibd.kmer.icrohn = ibd.kmer[ibd.meta$ibd_subtype %in% c("ileal Crohn's disease (non-resection)", "ileal Crohn's disease (resection)"),]
    ibd.kmer.icrohn.nr = ibd.kmer[ibd.meta$ibd_subtype == "ileal Crohn's disease (non-resection)",]
    ibd.kmer.icrohn.r = ibd.kmer[ibd.meta$ibd_subtype == "ileal Crohn's disease (resection)",]
    assign('ibd.kmer.L1', ibd.kmer[replace_na(ibd.meta$cd_location == "Ileal (L1)", replace=FALSE),], envir=.GlobalEnv)
    assign('ibd.kmer.L2', ibd.kmer[replace_na(ibd.meta$cd_location == "Colonic (L2)", replace=FALSE),], envir=.GlobalEnv)
    assign('ibd.kmer.L3', ibd.kmer[replace_na(ibd.meta$cd_location == "Ileocolonic (L3)", replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('ibd1939' %in% studies) {
    ibd1939.kmer.healthy.stool = ibd1939.kmer[ibd1939.meta$pheno == 'healthy' & ibd1939.meta$sample_type == 'feces',]
    ibd1939.kmer.cd.stool = ibd1939.kmer[ibd1939.meta$pheno == "Crohn's disease" & ibd1939.meta$sample_type == 'feces',]
    ibd1939.kmer.ic.stool = ibd1939.kmer[ibd1939.meta$pheno == 'ileum protect' & ibd1939.meta$sample_type == 'feces',]
    ibd1939.kmer.uc.stool = ibd1939.kmer[ibd1939.meta$pheno == 'ulcerative colitis' & ibd1939.meta$sample_type == 'feces',]
  }
  
  if ('ibs' %in% studies) {
    assign('ibs.kmer.healthy', ibs.kmer[!ibs.meta$ibs,], envir=.GlobalEnv)
    assign('ibs.kmer.ibs', ibs.kmer[ibs.meta$ibs,], envir=.GlobalEnv)
  }
  
  if ('kenya' %in% studies) {
    assign('kenya.kmer.human', kenya.kmer[kenya.meta$host == 'human',], envir=.GlobalEnv)
    assign('kenya.kmer.cow', kenya.kmer[kenya.meta$host == 'cow',], envir=.GlobalEnv)
    assign('kenya.kmer.chicken', kenya.kmer[kenya.meta$host == 'chicken',], envir=.GlobalEnv)
    assign('kenya.kmer.inorganic', kenya.kmer[kenya.meta$sample_type == 'swab',], envir=.GlobalEnv)
  }
  
  if ('knomics' %in% studies) {
    knomics.kmer.pre = knomics.kmer[replace_na(knomics.meta$diet == 'pre-diet', replace=FALSE),]
    knomics.kmer.post = knomics.kmer[replace_na(knomics.meta$diet == 'post-diet', replace=FALSE),]
  }
  
  if ('mice' %in% studies) {
    assign('mice.kmer.humandonor', mice.kmer[mice.meta$host == 'human' & mice.meta$pheno == "healthy (human FMT donor microbiome)",], envir=.GlobalEnv)
    assign('mice.kmer.humancdiff', mice.kmer[mice.meta$host == 'human' & mice.meta$pheno == 'cdiff',], envir=.GlobalEnv)
    assign('mice.kmer.humanfmtcdiff', mice.kmer[mice.meta$host == 'human' & mice.meta$pheno == 'cdiff (human FMT recipient)',], envir=.GlobalEnv)
    
    assign('mice.kmer.micedonor', mice.kmer[mice.meta$host == 'mouse' & mice.meta$pheno == "healthy (human FMT donor microbiome)",], envir=.GlobalEnv)
    assign('mice.kmer.micecdiff', mice.kmer[mice.meta$host == 'mouse' & mice.meta$pheno == 'cdiff',], envir=.GlobalEnv)
    assign('mice.kmer.micefmtcdiff', mice.kmer[mice.meta$host == 'mouse' & mice.meta$pheno == 'cdiff (human FMT recipient)',], envir=.GlobalEnv)
    
    assign('mice.kmer.health', mice.kmer[mice.meta$host == 'human' & mice.meta$pheno == 'healthy',], envir=.GlobalEnv)
    assign('mice.kmer.uc', mice.kmer[mice.meta$pheno == 'ulcerative colitis',], envir=.GlobalEnv)
    assign('mice.kmer.cd', mice.kmer[mice.meta$pheno == "Crohn's disease",], envir=.GlobalEnv)
  }
  
  if ('pascal' %in% studies) {
    assign('pascal.kmer.healthy', pascal.kmer[pascal.meta$ibd == 'healthy',], envir=.GlobalEnv)
    assign('pascal.kmer.cd', pascal.kmer[pascal.meta$ibd == "Crohn's disease",], envir=.GlobalEnv)
    assign('pascal.kmer.uc', pascal.kmer[pascal.meta$ibd == "ulcerative colitis",], envir=.GlobalEnv)
    
    assign('pascal.meta.healthy', pascal.meta[pascal.meta$ibd == 'healthy',], envir=.GlobalEnv)
    assign('pascal.meta.cd', pascal.meta[pascal.meta$ibd == "Crohn's disease",], envir=.GlobalEnv)
    assign('pascal.meta.uc', pascal.meta[pascal.meta$ibd == "ulcerative colitis",], envir=.GlobalEnv)
  }
  
  if ('protect' %in% studies) {
    assign('protect.kmer.inactive', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$pucai_c4_wkall == 'inactive', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.mild', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$pucai_c4_wkall == 'mild', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.moderate', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$pucai_c4_wkall == 'moderate', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.severe', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$pucai_c4_wkall == 'severe', replace=FALSE),], envir=.GlobalEnv)
    
    assign('protect.kmer.nr4', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$collection_week == 4 & protect.meta$remission_wk4 == 'No remission', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.nr12', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$collection_week == 12 & protect.meta$remission_wk12 == 'No remission', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.nr52', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$collection_week == 52 & protect.meta$remission_wk52 == 'No remission', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.nr', rbind(protect.kmer.nr4, protect.kmer.nr12, protect.kmer.nr52), envir=.GlobalEnv)
    
    assign('protect.kmer.r4', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$collection_week == 4 & protect.meta$remission_wk4 == 'Remission', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.r12', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$collection_week == 12 & protect.meta$remission_wk12 == 'Remission', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.r52', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$collection_week == 52 & protect.meta$remission_wk52 == 'Remission', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.r', rbind(protect.kmer.r4, protect.kmer.r12, protect.kmer.r52), envir=.GlobalEnv)
    
    assign('protect.kmer.epu', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$montreal_ord == 'Extensive or Pancolitis or Unasessable', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.left', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$montreal_ord == 'Left-sided colitis', replace=FALSE),], envir=.GlobalEnv)
    assign('protect.kmer.proc', protect.kmer[replace_na(protect.meta$sample_type == 'feces' & protect.meta$montreal_ord == 'Proctosigmoiditis', replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('risk' %in% studies) {
    assign('risk.kmer.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.cecum', risk.kmer[replace_na(risk.meta$sample_type == 'cecum', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.colon', risk.kmer[replace_na(risk.meta$sample_type %in% c('ascending colon', 'descending colon', 'sigmoid colon', 'transverse colon'), replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.healthy.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == 'healthy' & risk.meta$pheno_current == 'healthy', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.uc.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == 'ulcerative colitis' & risk.meta$pheno_current == 'ulcerative colitis', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.cd.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease", replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.healthy.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == 'healthy' & risk.meta$pheno_current == 'healthy', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.uc.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == 'ulcerative colitis' & risk.meta$pheno_current == 'ulcerative colitis', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.cd.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease", replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.healthy.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == 'healthy' & risk.meta$pheno_current == 'healthy', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.uc.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == 'ulcerative colitis' & risk.meta$pheno_current == 'ulcerative colitis', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.cd.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease", replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.healthy.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_current == 'healthy', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L1.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L1', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L2.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L2', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L3.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L3', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.healthy.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_current == 'healthy', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L1.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L1', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L2.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L2', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L3.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L3', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.healthy.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_current == 'healthy', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L1.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L1', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L2.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L2', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.L3.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$location == 'L3', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.B1.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B1', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B2.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B2', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B3.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B3', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B2B3.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B2+B3', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.B1.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B1', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B2.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B2', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B3.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B3', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B2B3.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B2+B3', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.B1.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B1', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B2.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B2', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B3.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B3', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.B2B3.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$pheno_enrollment == "Crohn's disease" & risk.meta$pheno_current == "Crohn's disease" & risk.meta$subtype_current == 'B2+B3', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.afro.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$race == 'African American', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.white.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$race == 'Caucasian', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.asian.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$race == 'South/East Asian', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.arab.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$race == 'West Asian/Arab', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.other.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$race == 'Other', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.afro.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$race == 'African American', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.white.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$race == 'Caucasian', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.asian.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$race == 'South/East Asian', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.arab.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$race == 'West Asian/Arab', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.other.rectum', risk.kmer[replace_na(risk.meta$sample_type == 'rectum' & risk.meta$race == 'Other', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.afro.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$race == 'African American', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.white.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$race == 'Caucasian', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.asian.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$race == 'South/East Asian', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.arab.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$race == 'West Asian/Arab', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.other.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$race == 'Other', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.noabdpain.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$history_of_abdominal_pain == 'false', replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.abdpain.ileum', risk.kmer[replace_na(risk.meta$sample_type == 'ileum' & risk.meta$history_of_abdominal_pain == 'true', replace=FALSE),], envir=.GlobalEnv)
    
    assign('risk.kmer.mild.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$severity < 0.5, replace=FALSE),], envir=.GlobalEnv)
    assign('risk.kmer.severe.stool', risk.kmer[replace_na(risk.meta$sample_type == 'stool' & risk.meta$severity >= 0.5, replace=FALSE),], envir=.GlobalEnv)
  }
  
  if ('stinki' %in% studies) {
    #stinki.kmer = cbind(stinki.kmer, stinki.meta$calprotectin)
    #stinki.kmer = stinki.kmer[!is.na(stinki.meta$calprotectin),]
    #stinki.meta = stinki.meta[!is.na(stinki.meta$calprotectin),]
    #colnames(stinki.kmer)[4097] = 'calprotectin'
    
    stinki.kmer.chrono = stinki.kmer[order(stinki.meta$collection_date),]
    stinki.meta.chrono = stinki.meta[order(stinki.meta$collection_date),]
    
    stinki.kmer.pretreat.resp = stinki.kmer[substr(stinki.meta$sample, nchar(stinki.meta$sample) - 1, nchar(stinki.meta$sample)) %in% c('00', '01') & stinki.meta$responder %in% c('UC responder', 'CD responder'),]
    stinki.kmer.pretreat.non = stinki.kmer[substr(stinki.meta$sample, nchar(stinki.meta$sample) - 1, nchar(stinki.meta$sample)) %in% c('00', '01') & stinki.meta$responder %in% c('UC non-responder', 'CD non-responder'),]
    
    stinki.kmer.pretreat = matrix(NA, nrow=0, ncol=4096)
    keep.indices = c()
    subjects = c()
    for (i in 1:nrow(stinki.meta.chrono)) {
      if (!(stinki.meta.chrono$subject[i] %in% subjects)) {
        stinki.kmer.pretreat = rbind(stinki.kmer.pretreat, stinki.kmer.chrono[i])
        subjects = c(subjects, stinki.meta.chrono$subject[i])
        keep.indices = c(keep.indices, i)
      }
    }
    stinki.meta.pretreat = stinki.meta.chrono[keep.indices,]
    
    stinki.kmer.resp = stinki.kmer[stinki.meta$responder %in% c('UC responder', 'CD responder'),]
    stinki.kmer.non = stinki.kmer[stinki.meta$responder %in% c('UC non-responder', 'CD non-responder'),]
  }
  
  if ('umich' %in% studies) {
    umich.kmer.normal = umich.kmer[umich.meta$pheno == 'normal',]
    umich.kmer.carcinoma = umich.kmer[umich.meta$pheno == 'carcinoma',]
  }
  
  if ('ubiome' %in% studies) {
    ubiome.kmer.stool = ubiome.kmer[ubiome.meta$sample_type == 'feces',]
  }
}