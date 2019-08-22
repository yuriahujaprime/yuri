### Script for loading kmer and metadata files for specific studies, optional kmer argument, e.g. k=6 (default) means 6mers.
### The matrix study.kmer and data frame study.meta are attached globally to the R session.
### The metadata data frame is trimmed as necessary so that its rows match the kmer matrix rows.
### Eric Proffitt

"amerigut------10317"
"anorex--------PRJEB11199"
"apcire--------apcireland"
"aus-----------PRJEB14739"
"cdiff307------PRJNA307992"
"cdiff311------PRJNA311224"
"cdiff379------PRJNA379979"
"cdiff419------PRJNA419097"
"cdiff427------PRJNA427597"
"celiac385-----PRJNA385740"
"celiac401-----PRJNA401920"
"chimp---------PRJNA82109"
"costello------449"
"crc258--------PRJNA258534"
"crc280--------PRJNA280026"
"crc284--------PRJNA284355"
"crc445--------PRJNA445346"
"crc607--------PRJEB6070"
"crc678--------PRJNA67873"
"emp147--------PRJEB14782"
"emp148--------PRJEB14801"
"emp150--------PRJEB15061"
"emp198--------PRJEB19825"
"emp219--------PRJEB21935"
"emp629--------PRJEB6292"
"finland_india-PRJNA327770"
"hazda---------PRJNA392012"
"hmp-----------2538"
"hmpbody-------1928"
"hmpbodyv4-----1928V4"
"ibd-----------1629"
"ibd1939-------1939"
"ibs-----------PRJNA268708"
"jungle--------PRJEB4489"
"kenya---------PRJNA345054"
"knomics-------PRJEB16344"
"mb1000--------EGAS00001002702"
"mice----------PRJNA413199"
"nafld---------PRJNA246121"
"pascal--------PRJNA422193"
"protect-------PRJNA436359"
"risk----------PRJNA237362"
"schizo--------PRJEB26004"
"stinki--------PRJNA324825"
"subra---------PRJNA237362"
"ubiome--------PRJEB20022"
"umich---------UMich2014"
"veg-----------PRJEB18693"

load_kmer = function(studies, k=6, readlen=150) {
  
  study_numbers = c('PRJEB11419', 'PRJEB11199', 'apcireland', 'PRJEB14739', 'PRJNA307992', 'PRJNA311224', 'PRJNA379979', 'PRJNA419097',
                    'PRJNA427597', 'PRJNA385740', 'PRJNA401920', 'PRJNA82109', '449', 'PRJNA258534', 'PRJNA280026', 'PRJNA284355', 'PRJNA445346', 'PRJEB6070', 'PRJNA67873', 'PRJEB14782', 'PRJEB14801', 'PRJEB15061', 'PRJEB19825', 'PRJEB21935', 'PRJEB6292',
                    'PRJNA327770', 'PRJNA392012', '2538', '1928', '1928V4', '1629', '1939', 'PRJNA268708', 'PRJEB4489', 'PRJNA345054',
                    'PRJEB16344', 'EGAS00001002702', 'PRJNA413199', 'PRJNA246121', 'PRJNA422193', 'PRJNA436359',
                    'PRJNA237362', 'PRJEB26004', 'PRJNA324825', 'PRJNA237362', 'PRJEB20022', 'UMich2014', 'PRJEB18693')
  
  names(study_numbers) = c('amerigut', 'anorex', 'apcire', 'aus', 'cdiff307', 'cdiff311', 'cdiff379', 'cdiff419',
                           'cdiff427', 'celiac385', 'celiac401', 'chimp', 'costello', 'crc258', 'crc280', 'crc284', 'crc445', 'crc607', 'crc678', 'emp147', 'emp148', 'emp150', 'emp198', 'emp219', 'emp629',
                           'finland_india', 'hazda', 'hmp', 'hmpbody', 'hmpbodyv4', 'ibd', 'ibd1939', 'ibs', 'jungle', 'kenya',
                           'knomics', 'mb1000', 'mice', 'nafld', 'pascal', 'protect', 'risk', 'schizo', 'stinki',
                           'subra', 'ubiome', 'umich', 'veg')
  
  if ('amerigut' %in% studies) {
    study_number = study_numbers['amerigut']
    amerigut.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    amerigut.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="")
    amerigut.meta$sample = as.character(amerigut.meta$sample)
  }
  
  if ('anorex' %in% studies) {
    study_number = study_numbers['anorex']
    anorex.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    anorex.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'factor', 'character'))
  }
  
  if ('apcire' %in% studies) {
    study_number = study_numbers['apcire']
    apcire.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    apcire.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', 'numeric', rep('factor', 3), 'integer', 'logical', rep('factor', 4), 'logical', rep('numeric', 2), 'logical', 'integer', 'factor', 'character'))
  }
  
  if ('aus' %in% studies) {
    study_number = study_numbers['aus']
    aus.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    aus.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                          colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'character'))
  }
  
  if ('cdiff307' %in% studies) {
    study_number = study_numbers['cdiff307']
    cdiff307.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    cdiff307.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                               colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), rep('logical', 3), rep('factor', 2), 'logical', rep('factor', 3), 'logical', 'character'))
  }
  
  if ('cdiff311' %in% studies) {
    study_number = study_numbers['cdiff311']
    cdiff311.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    cdiff311.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                               colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 7), 'factor', 'character'))
  }
  
  if ('cdiff379' %in% studies) {
    study_number = study_numbers['cdiff379']
    cdiff379.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    cdiff379.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                               colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), 'integer', 'numeric', 'factor', 'logical', rep('factor', 4), 'character'))
  }
  
  if ('cdiff419' %in% studies) {
    study_number = study_numbers['cdiff419']
    cdiff419.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    cdiff419.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                               colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', rep('factor', 2), 'integer', 'factor', 'integer', rep('factor', 2), rep('logical', 3), rep('factor', 6), 'character'))
  }
  
  if ('cdiff427' %in% studies) {
    study_number = study_numbers['cdiff427']
    cdiff427.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    cdiff427.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                               colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 5), 'character'))
  }
  
  if ('celiac385' %in% studies) {
    study_number = study_numbers['celiac385']
    celiac385.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_trimmed_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    celiac385.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                                colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', 'factor', 'numeric', rep('factor', 2), 'character'))
  }
  
  if ('celiac401' %in% studies) {
    study_number = study_numbers['celiac401']
    celiac401.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_trimmed_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    celiac401.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                                colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', 'numeric', rep('factor', 2), 'character'))
  }
  
  if ('chimp' %in% studies) {
    study_number = study_numbers['chimp']
    chimp.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    chimp.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                            colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), 'integer', rep('factor', 2), rep('logical', 5), 'factor', rep('logical', 8), 'integer', rep('logical', 12), 'character', rep('factor', 3), rep('logical', 6), rep('factor', 4), 'logical', 'factor', 'integer', 'numeric', rep('factor', 11), 'factor', 'numeric', 'factor', 'character'))
  }
  
  if ('costello' %in% studies) {
    study_number = study_numbers['costello']
    costello.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    costello.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                               colClasses=c(rep('character', 3), rep('integer', 2), rep('factor', 3), rep('integer', 2), rep('factor', 7), 'logical', 'factor', 'character'))
  }
  
  if ('crc258' %in% studies) {
    study_number = study_numbers['crc258']
    crc258.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    crc258.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), 'character'))
  }
  
  if ('crc280' %in% studies) {
    study_number = study_numbers['crc280']
    crc280.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    crc280.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), 'character'))
  }
  
  if ('crc284' %in% studies) {
    study_number = study_numbers['crc284']
    crc284.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    crc284.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', rep('factor', 3), 'integer', 'factor', 'character'))
  }
  
  if ('crc445' %in% studies) {
    study_number = study_numbers['crc445']
    crc445.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    crc445.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), 'numeric', rep('factor', 3), 'logical', rep('factor', 4), 'integer', 'character'))
  }
  
  if ('crc607' %in% studies) {
    study_number = study_numbers['crc607']
    crc607.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    crc607.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', 'numeric', rep('factor', 3), 'integer', 'factor', 'character'))
  }
  
  if ('crc678' %in% studies) {
    study_number = study_numbers['crc678']
    crc678.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    crc678.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), 'character'))
  }
  
  if ('emp147' %in% studies) {
    study_number = study_numbers['emp147']
    emp147.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    emp147.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', 'numeric', 'factor', 'character'))
  }
  
  if ('emp148' %in% studies) {
    study_number = study_numbers['emp148']
    emp148.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_100_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    emp148.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'numeric', rep('factor', 2), 'numeric', rep('factor', 6), 'character'))
  }
  
  if ('emp150' %in% studies) {
    study_number = study_numbers['emp150']
    emp150.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    emp150.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 3), rep('numeric', 4), 'factor', 'numeric', rep('factor', 2), 'integer', rep('factor', 3), 'character'))
  }
  
  if ('emp198' %in% studies) {
    study_number = study_numbers['emp198']
    emp198.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    emp198.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'numeric', rep('factor', 4), rep('logical', 3), 'character'))
  }
  
  if ('emp219' %in% studies) {
    study_number = study_numbers['emp219']
    emp219.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    emp219.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 5), rep('numeric', 8), rep('integer', 4), 'character'))
  }
  
  if ('emp629' %in% studies) {
    study_number = study_numbers['emp629']
    emp629.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    emp629.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), rep('integer', 2), rep('numeric', 2), rep('factor', 3), 'character'))
  }
  
  if ('finland_india' %in% studies) {
    study_number = study_numbers['finland_india']
    finland_india.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    finland_india.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                                    colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 7), 'character'))
  }
  
  if ('hazda' %in% studies) {
    study_number = study_numbers['hazda']
    hazda.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    hazda.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                            colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'character'))
  }
  
  if ('hmp' %in% studies) {
    study_number = study_numbers['hmp']
    hmp.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    hmp.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                          colClasses=c(rep('character', 4), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', 'factor', 'logical', 'factor',
                                       rep('logical', 2), 'numeric', 'factor', rep('integer', 2), 'logical', rep('factor', 3), 'logical', 'factor',
                                       'integer', 'factor', rep('logical', 2), rep('factor', 7), 'logical', 'factor', rep('logical', 4), 'numeric',
                                       'factor', rep('logical', 2), rep('factor', 2), rep('logical', 2), rep('factor', 3), rep('logical', 15), 'factor',
                                       'logical', 'factor', 'logical', 'factor', 'logical', 'integer', 'integer', rep('factor', 3), 'logical', 'factor',
                                       'logical', 'factor', rep('logical', 2), 'factor', 'logical', 'factor', rep('logical', 2), 'factor', 'logical', 'numeric', 'factor', 'character'))
  }
  
  if ('hmpbody' %in% studies) {
    study_number = study_numbers['hmpbody']
    hmpbody.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    hmpbody.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                              colClasses=c(rep('character', 3), rep('integer', 2), rep('factor', 3), 'integer', rep('factor', 3), 'character'))
  }
  
  if ('hmpbodyv4stool' %in% studies) {
    study_number = study_numbers['hmpbodyv4stool']
    hmpbodyv4stool.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    hmpbodyv4stool.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                                     colClasses=c(rep('character', 3), rep('integer', 2), rep('factor', 3), 'integer', rep('factor', 3), 'character'))
  }
  
  if ('ibd' %in% studies) {
    study_number = study_numbers['ibd']
    ibd.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer_50000.tsv"), sep='\t', header=TRUE)
    ibd.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                          colClasses=c(rep('character', 2), 'integer', 'character', 'integer', rep('factor', 3), 'integer', 'factor', 'character', 'numeric', rep('factor', 4), 'logical', 'numeric', 'character'))
  }
  
  if ('ibd1939' %in% studies) {
    study_number = study_numbers['ibd1939']
    ibd1939.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    ibd1939.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                              colClasses=c(rep('character', 3), rep('integer', 2), rep('factor', 3), 'integer', 'factor', 'numeric', rep('factor', 3), 'integer', 'factor', 'logical', 'factor', rep('logical', 4), 'factor', rep('logical', 2), 'character'))
  }
  
  if ('ibs' %in% studies) {
    study_number = study_numbers['ibs']
    ibs.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    ibs.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                          colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'logical', 'factor', rep('character', 3)))
  }
  
  if ('jungle' %in% studies) {
    study_number = study_numbers['jungle']
    jungle.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    jungle.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), 'character'))
  }
  
  if ('kenya' %in% studies) {
    study_number = study_numbers['kenya']
    kenya.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE) 
    kenya.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                            colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), 'numeric', 'factor', 'character'))
  }
  
  if ('knomics' %in% studies) {
    study_number = study_numbers['knomics']
    knomics.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE) 
    knomics.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                              colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 2), rep('logical', 3), 'factor', rep('logical', 3), 'numeric', rep('integer', 2), 'numeric', rep('integer', 48), 'factor', 'character'))
  }
  
  if ('mb1000' %in% studies) {
    study_number = study_numbers['mb1000']
    mb1000.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE) 
    mb1000.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', rep('numeric', 2), 'factor', rep('integer', 2), 'factor', rep('integer', 2), rep('factor', 4), rep('logical', 7), 'character'))
  }
  
  if ('mice' %in% studies) {
    study_number = study_numbers['mice']
    mice.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    mice.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                           colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 3), 'character'))
  }
  
  if ('nafld' %in% studies) {
    study_number = study_numbers['nafld']
    nafld.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    nafld.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                            colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'integer', 'factor', 'character'))
  }
  
  if ('pascal' %in% studies) {
    study_number = study_numbers['pascal']
    #pascal.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_trimmed_", k, "-mers_50000.tsv"), sep='\t', header=TRUE)
    pascal.kmer = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_50000.tsv"), sep='\t', header=TRUE)
    pascal.meta = read.table(paste0(getwd(), "/metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'character', 'character', 'factor', 'character'))
  }
  
  if ('protect' %in% studies) {
    study_number = study_numbers['protect']
    protect.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    protect.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                              colClasses=c(rep('character', 3), 'integer', rep('factor', 3), rep('integer', 3), rep('factor', 6), 'numeric', rep('factor', 2), rep('numeric', 6), 'factor', 'character'))
  }
  
  if ('risk' %in% studies) {
    study_number = study_numbers['risk']
    risk.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    risk.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                           colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 7), rep('logical', 2), rep('factor', 2), rep('logical', 3), rep('factor', 3), rep('character', 2)))
  }
  
  if ('schizo' %in% studies) {
    study_number = study_numbers['schizo']
    schizo.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    schizo.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), rep('integer', 2), rep('factor', 3), 'integer', 'factor', 'character'))
  }
  
  if ('stinki' %in% studies) {
    study_number = study_numbers['stinki']
    stinki.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    stinki.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 4), 'integer', 'numeric', 'character'))
  }
  
  if ('subra' %in% studies) {
    study_number = study_numbers['subra']
    subra.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    subra.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                            colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 4), 'character'))
  }
  
  if ('ubiome' %in% studies) {
    study_number = study_numbers['ubiome']
    ubiome.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    ubiome.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                             colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', 'factor', 'character'))
  }
  
  if ('umich' %in% studies) {
    study_number = study_numbers['umich']
    umich.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_chimera_removed_150_", k, "-mers_5000.tsv"), sep='\t', header=TRUE)
    umich.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                            colClasses=c(rep('character', 4), 'integer', rep('factor', 3), 'integer', rep('factor', 4), rep('integer', 3), rep('factor', 2), 'character', 'factor', 'character'))
  }
  
  if ('veg' %in% studies) {
    study_number = study_numbers['veg']
    veg.kmer = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/kmer/", study_number, "_", k, "mer.tsv"), sep='\t', header=TRUE)
    veg.meta = read.table(paste0(getwd(), "../metadata_curation/", study_number, "/final/", study_number, "_metadata.tsv"), sep='\t', header=TRUE, quote="",
                          colClasses=c(rep('character', 3), 'integer', rep('factor', 3), 'integer', rep('factor', 3), rep('numeric', 78), 'character'))
  }
  
  for (study in intersect(studies, names(study_numbers))) {
    study.kmer = eval(parse(text=paste0(study, '.kmer')))
    study.meta = eval(parse(text=paste0(study, '.meta')))
    
    study.common_samples = intersect(study.kmer$sample, study.meta$sample)
    study.kmer = as.matrix(study.kmer[,2:(4^k + 1)])
    study.meta = study.meta[study.meta$sample %in% study.common_samples,]
    
    assign(paste0(study, '.kmer'), study.kmer, envir=.GlobalEnv)
    assign(paste0(study, '.meta'), study.meta, envir=.GlobalEnv)
  }
}
