
human_maps <- c('duke', 'emory', 'mt_sinai_human', 'innsbruck', 'oxford',
								'fda', 'geneva', 'amc')
hamster_maps <- c('madison_pooled', 'madison_unpooled', 'galveston',
								  'emc_prnt', 'charite', 'madison_frnt', 'emc_calu', 'emc_vero')
mouse_maps <- c('st_louis', 'maryland')

pretty_labels <- c(
	# human maps
	`duke` = 'Duke',
	`emory` = 'Emory',
	`mt_sinai_human` = 'Mt. Sinai',
	`innsbruck` = 'Innsbruck',
	`oxford` = 'Oxford',
	`fda` = 'FDA',
	`geneva` = 'Geneva',
	`amc` = 'AMC',
	# hamster maps
	`madison_pooled` = 'Madison (pooled)',
	`madison_unpooled` = 'Madison (unpooled)',
	`galveston` = 'Galveston',
	`emc_prnt` = 'EMC (PRNT)',
	`charite` = 'Charité',
	`madison_frnt` = 'Madison (FRNT)',
	`emc_calu` = 'EMC (Calu-3)',
	`emc_vero` = 'EMC (VeroE6)',
	# mouse maps
	`st_louis` = 'WUSTL',
	`maryland` = 'Maryland (pooled)'
	)


map_colors <- c(
	# human maps
	`duke` = '#03569b',
	`emory` = '#F76A05',
	`mt_sinai_human` = '#ffc808',
	`innsbruck` = '#a2b324',
	`oxford` = '#808080',
	`fda` = '#742f32',
	`geneva` = '#333333',
	`amc` = 'magenta',
	# hamster maps
	`madison_pooled` = 'sandybrown',
	`madison_unpooled` = '#a020f0',
	`galveston` = '#e9a390',
	`emc_prnt` = '#93EDC3',
	`charite` = 'lightgoldenrod4',
	`madison_frnt` = '#75ada9',
	`emc_calu` = 'dodgerblue',
	`emc_vero` = 'firebrick',
	# mouse maps
	`st_louis` = '#ED93BD',
	`maryland` = '#9B84AD'
	)


animal_colors <- c(
	# human maps
	`duke` = '#d95f02',
	`emory` = '#d95f02',
	`mt_sinai_human` = '#d95f02',
	`innsbruck` = '#d95f02',
	`oxford` = '#d95f02',
	`fda` = '#d95f02',
	`geneva` = '#d95f02',
	`amc` = '#d95f02',
	# hamster maps
	`madison_pooled` = '#235c3f',
	`madison_unpooled` = '#235c3f',
	`galveston` = '#235c3f',
	`emc_prnt` = '#235c3f',
	`charite` = '#235c3f',
	`madison_frnt` = '#235c3f',
	`emc_calu` = '#235c3f',
	`emc_vero` = '#235c3f',
	# mouse maps
	`st_louis` = '#6383f2',
	`maryland` = '#6383f2'
	)


map_shapes <- c(
	# human maps
	`duke` = 16,
	`emory` = 16,
	`mt_sinai_human` = 16,
	`innsbruck` = 16,
	`oxford` = 16,
	`fda` = 16,
	`geneva` = 16,
	`amc` = 16,
	# hamster maps
	`madison_pooled` = 16,
	`madison_unpooled` = 16,
	`galveston` = 16,
	`emc_prnt` = 16,
	`charite` = 16,
	`madison_frnt` = 16,
	`emc_calu` = 16,
	`emc_vero` = 16,
	# mouse maps
	`st_louis` = 16,
	`maryland` = 16
	)

map_shapes_different <- c(
	# human maps
	`duke` = 15,
	`emory` = 16,
	`mt_sinai_human` = 17,
	`innsbruck` = 19,
	`oxford` = 18,
	`fda` = 14,
	`geneva` = 25,
	`amc` = 21,
	# hamster maps
	`madison_pooled` = 0,
	`madison_unpooled` = 1,
	`galveston` = 2,
	`emc_prnt` = 5,
	`charite` = 6,
	`madison_frnt` = 7,
	`emc_calu` = 10,
	`emc_vero` = 9,
	# mouse maps
	`st_louis` = 4,
	`maryland` = 8
	)

assay_colors <- c(
	# human maps
	`duke` = '#d95f02',
	`emory` = '#003366',
	`mt_sinai_human` = '#800080',
	`innsbruck` = '#003366',
	`oxford` = '#003366',
	`fda` = '#d95f02',
	`geneva` = '#336600',
	`amc` = '#d95f02',
	# hamster maps
	`madison_pooled` = '#008080',
	`madison_unpooled` = '#008080',
	`galveston` = '#003366',
	`emc_prnt` = '#336600',
	`charite` = '#336600',
	`madison_frnt` = '#003366',
	`emc_calu` = '#CCCC00',
	`emc_vero` = '#CCCC00',
	# mouse maps
	`st_louis` = '#003366',
	`maryland` = '#008080'
	)


sr_colors <- c(
	`mRNA-1273` = 'grey',
	`D614G convalescent` = '#333333',
	`B.1.1.7 convalescent` = '#637939',
	`B.1.351 convalescent` = '#e7ba52',
	`P.1 convalescent` = '#7b4173',
	`B.1.617.2 convalescent` = '#d18652',
	`BA.1 convalescent` = '#EF3737'
	)


animal_dataset_order <- c('duke', 'emory', 'fda', 'innsbruck',
												  'oxford', 'mt_sinai_human', 'amc', 'geneva',
                          'emc_prnt', 'charite', 'emc_vero', 'emc_calu',
                          'galveston', 'madison_frnt', 'madison_pooled', 'madison_unpooled',
                          'st_louis', 'maryland')

assay_dataset_order <- c('emory', 'innsbruck', 'oxford', 'galveston',
												 'st_louis', 'madison_frnt',
                         'duke', 'fda', 'amc',
                         'emc_vero', 'emc_calu',
                         'geneva', 'charite', 'emc_prnt',
                         'mt_sinai_human',
                         'madison_pooled', 'madison_unpooled', 'maryland')

duplicated_ags <- c('B.1.351', 'B.1.617.2', 'D614G', 'B.1.1.7', 'BA.1',
									  'P.1', '614D', 'B.1.621', 'B.1.1.7+E484K', 'BA.2',
									  'B.1.617.1', 'B.1.526+E484K', 'C.37', 'B.1.429', 'P.2',
									  'B.1.526+S477N', 'BA.5', 'R.1', 'B.1.617.2+K417N', 'C.36.3',
										'B.1.526', 'BA.1.1', 'BA.2.12.1')


assay <- data.frame(
  row.names = c('maryland', 'duke', 'emory',
                'mt_sinai_human', 'galveston', 'madison_pooled',
                'madison_unpooled', 'st_louis', 'oxford', 'emc_prnt',
                'innsbruck', 'madison_frnt', 'charite', 'fda',
                'geneva', 'amc', 'emc_calu', 'emc_vero'),
  val = c('CPE', 'LV-PV-neut', 'FRNT', 'Microneut', 'FRNT', 'CPE', 'CPE', 'FRNT', 'FRNT',
          'PRNT', 'FRNT', 'FRNT', 'PRNT', 'LV-PV-neut',
          'PRNT', 'LV-PV-neut', 'VSV-PV-neut', 'VSV-PV-neut')
)

model_animal <- data.frame(
  row.names = c('maryland', 'duke', 'emory',
                'mt_sinai_human', 'galveston', 'madison_pooled',
                'madison_unpooled', 'st_louis', 'oxford', 'emc_prnt',
                'innsbruck', 'madison_frnt', 'charite', 'fda',
                'geneva', 'amc', 'emc_calu', 'emc_vero'),
  val = c('Mouse', 'Human', 'Human', 'Human', 'Hamster', 'Hamster', 'Hamster', 'Mouse', 'Human',
          'Hamster', 'Human', 'Hamster', 'Hamster', 'Human',
          'Human', 'Human', 'Hamster', 'Hamster')
)

cell_type <- data.frame(
  row.names = c('maryland', 'duke', 'emory',
                'mt_sinai_human', 'galveston', 'madison_pooled',
                'madison_unpooled', 'st_louis', 'oxford', 'emc_prnt',
                'innsbruck', 'madison_frnt', 'charite', 'fda',
                'geneva', 'amc', 'emc_calu', 'emc_vero'),
  val = c('VeroE6-TMPRSS2', 'HEK293T-ACE2', 'VeroE6-TMPRSS2',
  	      'VeroE6', 'VeroE6', 'VeroE6-TMPRSS2',
  	      'VeroE6-TMPRSS2', 'VeroE6-TMPRSS2', 'VeroE6', 'calu-3',
  	      'Vero-TMPRSS2-ACE2', 'VeroE6-TMPRSS2', 'VeroE6', 'HEK293T-TMPRSS2-ACE2',
  	      'VeroE6', 'HEK293T-ACE2', 'calu-3', 'VeroE6')
)

cell_type_2 <- data.frame(
  row.names = c('maryland', 'duke', 'emory',
                'mt_sinai_human', 'galveston', 'madison_pooled',
                'madison_unpooled', 'st_louis', 'oxford', 'emc_prnt',
                'innsbruck', 'madison_frnt', 'charite', 'fda',
                'geneva', 'amc', 'emc_calu', 'emc_vero'),
  val = c('TMPRSS2', '-', 'TMPRSS2',
  	      '-', '-', 'TMPRSS2',
  	      'TMPRSS2', 'TMPRSS2', '-', '-',
  	      'TMPRSS2', 'TMPRSS2', '-', 'TMPRSS2',
  	      '-', '-', '-', '-')
)

dataset_labels_sub <- data.frame(
  row.names = c('fda', 'duke', 'innsbruck', 'amc', 'emory', 'geneva',
                'charite', 'emc_prnt', 'madison_frnt', 'emc_vero', 'emc_calu',
                'oxford', 'mt_sinai_human', 'st_louis', 'galveston',
                'madison_pooled', 'madison_unpooled', 'maryland'),
  val = c('Human, LV-PV-neut', 'Human, LV-PV-neut', 'Human, FRNT', 'Human, LV-PV-neut', 'Human, FRNT', 'Human, PRNT',
          'Hamster, PRNT', 'Hamster, PRNT', 'Hamster, FRNT', 'Hamster, VSV-PV-neut', 'Hamster, VSV-PV-neut',
          'Human, FRNT', 'Human, Microneut', 'Mouse, FRNT', 'Hamster, FRNT',
          'Hamster, CPE', 'Hamster, CPE', 'Mouse, CPE')
)

ag_colors <- list(
	'D614G' = '#393b79',
	'614D' = '#393b79',
	'B.1.1.7' = '#637939',
	'B.1.526' = '#b8bf58',
	'B.1.526+S477N' = '#72b370',
	'B.1.429' = '#9B9FD9',
	'C.36.3' = '#ed928c',
	'C.37' = '#79c9c9',
	'B.1.617.2' = '#d18652',
  'B.1.617.2+K417N' = '#d18652',
  'B.1.617.1' = '#ad494a',
  'B.1.351' = '#e7ba52',
  'P.1' = '#7b4173',
  'B.1.1.7+E484K' = '#637939',
  'B.1.621' = '#0096ad',
  'B.1.526+E484K' = '#9ab370',
  'R.1' = '#f2f53b',
  'P.2' = '#be7cf8',
  'BA.2' = '#d10fa2',
  'BA.1' = '#EF3737',
  'BA.5' = '#F08DA5',
  'BA.1.1' = '#EF3737',
  'BA.2.12.1' = '#d10fa2',
  'A.23.1' = '#3bf5f5',
  'BA.3' = '#EF3737',
  'AY.4.2' = '#d18652',
  'BA.2-12' = '#d10fa2',
  'BA.4' = '#F08DA5',
  'B.1+E484K' = 'cyan',
  'N501Y+D614G' = '#637939',
  'B.1.525' = '#c03bf5',
  'B.1.617.2 (AY.1)+K417N' = '#d18652',
  'B.1.617.2 (AY.2)+K417N' = '#d18652',
  'B.1.617.2 (AY.3)+E484Q' = '#d18652',
  'B.1.630' = '#9dd455'
)


ag_shapes <- list(
	'D614G' = 16,
	'614D' = 20,
	'B.1.1.7' = 16,
	'B.1.526' = 16,
	'B.1.526+S477N' = 16,
	'B.1.429' = 16,
	'C.36.3' = 16,
	'C.37' = 16,
	'B.1.617.2' = 16,
  'B.1.617.2+K417N' = 20,
  'B.1.617.1' = 16,
  'B.1.351' = 16,
  'P.1' = 16,
  'B.1.1.7+E484K' = 20,
  'B.1.621' = 16,
  'B.1.526+E484K' = 16,
  'R.1' = 16,
  'P.2' = 16,
  'BA.2' = 16,
  'BA.1' = 16,
  'BA.5' = 16,
  'BA.1.1' = 20,
  'BA.2.12.1' = 20,
  'A.23.1' = 20,
  'BA.3' = 20,
  'AY.4.2' = 20,
  'BA.2-12' = 16,
  'BA.4' = 20,
  'B.1+E484K' = 16,
  'N501Y+D614G' = 20,
  'B.1.525' = 20,
  'B.1.617.2 (AY.1)+K417N' = 20,
  'B.1.617.2 (AY.2)+K417N' = 20,
  'B.1.617.2 (AY.3)+E484Q' = 20,
  'B.1.630' = 16
)

# Antigen colors by substitutions
wt <- '#000080'
tr <- '#FFA500'
bl <- '#d95f02'
cr <- '#8aa856'
br <- '#800000'
other <- 'grey'
ba1 <- '#800080'
ba2 <- '#d10fa2'
ba45 <- '#BC8F8F'
b117 <- '#008B8B'

ag_colors_info <- data.frame(
	ags = c('D614G', 'B.1.1.7', 'B.1.1.7+E484K', 'P.1', 'B.1.351', 'B.1.429',
			    'B.1.526+E484K', 'B.1.617.1', 'B.1.617.2', 'B.1.617.2+K417N',
			    'B.1.617.2 (AY.1)+K417N', 'B.1.617.2 (AY.2)+K417N',
			    'B.1.617.2 (AY.3)+E484Q', 'B.1.621', 'C.37', 'BA.1', 'BA.1.1', 'BA.2',
			    'BA.2.12.1', 'BA.3', 'BA.5', '614D', 'R.1', 'A.23.1', 'B.1.525',
			    'B.1.526', 'B.1.630', 'B.1.526+S477N', 'P.2', 'C.36.3',
			    'N501Y+D614G', 'AY.4.2', 'BA.2-12', 'BA.4', 'B.1+E484K'),
	ag_cols = c(wt, b117, tr, tr, tr, bl, cr, br,
             bl, bl, bl, bl, br,
             tr, bl, ba1, ba1, ba2, ba2, ba2, ba45,
             wt, cr, wt, cr, wt, br, wt,
             cr, wt, b117, bl, ba2, ba45, cr)
)

ag_who_labels <- c(
  'D614G'                  = 'D614G',
  'B.1.1.7'                = 'B.1.1.7 (__Alpha__)',
  'B.1.1.7+E484K'          = 'B.1.1.7+E484K',
  'P.1'                    = 'P.1 (__Gamma__)',
  'B.1.351'                = 'B.1.351 (__Beta__)',
  'B.1.429'                = 'B.1.429 (__Epsilon__)',
  'B.1.526+E484K'          = 'B.1.526+E484K (__Iota__)',
  'B.1.617.1'              = 'B.1.617.1 (__Kappa__)',
  'B.1.617.2'              = 'B.1.617.2 (__Delta__)',
  'B.1.617.2+K417N'        = 'B.1.617.2+K417N',
  'B.1.617.2 (AY.1)+K417N' = 'B.1.617.2 (AY.1)+K417N',
  'B.1.617.2 (AY.2)+K417N' = 'B.1.617.2 (AY.2)+K417N',
  'B.1.617.2 (AY.3)+E484Q' = 'B.1.617.2 (AY.3)+E484Q',
  'B.1.621'                = 'B.1.621 (__Mu__)',
  'C.37'                   = 'C.37 (__Lambda__)',
  'BA.1'                   = 'BA.1 (__Omicron__)',
  'BA.1.1'                 = 'BA.1.1 (__Omicron__)',
  'BA.2'                   = 'BA.2 (__Omicron__)',
  'BA.2.12.1'              = 'BA.2.12.1 (__Omicron__)',
  'BA.3'                   = 'BA.3 (__Omicron__)',
  'BA.4/BA.5'              = 'BA.4/BA.5 (__Omicron__)'
)

wt_ags <- c(
  'duke' = 'D614G',
  'innsbruck' = 'D614G',
  'emory' = 'D614G',
  'emc_prnt' = 'D614G',
  'oxford' = '614D',
  'mt_sinai_human' = '614D',
  'st_louis' = 'D614G',
  'kcl' = '614D',
  'galveston' = '614D',
  'madison_pooled' = 'D614G',
  'madison_unpooled' = 'D614G',
  'maryland' = '614D',
  'charite' = 'D614G',
  'madison_frnt' = '614D',
  'fda' = 'D614G',
  'geneva' = 'D614G',
  'amc' = 'D614G',
  'emc_calu' = 'D614G',
  'emc_vero' = 'D614G'
)
