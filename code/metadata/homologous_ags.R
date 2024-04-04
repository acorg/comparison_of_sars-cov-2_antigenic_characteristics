
homologous_ags_mt_sinai_human <- c(
  'B.1.351 convalescent' = 'B.1.351',
  'mRNA-1273' = '614D'
)

homologous_ags_emory <- c(
  'B.1.351 convalescent' = 'B.1.351',
  'mRNA-1273' = 'D614G',
  'WT convalescent' = 'D614G',
  'B.1.617.2 convalescent' = 'B.1.617.2'
)

homologous_ags_oxford <- c(
  'B.1.351 convalescent' = 'B.1.351',
  'P.1 convalescent' = 'P.1',
  'B.1.1.7 convalescent' = 'B.1.1.7',
  'mRNA-1273' = '614D',
  'AstraZeneca' = '614D',
  'WT convalescent' = '614D',
  'B.1.617.2 convalescent' = 'B.1.617.2'
)

homologous_ags_innsbruck <- c(
  'B.1.617.2 convalescent' = 'B.1.617.2',
  'B.1.351 convalescent' = 'B.1.351',
  'BA.1 convalescent' = 'BA.1',
  'BA.2 convalescent' = 'BA.2',
  'B.1.1.7 convalescent' = 'B.1.1.7',
  'mRNA-1273' = 'D614G',
  'AstraZeneca' = 'D614G',
  'D614G convalescent' = 'D614G',
  'Pfizer' = 'D614G',
  'AstraZeneca-Pfizer' = 'D614G'
)

homologous_ags_amc <- c(
  'B.1.617.2 convalescent' = 'B.1.617.2',
  'B.1.351 convalescent' = 'B.1.351',
  'BA.1 convalescent' = 'BA.1',
  'P.1 convalescent' = 'P.1',
  'BA.2 convalescent' = 'BA.2',
  'B.1.1.7 convalescent' = 'B.1.1.7',
  'D614G convalescent' = 'D614G'
)

homologous_ags_geneva <- c(
  'B.1.617.2 convalescent' = 'B.1.617.2',
  'B.1.351 convalescent' = 'B.1.351',
  'P.1 convalescent' = 'P.1',
  'B.1.1.7 convalescent' = 'B.1.1.7',
  'D614G convalescent' = 'D614G',
  'mRNA-1273' = 'D614G'
)

homologous_ags_fda <- c(
  'D614G convalescent' = 'D614G',
  'B.1.429 convalescent' = 'B.1.429',
  'B.1.1.7 convalescent' = 'B.1.1.7',
  'B.1.617.2 convalescent' = 'B.1.617.2',
  'C.37 convalescent' = 'C.37',
  'B.1.351 convalescent' = 'B.1.351',
  'P.1 convalescent' = 'P.1',
  'BA.1 convalescent' = 'BA.1',
  'D614G convalescent' = 'D614G',
  'B.1.526+E484K convalescent' = 'B.1.526+E484K',
  'other' = 'D614G'
)


srHomologousLogTiters <- function(map) {
  logtiters <- logtiterTable(map)
  homologous_ags <- srHomologousAgs(map)
  lapply(seq_len(numSera(map)), function(i) {
    logtiters[homologous_ags[[i]], i]
  })
}
