Endothelial = c('Cldn5', 'Pecam1', 'Ocln', 'Slc2a1', 'Abcb1', 'Sema3g', 'Bmx', 'Mfsd2a', 'Tfrc', 'Slc38a5')
Mural = c('Pdgfrb', 'Kcnj8', 'Acta2', 'Myl9')
Astrocyte = c('Gfap', 'Aqp4', 'Clu', 'Aldh1l1')
Oligodendrocyte = c('Mog', 'Mbp', 'Plp1', 'Apod')
Oligo_Prog = c('Pdgfra','Vcan', 'Olig2', 'C1ql1')
Microglia = c('C1qb', 'C1qa', 'Csf1r', 'Ctss', 'Tmem119')
dSPN = c('Drd1', 'Pdyn', 'Chrm4', 'Foxp2', 'Chrm4')
iSPN = c('Drd2', 'Adora2a', 'Gpr6', 'Penk', 'Sp9')
Cholinergic_IN = c('Chat', 'Slc5a7')
GABAergic_IN = c('Npy', 'Sst', 'Nos1')
PV_IN = c('Pvalb', 'Cox6a2', 'Kit', 'Crtac1', 'Adamts5')
T_cell = c('Cd4', 'Cd8a')
Cil_ependymal = c('Hydin', 'Armc4', 'Dnali1', 'Spag17')
Sec_ependymal = c('Npr3', 'Prlr', 'Slc4a5')

marker_genes = list('Endothelial' = Endothelial,
                      'Mural' = Mural,
                      'Astrocyte' = Astrocyte,
                      'Oligodendrocyte' = Oligodendrocyte,
                      'OPC' = Oligo_Prog,
                      'Microglia' = Microglia,
                      'dSPN' = dSPN,
                      'iSPN' = iSPN,
                      'Cholinergic_Interneuron' = Cholinergic_IN,
                      'GABAergic_Interneuron' = GABAergic_IN,
                      'T-Cell' = T_cell,
                      'PV_Interneuron' = PV_IN,
                      'Cil_Ependymal' = Cil_ependymal,
                      'Sec_Ependymal' = Sec_ependymal
                    )

marker_genes <- lapply(marker_genes[order(names(marker_genes))], toupper)
