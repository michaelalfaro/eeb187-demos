EEB 187 — Instructor Walkthrough Bundle: Chaetodontidae
=======================================================
Exemplars: 130 species (curated for the divergence analysis)
Tree tips: 109 (chronogram, crown ~43 Ma, ladderized)
Tips matching exemplars: 109

Files (mirroring the student FishView export):
  images/<Genus_species>.png       — one exemplar per species
  chaetodontidae-tree.tre          — time-calibrated phylogeny
  picked_species.txt               — exemplar species list
  picked_tip_labels.txt            — tree tips matching exemplars

Cached intermediates (recipe runs these once, then comments the call):
  pavo_classify.rds                — pavo::classify(..., kcols=k) per species
  phylo_pca.rds                    — phytools::phyl.pca on pavo metrics

Source of original analysis:
  /Users/michaelalfaro/Dropbox/git/chaets-divergence-2026/analysis
