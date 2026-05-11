# ============================================================================
# Demo 7 — Comparative Phylogenetics, Part 1
# EEB 187: Ecology and Evolution of Color
# Week 8, Wednesday May 20, 2026
#
# WHAT WE'RE DOING:
#   Learning to map COLOR PATTERN traits onto a phylogenetic tree and
#   test whether closely related species are more similar in their
#   color patterns than expected by chance.
#
#   We'll work with surgeonfishes (Acanthuridae) — a reef fish family
#   with dramatic ontogenetic color shifts and diverse pattern elements
#   (stripes, spots, marbling, eyespots, saddles).
#
#   Color pattern data: Frederich et al. (2026) — 59 pattern traits
#   coded across head, trunk, and tail for 83 species.
#
#   Phylogeny: Frederich et al. — 65-tip time-calibrated tree.
#
# WHAT YOU NEED:
#   R (version 4.2+), with packages: ape, phytools
#   If you don't have them: install.packages(c("ape", "phytools"))
#
# HOW TO USE THIS SCRIPT:
#   Don't run the whole thing at once!
#   Go section by section. Read the comments. Run each block.
#   Look at what it produces. Talk to your neighbor about what you see.
# ============================================================================


# ---- SECTION 1: Load packages ----

library(ape)       # reading and manipulating phylogenetic trees
library(phytools)  # visualization + phylogenetic signal tests


# ---- SECTION 2: Load the surgeonfish tree ----
# This is a time-calibrated tree — branch lengths are in MILLIONS OF YEARS.
# "Time-calibrated" means we can read off when lineages diverged.

tree <- read.nexus("data/frederich/acanthuridae-tree.tree")

# What did we just load?
tree

# How many species?
cat("Number of tips:", length(tree$tip.label), "\n")

# What are some of the species?
head(tree$tip.label, 10)

# What genera are in the tree?
genera <- unique(sub("_.*", "", tree$tip.label))
cat("Genera:", paste(genera, collapse = ", "), "\n")

# You should see: Acanthurus, Ctenochaetus, Naso, Zebrasoma,
# Paracanthurus (Dory!), and Prionurus.


# ---- SECTION 3: Plot the tree ----
# Let's see what it looks like.

# Simple rectangular tree:
plot(tree, cex = 0.45, no.margin = TRUE)

# TAKE A MOMENT: Can you spot the major genera?
# Acanthurus is the biggest genus. Naso (unicornfishes) is distinctive.
# Zebrasoma (tangs) should cluster together.

# Fan-shaped version (fits better on screen):
plotTree(tree, type = "fan", fsize = 0.4, lwd = 1.5)

# How old is this family?
cat("Root age:", round(max(branching.times(tree)), 1), "million years\n")

# QUESTION FOR YOUR NEIGHBOR:
# Surgeonfishes first diverged ~65 million years ago.
# That's right around the end-Cretaceous mass extinction.
# Coincidence? (We'll come back to this in Lec 17.)


# ---- SECTION 4: Load the color pattern data ----
# Frederich et al. coded 52 binary color-pattern traits for each species,
# organized by body region:
#   H_  = Head
#   Tr_ = Trunk (body)
#   Ta_ = Tail
#
# Each trait is 0 (absent) or 1 (present).
# Examples: H_marbling, Tr_1_h_strip, Ta_eyespot, Ta_saddle

traits <- read.csv("data/frederich/acanthuridae-color-patterns.csv",
                   stringsAsFactors = FALSE, check.names = FALSE)

# How many species? How many traits?
cat("Species:", nrow(traits), "\n")
cat("Columns:", ncol(traits), "\n")

# Look at the trait names:
names(traits)

# The first 7 columns are metadata (species name + biogeographic regions).
# Columns 8 onward are the color pattern traits.


# ---- SECTION 5: Match tree tips to trait data ----
# THE GOLDEN RULE OF COMPARATIVE METHODS:
# The species in your tree and your data MUST match exactly.

overlap <- intersect(tree$tip.label, traits$species)
cat("Species in BOTH tree and data:", length(overlap), "\n")
cat("In tree but NOT in data:", length(setdiff(tree$tip.label, traits$species)), "\n")
cat("In data but NOT in tree:", length(setdiff(traits$species, tree$tip.label)), "\n")

# Prune the tree to only the species we have data for:
tree_matched <- drop.tip(tree, setdiff(tree$tip.label, overlap))
cat("Matched tree:", length(tree_matched$tip.label), "tips\n")

# Also filter the trait table:
tr <- traits[traits$species %in% overlap, ]

# LESSON: Real datasets NEVER have perfect overlap.
# You always lose some species. That's normal.
# The important thing is to be explicit about what you kept and what you dropped.


# ---- SECTION 6: Create a "pattern complexity" score ----
# Before mapping individual traits, let's create a SUMMARY variable:
# how many total color pattern elements does each species have?
#
# This is like asking: "How visually complex is this fish?"

# Sum all the binary pattern traits (columns 8 onward):
pattern_cols <- names(tr)[8:(ncol(tr))]
tr$pattern_complexity <- rowSums(tr[, pattern_cols], na.rm = TRUE)

# Look at the distribution:
hist(tr$pattern_complexity, breaks = 12, col = "#6fa8dc",
     main = "Color pattern complexity across surgeonfishes",
     xlab = "Number of pattern elements", ylab = "Number of species")

# Who's the most complex? Who's the simplest?
cat("Most complex:", tr$species[which.max(tr$pattern_complexity)],
    "—", max(tr$pattern_complexity), "elements\n")
cat("Simplest:    ", tr$species[which.min(tr$pattern_complexity)],
    "—", min(tr$pattern_complexity), "elements\n")

# Compare genera:
cat("\n=== Average complexity by genus ===\n")
genus <- sub("_.*", "", tr$species)
for (g in unique(genus)) {
  vals <- tr$pattern_complexity[genus == g]
  cat(sprintf("  %-15s  n = %2d   mean = %.1f\n", g, length(vals), mean(vals)))
}

# QUESTION: Which genus has the most complex color patterns on average?
# Does that surprise you?


# ---- SECTION 7: Map pattern complexity onto the tree (contMap) ----
# contMap "paints" a continuous trait along the branches of the tree.
# It estimates ancestral values by interpolation and shows you
# WHERE on the tree the trait changed.

# Create a named vector (phytools requires this format):
pc <- tr$pattern_complexity
names(pc) <- tr$species

# Generate the contMap:
cm <- contMap(tree_matched, pc, plot = FALSE)

# Plot it:
plot(cm, type = "fan", fsize = 0.4, lwd = 3,
     legend = 0.7 * max(nodeHeights(tree_matched)))
title("Color pattern complexity mapped on Acanthuridae phylogeny")

# WHAT TO LOOK FOR:
# - Cool colors (blue) = SIMPLE patterns (few elements)
# - Warm colors (red)  = COMPLEX patterns (many elements)
# - Are closely related species similar in complexity?
# - Can you spot any SHIFTS — branches where complexity changes abruptly?
# - Do any genera stand out as consistently complex or simple?


# ---- SECTION 8: Test phylogenetic signal — Blomberg's K ----
# The key question:
#   "Is color pattern complexity MORE similar among close relatives
#    than expected by chance?"
#
# Blomberg's K tells us:
#   K = 1  → trait matches Brownian motion on the tree exactly
#   K > 1  → relatives are MORE similar than BM predicts (strong conservation)
#   K < 1  → relatives are LESS similar than BM predicts (more labile)
#   K ≈ 0  → no relationship between tree and trait
#
# The p-value comes from shuffling the trait values across tips 1000 times.

K_result <- phylosig(tree_matched, pc, method = "K", test = TRUE)
K_result

cat("\n--- YOUR INTERPRETATION ---\n")
cat("Blomberg's K =", round(K_result$K, 3), "\n")
cat("p-value      =", round(K_result$P, 4), "\n")

if (K_result$P < 0.05) {
  cat("→ SIGNIFICANT phylogenetic signal in color pattern complexity!\n")
  cat("  Closely related surgeonfishes DO have more similar patterns.\n")
} else {
  cat("→ No significant signal. Patterns are NOT predicted by relatedness.\n")
}

# THINK ABOUT THIS:
# K < 1 means there IS signal, but it's weaker than pure Brownian motion.
# What could cause that? (Hint: think about convergent evolution,
# or natural selection pushing distantly related species toward
# similar patterns in similar habitats.)


# ---- SECTION 9: Pagel's lambda — a second test ----
# Lambda asks the same question a different way:
#   λ = 1  → the tree perfectly predicts trait similarity
#   λ = 0  → the tree is irrelevant (star phylogeny)

lambda_result <- phylosig(tree_matched, pc, method = "lambda", test = TRUE)
lambda_result

cat("\nPagel's λ =", round(lambda_result$lambda, 3), "\n")
cat("p-value   =", round(lambda_result$P, 4), "\n")

# Do K and λ agree? They often give slightly different answers.
# λ is more sensitive to deep phylogenetic signal;
# K is more sensitive to tip-level similarity.


# ---- SECTION 10: Map a binary trait — trunk marbling ----
# Now let's look at a SINGLE, specific color pattern element.
# "Trunk marbling" = irregular, worm-like markings on the body.
# This connects directly to Lecture 12's Turing patterns!

marbling <- tr$Tr_marbling
names(marbling) <- tr$species

# How many species have trunk marbling?
cat("Trunk marbling present:", sum(marbling), "of", length(marbling), "species\n")

# Color the tree tips by marbling presence:
cols <- c("0" = "#6fa8dc", "1" = "#e07a5f")
plot(tree_matched, type = "fan", fsize = 0.35, lwd = 1.5,
     tip.color = cols[as.character(marbling[tree_matched$tip.label])])
legend("topleft",
       legend = c("No trunk marbling", "Trunk marbling present"),
       fill = c("#6fa8dc", "#e07a5f"), cex = 0.9, bty = "n")
title("Trunk marbling in Acanthuridae")

# QUESTIONS:
# 1. Does marbling cluster on the tree, or is it scattered randomly?
# 2. How many times do you think marbling evolved independently?
# 3. Could this pattern be a Turing-generated trait? (Think back to Lec 12)


# ---- SECTION 11: Try another trait on your own ----
# Pick ONE of these traits and repeat Sections 7-9:
#
# SUGGESTIONS (pick one):
#   tr$H_marbling          — marbling on the HEAD
#   tr$Ta_saddle           — saddle markings on the TAIL
#   tr$H_1_v_strip         — a single vertical stripe on the HEAD
#   tr$Ta_colored_scalpel  — colored caudal spine (surgeonfish weapon!)
#   tr$Ta_eyespot          — eyespot on the TAIL
#
# STEPS:
#   1. Create a named vector
#   2. Run contMap (if continuous) or plot tip colors (if binary)
#   3. Run phylosig with method="K"
#   4. Write down: What is K? Is p < 0.05? What does it mean?

# YOUR CODE HERE:
# my_trait <- tr$___________
# names(my_trait) <- tr$species
# ...


# ---- SECTION 12: Wrap-up discussion ----
#
# TODAY YOU LEARNED:
#   1. How to LOAD a phylogeny and trait data into R
#   2. How to MATCH tips to data (and handle mismatches)
#   3. How to VISUALIZE traits on a tree (contMap, tip colors)
#   4. How to TEST for phylogenetic signal (Blomberg's K, Pagel's λ)
#   5. How to INTERPRET what signal means for color pattern evolution
#
# THE BIG QUESTION:
#   Are color patterns phylogenetically conserved in surgeonfishes?
#   Your answer should be nuanced — some traits show signal, others don't.
#   WHY might that be?
#
# NEXT WEEK (Part 2):
#   You'll apply these exact methods to YOUR group's fish lineage,
#   using images from FishView and color metrics from pavo.
#   The question becomes: "Is color diversity conserved in MY clade?"
#
# ============================================================================
# DATA CITATIONS:
#
# Tree + color patterns:
#   Frederich, B., et al. (2026). Rapid and repeated evolution of pigmentation
#   patterns in reef fishes. BMC Evolutionary Biology.
#
# See also (labrid phylogeny for your final projects):
#   Brownstein, C. D., et al. (2025). Phylogenomics establishes an Early
#   Miocene reconstruction of reef vertebrate diversity. Science Advances,
#   11, eadu6149. https://doi.org/10.1126/sciadv.adu6149
# ============================================================================
