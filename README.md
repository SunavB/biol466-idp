# biol466-idp
Project Summary
Intrinsically disordered proteins (IDPs) lack a stable three-dimensional structure yet carry out
essential biological roles through context-dependent mechanisms. The most studied of these is
coupled folding-and-binding — disorder-to-order transition — in which a disordered region
adopts a defined structure on engaging a partner. This project asks not simply whether Gene
Ontology (GO) annotations — Biological Process (BP), Molecular Function (MF) and Cellular
Component (CC) — help predict this canonical IDP functional behaviour, but how much they add
beyond a strong modern sequence representation, which aspect contributes most, whether the
aspects are redundant or synergistic, and whether their benefit is concentrated in proteins the
sequence model predicts poorly.
Using the curated DisProt database, a single fixed machine-learning classifier (a random forest)
will be trained under a 2×2×2 factorial design that switches each GO aspect on or off, producing
eight feature configurations including a sequence-only control. The prediction target is a binary
IDPO label — does the disordered region undergo a disorder-to-order transition — drawn from
DisProt's Structural transition layer and deliberately kept outside the GO ontology, so the design
carries no GO-to-GO circularity. The sequence-only baseline is deliberately strong: amino acid
sequences are represented using embeddings from the pre-trained ESM-2 protein language
model, so that any measured contribution of GO context is a contribution over and above a
competitive model rather than over a weak strawman. Performance is compared with stratified
cross-validation, paired statistical testing with multiple-comparison correction, a label-shuffling
negative control, and analysis of main and interaction effects.
