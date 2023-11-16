# ModrRNA1: Potential Bacterial rRNA Modification Site Identification Tool

## Overview

Ribosomal RNA modifications play a crucial role in bacteria, influencing function, ribosome formation, and antibiotic resistance. ModrRNA1 is a user-friendly, streamlined bioinformatics tool designed to identify potential rRNA modification sites and associated enzymes across various bacteria. It replaces manual methods, utilizing a curated database of known bacterial rRNA modifications and enzymes. The tool employs sequence alignment with a customized scoring system for precise site identification, generating annotated RNA sequence plots for high-confidence alignments. Enzyme detection is facilitated through structural comparisons, phylogenetics, and BLAST.

## Paper

Thilanka, U.; Warushavithana, P. S. ModrRNA1: A Tool for Rapid General Identification
of Possible RNA Modification Sites in Bacterial rRNA and Prediction of the Associated
Putative Enzymes. bioRxiv 2023. [https://doi.org/10.1101/2023.11.01.564906](https://doi.org/10.1101/2023.11.01.564906)
(Currently in Peer-review)

## Features

- **Site Annotation:**
  - Redundant sequences are removed based on a threshold accuracy score of 60.
  - Annotated RNA sequence plots are generated using the DnaFeaturesViewer library.
  - Annotations include modification type, nucleotide location, and unique identifiers for enzyme retrieval.

- **Enzyme Detection:**
  - Users can identify potential enzymes involved in a modification if the source bacterium is known.
  - Structural comparisons use the Protein Data Bank (PDB) and AlphaFold.
  - Evolutionary distance calculations aid in selecting phylogenetically related enzymes.
  - Results are outputted as entries to the Protein Data Bank and the AlphaFold Protein Structure Database.

- **Web Application:**
  - Developed using Streamlit, an open-source Python-based application framework.
  - User-friendly interface hosted on an Azure app service instance at [https://modrrna.biotechie.org](https://modrrna.biotechie.org).


## Acknowledgments

We thank RMBase, MODOMICS, ESMFold, Streamlit, and all contributors to the databases and libraries used in the development of ModrRNA1. Your work has been instrumental in making this tool possible.

**Figure 1: Workflow Diagram**
![Workflow Diagram](https://i.ibb.co/3CCNdFM/csd.png)


## Notes
The open-source codebase, available here without restrictions, can be utilized to develop novel tools built upon this pipeline. Especially focusing on enzymes that operate on conserved sequences, similar to the methodology employed in this study.
