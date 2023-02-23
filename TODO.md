## TO DO

### Now

- Verify backwards compatibility using experiment 1.0_0. Move active work to new experiments repo.
    - Remove hardcoded accessory_data paths from evaluator.py
- Test ggrn_backend3 fully, including demo bug and tests outlined below
    - Implement closest-control or OT matching in ggrn_backend3
    - Clean up prediction code?

- Fix pLI + degree plots etc: build better bins than true versus false and fix the scale on the Y axis and consider some more ESC specific networks and implement a minimum expression filter for Nakatake
- Then try out GGRN backend 3 on real data. 
- Add GEARS to GGRN: https://www.biorxiv.org/content/10.1101/2022.07.12.499735v2.full.pdf
- Add BINGO to GGRN: https://www.nature.com/articles/s41467-020-17217-1#Sec2

### Misc recent

- Simulated data: unfuck the GGRN backend1 simulations
- Clean up do_one_experiment in terms of discrete steps. Save the data splits instead of keeping them all in RAM.
- Fix 1.1.2_1 MAE benefit plot. Fix the system for specifying baseline comparisons so that it properly handles factors like data split, for which the baseline should be vary, and factors like regression model, for which the baseline should be held constant at a specific level. 
- Join Ji & Hansen lab meetings where possible
- GGRN: separate the error model from the biological stochasticity
- Make a schematic for the benchmarking experiments

### All tasks for "Subproject 1"

- Read examples of nature methods ["Analysis" format](https://www.nature.com/nmeth/content)
- Read papers citing stuff besides CO. Keep track of what they are being used for. 
- Backend 3 prototyping: 
    - `torch.linalg.norm(x)` versus `torch.sum(x**2)`???
    - Include sparsity penalty for Q & R???
- Add new datasets and networks
    - ingest perturbations: 
        - Joung atlas https://www.sciencedirect.com/science/article/pii/S0092867422014702
        - ENCODE perturbation data used by SCENIC+ preprint (URL unknown) 
        - KnockTF data? https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6943067/
        - scperturb datasets? https://www.biorxiv.org/content/10.1101/2022.08.20.504663v2.full.pdf
        - FANTOM4 data THP1 data
        - CarPOOL seq THP1 https://www.nature.com/articles/s41592-022-01705-x
        - 2016 K562 perturb-seq data: is it consistent with big modern K562 perturb-seq?
        - Could get [metabolic labeling velocity data in K562 and RPE-1](https://pubmed.ncbi.nlm.nih.gov/32139547/) to complement perturbation data
        - CMAP: downloaded and mostly ingested (look for the notebook) but it has problems
        - Treutlein group brain organoids [here](https://www.nature.com/articles/s41586-022-05279-8)
    - ingest networks:
        - [DoRothEA](https://saezlab.github.io/dorothea/articles/single_cell_vignette.html) (BiocManager::install("dorothea"))
        - [TRRUST](https://www.grnpedia.org/trrust/downloadnetwork.php) ("Transcriptional Regulatory Relationships Unraveled by Sentence-based Text mining")
        - chip-atlas?
        - remap?
        - RegNetwork http://www.regnetworkweb.org/home.jsp (transcription factor and microRNA mediated gene regulations)
        - mouse atlas regulons: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6281296/ and there may be other examples like this; could look into it. 
        - I had previously downloaded a collection used by IRENE, but never finished ingesting it.
        - gTEX has more modern options, including Ashis and Englehardt lab
- Dataset Due Diligence: 
    - Revise QC for it to be more uniform across all datasets?
