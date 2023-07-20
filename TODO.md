## TO DO

### Personal

- Train tickets: Portland bike tour rezi
- Lab captains: Monika and Elise??
- Make a list of MRB organizers (add to BroadStripes?)
    - Charlie Walters
    - Jacob Mitchell
    - Noah Workman
    - Mara Lanis
- Onboarding doc

### Now

- Add suggested reviewers
    - Meinshausen & Buhlmann?
    - Sushmita Roy
    - TM Murali or other BEELINE authors
    - SCENIC-related people: Stein Aerts, Pierre Guerts, Sara Aibar, VA Huynh-Thu, or Carmen Bravo Gonzales-Blas?
    - Geoffrey Schiebinger? (the WOT guy)
    - Julio Saez-Rodrigues (wrote recent GRN review https://www.nature.com/articles/s41576-023-00618-5)
- contact GEARS and let them know what's going on
    > Dear Dr. Leskovic,
    I’m a PhD student at JHU working with Patrick Cahan and Alexis Battle. We have been working with GEARS recently. We wanted to give you a heads up that we haven’t been able to reproduce figure 2b (as discussed with Yusuf Roohani on github issue #5), and are seeing overall performance that doesn’t exceed baseline levels on some of our evaluations. We would be happy to hear feedback from your team. Our code will be made public soon, so if there any errors or other comments on how we have been running GEARS, we can discuss further and update our analysis.
    Best,
    Eric Kernfeld

- invite Josh W to view the code
- Show MAE instead of MAE benefit, or show MAE benefit as a percentage of baseline MAE?
- REMOVE autoregressive models and separate software
- Update simulations to include more base networks??
- **Goal: post benchmarking preprint week of July 24th**

### Eventually maybe

- Handle indirect effects:
    - Treat perturbed TF as unknown / random
    - Try CO-style iteration
    - try alternate matching scheme

- Add scFormer (take over from Zexi)
- Try out the "simple" split again
- PRESCIENT:
    - NN architecture?
    - Project then simulate for KO's?
    - Beta-cell figure bottom panel is real data or all computational? 
- Implement & test [OT matching](https://pythonot.github.io/auto_examples/plot_OT_1D.html) in ggrn_backend3
- Read examples of nature methods ["Analysis" format](https://www.nature.com/nmeth/content)
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
- Briefly review supervised GRN inference and summarize the shared and unique elements of the following projects. 
    - [GIANT](https://giant.princeton.edu/)
    - [GRGNN](https://pubmed.ncbi.nlm.nih.gov/33294129/)
    - [GRADIS](https://www.nature.com/articles/s41540-020-0140-1)
    - [SVM](https://academic.oup.com/bib/article/15/2/195/212108)
    - [DeepTFni](https://www.nature.com/articles/s42256-022-00469-5#Abs1)
- LLM's:
    - GeneFormer: **later**: There is a table naming all included data.Training data do not include Frangieh, adamson, dixit, or replogle perturb-seqs, and actually are mostly huge whole-organism atlases with not a large total number of studies covered. They specifically exclude genetically weird cells, e.g. cancer, so melanoma (Frangieh) or K562 (replogle, adamson) would be excluded. 30M cells total. 
    - scGPT: **later**. No perturbation code available yet. Training data are "10M blood and bone marrow cells from CELLxGENE". There is no table explicitly naming all included data. Best guess is that it does not include any perturb-seq. 
    - scFoundation: **later** Their [github link](https://github.com/biomap-research/scFoundation) is dead! Training data include 50M cells of diverse types, scraped from diverse sources e.g. GEO. There is no table naming all included data.
    - scFormer: **Zexi** They have a [perturbation demo](https://github.com/bowang-lab/scFormer/blob/main/examples/perturbation/dev_perturb.py), but it suffers from grievous lack of encapsulation. This model seems to do no pretraining, so no worries about training set contamination.
    - scBERT: **never** This is another scRNA Transformer, pretrained on PanglaoDB, but it does not offer zero-shot prediction of genetic perturbation outcomes.
    - Exceiver: **never** Trained on Tabula Sapiens, 500k cells. Does not offer zero-shot prediction of genetic perturbation outcomes.