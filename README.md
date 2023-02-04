# Tcell.5mC.ID

# <img src="figures/DALLE2.png" align="center" width="80" />

## Targeted Epigenetic T cell Subtyping Using Long Reads

Code required to reproduce the analyses published in Goldsmith et al. BiorXiv 2023, *"Single molecule DNA methylation and hydroxymethylation reveal unique epigenetic identity profiles of T helper cells"*, by Chloe Goldsmith, Olivier Fesneau, Valentin Thevin, Maria I. Matias, Julie Perrault, Ali Hani Abid, Naomi Taylor, Valérie Dardalhon, Julien C. Marie, and Hector Hernandez-Vargas.

# Rationale

Both identity and plasticity of CD4 T helper Th cells are regulated in part by epigenetic mechanisms2,4 However, a method that reliably and readily profiles DNA base modifications is still needed to finely study Th cell differentiation. Cytosine methylation (5mC) and, as well as cytosine hydroxymethylation (5hmC) are  DNA modifications that identify stable cell phenotypes but their potential to  characterize intermediate cell transitions has not yet been evaluated. To assess transition states in Th cells, we developed a new method to profile Th cell identity using cas9-targeted single molecule nanopore sequencing and find that 5mC and 5hmC can be used as markers of cellular identity. Targeting as few as 10 selected genomic loci, we were able to distinguish major differentiated T cell subtypes as well as intermediate phenotypes by their native DNA 5mC/5hmC patterns. Moreover, by using off-target sequences we were able to infer transcription factor activities relevant to each cell subtype. Our analysis demonstrates the importance of epigenetic regulation by 5mC and 5mhC modifications in the establishment of Th cell identity. Furthermore, our data highlight the potential to exploit this immune profiling application to elucidate the pathogenic role of Th transition states in autoimmune diseases.

# Methods

## In vitro differentiation 
Naive T cells were isolated from male Foxp3-GFP reporter mice and were cell sorted and activated with anti-CD3/CD28 under different polarizing conditions: Th0 (neutral), Th1, Th2, Th17s, Th1/17s and Tregs. 

## Library prep and targeted nanopore sequencing 
DNA was extracted from 1x10^6 CD4+ T cells/ replicate using the Qiagen DNA mini kit according to manufacturer's instructions (Qiagen Hilden, Germany).
A double cutting approach was used whereby 2x S. pyogenes Cas9 Alt-R™ trans-activating crispr RNAs (crRNAs) were designed for both upstream and downstream of each target loci. crRNAs were pooled in equimolar ratios (100 µM) and annealed to S. pyogenes Cas9 Alt-R™ tracrRNA (100 µM) (Integrated DNA Technologies, Iowa United States). crRNA•tracrRNA pool (10 µM) incubated with Cut smart buffer (NEB Cat # B7204) and Alt-R® S. pyogenes HiFi Cas9 nuclease V3 (62 µM) forming crRNA-tracrRNA-Cas9 ribonucleoprotein complexes (RNPs). 3000ng of high molecular weight DNA was prepared by blocking available DNA ends with calf intestinal phosphatase (NEB cat #M0525). DNA and dA-tails on all available DNA ends were cleaved by RNPs and Taq polymerase (NEB Cat # M0273) and dATP (NEB Cat # N0440), activating the Cas9 cut sites for ligation. Nanopore Native Barcodes (Oxford Nanopore Technology, #EXP-NBD104) were ligated to available DNA ends with NEB Blunt/TA Ligase Master Mix (NEB #M0367) and up to 4 samples were pooled before ligation of Nanopore adapters with NEBNext® Quick T4 DNA Ligase (NEB #E6057) and loading onto a primed MinION flow cell (pore chemistry 9.4.1) and sequencing for 72h. Raw read quality was determined with pycoQC. 

## Detection of 5mC and 5hmC from Nanopore signal data 
Raw fast5 files were basecalled with Megalodon (version 2.5.0, ONT) allowing the simultaneous delineation of 5mC and 5hmC by remora models (ONT) for the detection of modified bases. To this end, the following general megalodon parameters were used:

*megalodon $FAST5PATH/$barcode/  --guppy-server-path $GUPPY --guppy-config $CONFIG --remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0 --outputs basecalls mods per_read_mods mod_basecalls mappings mod_mappings mods --reference $REF --devices 0 --processes 8 --output-directory ${OUTDIR}/${barcode}*

## Differential methylation
Differential methylation was detected by DSS as previously described. Data processing and statistical analyses were performed using R/Bioconductor (R version 4.0.3). Bed files containing modified cytosine information were transformed into a BSseq object for differential methylation analysis with dispersion shrinkage for sequencing data (DSS) as previously described. DSS tests for differential methylation at single CpG-sites were evaluated using a Wald test on the coefficients of a beta-binomial regression of count data with an ‘arcsine’ link function. For DSS, a multifactor model of experimental conditions (i.e., Th0, Th1, Th2, Th17 ,Th1/17 and Tregs) was used to account for the biological replicates as well as on- and off-target regions, regardless of coverage. Differentially methylated regions (DMRs) were defined as those loci with at least 3 CpG sites within a distance of less than 50 bp, and with changes in > 50% of all CpG sites exhibiting p < 0.05. DMRs were plotted using the plotDMRs function of the dmrseq package. 

Off-targets were analysed for transcription factor (TF) activity using the cistrome information associated to the MIRA R package. After extracting all genomic binding regions corresponding to selected TFs, MIRA functions were used for aggregation of DNA methylation data and visualization of TF activity.

# Associated website

An interactive shiny application can be accessed **<a href="http://20.56.136.251:3838/test.app/">here</a>.**

# GitHub Contributors

- Chloe Goldsmith.
- Ali Hani Abid.
- Hector Hernandez-Vargas.

---
