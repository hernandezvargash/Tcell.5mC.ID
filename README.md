# Tcell.5mC.ID

## Targeted Epigenetic T cell Subtyping Using Long Reads

Code required to reproduce the analyses published in Goldsmith et al. BiorXiv 2022, "T-cell epigenetic identity profile sequencing (TEIP-seq) for single molecule methylation (5mC) and hydroxymethylation (5hmC) in native DNA long reads".

# Rationale

CD4+ T cells are key mediators of immunity and immunologic memory. Naive CD4+ T cells differentiate into effector T cells with different effector phenotypes which regulate the adaptive immune response. Importantly, T cells are often dysregulated in chronic inflammatory conditions, in particular autoimmune disease and cancer.

Differentiated CD4+ T cells have classically been divided into distinct stable subtypes, including Th1, Th2, Th9, Th17, Th22, Tfh, Treg and Tr1 cells. However, plasticity is an important feature of T cells, whereby they adapt their functions in response to changing circumstances forming a continuum of polarized phenotypes. Plasticity provides advantages to the host for protective immunity and immune control. A loss of effector plasticity is a hallmark of T cell ageing and is linked with reduced immunity and immune function. 

Importantly, the identity and plasticity of T cells is regulated by mechanisms such as DNA methylation. The so-called "epigenetic"" landscape defined by DNA methylation makes each cell type unique. Therefore, mapping epigenetics represents an ideal way to discriminate cell types by ‘epigenotype’. DNA methylation (5mC) is the most stable and most widely studied epigenetic marks and is involved in gene silencing. DNA hydroxymethylation (5hmC) is largely understudied, but recently has been linked to gene activity during differentiation.

Despite the evidence for 5mC/5hmC as a marker of cellular identity, in particular in T cells, technical limitations have prevented this knowledge from reaching clinical use. Simultaneous detection of cytosine modifications (including 5mC and 5hmC) has been made recently available by the next generation of Nanopore long-read sequencers, that are able to provide such information in the native DNA context. Here, we take advantage of Nanopore sequencing for simultaneously profiling 5mC and 5hmC in native DNA on naive murine T cells polarized under Th1, Th2, Th17, Th1/17 and Treg conditions.

Our work is a first step towards the implementation of 5mC/5hmC in clinical immune profiling and identification of pathogenic immune cell subtypes. 


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

...


---