# README

Description of columns in `annotation.txt`:

*  **id** - unique ID for sample. Composed of sample number, cell line, and dosage information.
*  **sample** - sample number. Range from 1-230. Each set of 5 are the treatments for one individual, in increasing concentrations of dox.
*  **cell_line** - the cell line number. Each cell line corresponds one individual.
*  **dbgap** - anonymized ID for matching cell_line to genotypes.
*  **dosage** - the concentration of dox used in the treatment.
*  **rin** - the RNA Integrity Number (RIN). Range from 1-10.
*  **rna_conc** - the RNA concentration (ng/uL) after RNA extraction.
*  **lib_prep_batch** - the batches the libraries were prepared in.
*  **index** - the TruSeq adapter index used for multiplexing.
*  **library_conc** - the concentration (ng/uL) of the library measured on the Bioanalyzer.
*  **fragment_size** - the mean fragment size as measured on the Bioanalyzer.
*  **qpcr** - the library concentration (nmol) as measured via qPCR.
*  **qpcr_dilute** - the concentration (nmol) of the diluted library.
*  **master_mix** - the master mix for pooling libraries. 10 samples (2 individuals) per master mix.
*  **lane_perc** - The percentage of sequences in one lane that were assigned to a sample.
