# Extract COG Sequences



This script extracts nucleotide and amino acid sequences corresponding to a list of COGs (Clusters of Orthologous Groups) from annotated microbial genomes. It works with Prokka-annotated genomes, eggNOG-mapper output, and genome FASTA files, and organizes the extracted sequences into per-COG folders.

---

## Features

- Parses `.emapper.annotations` files from eggNOG-mapper to identify gene-to-COG mappings.
- Matches genes to coordinates in Prokka GFF files and extracts nucleotide and amino acid sequences.
- Organizes outputs into folders per COG (with gene name included).
- Supports different overwrite policies:
  - `ask`: prompts for each existing file
  - `always`: always overwrite existing files
  - `never`: never overwrite
  - `dry-run`: simulate what would be done, without writing anything
- Logs all steps to `cog_extraction_log.txt` and prints progress to terminal.

---

## Required Input Files

1. **Genome annotations (Prokka):**
   ```
   ./prokka/<genome_id>/<genome_id>.gff
   ./prokka/<genome_id>/<genome_id>.faa
   ```

2. **eggNOG-mapper outputs:**
   ```
   ./eggnog/*.emapper.annotations
   ```

3. **Genome nucleotide FASTA files (.fna):**
   ```
   ./genomes_selected/*.fna
   ```

4. **COG list (`cog_gene_list.txt`)**  
   A **tab-delimited file** with:
   - a required column `cog`
   - and an optional column `gene` (used for folder naming)

   Example:
   ```
   cog	gene
   COG0605	Superoxide_dismutase
   COG2032	CU/Zn_SOD
   COG2077
   ```

   - If `gene` is missing or empty, the COG ID will be used as the gene name.
   - This file **must exist** in the same directory as the script.

5. *(Optional)* **Genome list (`genome_list.txt`)**  
   If present, this file restricts the analysis to the listed genomes.

   - Format: plain text, one genome ID per line.
   - Example:
     ```
     GCA_00012345.1
     GCA_00067890.1
     ```

   - If this file is **not present**, the script will use **all genome folders** in the Prokka directory.

---

## Output

- Nucleotide sequences:
  ```
  ./cog_sequences_nt/COG0605_Superoxide_dismutase/GCA_00012345_COG0605.fna
  ```

- Amino acid sequences:
  ```
  ./cog_sequences_aa/COG0605_Superoxide_dismutase/GCA_00012345_COG0605.faa
  ```

- Log file:
  ```
  cog_extraction_log.txt
  ```

---

## Usage

```bash
python Extract_COG_sequences.py [ask|always|never|dry-run]
```

- If no argument is passed, the default is `ask`.

---

## Examples

**Run using all genome folders (default behavior):**
```bash
python Extract_COG_sequences.py always
```

**Run using a genome ID list file:**
1. Create a `genome_list.txt` file in the same directory as the script:
   ```
   GCA_00012345.1
   GCA_00067890.1
   ```

2. Then run the script:
   ```bash
   python Extract_COG_sequences.py ask
   ```

**Perform a dry run to preview actions:**
```bash
python Extract_COG_sequences.py dry-run
```

---

## Dependencies

- Python 3.6+
- Biopython (`pip install biopython`)
- A properly structured input directory tree as described above

- **Date:** July 2025  
**Institution:** Okinawa Institute of Science and Technology (OIST)  
**Author:** Fatima Li-Hau (with assistance from ChatGPT)
##  Citation

If you use this code in your research, please cite it as:

Li-Hau, F. (2024). Extract COG Sequences Tool (Version 1.0) [Computer software]. Okinawa Institute of Science and Technology. https://github.com/microfafa-gh/Utilities/

Please also consider citing the tools used in this workflow, such as:

- Cock et al., 2009. [Biopython: freely available Python tools for computational molecular biology and bioinformatics](https://doi.org/10.1093/bioinformatics/btp163). Bioinformatics.
- Seemann, T. 2014. [Prokka: rapid prokaryotic genome annotation](https://doi.org/10.1093/bioinformatics/btu153). Bioinformatics.
- Huerta-Cepas et al., 2017. [eggNOG-Mapper: functional annotation of orthologs](https://doi.org/10.1093/molbev/msx148). Molecular Biology and Evolution.



---
