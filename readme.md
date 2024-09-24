
# Explicación de las columnas de `GSE174255_sgRNA-Read-Counts.xls`

1.	sgRNA: The specific single guide RNA (sgRNA) used for each perturbation. sgRNAs guide the CRISPR system to target specific  genes for activation (CRISPRa) or inhibition (CRISPRi).
2.	Gene: The target gene associated with each sgRNA. This indicates which gene is being perturbed by the CRISPR system.
3.	01_CalabreseSetA_Plasmid: Control or baseline condition, representing sgRNA abundance from the plasmid library before cells were sorted based on cytokine expression. This column helps to normalize the results of the screen.
4.	03_CalabreseSetA_Donor1_IFNG_low: The sgRNA abundance for Donor 1 in the subset of cells that express low levels of interferon gamma (IFN-γ). This indicates how sgRNAs are distributed among cells that produce low IFN-γ, suggesting genes that might be downregulating IFN-γ production.
5.	04_CalabreseSetA_Donor1_IFNG_high: The sgRNA abundance for Donor 1 in the subset of cells that express high levels of IFN-γ. This shows how sgRNAs are distributed among cells with high IFN-γ production, indicating potential genes involved in upregulating IFN-γ.
6.	05_CalabreseSetA_Donor1_IL2_low: The sgRNA abundance for Donor 1 in the subset of cells expressing low levels of interleukin 2 (IL-2), identifying genes that might negatively regulate IL-2 production.
7.	06_CalabreseSetA_Donor1_IL2_high: The sgRNA abundance for Donor 1 in the subset of cells with high levels of IL-2. Genes associated with these sgRNAs might positively regulate IL-2 production.
8.	07_CalabreseSetA_Donor1_IFNG_unsorted: The sgRNA abundance for Donor 1 in an unsorted population of cells for IFN-γ. This acts as a reference to compare sgRNA distribution before sorting based on IFN-γ expression.
9.	08_CalabreseSetA_Donor1_IL2_unsorted: The sgRNA abundance for Donor 1 in the unsorted population for IL-2. Similar to IFN-γ unsorted, this provides a baseline to see sgRNA distribution before sorting by IL-2 expression.
10.	15_CalabreseSetA_Donor2_IFNG_low: Same as column 3, but for Donor 2, showing sgRNA distribution in cells with low IFN-γ expression.
11.	16_CalabreseSetA_Donor2_IFNG_high: Same as column 4, but for Donor 2, showing sgRNA distribution in cells with high IFN-γ expression.
12.	17_CalabreseSetA_Donor2_IL2_low: Same as column 5, but for Donor 2, showing sgRNA distribution in cells with low IL-2 expression.
13.	18_CalabreseSetA_Donor2_IL2_high: Same as column 6, but for Donor 2, showing sgRNA distribution in cells with high IL-2 expression.
14.	19_CalabreseSetA_Donor2_IFNG_unsorted: Same as column 7, but for Donor 2, showing sgRNA distribution in the unsorted population for IFN-γ.
15.	20_CalabreseSetA_Donor2_IL2_unsorted: Same as column 8, but for Donor 2, showing sgRNA distribution in the unsorted population for IL-2.

Key Points:

	•	The columns represent different conditions and experimental groups based on two donors and the cytokines IL-2 and IFN-γ.
	•	The “low” and “high” columns correspond to bins of cells sorted by the amount of cytokine they produce.
	•	Unsorted columns provide a baseline sgRNA distribution before sorting the cells based on their cytokine production levels.

This structure allows researchers to assess which genes, when activated or repressed, lead to increased or decreased IL-2 and IFN-γ production, thereby identifying key regulators of immune function.

# Explicación GSE174284_gene_counts_raw.txt

1.	‘Geneid’: This column contains the unique identifiers for the genes. These IDs are usually gene symbols or IDs from a reference genome, allowing you to track which gene each row refers to.

The remaining columns describe gene expression data for different donors and experimental conditions. The labels are structured as:

	•	DonorX: This indicates the biological donor (e.g., Donor1, Donor2, etc.). There are four donors (Donor1 to Donor4).
	•	sgControl vs sgFOXQ1:
	•	sgControl refers to a control group where no specific gene is targeted.
	•	sgFOXQ1 refers to an experimental group where the FOXQ1 gene has been knocked out or altered using a CRISPR/Cas9-based technique.
	•	Stim vs NoStim:
	•	Stim refers to cells that have been stimulated with some treatment or condition.
	•	NoStim refers to unstimulated cells.

Thus, each sample in the dataset is a combination of:

	•	A biological donor,
	•	A control (sgControl) or FOXQ1-targeted group (sgFOXQ1),
	•	Stimulated or unstimulated condition.

Example:

	•	‘1_Donor1_sgControl_NoStim’: Gene expression data from Donor 1, in the control group (no gene knockout), under unstimulated conditions.
	•	‘4_Donor1_sgFOXQ1_Stim’: Gene expression data from Donor 1, with the FOXQ1 gene knocked out, under stimulated conditions.

This structure is repeated for the four donors, giving a total of 16 samples across the various conditions.

# Estructura de la Matriz

Structure of the Matrix

	1.	Rows:
	•	Each row corresponds to a unique gene, identified by either an Ensembl Gene ID (e.g., 'ENSG00000243485') or a common gene name with an isoform suffix (e.g., 'TNFRSF9-1', 'TNFRSF9-2').
	2.	Columns:
	•	Each column corresponds to a unique cell, identified by a cell barcode (e.g., 'AAACCCACAACAAGAT-1').
	3.	Values:
	•	The entries in the matrix represent the expression levels of the corresponding gene in the respective cell. These values are typically counts of RNA transcripts detected in each cell, indicating how much of each gene is expressed.

Summary of Contents

	•	Matrix Dimensions: The matrix is likely structured as follows:
	•	Number of Rows: Represents the total number of genes.
	•	Number of Columns: Represents the total number of cells.
	•	Matrix Entries: Each entry in the matrix can be described as:
	•	Non-zero Entries: Indicate that the corresponding gene is expressed in that particular cell, with the value reflecting the expression level (e.g., the count of transcripts).
	•	Zero Entries: Indicate that the gene is not expressed in that cell, and these entries may be omitted in sparse representations.

Example Interpretation

If you have an entry at row ENSG00000243485 (gene A) and column AAACCCACAACAAGAT-1 (cell 1) with a value of 10, this would indicate that in cell 1, gene A is expressed at a level of 10 counts.

Use in Analysis

This matrix allows researchers to perform various analyses, including:

	•	Clustering: Identifying groups of similar cells based on their expression profiles.
	•	Differential Expression Analysis: Comparing gene expression levels between different cell types or conditions.
	•	Cell Type Identification: Using gene expression patterns to classify cells into specific types or states.

In summary, the matrix encapsulates the gene expression data across all cells in the study, providing a foundation for understanding cellular behavior and biology at a granular level. If you have specific analyses in mind or need assistance with interpretation, let me know!

# Para la matriz Raw_Read_Count del paper original
# 
In a matrix where the rows represent gene names (like ‘Gtf2a2’, ‘Wwp1’, ‘Pak2’, etc.) and the columns represent conditions or experimental groups (like ‘WT’ for wild-type and ‘Lmo4OE’ for Lmo4 overexpression), the values in the matrix typically represent the expression levels of the corresponding genes under each condition. Here’s how to interpret the matrix:

	1.	Rows: Each row corresponds to a specific gene. For example, the first row might represent the expression level of the gene ‘Gtf2a2’.
	2.	Columns: Each column represents a condition or treatment. In this case, ‘WT’ refers to the wild-type condition, and ‘Lmo4OE’ refers to the condition where the Lmo4 gene is overexpressed.
	3.	Values:
	•	The values in the matrix can represent various types of measurements related to gene expression, including:
	•	Raw expression counts: The number of RNA transcripts detected for each gene.
	•	Normalized expression levels: Adjusted counts that account for differences in sequencing depth or other experimental variables, often expressed in counts per million (CPM) or fragments per kilobase of transcript per million mapped reads (FPKM).
	•	Fold changes: If the values are calculated as the ratio of expression levels between the two conditions, they may indicate how much more or less a gene is expressed in the ‘Lmo4OE’ condition compared to ‘WT’.
	•	Statistical significance: In some cases, the values might represent p-values or other statistical metrics indicating the significance of the differences in expression between the two conditions.

Example Interpretation

	•	Gene Gtf2a2:
	•	WT: 100
	•	Lmo4OE: 250
In this example, the expression of Gtf2a2 is significantly higher in the Lmo4OE condition compared to the WT condition. If you were to calculate the fold change, it would be 2.5, indicating that Gtf2a2 is 2.5 times more expressed in the Lmo4OE condition.

Summary

The matrix provides a comparative view of gene expression across different conditions, allowing researchers to analyze the effects of the Lmo4 overexpression on the expression of specific genes. This can help in understanding the biological processes and pathways influenced by Lmo4.