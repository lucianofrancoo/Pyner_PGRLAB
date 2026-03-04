# 📊 PROPOSED EXPANDED COLUMNS (CSV)

## Category: EXPERIMENTAL DETAILS

### 1. **Sample Collection Conditions**
- `Sample_Collection_Conditions`: Temperature, light, humidity, pressure, pH, etc.
- `Temperature_Range`: Temperature range (e.g., "22-25°C", "37°C", "Room temperature")
- `Light_Conditions`: Photoperiod and intensity (e.g., "16h light/8h dark, 350 μmol·m−2·s−1")
- `Growth_Medium`: Substrate/medium (e.g., "soil mixture 3:1", "Murashige and Skoog medium", "LB agar")
- `Environmental_Stress`: Special environmental stress conditions

### 2. **Age / Developmental Stage**
- `Organism_Age`: Age of the organism (e.g., "14 days post-germination", "8-week-old")
- `Developmental_Stage`: Specific stage (e.g., "flowering", "seedling", "vegetative", "larval stage 3")
- `Organism_Stage_Details`: Additional details (e.g., "3-leaf stage", "boot stage", "adult")

### 3. **Molecules Extracted / Analyzed**
- `Molecules_Extracted`: What was extracted (can be multiple: DNA; RNA; Protein; Metabolites; Lipids)
- `DNA_Type`: If applicable (e.g., "genomic DNA", "mtDNA", "plasmid DNA")
- `RNA_Type`: If applicable (e.g., "total RNA", "mRNA", "small RNA", "miRNA")
- `Protein_Fraction`: If applicable (e.g., "total protein", "soluble protein", "membrane protein")
- `Other_Molecules`: If applicable (e.g., "metabolites", "lipids", "secondary metabolites", "hormones")

### 4. **Temporal Design (Time Course)**
- `Time_Course_Design`: Is it a temporal study? (Yes/No)
- `Time_Points`: Number of time points
- `Time_Intervals`: Sampling interval (e.g., "every 24h", "every 6h", "0, 2, 4, 6, 12, 24 hours")
- `Time_Duration`: Total duration of the experiment (e.g., "7 days", "48 hours")

### 5. **Replication Design**
- `Sample_Size`: Total number of samples or animals
- `Biological_Replicates`: Number (e.g., "3")
- `Technical_Replicates`: Number (e.g., "2")
- `Experimental_Design`: Type of design (e.g., "factorial", "blocked", "randomized")

### 6. **Treatments / Comparisons**
- `Treatment_Groups`: How many treatment groups? (e.g., "4 groups")
- `Control_Type`: Type of control (e.g., "untreated control", "mock treatment", "wild-type")
- `Dose_Range`: Dose range if applicable (e.g., "0-100 μM", "0.5-5 mg/kg")

### 7. **Measurement / Analysis Tools**
- `Measurement_Tools`: Specific instruments used (e.g., "Illumina NovaSeq 6000", "MALDI-TOF", "HPLC")
- `Detection_Method`: Detection method (e.g., "qRT-PCR on ABI 7500", "LC-MS/MS")
- `Sequencing_Depth`: If applicable (e.g., "100x coverage", "50 million reads")

### 8. **Organism Specificity (beyond species)**
- `Species`: Species name (e.g., "Solanum lycopersicum", "Homo sapiens", "Escherichia coli")
- `Strain_Variety`: Specific strain/variety/cultivar (e.g., "Micro-Tom cultivar", "A549 cell line", "K12 strain")
- `Genotype`: Genotype if relevant (e.g., "wild-type", "SlAREB knockout", "CRISPR-edited")
- `Source_Tissue_Origin`: Tissue origin if it is a cell culture (e.g., "leaf mesophyll", "root apex")

### 9. **Quality Parameters**
- `Quality_Metrics`: Reported quality metrics (e.g., "RNA integrity number >8", "Protein purity >95%")
- `Contamination_Check`: Was contamination verified? (e.g., "mycoplasma screening negative")

### 10. **Additional Biological Information**
- `Pathway_Focus`: What pathways/processes were studied? (e.g., "stress response", "photosynthesis")
- `Biomarkers_Measured`: Specific biomarkers (e.g., "chlorophyll content", "enzyme activity")
- `Differential_Expression_Threshold`: What threshold was used (e.g., "log2FC > 1.5, p < 0.05")

---

## 📋 COLUMN ORDER PROPOSAL (CSV)

```
PMID | PMCID | Title | 
Year | Journal | DOI | 
Relevance_Score | Is_Relevant | 
Organisms | Species | Strain_Variety | Genotype | 
Tissues_Organs | Source_Tissue_Origin | 
Developmental_Stage | Organism_Age | 
Conditions | Environmental_Stress | 
Temperature_Range | Light_Conditions | Growth_Medium | 
Sample_Collection_Conditions | 
Molecules_Extracted | DNA_Type | RNA_Type | Protein_Fraction | Other_Molecules | 
Strategies | Measurement_Tools | Detection_Method | 
Time_Course_Design | Time_Points | Time_Intervals | Time_Duration | 
Sample_Size | Biological_Replicates | Technical_Replicates | 
Treatment_Groups | Control_Type | Dose_Range | 
Quality_Metrics | Contamination_Check | 
Pathway_Focus | Biomarkers_Measured | Differential_Expression_Threshold | 
Abstract_Preview
```
