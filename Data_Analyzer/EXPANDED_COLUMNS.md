# 📊 COLUMNAS EXPANDIDAS PROPUESTAS

## Categoría: DETALLES EXPERIMENTALES

### 1. **Condiciones de Toma de Muestras**
- `Sample_Collection_Conditions` - Temperatura, luz, humedad, presión, pH, etc.
- `Temperature_Range` - Rango de temperatura (ej: "22-25°C", "37°C", "Room temperature")
- `Light_Conditions` - Fotoperiodo y intensidad (ej: "16h light/8h dark, 350 μmol·m−2·s−1")
- `Growth_Medium` - Sustrato/medio (ej: "soil mixture 3:1", "Murashige and Skoog medium", "LB agar")
- `Environmental_Stress` - Condiciones especiales de estrés ambiental

### 2. **Edad / Estadio de Desarrollo**
- `Organism_Age` - Edad del organismo (ej: "14 days post-germination", "8-week-old")
- `Developmental_Stage` - Estadio específico (ej: "flowering", "seedling", "vegetative", "larval stage 3")
- `Organism_Stage_Details` - Detalles adicionales (ej: "3-leaf stage", "boot stage", "adult")

### 3. **Moléculas Extraídas / Analizadas**
- `Molecules_Extracted` - Qué se extrajo (puede ser múltiple: DNA;RNA;Protein;Metabolites;Lipids)
- `DNA_Type` - Si aplica (ej: "genomic DNA", "mtDNA", "plasmid DNA")
- `RNA_Type` - Si aplica (ej: "total RNA", "mRNA", "small RNA", "miRNA")
- `Protein_Fraction` - Si aplica (ej: "total protein", "soluble protein", "membrane protein")
- `Other_Molecules` - Si aplica (ej: "metabolites", "lipids", "secondary metabolites", "hormones")

### 4. **Diseño Temporal (Time Course)**
- `Time_Course_Design` - ¿Es un estudio temporal? (Yes/No)
- `Time_Points` - Número de puntos temporales
- `Time_Intervals` - Intervalo de muestreo (ej: "every 24h", "every 6h", "0, 2, 4, 6, 12, 24 hours")
- `Time_Duration` - Duración total del experimento (ej: "7 days", "48 hours")

### 5. **Diseño Replicacional**
- `Sample_Size` - Número de replicados (ej: "3 biological replicates")
- `Biological_Replicates` - Número (ej: "3")
- `Technical_Replicates` - Número (ej: "2")
- `Experimental_Design` - Tipo de diseño (ej: "factorial", "blocked", "randomized")

### 6. **Tratamientos / Comparaciones**
- `Treatment_Groups` - ¿Cuántos grupos de tratamiento? (ej: "4 groups")
- `Control_Type` - Tipo de control (ej: "untreated control", "mock treatment", "wild-type")
- `Dose_Range` - Rango de dosis si aplica (ej: "0-100 μM", "0.5-5 mg/kg")

### 7. **Herramientas de Análisis / Medición**
- `Measurement_Tools` - Instrumentos específicos usados (ej: "Illumina TrueSeq", "MALDI-TOF", "HPLC")
- `Detection_Method` - Método de detección (ej: "qRT-PCR on ABI 7500", "LC-MS/MS")
- `Sequencing_Depth` - Si aplica (ej: "100x coverage", "50 million reads")
- `Resolution` - Si aplica (ej: "0.1 μm", "1 second intervals")

### 8. **Especificidad del Organismo (más allá de especie)**
- `Species` - Especie (ej: "Solanum lycopersicum", "Homo sapiens", "Escherichia coli")
- `Strain_Variety` - Cepa/variedad específica (ej: "Micro-Tom cultivar", "A549 cell line", "K12 strain")
- `Genotype` - Genotipo si es relevante (ej: "wild-type", "SlAREB knockout", "CRISPR-edited")
- `Source_Tissue_Origin` - Origen del tejido si es cultivo celular (ej: "leaf mesophyll", "root apex", "cortex tissue")

### 9. **Parámetros de Calidad**
- `Quality_Metrics` - Métricas de calidad reportadas (ej: "RNA integrity number >8", "Protein purity >95%")
- `Contamination_Check` - ¿Se verificó contaminación? (ej: "mycoplasma screening negative")

### 10. **Información Biológica Adicional**
- `Pathway_Focus` - ¿Qué vías/procesos fueron estudiados? (ej: "stress response", "photosynthesis", "immune response")
- `Biomarkers_Measured` - Biomarcadores específicos (ej: "chlorophyll content", "enzyme activity", "cytokine levels")
- `Differential_Expression_Threshold` - Qué umbral se usó (ej: "log2FC > 1.5, p < 0.05")

---

## 📋 PROPUESTA DE ORDEN DE COLUMNAS (CSV)

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

---

## 🎯 INSTRUCCIONES PARA EL LLM

Agregado al prompt que le enviamos a Ollama:

```
ADDITIONAL EXPERIMENTAL DETAILS EXTRACTION:

6. Sample Collection Conditions: How were samples collected?
   - Temperature during collection (e.g., "22-25°C", "room temperature")
   - Light conditions (fotoperiod, intensity)
   - Growth medium or culture conditions
   - Time of day collected, if specified

7. Organism Age / Developmental Stage: At what stage/age were samples taken?
   - Specific age (e.g., "14 days post-germination", "8-week-old", "log phase bacteria")
   - Developmental stage (e.g., "flowering", "seedling stage 3", "exponential phase")
   - Include botanical stage names if plants (e.g., "boot stage", "pre-anthesis")
   - Include cell cycle phase if applicable
   - Include life stage if animal (larval, adult, post-metamorphosis, etc)

8. Molecules Extracted/Analyzed: What biological molecules were extracted?
   - Can be multiple: DNA (type?), RNA (type?), Protein (type?), Metabolites, Lipids, etc.
   - Be specific: "total RNA OR mRNA OR small RNA OR miRNA"
   - "total protein OR soluble protein OR membrane protein"
   - "genomic DNA OR mtDNA OR plasmid DNA"

9. Temporal Design: Was this a time-course study?
   - Yes/No for time-course
   - Number of time points
   - Sampling intervals (e.g., "every 24 hours", "0, 2, 4, 6, 12, 24 hours")
   - Total duration

10. Replication Design: How many replicates?
    - Biological replicates (n=?)
    - Technical replicates
    - Type of experimental design

11. Treatments/Comparisons:
    - Number of treatment groups
    - Control type (untreated, mock, wild-type, etc)
    - Dose ranges if applicable

12. Measurement Tools:
    - Specific instruments (e.g., "Illumina TrueSeq", "qRT-PCR on ABI 7500")
    - Sequencing depth if applicable
    - Detection resolution

13. Organism Specificity:
    - Strain/cultivar/cell line (e.g., "Micro-Tom cultivar", "A549 cell line")
    - Genotype if relevant (wild-type, knockout, transgenic)

RESPONSE FORMAT (JSON):
{
  "relevance_score": <0-10>,
  "organisms": ["organism1", ...],
  "species": "Solanum lycopersicum",
  "strain_variety": "Micro-Tom cultivar",
  "genotype": "wild-type OR transgenic OR knockout",
  "developmental_stage": "3-leaf stage",
  "organism_age": "14 days post-germination",
  "tissues_organs": ["root", "leaf", ...],
  "source_tissue_origin": "leaf mesophyll",
  "conditions": ["drought stress", ...],
  "temperature_range": "22-25°C",
  "light_conditions": "16h light/8h dark, 350 μmol·m−2·s−1",
  "growth_medium": "soil mixture 3:1",
  "sample_collection_conditions": "collected at noon, ambient conditions",
  "molecules_extracted": ["RNA", "Protein"],
  "rna_type": "total RNA ; mRNA",
  "protein_type": "total protein",
  "dna_type": "not mentioned",
  "other_molecules": "metabolites",
  "time_course_design": "Yes",
  "time_points": 6,
  "time_intervals": "every 24 hours",
  "time_duration": "144 hours",
  "sample_size": "3 biological replicates",
  "biological_replicates": 3,
  "technical_replicates": 2,
  "treatment_groups": 4,
  "control_type": "untreated control",
  "strategies": ["RNA-Seq", "qRT-PCR"],
  "measurement_tools": "Illumina TrueSeq",
  "detection_method": "next-generation sequencing",
  "quality_metrics": "RNA integrity number > 8",
  "pathway_focus": "stress response pathway",
  "biomarkers_measured": "gene expression levels",
  "not_mentioned_fields": ["dose_range", "contamination_check", "sequencing_depth"]
}

IMPORTANT:
- Return ONLY valid JSON
- If information is NOT FOUND in the paper, use "not described" or "not mentioned"
- NEVER INVENT data - if you don't find it, say so
- Use exact quotes from Methods/Results when available
- Multiple values separated by " ; " (semicolon with spaces)
- Be precise with units (e.g., "μmol·m−2·s−1", "°C", "hours")
- For animals: include life stage (larval, adult, embryonic stage X, etc)
- For microorganisms: include growth phase (exponential, stationary, lag phase)
- For cell cultures: include cell line name and origin tissue if known
- For tissues/organs: be specific (not just "root" but "root meristem" if described)
```

---

## 💡 COLUMNAS ADICIONALES SUGERIDAS (Bonus)

Si queremos ir aún más lejos, podríamos agregar:

- `NCBI_Accession` - Acceso en NCBI (ej: "PRJNA123456")
- `GEO_Accession` - Acceso en GEO (ej: "GSE12345")
- `Raw_Data_Available` - ¿Datos crudos disponibles? (Yes/No)
- `Normalization_Method` - Método de normalización usado
- `Statistical_Method` - TestS estadísticos usado
- `Co_Authors_Count` - Número de autores
- `PubMed_Central_Link` - URL directo a PMC
- `Tissue_Specificity` - ¿Es tissue-specific? (Yes/No)
- `Disease_Model` - Si aplica (ej: "Type 2 Diabetes", "Alzheimer's model")
- `Animal_Model_Type` - Si es animal (ej: "mouse", "zebrafish", "C. elegans")
- `Cell_Type` - Si es cultivo celular (ej: "HEK293 cells", "primary neurons")

---

## 📈 COMPARATIVA: ANTES vs DESPUÉS

### ANTES (Actual)
```
PMID,PMCID,Title,Relevance_Score,Is_Relevant,Organisms,Tissues,Conditions,Strategies,Year,Journal,DOI,Abstract_Preview

41068586,PMC12513158,"Transcriptome-based...",10,Yes,"Solanum lycopersicum ; Solanum tuberosum","root","drought stress ; water",,"RNA-Seq ; qRT-PCR",2025,"BMC plant biology",...
```

### DESPUÉS (Expandido)
```
PMID,PMCID,Title,Relevance_Score,Is_Relevant,Species,Strain,Genotype,Developmental_Stage,Organism_Age,Tissues,Conditions,Temperature,Light,Growth_Medium,Molecules_Extracted,RNA_Type,Protein_Type,Time_Course,Time_Points,Time_Intervals,Sample_Size,Bio_Replicates,Treatment_Groups,Control_Type,Strategies,Measurement_Tools,Quality_Metrics,Pathway_Focus,Year,Journal,DOI,Abstract_Preview

41068586,PMC12513158,"Transcriptome-based...",10,Yes,"Solanum lycopersicum ; Solanum tuberosum","not described","not described","seedling ; vegetative","14 days post-germination","root","drought stress ; water shortage ; heat","25°C","16h light/8h dark, 350 μmol","soil mixture 3:1","RNA ; Protein","total RNA ; mRNA","total protein","Yes",6,"every 24-48 hours","not described",3,"untreated control ; drought 50% ; drought 75%","RNA-Seq ; qRT-PCR ; microarray","Illumina TrueSeq ; ABI 7500","RNA integrity number > 8","drought response ; water stress tolerance",2025,"BMC plant biology",...
```

---

## 🔍 EJEMPLOS DE EXTRACCIÓN POR TIPO DE ORGANISMO

### 📌 Ejemplo 1: PLANTA (Tomato)
```
Species: Solanum lycopersicum
Strain_Variety: Micro-Tom cultivar
Genotype: wild-type
Developmental_Stage: 3-leaf stage ; flowering stage
Organism_Age: 14 days post-germination ; 8 weeks old
Tissues_Organs: root ; leaf ; stem
Source_Tissue_Origin: not described
Conditions: drought stress ; 50% soil water deficit
Temperature_Range: 22-25°C
Light_Conditions: 16h light/8h dark, 350 μmol·m−2·s−1
Growth_Medium: soil mixture (peat:perlite 3:1)
Molecules_Extracted: RNA ; protein
RNA_Type: total RNA ; mRNA
Protein_Type: soluble protein
Time_Course_Design: Yes
Time_Points: 6
Time_Intervals: every 24 hours
Time_Duration: 144 hours (6 days)
Sample_Size: 3 biological replicates ; 2 technical replicates
Treatment_Groups: 3 (control ; mild drought ; severe drought)
Control_Type: well-watered control (100% field capacity)
Strategies: RNA-Seq ; qRT-PCR validation ; enzyme assays
Measurement_Tools: Illumina NovaSeq 6000 ; qRT-PCR on ABI 7500
```

### 📌 Ejemplo 2: ANIMAL (Mouse)
```
Species: Mus musculus
Strain_Variety: C57BL/6 mice
Genotype: wild-type
Developmental_Stage: adult
Organism_Age: 8-12 weeks old
Tissues_Organs: brain ; liver ; kidney
Source_Tissue_Origin: hippocampus ; hepatic tissue
Conditions: Alzheimer's disease model ; aging
Temperature_Range: ~37°C (body temperature)
Light_Conditions: 12h light/12h dark cycle
Growth_Medium: in vivo (whole organism)
Molecules_Extracted: RNA ; protein ; DNA
RNA_Type: total RNA ; mRNA
Protein_Type: membrane protein ; cytoplasmic fraction
DNA_Type: genomic DNA
Time_Course_Design: Yes
Time_Points: 4 (early, mid, late, terminal stages)
Time_Intervals: every 3 months
Time_Duration: 12 months
Sample_Size: n=6-8 per group
Biological_Replicates: 6
Technical_Replicates: 2
Treatment_Groups: 2 (control ; AD model)
Control_Type: age-matched wild-type littermates
Strategies: RNA-Seq ; Western blot ; mass spectrometry ; behavioral tests
Measurement_Tools: Illumina HiSeq ; MALDI-TOF-MS
Tissue_Specificity: Yes (brain-specific analysis)
Disease_Model: Transgenic Alzheimer's model (5xFAD)
```

### 📌 Ejemplo 3: BACTERIA
```
Species: Escherichia coli
Strain_Variety: K12 strain ; MG1655
Genotype: wild-type ; lacZ knockout
Developmental_Stage: exponential phase ; stationary phase
Organism_Age: 6-8 hours from inoculation (exponential) ; 24 hours (stationary)
Tissues_Organs: whole cell
Growth_Medium: M9 minimal medium ; LB broth
Conditions: iron limitation ; heat shock 42°C
Temperature_Range: 37°C (control) ; 42°C (heat shock)
Light_Conditions: anaerobic conditions
Time_Course_Design: Yes
Time_Points: 5
Time_Intervals: every 30 minutes
Time_Duration: 150 minutes
Sample_Size: 3 biological replicates
Biological_Replicates: 3
Technical_Replicates: 2
Treatment_Groups: 2 (control ; heat shock)
Control_Type: no treatment, 37°C
Strategies: RNA-Seq ; qRT-PCR ; ChIP-Seq ; metabolomics
Measurement_Tools: Illumina NextSeq 500 ; Agilent HPLC-MS
Quality_Metrics: >95% viability (CFU counting)
Contamination_Check: PCR screening for plasmid contamination (negative)
Growth_Phase: exponential phase (OD600 0.4-0.6)
```

### 📌 Ejemplo 4: CULTIVO CELULAR
```
Species: Homo sapiens
Strain_Variety: HEK293T cells (human embryonic kidney)
Genotype: not described
Tissues_Organs: not applicable (cultured cells)
Source_Tissue_Origin: embryonic kidney tissue
Cell_Type: HEK293T (human embryonic kidney cells)
Conditions: serum starvation ; growth factor stimulation (EGF)
Temperature_Range: 37°C
Light_Conditions: standard incubator (no light requirement)
Growth_Medium: DMEM supplemented with 10% FBS
Molecules_Extracted: RNA ; protein ; metabolites
RNA_Type: total RNA ; mRNA
Protein_Type: nuclear protein ; cytoplasmic protein
Other_Molecules: phosphorylated metabolites
Time_Course_Design: Yes
Time_Points: 6
Time_Intervals: every 15 minutes
Time_Duration: 90 minutes
Sample_Size: 3 biological replicates
Biological_Replicates: 3 (from different culture passages)
Technical_Replicates: 2
Treatment_Groups: 2 (control serum-free ; EGF-stimulated)
Control_Type: DMEM without serum, no growth factors
Dose_Range: EGF 10-100 ng/mL
Strategies: RNA-Seq ; phosphoproteomics ; ELISA ; Western blot
Measurement_Tools: Illumina MiSeq ; ESI-MS/MS
Quality_Metrics: >99% cell purity (flow cytometry)
Contamination_Check: mycoplasma screening (negative), bacterial culture (negative)
```

---

## ✅ RESUMEN: QUÉ SE EXTRAE AHORA

| Aspecto | Columnas | Ejemplos |
|---------|----------|----------|
| **Identificación** | PMID, PMCID, Title, Year, Journal, DOI | 41068586, PMC12513158, ... |
| **Relevancia** | Relevance_Score, Is_Relevant | 10, Yes |
| **Organismo Principal** | Species, Strain, Genotype, Tissue, Cell_Type | Solanum lycopersicum, Micro-Tom, wild-type, root |
| **Timing Biológico** | Developmental_Stage, Organism_Age | seedling, 14 days post-germination |
| **Condiciones** | Temperature, Light, Growth_Medium, Environmental_Stress | 25°C, 16h light, soil, drought |
| **Moléculas** | Molecules_Extracted, RNA_Type, Protein_Type, DNA_Type | RNA, mRNA, total protein, genomic DNA |
| **Tiempo (Course)** | Time_Course_Design, Time_Points, Time_Intervals, Duration | Yes, 6, every 24h, 144h |
| **Replicación** | Sample_Size, Bio_Replicates, Tech_Replicates, Replication_Design | 3, 3, 2, factorial |
| **Tratamientos** | Treatment_Groups, Control_Type, Dose_Range | 3, untreated control, 0-100 μM |
| **Herramientas** | Strategies, Measurement_Tools, Detection_Method | RNA-Seq, Illumina TrueSeq, LC-MS/MS |
| **Calidad** | Quality_Metrics, Contamination_Check | RNA integrity > 8, mycoplasma negative |
| **Enfoque Biológico** | Pathway_Focus, Biomarkers_Measured, Disease_Model | stress response, gene expression, Alzheimer's |
| **Preview** | Abstract_Preview | "Plants possess..." |

