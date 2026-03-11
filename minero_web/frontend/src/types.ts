export type SourceMode = 'pubmed' | 'pmc' | 'bioproject';

export type RelevanceLabel = 'alta' | 'media' | 'baja';
export type EvidenceLevel = 'directa' | 'indirecta' | 'débil';
export type ModelSource = 'ollama' | 'heuristic';

export interface Classification {
  relevance_label: RelevanceLabel;
  relevance_score: number;
  reason_short: string;
  tags: string[];
  evidence_level: EvidenceLevel;
  model_source: ModelSource;
}

export interface QueryGeneration {
  user_input: string;
  extracted: {
    organism?: string | null;
    organism_variants?: string[];
    strategies?: string[];
    tissues?: string[];
    conditions?: string[];
    free_terms?: string[];
  };
  synonyms?: {
    organism?: string[];
    strategies?: string[];
    tissues?: string[];
    conditions?: string[];
  };
  ncbi_query: string;
  ready_to_use: boolean;
  clarification_needed: boolean;
  clarification_message: string;
}

export interface MineroMetadata {
  status: 'success' | 'partial-success' | 'empty';
  query: string;
  source: SourceMode;
  total_results: number;
  classification_version: string;
  classification_timestamp: string;
  llm_runtime_available: boolean;
  model_default: ModelSource;
}

export interface PubmedResult {
  pmid: string;
  title: string;
  abstract: string;
  authors?: string[];
  year?: string;
  journal?: string;
  publication_type?: string;
  doi?: string;
  pmcid?: string;
  url?: string;
  fetched_at?: string;
  classification: Classification;
}

export interface BioprojectResult {
  bioproject: string;
  title: string;
  submission_date?: string;
  organism?: string;
  project_type?: string;
  description?: string;
  sra_experiments_count?: number;
  biosamples_count?: number;
  sra_runs_count?: number;
  sra_hierarchy?: Record<string, unknown>;
  publications_found?: number;
  search_method?: string;
  papers_summary?: string;
  error?: string;
  classification: Classification;
}

export type MineroResult = PubmedResult | BioprojectResult;

export interface MineroResponse {
  metadata: MineroMetadata;
  query_generation: QueryGeneration;
  results: MineroResult[];
}

export interface SearchHistoryEntry {
  id: string;
  created_at: string;
  source: SourceMode;
  natural_query: string;
  ncbi_query: string;
  status: MineroMetadata['status'];
  total_results: number;
  response: MineroResponse;
}

export interface SearchPayload {
  natural_query: string;
  source: SourceMode;
  max_results: number;
  use_llm: boolean;
}

export type AppView = 'buscar' | 'resultados' | 'analisis' | 'ayuda';
export type AppStatus = 'idle' | 'loading' | 'success' | 'partial-success' | 'empty' | 'error';

export type ProgressState = 'pending' | 'in_progress' | 'completed' | 'failed';

export interface ProgressStep {
  id: number;
  label: string;
  detail: string;
  state: ProgressState;
}

// ─── Modo Pro — PaperAnalyzer (53 columnas) ─────────────────────────────────

export interface ProAnalysisResult {
  PMID: string;
  PMCID?: string;
  Title: string;
  Year?: string;
  Journal?: string;
  DOI?: string;
  Relevance_Score: number;
  Is_Relevant: 'Yes' | 'No' | 'Error';
  Summary: string;
  Relevance_Explanation: string;
  Organisms: string;
  Species: string;
  Strain_Variety: string;
  Genotype: string;
  Tissues_Organs: string;
  Source_Tissue_Origin: string;
  Cell_Type: string;
  Developmental_Stage: string;
  Organism_Age: string;
  Growth_Phase: string;
  Conditions: string;
  Environmental_Stress: string;
  Temperature_Range: string;
  Light_Conditions: string;
  Growth_Medium: string;
  Sample_Collection_Conditions: string;
  Molecules_Extracted: string;
  RNA_Type: string;
  DNA_Type: string;
  Protein_Type: string;
  Other_Molecules: string;
  Strategies: string;
  Measurement_Tools: string;
  Detection_Method: string;
  Time_Course_Design: string;
  Time_Points: string;
  Time_Intervals: string;
  Time_Duration: string;
  Sample_Size: string;
  Biological_Replicates: string;
  Technical_Replicates: string;
  Replication_Design: string;
  Treatment_Groups: string;
  Control_Type: string;
  Dose_Range: string;
  Quality_Metrics: string;
  Contamination_Check: string;
  Pathway_Focus: string;
  Biomarkers_Measured: string;
  Disease_Model: string;
  Normalization_Method: string;
  Statistical_Method: string;
  Differential_Expression_Threshold: string;
  Raw_Data_Available: string;
  Abstract_Preview: string;
}

export interface ProAnalysisMetadata {
  total_analyzed: number;
  total_relevant: number;
  total_errors: number;
  pmc_full_text_used: number;
  abstract_only: number;
  techniques_enhanced: number;
  model: string;
  analyzed_at: string;
}

export interface ProAnalysisResponse {
  metadata: ProAnalysisMetadata;
  results: ProAnalysisResult[];
  network_html?: string;
}
