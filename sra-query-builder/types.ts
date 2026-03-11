
export interface SRAInputs {
  naturalQuery: string;
}

export interface GeneStats {
  condition: string;
  up: number;
  down: number;
}

export interface SRAResultEntry {
  bioproject: string;
  title: string;
  organism: string;
  tissue: string;
  conditions: string;
  isTimeSeries: boolean;
  strategy: string;
  libraryCount: number;
  experimentCount: number;
  timePoints: number[];
  cultivar: string;
  geneStats: GeneStats[];
}

export interface GeneratedQuery {
  natural_query: string;
  esearch_query: string;
  results: SRAResultEntry[];
}

export enum AppStatus {
  IDLE = 'IDLE',
  GENERATING = 'GENERATING',
  SUCCESS = 'SUCCESS',
  ERROR = 'ERROR'
}

export type AppView = 'SEARCH' | 'RESULTS' | 'ANALYSIS' | 'DOCS';
