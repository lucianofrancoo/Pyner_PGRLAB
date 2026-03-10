import type { BioprojectResult, MineroResponse, MineroResult, PubmedResult, SourceMode } from '../types';

function downloadBlob(data: BlobPart, filename: string, mimeType: string): void {
  const blob = new Blob([data], { type: mimeType });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;
  anchor.click();
  URL.revokeObjectURL(url);
}

export function exportJson(response: MineroResponse): void {
  const stamp = new Date().toISOString().replace(/[:.]/g, '-');
  const source = response.metadata.source;
  downloadBlob(JSON.stringify(response, null, 2), `minero_${source}_${stamp}.json`, 'application/json');
}

function csvEscape(value: unknown): string {
  const text = String(value ?? '');
  if (text.includes(',') || text.includes('"') || text.includes('\n')) {
    return `"${text.replace(/"/g, '""')}"`;
  }
  return text;
}

function buildPubmedRow(item: PubmedResult): Record<string, string | number> {
  const classification = item.classification;

  return {
    pmid: item.pmid,
    title: item.title,
    year: item.year ?? '',
    journal: item.journal ?? '',
    publication_type: item.publication_type ?? '',
    relevance_label: classification.relevance_label,
    relevance_score: classification.relevance_score,
    tags: Array.isArray(classification.tags) ? classification.tags.join(';') : '',
    evidence_level: classification.evidence_level,
    model_source: classification.model_source,
    reason_short: classification.reason_short,
  };
}

function buildBioprojectRow(item: BioprojectResult): Record<string, string | number> {
  const classification = item.classification;

  return {
    bioproject: item.bioproject,
    title: item.title,
    organism: item.organism ?? '',
    sra_experiments_count: item.sra_experiments_count ?? '',
    biosamples_count: item.biosamples_count ?? '',
    relevance_label: classification.relevance_label,
    relevance_score: classification.relevance_score,
    tags: Array.isArray(classification.tags) ? classification.tags.join(';') : '',
    evidence_level: classification.evidence_level,
    model_source: classification.model_source,
    reason_short: classification.reason_short,
  };
}

function getTagValue(tags: string[] | undefined, prefix: string): string {
  const tag = (tags ?? []).find((entry) => entry.startsWith(prefix));
  if (!tag) return '';
  return tag.slice(prefix.length).trim();
}

type SraRunRow = {
  bioproject: string;
  project_title: string;
  organism: string;
  biosample: string;
  experiment_id: string;
  experiment_title: string;
  run_accession: string;
  spots: number | string;
  size: string;
  library_name: string;
  library_strategy: string;
  library_source: string;
  library_selection: string;
  library_layout: string;
  instrument: string;
  tissue: string;
  condition: string;
};

function flattenBioprojectSraRows(item: BioprojectResult): SraRunRow[] {
  const tags = item.classification?.tags ?? [];
  const tissue = getTagValue(tags, 'tejido:');
  const condition = getTagValue(tags, 'condicion:');
  const hierarchy = item.sra_hierarchy;

  if (!hierarchy || typeof hierarchy !== 'object') {
    return [];
  }

  const rows: SraRunRow[] = [];
  const biosampleEntries = Object.entries(hierarchy as Record<string, unknown>);

  biosampleEntries.forEach(([biosampleId, biosampleValue]) => {
    if (!biosampleValue || typeof biosampleValue !== 'object') return;

    const biosampleObj = biosampleValue as Record<string, unknown>;
    const experiments = Array.isArray(biosampleObj.experiments) ? biosampleObj.experiments : [];

    experiments.forEach((expValue) => {
      if (!expValue || typeof expValue !== 'object') return;

      const experiment = expValue as Record<string, unknown>;
      const metadata =
        experiment.metadata && typeof experiment.metadata === 'object'
          ? (experiment.metadata as Record<string, unknown>)
          : {};
      const runs = Array.isArray(experiment.runs) ? experiment.runs : [];

      if (!runs.length) {
        rows.push({
          bioproject: item.bioproject,
          project_title: item.title ?? '',
          organism: item.organism ?? '',
          biosample: biosampleId,
          experiment_id: String(experiment.experiment_id ?? ''),
          experiment_title: String(experiment.title ?? ''),
          run_accession: '',
          spots: '',
          size: '',
          library_name: String(metadata.library_name ?? ''),
          library_strategy: String(metadata.library_strategy ?? ''),
          library_source: String(metadata.library_source ?? ''),
          library_selection: String(metadata.library_selection ?? ''),
          library_layout: String(metadata.library_layout ?? ''),
          instrument: String(metadata.instrument ?? ''),
          tissue,
          condition,
        });
        return;
      }

      runs.forEach((runValue) => {
        if (!runValue || typeof runValue !== 'object') return;
        const run = runValue as Record<string, unknown>;
        rows.push({
          bioproject: item.bioproject,
          project_title: item.title ?? '',
          organism: item.organism ?? '',
          biosample: biosampleId,
          experiment_id: String(experiment.experiment_id ?? ''),
          experiment_title: String(experiment.title ?? ''),
          run_accession: String(run.accession ?? ''),
          spots: Number(run.spots ?? 0),
          size: String(run.size ?? ''),
          library_name: String(metadata.library_name ?? ''),
          library_strategy: String(metadata.library_strategy ?? ''),
          library_source: String(metadata.library_source ?? ''),
          library_selection: String(metadata.library_selection ?? ''),
          library_layout: String(metadata.library_layout ?? ''),
          instrument: String(metadata.instrument ?? ''),
          tissue,
          condition,
        });
      });
    });
  });

  return rows;
}

export function exportCsv(results: MineroResult[], source: SourceMode): void {
  if (!results.length) {
    return;
  }

  const rows =
    source === 'pubmed'
      ? results.filter((item): item is PubmedResult => 'pmid' in item).map(buildPubmedRow)
      : results.filter((item): item is BioprojectResult => 'bioproject' in item).map(buildBioprojectRow);
  if (!rows.length) {
    return;
  }
  const headers = Object.keys(rows[0]);
  const csv = [headers.join(','), ...rows.map((row) => headers.map((h) => csvEscape(row[h])).join(','))].join('\n');

  const stamp = new Date().toISOString().replace(/[:.]/g, '-');
  downloadBlob(csv, `minero_${source}_${stamp}.csv`, 'text/csv;charset=utf-8');
}

export function exportBioprojectSraTsv(results: MineroResult[]): void {
  const bioprojectRows = results.filter((item): item is BioprojectResult => 'bioproject' in item);
  if (!bioprojectRows.length) return;

  const sraRows = bioprojectRows.flatMap(flattenBioprojectSraRows);
  if (!sraRows.length) return;

  const headers = Object.keys(sraRows[0]);
  const tsv = [headers.join('\t'), ...sraRows.map((row) => headers.map((h) => String(row[h as keyof SraRunRow] ?? '')).join('\t'))].join('\n');

  const stamp = new Date().toISOString().replace(/[:.]/g, '-');
  downloadBlob(tsv, `minero_bioproject_sra_libraries_${stamp}.tsv`, 'text/tab-separated-values;charset=utf-8');
}

export function exportBioprojectMetadataTsv(response: MineroResponse): void {
  if (response.metadata.source !== 'bioproject') return;

  const projects = response.results.filter((item): item is BioprojectResult => 'bioproject' in item);
  if (!projects.length) return;

  const payload = {
    exported_at: new Date().toISOString(),
    query: response.metadata.query,
    source: response.metadata.source,
    total_results: projects.length,
    projects: projects.map((item) => ({
      bioproject: item.bioproject,
      title: item.title,
      description: item.description ?? '',
      organism: item.organism ?? '',
      submission_date: item.submission_date ?? '',
      project_type: item.project_type ?? '',
      sra_experiments_count: item.sra_experiments_count ?? 0,
      biosamples_count: item.biosamples_count ?? 0,
      sra_runs_count: item.sra_runs_count ?? 0,
      publications_found: item.publications_found ?? 0,
      search_method: item.search_method ?? '',
      papers_summary: item.papers_summary ?? '',
      classification: item.classification,
      sra_hierarchy: item.sra_hierarchy ?? {},
    })),
  };

  const stamp = new Date().toISOString().replace(/[:.]/g, '-');
  const rows = payload.projects.map((project) => ({
    bioproject: project.bioproject,
    title: project.title,
    description: project.description,
    organism: project.organism,
    submission_date: project.submission_date,
    project_type: project.project_type,
    sra_experiments_count: project.sra_experiments_count,
    biosamples_count: project.biosamples_count,
    sra_runs_count: project.sra_runs_count,
    publications_found: project.publications_found,
    search_method: project.search_method,
    papers_summary: project.papers_summary,
    relevance_label: project.classification?.relevance_label ?? '',
    relevance_score: project.classification?.relevance_score ?? '',
    evidence_level: project.classification?.evidence_level ?? '',
    model_source: project.classification?.model_source ?? '',
    reason_short: project.classification?.reason_short ?? '',
    tags: Array.isArray(project.classification?.tags) ? project.classification.tags.join(';') : '',
  }));
  if (!rows.length) return;

  const headers = Object.keys(rows[0]);
  const tsv = [headers.join('\t'), ...rows.map((row) => headers.map((h) => String(row[h as keyof typeof row] ?? '')).join('\t'))].join('\n');
  downloadBlob(tsv, `minero_bioproject_metadata_${stamp}.tsv`, 'text/tab-separated-values;charset=utf-8');
}
