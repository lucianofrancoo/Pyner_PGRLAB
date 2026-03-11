export function HelpView() {
  return (
    <section className="panel">
      <header className="panel-header">
        <h2>Documentation & Help</h2>
        <p>Quick guide for scientific users with minimal bioinformatics background.</p>
      </header>

      <div className="help-grid">
        <article>
          <h3>Recommended workflow</h3>
          <ol>
            <li>Write your biological question in natural language.</li>
            <li>Choose PubMed or BioProject based on your goal.</li>
            <li>Review the generated NCBI query.</li>
            <li>Inspect relevance-ranked classified results.</li>
            <li>Export CSV/JSON for sharing and reproducibility.</li>
          </ol>
        </article>

        <article>
          <h3>Relevance meaning</h3>
          <ul>
            <li><strong>high</strong>: strong query alignment.</li>
            <li><strong>medium</strong>: partial alignment, needs review.</li>
            <li><strong>low</strong>: weak alignment with the objective.</li>
          </ul>
        </article>

        <article>
          <h3>Minimal glossary</h3>
          <ul>
            <li><strong>BioProject</strong>: umbrella project record for sequencing studies.</li>
            <li><strong>SRA</strong>: repository of sequencing runs and metadata.</li>
            <li><strong>PubMed</strong>: database of scientific publications.</li>
            <li><strong>Direct evidence</strong>: record tightly aligned with your query.</li>
          </ul>
        </article>
      </div>
    </section>
  );
}
