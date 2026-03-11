import { CheckCircle2, Circle, LoaderCircle, XCircle } from 'lucide-react';
import type { ProgressStep } from '../types';

interface ProgressTrackerProps {
  steps: ProgressStep[];
}

function StepIcon({ state }: { state: ProgressStep['state'] }) {
  if (state === 'completed') return <CheckCircle2 size={16} />;
  if (state === 'in_progress') return <LoaderCircle size={16} className="spin" />;
  if (state === 'failed') return <XCircle size={16} />;
  return <Circle size={16} />;
}

export function ProgressTracker({ steps }: ProgressTrackerProps) {
  return (
    <section className="panel compact">
      <header className="panel-header">
        <h2>Query Progress</h2>
        <p>Workflow tracking inspired by `test_fetcher_integrator.sh`.</p>
      </header>

      <ol className="progress-list">
        {steps.map((step) => (
          <li key={step.id} className={`progress-item ${step.state}`}>
            <span className="progress-icon">
              <StepIcon state={step.state} />
            </span>
            <div>
              <strong>{step.label}</strong>
              <p>{step.detail}</p>
            </div>
          </li>
        ))}
      </ol>
    </section>
  );
}
