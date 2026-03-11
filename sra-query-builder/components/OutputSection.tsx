
import React, { useState, useEffect, useRef } from 'react';
import { GeneratedQuery, AppStatus } from '../types';

interface Props {
  query: GeneratedQuery | null;
  status: AppStatus;
  error: string | null;
}

export const OutputSection: React.FC<Props> = ({ query, status, error }) => {
  const [copied, setCopied] = useState(false);
  const terminalRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    if (terminalRef.current) {
      terminalRef.current.scrollTop = terminalRef.current.scrollHeight;
    }
  }, [status, query]);

  const handleCopy = () => {
    if (query?.esearch_query) {
      navigator.clipboard.writeText(query.esearch_query);
      setCopied(true);
      setTimeout(() => setCopied(false), 2000);
    }
  };

  const getTime = () => {
    const d = new Date();
    return d.toLocaleTimeString([], { hour12: false, hour: '2-digit', minute: '2-digit', second: '2-digit' });
  };

  return (
    <div className="h-full flex flex-col p-6">
      <div className="flex items-center justify-between mb-4">
        <div className="flex items-center gap-2 text-white/40">
          <svg className="w-4 h-4" fill="none" viewBox="0 0 24 24" stroke="currentColor"><path strokeLinecap="round" strokeLinejoin="round" strokeWidth={2} d="M10 20l4-16m4 4l4 4-4 4M6 16l-4-4 4-4" /></svg>
          <h2 className="text-[11px] font-bold uppercase tracking-[0.2em]">Console Output</h2>
        </div>
        {query && (
          <button
            onClick={handleCopy}
            className="px-3 py-1 bg-white/5 border border-white/10 hover:border-emerald-500/50 text-emerald-500 rounded text-[10px] uppercase font-bold tracking-widest transition-all"
          >
            {copied ? 'Copied' : 'Export Query'}
          </button>
        )}
      </div>

      <div 
        ref={terminalRef}
        className="flex-1 bg-[#05070a] border border-white/5 rounded-lg p-6 font-mono text-xs overflow-y-auto scroll-smooth shadow-inner"
      >
        {status === AppStatus.IDLE && (
          <div className="text-neutral-600">
            <span className="text-emerald-500/50 mr-2">[{getTime()}] [INFO]</span> System initialized. Waiting for input...
          </div>
        )}

        {status === AppStatus.GENERATING && (
          <div className="space-y-1">
            <div className="text-emerald-500/80">
              <span className="text-neutral-600 mr-2">[{getTime()}] [INFO]</span> Initiating LLM context for sequence data mining...
            </div>
            <div className="text-blue-400/80 animate-pulse">
              <span className="text-neutral-600 mr-2">[{getTime()}] [PROC]</span> Mapping phylogenetic synonyms...
            </div>
          </div>
        )}

        {status === AppStatus.SUCCESS && query && (
          <div className="space-y-6">
            <div className="text-emerald-500">
              <span className="text-neutral-600 mr-2">[{getTime()}] [SUCCESS]</span> Generation successful. Search agent synchronized.
            </div>
            
            <div className="bg-[#0d121d] p-4 rounded border border-white/5 border-l-2 border-l-emerald-500">
              <div className="text-[10px] text-neutral-500 uppercase font-bold mb-2 tracking-widest">Natural Intelligence Summary</div>
              <div className="text-neutral-200 leading-relaxed italic">
                "{query.natural_query}"
              </div>
            </div>

            <div className="bg-[#0d121d] p-4 rounded border border-white/5 border-l-2 border-l-blue-500">
              <div className="text-[10px] text-neutral-500 uppercase font-bold mb-2 tracking-widest">Optimized E-Search String</div>
              <div className="text-blue-400 leading-relaxed break-words font-medium">
                {query.esearch_query}
              </div>
            </div>

            <div className="text-[10px] text-emerald-500/50 uppercase tracking-widest animate-pulse">
              {">"} Ready for new search parameters.
            </div>
          </div>
        )}

        {status === AppStatus.ERROR && (
          <div className="text-red-500 bg-red-500/5 p-4 rounded border border-red-500/20">
            <span className="font-bold mr-2">[{getTime()}] [ERROR]</span> {error}
          </div>
        )}
      </div>

      <div className="mt-4 flex items-center justify-between px-2">
         <div className="flex gap-4">
            <span className="text-[9px] text-neutral-600 font-bold uppercase tracking-widest flex items-center gap-1.5">
               <span className="w-1.5 h-1.5 rounded-full bg-emerald-500"></span> NCBI SRA
            </span>
            <span className="text-[9px] text-neutral-600 font-bold uppercase tracking-widest flex items-center gap-1.5">
               <span className="w-1.5 h-1.5 rounded-full bg-blue-500"></span> Gemini 2.0
            </span>
         </div>
         <div className="text-[9px] text-neutral-700 font-mono">
            RT_LATENCY: 242ms // THREADS: 12
         </div>
      </div>
    </div>
  );
};
