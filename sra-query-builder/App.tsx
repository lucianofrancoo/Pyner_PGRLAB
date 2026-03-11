
import React, { useState } from 'react';
import { SRAInputs, GeneratedQuery, AppStatus, AppView } from './types';
import { InputSection } from './components/InputSection';
import { OutputSection } from './components/OutputSection';
import { ResultsTable } from './components/ResultsTable';
import { AnalysisDashboard } from './components/AnalysisDashboard';
import { Documentation } from './components/Documentation';
import { 
  Search, 
  Database, 
  BarChart3, 
  FileText, 
  Github, 
  Cpu, 
  Box,
  AlertCircle
} from 'lucide-react';

const LabLogo = () => (
  <svg width="42" height="42" viewBox="0 0 100 100" fill="none" xmlns="http://www.w3.org/2000/svg" className="drop-shadow-[0_0_8px_rgba(168,85,247,0.4)]">
    <path d="M50 5L95 85H5L50 5Z" stroke="#7e22ce" strokeWidth="5" />
    <path d="M50 15C55 25 45 35 50 45C55 55 45 65 50 75" stroke="#a855f7" strokeWidth="4" strokeLinecap="round" />
    <circle cx="50" cy="50" r="8" fill="#10b981" />
  </svg>
);

const App: React.FC = () => {
  const [view, setView] = useState<AppView>('SEARCH');
  const [inputs, setInputs] = useState<SRAInputs>({
    naturalQuery: ''
  });

  const [status, setStatus] = useState<AppStatus>(AppStatus.IDLE);
  const [data, setData] = useState<GeneratedQuery | null>(null);
  const [error, setError] = useState<string | null>(null);

  // Simulación de datos biológicos para la interfaz
  const mockResult: GeneratedQuery = {
    natural_query: "Búsqueda expandida para Arabidopsis thaliana en experimentos de RNA-Seq relacionados con estrés por nitrógeno y sequía.",
    esearch_query: "(\"Arabidopsis thaliana\"[Organism] AND \"RNA-Seq\"[Strategy]) AND (\"Nitrogen\"[All Fields] OR \"Drought\"[All Fields])",
    results: Array.from({ length: 12 }, (_, i) => ({
      bioproject: `PRJNA${100000 + i}`,
      title: i % 2 === 0 ? "Transcriptomic response to low nitrogen availability" : "Drought stress dynamics in Arabidopsis seedlings",
      organism: "Arabidopsis thaliana",
      tissue: i % 3 === 0 ? "Leaf" : i % 3 === 1 ? "Root" : "Seedling",
      conditions: i % 2 === 0 ? "Nitrogen Starvation" : "Water Deficit",
      isTimeSeries: i % 4 === 0,
      strategy: "RNA-Seq",
      libraryCount: Math.floor(Math.random() * 100) + 10,
      experimentCount: Math.floor(Math.random() * 20) + 5,
      timePoints: [0, 12, 24, 48],
      cultivar: i % 2 === 0 ? "Col-0" : "Ler",
      geneStats: [
        { condition: "Stress", up: 450 + i * 10, down: 320 + i * 5 },
        { condition: "Control", up: 10, down: 15 }
      ]
    }))
  };

  const handleGenerate = async () => {
    if (inputs.naturalQuery.length < 5) return;
    setStatus(AppStatus.GENERATING);
    setError(null);
    
    // Simulación de tiempo de respuesta de la UI
    setTimeout(() => {
      setData(mockResult);
      setStatus(AppStatus.SUCCESS);
    }, 1200);
  };

  const NavItem = ({ id, icon: Icon, label }: { id: AppView, icon: any, label: string }) => (
    <button
      onClick={() => setView(id)}
      className={`group relative flex items-center gap-4 px-5 py-4 rounded-2xl transition-all duration-300 ${
        view === id 
          ? 'bg-emerald-500 text-white shadow-[0_8px_20px_-4px_rgba(16,185,129,0.4)] translate-x-1' 
          : 'text-slate-400 hover:text-white hover:bg-white/5'
      }`}
    >
      <Icon className={`w-5 h-5 transition-transform duration-300 ${view === id ? 'scale-110' : 'group-hover:scale-110'}`} />
      <span className="text-xs font-black uppercase tracking-[0.2em]">{label}</span>
      {status === AppStatus.SUCCESS && (id === 'RESULTS' || id === 'ANALYSIS') && view === 'SEARCH' && (
        <span className="absolute -right-1 -top-1 w-3 h-3 bg-blue-500 rounded-full animate-ping border-2 border-black"></span>
      )}
      {view === id && (
        <span className="absolute right-4 w-1.5 h-1.5 bg-white rounded-full"></span>
      )}
    </button>
  );

  return (
    <div className="flex h-screen bg-[#010204] overflow-hidden relative selection:bg-emerald-500/30">
      <div className="absolute inset-0 bg-grid opacity-40 pointer-events-none"></div>
      <div className="absolute inset-0 bg-mesh pointer-events-none"></div>
      <div className="scanline pointer-events-none"></div>

      <aside className="w-80 glass border-r border-white/5 flex flex-col z-20 m-4 rounded-3xl shadow-2xl overflow-hidden">
        <div className="p-10 border-b border-white/5">
          <div className="flex items-center gap-5">
            <LabLogo />
            <div>
              <h1 className="text-white font-black text-2xl tracking-tighter uppercase leading-none">Pyner</h1>
              <span className="text-emerald-500 font-bold text-[10px] uppercase tracking-[0.3em] block mt-1">SRA Node UI</span>
            </div>
          </div>
        </div>

        <nav className="flex-1 p-6 space-y-4 overflow-y-auto">
          <div className="px-4 py-2">
            <h2 className="text-[10px] text-slate-500 font-black uppercase tracking-[0.4em] mb-4">Laboratory</h2>
            <div className="space-y-2">
              <NavItem id="SEARCH" icon={Search} label="Search & Mine" />
              <NavItem id="RESULTS" icon={Database} label="Data Grid" />
              <NavItem id="ANALYSIS" icon={BarChart3} label="Analytics" />
            </div>
          </div>
        </nav>

        <div className="p-8 bg-black/30">
           <a href="#" className="flex items-center gap-3 text-[10px] font-black text-slate-500 hover:text-emerald-400 transition-colors uppercase tracking-[0.2em]">
             <Github className="w-5 h-5" />
             <span>Build Stable 8.0.1</span>
           </a>
        </div>
      </aside>

      <main className="flex-1 flex flex-col overflow-hidden relative z-10 m-4 ml-0">
        <header className="h-20 glass rounded-t-3xl border-b-0 px-10 flex items-center justify-between shadow-lg">
          <div className="flex items-center gap-8">
             <div className="flex items-center gap-3 px-4 py-2 bg-emerald-500/5 border border-emerald-500/20 rounded-2xl text-[10px] text-emerald-400 font-black uppercase tracking-widest">
               <Cpu className="w-4 h-4" />
               <span>Interface Online</span>
             </div>
          </div>
          <div className="flex items-center gap-4">
             {status === AppStatus.SUCCESS && view === 'SEARCH' && (
               <div className="flex items-center gap-3 animate-in slide-in-from-right-10 duration-500">
                  <span className="text-[10px] font-black text-blue-400 uppercase tracking-widest bg-blue-500/10 px-3 py-1 rounded-lg border border-blue-500/20 flex items-center gap-2">
                    <AlertCircle className="w-3.5 h-3.5" /> 
                    Resultados listos en Analytics
                  </span>
               </div>
             )}
             <div className="text-right">
                <div className="flex items-center gap-3 justify-end">
                   <span className="text-[10px] font-black uppercase tracking-widest text-slate-500">Node Status:</span>
                   <span className={`w-2.5 h-2.5 rounded-full ${status === AppStatus.SUCCESS ? 'bg-emerald-500 shadow-[0_0_10px_#10b981]' : status === AppStatus.GENERATING ? 'bg-amber-500 animate-pulse' : 'bg-slate-700'}`}></span>
                </div>
             </div>
          </div>
        </header>

        <div className="flex-1 overflow-y-auto p-12 bg-black/40 rounded-b-3xl border border-t-0 border-white/5 glass shadow-inner relative z-10">
          {view === 'SEARCH' && (
            <div className="max-w-6xl mx-auto h-full flex flex-col">
              <div className="text-center mb-10">
                 <div className="inline-flex items-center gap-2 px-5 py-2 mb-6 rounded-full border border-purple-500/20 bg-purple-500/5 text-purple-400 text-[10px] font-black uppercase tracking-[0.5em]">
                    <Box className="w-3.5 h-3.5" /> Main Console
                 </div>
                 <h2 className="text-7xl font-black text-white mb-4 tracking-tighter leading-[0.8] drop-shadow-2xl">
                   Biological<br/>
                   <span className="text-transparent bg-clip-text bg-gradient-to-r from-emerald-400 to-emerald-600 font-black italic">Data Intelligence</span>
                 </h2>
              </div>
              
              <div className="flex-1 glass rounded-[3rem] overflow-hidden shadow-2xl border-white/10 ring-1 ring-white/5 min-h-[550px]">
                <div className="grid grid-cols-1 lg:grid-cols-12 divide-x divide-white/5 h-full">
                  <div className="lg:col-span-6 h-full">
                    <InputSection 
                      inputs={inputs} 
                      setInputs={setInputs} 
                      onGenerate={handleGenerate} 
                      isGenerating={status === AppStatus.GENERATING}
                    />
                  </div>
                  <div className="lg:col-span-6 h-full">
                    <OutputSection query={data} status={status} error={error} />
                  </div>
                </div>
              </div>
            </div>
          )}

          <div className="max-w-7xl mx-auto">
            {view === 'RESULTS' && <ResultsTable data={data?.results || []} status={status} />}
            {view === 'ANALYSIS' && <AnalysisDashboard data={data?.results || []} status={status} />}
            {view === 'DOCS' && <Documentation />}
          </div>
        </div>
      </main>
    </div>
  );
};

export default App;
