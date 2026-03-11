
import React, { useRef, useEffect } from 'react';
import { SRAInputs } from '../types';
import { Search, Sparkles, Terminal, Command, Beaker } from 'lucide-react';

interface Props {
  inputs: SRAInputs;
  setInputs: (val: SRAInputs) => void;
  onGenerate: () => void;
  isGenerating: boolean;
}

export const InputSection: React.FC<Props> = ({ inputs, setInputs, onGenerate, isGenerating }) => {
  const textareaRef = useRef<HTMLTextAreaElement>(null);

  useEffect(() => {
    if (textareaRef.current) {
      textareaRef.current.focus();
    }
  }, []);

  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter' && e.ctrlKey) {
      onGenerate();
    }
  };

  const loadExample = () => {
    setInputs({ 
      naturalQuery: "Toda información de Arabidopsis thaliana en RNA-Seq, sin importar el tejido, en condiciones de nitrógeno o sequía" 
    });
    if (textareaRef.current) {
      textareaRef.current.focus();
    }
  };

  return (
    <div className="p-8 h-full flex flex-col bg-black/40 relative z-30">
      <div className="flex items-center justify-between mb-6">
        <div className="flex items-center gap-3">
          <div className="w-1 h-8 bg-emerald-500 rounded-full"></div>
          <div>
            <h2 className="text-sm font-black text-white uppercase tracking-[0.3em]">Query Interface</h2>
            <p className="text-[10px] text-slate-500 font-bold uppercase tracking-widest">Natural Language Processing Mode</p>
          </div>
        </div>
        <button 
          onClick={loadExample}
          className="px-4 py-2 bg-emerald-500/10 border border-emerald-500/20 hover:bg-emerald-500/20 text-emerald-400 rounded-xl flex items-center gap-2 text-[10px] font-black uppercase tracking-widest transition-all group active:scale-95"
        >
          <Beaker className="w-3.5 h-3.5 group-hover:rotate-12 transition-transform" />
          Cargar Ejemplo
        </button>
      </div>

      <div className="flex-1 flex flex-col gap-4 min-h-0">
        <div className="relative flex-1 group">
          <textarea
            ref={textareaRef}
            value={inputs.naturalQuery}
            onChange={(e) => setInputs({ naturalQuery: e.target.value })}
            onKeyDown={handleKeyDown}
            placeholder="Escriba aquí su consulta biológica (Ej: 'Arabidopsis thaliana RNA-Seq en sequía')..."
            className="w-full h-full bg-[#030508] border-2 border-white/5 text-emerald-400 text-lg rounded-[2.5rem] p-10 focus:border-emerald-500/40 focus:ring-8 focus:ring-emerald-500/5 focus:outline-none transition-all duration-500 placeholder-slate-800 font-medium resize-none leading-relaxed shadow-2xl relative z-40 selection:bg-emerald-500/30"
          />
          
          <div className="absolute bottom-8 left-10 flex items-center gap-6 pointer-events-none opacity-40">
             <div className="flex items-center gap-2">
                <Command className="w-3.5 h-3.5" />
                <span className="text-[10px] font-black uppercase tracking-widest">Ctrl + Enter para Ejecutar</span>
             </div>
          </div>

          <div className="absolute bottom-8 right-10 pointer-events-none">
            <Sparkles className="w-5 h-5 text-emerald-900 group-focus-within:text-emerald-500 transition-colors duration-700" />
          </div>
        </div>

        <div className="bg-[#05070a] border border-white/5 p-6 rounded-3xl flex items-center gap-5">
           <div className="p-3 bg-emerald-500/10 rounded-2xl">
              <Search className="w-5 h-5 text-emerald-500" />
           </div>
           <div>
              <h4 className="text-[11px] font-black text-slate-300 uppercase tracking-widest mb-1">Syntactic Analysis</h4>
              <p className="text-[11px] text-slate-500 font-medium italic">
                "Pyner identificará automáticamente taxones, tejidos y condiciones experimentales."
              </p>
           </div>
        </div>
      </div>

      <div className="mt-8">
        <button
          onClick={onGenerate}
          disabled={isGenerating || inputs.naturalQuery.length < 5}
          className={`w-full py-6 rounded-[2rem] text-xs font-black uppercase tracking-[0.5em] transition-all duration-500 flex items-center justify-center gap-4 border ${
            isGenerating || inputs.naturalQuery.length < 5
              ? 'bg-transparent border-white/5 text-slate-700 cursor-not-allowed'
              : 'bg-emerald-600 border-emerald-400/20 hover:bg-emerald-500 text-white shadow-[0_20px_40px_-12px_rgba(16,185,129,0.5)] active:scale-[0.98]'
          }`}
        >
          {isGenerating ? (
            <div className="flex items-center gap-3">
              <div className="w-4 h-4 border-2 border-white/20 border-t-white rounded-full animate-spin"></div>
              <span>Procesando Petición</span>
            </div>
          ) : (
            <>
              <Search className="w-4 h-4" />
              Analizar y Ejecutar Consulta
            </>
          )}
        </button>
      </div>
    </div>
  );
};
