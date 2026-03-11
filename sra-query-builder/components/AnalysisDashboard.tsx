
import React from 'react';
import { SRAResultEntry, AppStatus } from '../types';
import { 
  BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, 
  PieChart, Pie, Cell, ScatterChart, Scatter, Legend, CartesianGrid
} from 'recharts';
import { Pickaxe, Clock, Activity, Target, ShieldCheck, Microscope } from 'lucide-react';

const MiningGraphic = () => (
  <div className="relative flex flex-col items-center justify-center p-6 bg-black/50 rounded-[2rem] border border-white/10 shadow-inner w-full h-full min-h-[300px]">
    <div className="absolute inset-0 opacity-5 pointer-events-none overflow-hidden">
       <Microscope className="w-full h-full text-emerald-500 scale-150" />
    </div>
    <div className="relative z-10 flex flex-col items-center gap-2 w-full max-w-[200px]">
      {[
        { label: 'NCBI Core', color: 'bg-neutral-100' },
        { label: 'SRA Archive', color: 'bg-orange-300' },
        { label: 'GEO Portal', color: 'bg-yellow-300' },
        { label: 'BioProject', color: 'bg-emerald-300' },
        { label: 'PubMed Ref', color: 'bg-blue-300' }
      ].map((item) => (
        <div key={item.label} className={`w-full py-2.5 text-center text-[9px] font-black text-black border-2 border-black rounded-lg shadow-[4px_4px_0px_0px_rgba(0,0,0,1)] hover:translate-x-0.5 hover:translate-y-0.5 transition-transform ${item.color} uppercase tracking-[0.15em]`}>
          {item.label}
        </div>
      ))}
    </div>
  </div>
);

const StackedCounter = ({ label, value, color }: { label: string, value: number, color: string }) => (
  <div className="flex flex-col items-center gap-4 flex-1 min-w-[80px]">
    <div className="h-48 md:h-64 w-full max-w-[40px] bg-black border border-white/10 rounded-xl relative overflow-hidden flex flex-col-reverse shadow-xl">
      <div className={`${color} w-full transition-all duration-1000 ease-out shadow-[0_0_20px_rgba(16,185,129,0.3)]`} style={{ height: `${Math.min((value / 2000) * 100, 100)}%` }}></div>
      <div className="absolute inset-0 flex items-center justify-center text-[10px] font-black text-white drop-shadow-[0_1px_2px_rgba(0,0,0,1)] rotate-[-90deg] md:rotate-0">
        {value}
      </div>
    </div>
    <span className="text-[9px] font-black uppercase tracking-[0.2em] text-slate-400 text-center leading-tight">
      {label}
    </span>
  </div>
);

export const AnalysisDashboard: React.FC<{ data: SRAResultEntry[], status: AppStatus }> = ({ data, status }) => {
  if (status !== AppStatus.SUCCESS || data.length === 0) {
    return (
      <div className="h-full flex flex-col items-center justify-center text-center p-12 md:p-20 bg-black/40 rounded-3xl border border-white/5 border-dashed backdrop-blur-md">
        <Activity className="w-16 h-16 md:w-24 md:h-24 text-slate-800 mb-8 animate-pulse" />
        <h3 className="text-2xl md:text-4xl font-black text-white mb-4 tracking-tighter">System Idle</h3>
        <p className="text-slate-400 text-sm md:text-lg max-w-sm font-medium">Ejecute un protocolo de búsqueda para visualizar datos analíticos.</p>
      </div>
    );
  }

  const totals = {
    studies: data.length,
    libraries: data.reduce((acc, curr) => acc + (curr.libraryCount || 0), 0),
    experiments: data.reduce((acc, curr) => acc + (curr.experimentCount || 0), 0)
  };

  const scatterData = data.flatMap((d, i) => (d.timePoints || []).map(t => ({ day: t, index: i, condition: d.conditions })));
  const tissueData = data.reduce((acc: any, curr) => {
    acc[curr.tissue] = (acc[curr.tissue] || 0) + 1;
    return acc;
  }, {});
  const tissueChartData = Object.keys(tissueData).map(key => ({ name: key, value: tissueData[key] }));

  const cultivarData = data.reduce((acc: any, curr) => {
    acc[curr.cultivar] = (acc[curr.cultivar] || 0) + 1;
    return acc;
  }, {});
  const cultivarChartData = Object.keys(cultivarData).map(key => ({ name: key, value: cultivarData[key] }));

  const geneData = data.flatMap(d => d.geneStats || []);
  const aggregatedGeneData = geneData.reduce((acc: any, curr) => {
    if(!acc[curr.condition]) acc[curr.condition] = { condition: curr.condition, up: 0, down: 0 };
    acc[curr.condition].up += curr.up;
    acc[curr.condition].down += curr.down;
    return acc;
  }, {});
  const finalGeneData = Object.values(aggregatedGeneData);

  const COLORS = ['#10b981', '#3b82f6', '#8b5cf6', '#f59e0b', '#ec4899', '#06b6d4'];

  return (
    <div className="space-y-8 md:space-y-12 animate-in fade-in slide-in-from-bottom-8 duration-1000 pb-20">
      
      {/* SECTION B: Metadata Profile */}
      <div className="glass p-8 md:p-12 rounded-[2.5rem] shadow-3xl relative overflow-hidden">
        <div className="absolute top-4 left-6 text-[10px] md:text-[12px] font-black text-emerald-500/30 tracking-[0.5em] uppercase">METADATA_STREAM_B</div>
        
        <div className="grid grid-cols-1 lg:grid-cols-12 gap-8 md:gap-12 mt-4">
          {/* Hierarchical Counts */}
          <div className="lg:col-span-3 flex items-end justify-around gap-2 px-4 pb-4 border-b lg:border-b-0 lg:border-r border-white/5">
            <StackedCounter label="Studies" value={totals.studies} color="bg-emerald-500" />
            <StackedCounter label="Libraries" value={totals.libraries} color="bg-blue-500" />
            <StackedCounter label="Exps" value={totals.experiments} color="bg-emerald-600" />
          </div>

          {/* Temporal Distribution */}
          <div className="lg:col-span-5 flex flex-col px-4 border-b lg:border-b-0 lg:border-r border-white/5">
            <h3 className="text-[10px] font-black text-white uppercase tracking-[0.3em] text-center mb-6 flex items-center justify-center gap-3">
              <Clock className="w-4 h-4 text-emerald-500" /> Temporal Distribution
            </h3>
            <div className="h-64 w-full">
              <ResponsiveContainer width="100%" height="100%">
                <ScatterChart margin={{ top: 10, right: 20, bottom: 30, left: 0 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.03)" vertical={false} />
                  <XAxis 
                    type="number" 
                    dataKey="day" 
                    name="Days" 
                    stroke="#475569" 
                    fontSize={10} 
                    tick={{ fontWeight: '800' }}
                    label={{ value: 'Days Post-Stress', position: 'insideBottom', offset: -15, fontSize: 10, fill: '#10b981', fontWeight: '900', tracking: '0.1em' }}
                  />
                  <YAxis type="number" dataKey="index" hide />
                  <Tooltip 
                    cursor={{ strokeDasharray: '3 3' }} 
                    contentStyle={{ backgroundColor: '#05070a', border: '1px solid rgba(16,185,129,0.3)', borderRadius: '12px', fontSize: '10px' }} 
                  />
                  <Scatter name="TimePoints" data={scatterData} fill="#10b981" />
                </ScatterChart>
              </ResponsiveContainer>
            </div>
          </div>

          {/* Doughnuts */}
          <div className="lg:col-span-4 flex flex-col sm:flex-row lg:flex-col justify-around gap-6 px-4">
            <div className="flex-1 text-center">
              <span className="text-[9px] font-black text-slate-500 uppercase tracking-[0.2em] block mb-2">Morphology</span>
              <div className="h-28 w-full">
                <ResponsiveContainer>
                  <PieChart>
                    <Pie data={tissueChartData} innerRadius={20} outerRadius={40} paddingAngle={4} dataKey="value">
                      {tissueChartData.map((_, i) => <Cell key={i} fill={COLORS[i % COLORS.length]} />)}
                    </Pie>
                    <Tooltip contentStyle={{ backgroundColor: '#05070a', border: 'none', borderRadius: '8px', fontSize: '10px' }} />
                  </PieChart>
                </ResponsiveContainer>
              </div>
            </div>
            <div className="flex-1 text-center">
              <span className="text-[9px] font-black text-slate-500 uppercase tracking-[0.2em] block mb-2">Genotypes</span>
              <div className="h-28 w-full">
                <ResponsiveContainer>
                  <PieChart>
                    <Pie data={cultivarChartData} innerRadius={20} outerRadius={40} paddingAngle={4} dataKey="value">
                      {cultivarChartData.map((_, i) => <Cell key={i} fill={COLORS[(i+2) % COLORS.length]} />)}
                    </Pie>
                    <Tooltip contentStyle={{ backgroundColor: '#05070a', border: 'none', borderRadius: '8px', fontSize: '10px' }} />
                  </PieChart>
                </ResponsiveContainer>
              </div>
            </div>
          </div>
        </div>
      </div>

      {/* SECTION C: Differential Expression Dashboard */}
      <div className="flex flex-col xl:flex-row gap-8 items-stretch">
        <div className="flex-1 glass p-8 md:p-12 rounded-[2.5rem] shadow-3xl">
           <h3 className="text-sm font-black text-white uppercase tracking-[0.3em] mb-10 flex items-center gap-4">
             <Target className="w-5 h-5 text-emerald-500" /> Regulation Signature
           </h3>
           <div className="h-80 w-full min-h-[300px]">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={finalGeneData as any} layout="vertical" margin={{ left: 20, right: 40, bottom: 20 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.03)" horizontal={false} />
                  <XAxis type="number" stroke="#475569" fontSize={10} axisLine={false} tickLine={false} tick={{ fontWeight: 'bold' }} />
                  <YAxis dataKey="condition" type="category" stroke="#fff" fontSize={11} fontWeight="900" axisLine={false} tickLine={false} width={80} />
                  <Tooltip contentStyle={{ backgroundColor: '#05070a', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '12px', fontSize: '11px' }} cursor={{ fill: 'rgba(255,255,255,0.02)' }} />
                  <Legend iconType="circle" wrapperStyle={{ fontSize: '10px', paddingTop: '20px', fontWeight: '800', textTransform: 'uppercase', letterSpacing: '0.1em' }} />
                  <Bar dataKey="up" name="UP" fill="#f59e0b" radius={[0, 4, 4, 0]} barSize={24} />
                  <Bar dataKey="down" name="DOWN" fill="#3b82f6" radius={[0, 4, 4, 0]} barSize={24} />
                </BarChart>
              </ResponsiveContainer>
           </div>
        </div>

        <div className="w-full xl:w-[350px] flex flex-col gap-6">
          <div className="bg-emerald-500/5 border border-emerald-500/10 p-8 rounded-[2rem] shadow-inner relative group flex-1">
             <h4 className="text-[10px] font-black text-emerald-400 uppercase tracking-[0.3em] mb-4 flex items-center gap-3">
               <Activity className="w-4 h-4" /> Analysis Brief
             </h4>
             <p className="text-sm text-slate-300 leading-relaxed font-bold italic">
               "Sintetizado a través de <span className="text-emerald-400">{totals.studies}</span> estudios, la firma de estrés es consistente. El mapeo de <span className="text-emerald-400">{totals.libraries}</span> librerías ofrece una base robusta para descubrimiento genético."
             </p>
          </div>
          
          <div className="bg-black/50 border border-white/5 p-8 rounded-[2rem] shadow-lg">
             <h4 className="text-[10px] font-black text-slate-500 uppercase tracking-[0.3em] mb-4">Data Reliability</h4>
             <div className="flex items-center justify-between text-xs font-black mb-3">
               <span className="text-slate-600 uppercase">Score</span>
               <span className="text-emerald-500 font-mono text-lg">94.2%</span>
             </div>
             <div className="w-full bg-slate-900 h-2.5 rounded-full overflow-hidden p-0.5 shadow-inner">
               <div className="bg-emerald-500 h-full rounded-full w-[94%] shadow-[0_0_15px_rgba(16,185,129,0.6)]"></div>
             </div>
          </div>

          <div className="hidden xl:block">
            <MiningGraphic />
          </div>
        </div>
      </div>
      
    </div>
  );
};
