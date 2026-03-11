
import React, { useState } from 'react';
import { SRAResultEntry, AppStatus } from '../types';
// Fixed misplaced import by moving Database to the top-level imports
import { Download, Filter, Search, ChevronRight, Database } from 'lucide-react';

export const ResultsTable: React.FC<{ data: SRAResultEntry[], status: AppStatus }> = ({ data, status }) => {
  const [filter, setFilter] = useState('');

  const filteredData = data.filter(item => 
    item.title.toLowerCase().includes(filter.toLowerCase()) || 
    item.organism.toLowerCase().includes(filter.toLowerCase())
  );

  const handleExport = () => {
    const csv = [
      ['BioProject', 'Title', 'Organism', 'Tissue', 'Conditions', 'Strategy'],
      ...filteredData.map(r => [r.bioproject, r.title, r.organism, r.tissue, r.conditions, r.strategy])
    ].map(e => e.join(",")).join("\n");
    
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `sra_export_${new Date().getTime()}.csv`;
    a.click();
  };

  if (status === AppStatus.IDLE) {
    return (
      <div className="h-full flex flex-col items-center justify-center text-center p-12 bg-white/5 rounded-2xl border border-white/10 border-dashed">
        <Database className="w-12 h-12 text-neutral-700 mb-4" />
        <h3 className="text-xl font-bold text-white mb-2">No Data Mined</h3>
        <p className="text-neutral-500 text-sm max-w-sm">Please return to the Search tab and execute a mining sequence to populate this repository.</p>
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-2xl font-black text-white tracking-tighter">Mined Repository</h2>
          <p className="text-xs text-neutral-500 uppercase tracking-widest font-bold">Synchronized SRA Metadata</p>
        </div>
        <div className="flex gap-3">
          <div className="relative">
            <Search className="w-4 h-4 absolute left-3 top-1/2 -translate-y-1/2 text-neutral-500" />
            <input 
              type="text" 
              placeholder="Quick Filter..." 
              value={filter}
              onChange={(e) => setFilter(e.target.value)}
              className="bg-[#0a0f18] border border-white/10 rounded-lg pl-10 pr-4 py-2 text-xs focus:border-emerald-500 outline-none w-64 text-neutral-200"
            />
          </div>
          <button 
            onClick={handleExport}
            className="flex items-center gap-2 px-4 py-2 bg-emerald-600 hover:bg-emerald-500 text-white rounded-lg text-xs font-bold transition-all shadow-lg shadow-emerald-500/20"
          >
            <Download className="w-4 h-4" /> Export CSV
          </button>
        </div>
      </div>

      <div className="bg-[#0a0f18] border border-white/5 rounded-xl overflow-hidden">
        <table className="w-full text-left text-xs">
          <thead>
            <tr className="bg-white/5 text-neutral-400 border-b border-white/5 uppercase tracking-widest text-[10px] font-bold">
              <th className="px-6 py-4">BioProject ID</th>
              <th className="px-6 py-4">Title / Context</th>
              <th className="px-6 py-4">Organism</th>
              <th className="px-6 py-4">Tissue</th>
              <th className="px-6 py-4 text-center">T-Series</th>
              <th className="px-6 py-4 text-right">Strategy</th>
            </tr>
          </thead>
          <tbody className="divide-y divide-white/5">
            {filteredData.map((row, i) => (
              <tr key={i} className="hover:bg-white/5 transition-colors group">
                <td className="px-6 py-4 font-mono text-emerald-500 font-bold">{row.bioproject}</td>
                <td className="px-6 py-4 max-w-xs truncate text-neutral-300 font-medium">{row.title}</td>
                <td className="px-6 py-4 text-neutral-400 italic">{row.organism}</td>
                <td className="px-6 py-4">
                  <span className="px-2 py-1 bg-white/5 rounded border border-white/10 text-[9px] uppercase font-bold text-neutral-500">
                    {row.tissue}
                  </span>
                </td>
                <td className="px-6 py-4 text-center">
                  {row.isTimeSeries ? (
                    <span className="text-blue-500 font-bold">TRUE</span>
                  ) : (
                    <span className="text-neutral-700">FALSE</span>
                  )}
                </td>
                <td className="px-6 py-4 text-right">
                   <div className="flex items-center justify-end gap-2 text-neutral-500 group-hover:text-emerald-500 transition-colors">
                      {row.strategy}
                      <ChevronRight className="w-3 h-3" />
                   </div>
                </td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
    </div>
  );
};
