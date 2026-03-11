
import React from 'react';
import { ExternalLink, Book, Map, Users, ChevronRight } from 'lucide-react';

const DocCard = ({ icon: Icon, title, desc, link }: any) => (
  <a 
    href={link} 
    className="block p-6 bg-[#0a0f18] border border-white/5 rounded-xl hover:border-emerald-500/30 transition-all group"
  >
    <div className="flex items-center gap-4 mb-3">
      <div className="w-10 h-10 bg-white/5 rounded-lg flex items-center justify-center text-neutral-400 group-hover:text-emerald-500 transition-colors">
        <Icon className="w-5 h-5" />
      </div>
      <h3 className="text-white font-bold group-hover:translate-x-1 transition-transform">{title}</h3>
    </div>
    <p className="text-xs text-neutral-500 leading-relaxed mb-4">{desc}</p>
    <div className="flex items-center gap-2 text-[10px] font-bold uppercase tracking-widest text-emerald-500 opacity-0 group-hover:opacity-100 transition-opacity">
      View Reference <ExternalLink className="w-3 h-3" />
    </div>
  </a>
);

export const Documentation: React.FC = () => {
  return (
    <div className="max-w-4xl mx-auto space-y-12 animate-in slide-in-from-bottom duration-500">
      <div className="text-center">
        <h2 className="text-3xl font-black text-white mb-2 tracking-tighter">Resource Center</h2>
        <p className="text-neutral-500 text-sm">Official documentation and developer guidelines for the Pyner Ecosystem.</p>
      </div>

      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        <DocCard 
          icon={Book} 
          title="README.md" 
          desc="Quickstart guide, installation instructions, and basic usage of the Pyner CLI and API."
          link="#"
        />
        <DocCard 
          icon={Map} 
          title="ROADMAP.md" 
          desc="The technical trajectory of the MAInr agent, including planned GPU cluster optimizations."
          link="#"
        />
        <DocCard 
          icon={Users} 
          title="CONTRIBUTING.md" 
          desc="Join the community. Guidelines for code standards, PRs, and biological metadata mapping."
          link="#"
        />
        <DocCard 
          icon={ExternalLink} 
          title="NCBI E-Utilities API" 
          desc="Deep dive into the Entrez Programming Utilities used for underlying data retrieval."
          link="https://www.ncbi.nlm.nih.gov/books/NBK25501/"
        />
      </div>

      <div className="bg-[#0a0f18] border border-white/5 p-8 rounded-2xl prose prose-invert max-w-none text-neutral-400 text-xs">
         <h4 className="text-white font-bold mb-4 uppercase tracking-widest">Search Philosophy</h4>
         <p className="mb-4">
           The Pyner engine utilizes a two-stage boolean expansion. First, natural language biological concepts (e.g., "drought") are cross-referenced with taxonomic 
           dictionaries to include synonyms like "water stress" or "osmotic pressure."
         </p>
         <ul className="space-y-2 list-none p-0">
           <li className="flex items-center gap-2"><ChevronRight className="w-3 h-3 text-emerald-500" /> [Organism] tags are strictly prioritized for search performance.</li>
           <li className="flex items-center gap-2"><ChevronRight className="w-3 h-3 text-emerald-500" /> [All Fields] is used for broad tissue and condition discovery.</li>
           <li className="flex items-center gap-2"><ChevronRight className="w-3 h-3 text-emerald-500" /> Boolean logic avoids NCBI-specific tags not valid in the SRA domain.</li>
         </ul>
      </div>
    </div>
  );
};
