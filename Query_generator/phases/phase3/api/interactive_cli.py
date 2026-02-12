"""
Interactive CLI for Query Generator
====================================
User interaction flow: generate → validate → correct → learn
"""

import sys
import logging
from pathlib import Path
from typing import Dict, Optional

from phase3.api.query_generator import QueryGeneratorService

logger = logging.getLogger(__name__)


class InteractiveQueryGenerator:
    """Generador de queries interactivo con aprendizaje"""
    
    def __init__(self, kb_path: Path, ollama_client=None, query_cache_path: Optional[Path] = None, cache_dir: Optional[Path] = None):
        self.cache_dir = cache_dir or Path(__file__).parent.parent / "support_dictionary"
        self.service = QueryGeneratorService(kb_path, ollama_client, query_cache_path, cache_dir=self.cache_dir)
    
    def run_interactive(self, user_input: str) -> bool:
        """
        Run query generation in interactive mode
        
        Flow:
        1. Generate query
        2. Display results and return query (no satisfaction prompt)
        
        Returns:
            True if user is satisfied with final query
        """
        
        print("\n" + "="*80)
        print("🔍 PYNER QUERY GENERATOR")
        print("="*80)
        print("I turn natural language into NCBI SRA boolean queries.")
        print("Try adding organism + strategy + tissue/condition for best results.")
        print("Suggestions:")
        print("  - Arabidopsis thaliana drought stress RNA-Seq")
        print("  - Mus musculus liver transcriptome")
        print("  - rice root heat stress")
        print("="*80)
        
        # Paso 1: Generar query inicial
        result = self.service.generate_query(user_input, use_llm=True)
        
        # Show if clarification is needed
        if result.get('clarification_needed'):
            print(f"\n⚠️ Need more details:")
            print(f"   {result.get('clarification_message', '')}")
            return False
        
        # Mostrar extracción
        self._print_extraction(result)

        # Broad requests: "quiero saber todo sobre ..."
        if self._is_broad_request(user_input) and result['extracted'].get('organism_variants'):
            return self._handle_broad_request(result)

        print(f"\n✅ NCBI Query:")
        print(f"   {result['ncbi_query']}")
        print("\n" + "="*80)
        print("📋 Copy the query above and paste it into NCBI SRA search")
        print("="*80 + "\n")
        return True
    
    def _print_extraction(self, result: Dict):
        """Display extracted terms"""
        print(f"\n📝 Input: {result['user_input']}")
        print(f"\n🔎 Extracted Terms:")
        print(f"   Organism:    {result['extracted']['organism'] or '(none)'}")
        print(f"   Strategies:  {', '.join(result['extracted']['strategies']) or '(none)'}")
        print(f"   Tissues:     {', '.join(result['extracted']['tissues']) or '(none)'}")
        print(f"   Conditions:  {', '.join(result['extracted']['conditions']) or '(none)'}")
        print(f"   Keywords:    {', '.join(result['extracted']['free_terms']) or '(none)'}")

        syn = result.get('synonyms', {})
        org_syn = ', '.join(syn.get('organism', []) or []) or '(none)'
        strat_syn = ', '.join(syn.get('strategies', []) or []) or '(none)'
        tissue_syn = ', '.join(syn.get('tissues', []) or []) or '(none)'
        cond_syn = ', '.join(syn.get('conditions', []) or []) or '(none)'
        print(f"\n🔗 Synonyms:")
        print(f"   Organism:    {org_syn}")
        print(f"   Strategies:  {strat_syn}")
        print(f"   Tissues:     {tissue_syn}")
        print(f"   Conditions:  {cond_syn}")
        
        print(f"\n📊 Generated Query:")
        print(f"   {result['ncbi_query']}")
    
    def _ask_satisfaction(self) -> bool:
        """Ask if user is satisfied with the query"""
        while True:
            response = input("\n✓ Are you satisfied with this query? (y/n): ").strip().lower()
            if response in ['s', 'si', 'yes', 'y']:
                return True
            elif response in ['n', 'no']:
                return False
            else:
                print("Please answer 'y' or 'n'")

    def _is_broad_request(self, user_input: str) -> bool:
        text = user_input.lower()
        return any(
            phrase in text
            for phrase in [
                "todo sobre", "todo acerca", "todo de", "quiero saber todo",
                "all about", "everything about"
            ]
        )

    def _handle_broad_request(self, result: Dict) -> bool:
        org_variants = result['extracted'].get('organism_variants', [])
        organism = org_variants[0] if org_variants else result['extracted'].get('organism')
        if not organism:
            return False

        print("\nI can search everything for this organism or narrow by strategy.")
        print("Do you want EVERYTHING? (y/n)")
        response = input("> ").strip().lower()

        if response in ['y', 'yes', 's', 'si']:
            query = self.service.query_builder.build_query(organisms=org_variants)
            print(f"\n✅ NCBI Query (all data):")
            print(f"   {query}")
            print("\nTip: add tissue/condition for a focused search.")
            return True

        # Show strategy options
        strategies = sorted(self.service.validator.strategies)
        print("\nAvailable strategies:")
        for i, strat in enumerate(strategies, 1):
            print(f"  {i:2}. {strat}")

        choice = input("\nPick a strategy by number (or press Enter to skip): ").strip()
        if not choice:
            print("\nTip: add tissue/condition for a focused search.")
            return False

        selected = None
        if choice.isdigit():
            idx = int(choice)
            if 1 <= idx <= len(strategies):
                selected = strategies[idx - 1]

        if not selected:
            print("Invalid selection.")
            return False

        query = self.service.query_builder.build_query(
            organisms=org_variants,
            strategies=[selected]
        )
        print(f"\n✅ NCBI Query:")
        print(f"   {query}")
        print("\nTip: add tissue/condition for a focused search.")
        return True
    
    def show_statistics(self):
        """Display technical vocabulary statistics"""
        vocab = self.service.validator.vocab or {}
        tech_path = self.cache_dir / "technical_vocabulary.json"
        plant_tissues = len(vocab.get("plant_tissues", []))
        animal_tissues = len(vocab.get("animal_tissues", []))
        generic_tissues = len(vocab.get("generic_tissues", []))
        process_keywords = len(vocab.get("process_keywords", []))
        organism_aliases = len(vocab.get("organism_aliases", {}))
        strategy_keywords = len(vocab.get("strategy_keywords", {}))
        print("\n" + "="*80)
        print("📊 TECHNICAL VOCABULARY STATISTICS")
        print("="*80)
        print(f"  Plant tissues:       {plant_tissues} entries")
        print(f"  Animal tissues:      {animal_tissues} entries")
        print(f"  Generic tissues:     {generic_tissues} entries")
        print(f"  Process keywords:    {process_keywords} entries")
        print(f"  Organism aliases:    {organism_aliases} entries")
        print(f"  Strategy keywords:   {strategy_keywords} entries")
        print(f"\n  📁 Technical vocab file: {tech_path}")
        print("")
