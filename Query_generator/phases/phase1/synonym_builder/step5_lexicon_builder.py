import os
import re
import json
import time
from collections import Counter
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
INPUT_FILE = BASE_DIR / "output" / "step2_mesh_consolidated_parallel_1000s.json"
OUTPUT_FILE = BASE_DIR / "output" / "experimental_conditions_lexicon.json"

# Núcleos de condición experimental que por sí solos indican una condición
CONDITION_CORES = {
    "stress", "stresses", "exposure", "treatment", "deficiency", 
    "deprivation", "starvation", "toxicity", "challenge", "infection", 
    "injury", "hypoxia", "anoxia", "irradiation", "salinity",
    "limitation", "depletion", "restriction", "damage", "perturbation", "inflammation",
    "drought", "heat", "cold", "salt", "freezing", "chilling"
}

# Modificadores que suelen acompañar a un núcleo (o indicar condición en combinación)
CONDITION_MODIFIERS = {
    "heat", "cold", "freezing", "chilling", "temperature", "temperatures",
    "thermal", "radiation", "uv", "ultraviolet", "drug", "drugs", "chemical",
    "water", "osmotic", "salt", "drought", "oxidative", "nutrient", 
    "metal", "light", "dark", "shade",
    "oxygen", "glucose", "iron", "phosphate", "mechanical", "pressure", "shear",
    "immune", "viral", "bacterial", "toxin", "compound"
}

# Palabras compatibles adicionales permitidas junto a un modificador para validar
COMPATIBLE_WORDS = CONDITION_CORES.union({
    "induced", "mediated", "response", "treated", "tolerance", 
    "resistant", "resistance", "damaged", "inflammatory"
})


# Simple regex for tokenization (words only, length >= 2)
TOKEN_RE = re.compile(r'\b[a-z]{2,}\b')

def get_ngrams(text, n_min=2, n_max=4):
    """Extract 2-4 word n-grams from lowercased text."""
    # Split text into sentences roughly, then tokenize to prevent n-grams across sentence boundaries
    sentences = re.split(r'[.!?]', text.lower())
    ngrams = []
    
    for sentence in sentences:
        tokens = TOKEN_RE.findall(sentence)
        if len(tokens) < n_min:
            continue
            
        for n in range(n_min, n_max + 1):
            for i in range(len(tokens) - n + 1):
                gram_tokens = set(tokens[i:i+n])
                
                # Rule 1: Aceptar si contiene un núcleo de condición
                if gram_tokens & CONDITION_CORES:
                    ngrams.append(" ".join(tokens[i:i+n]))
                # Rule 2: Aceptar si contiene un modificador AND una palabra compatible
                elif (gram_tokens & CONDITION_MODIFIERS) and (gram_tokens & COMPATIBLE_WORDS):
                    ngrams.append(" ".join(tokens[i:i+n]))
                    
    return ngrams

def process_large_json(file_path):
    """
    Reads the 10GB JSON without loading everything in memory.
    Since it is formatted consistently as:
      "sample_texts": [
        "Text 1",
        "Text 2"
      ]
    We can read line by line and extract the strings directly by Regex.
    """
    print(f"Abriendo archivo consolidado: {file_path.name}")
    
    ngram_counts = Counter()
    total_lines = 0
    t0 = time.time()
    
    # Regex to capture strings that look like elements of the "sample_texts" array
    # Usually they look like:       "The abstract text...",
    text_pattern = re.compile(r'^\s*"(.+)",?$')
    in_sample_texts = False
    
    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            total_lines += 1
            if total_lines % 5_000_000 == 0:
                print(f"  ... procesadas {total_lines:,} líneas ({len(ngram_counts):,} n-gramas únicos hasta ahora)")
            
            # Very basic state machine
            if '"sample_texts": [\n' in line or '"sample_texts": [' in line:
                in_sample_texts = True
                continue
            elif in_sample_texts and line.strip() == '],':
                in_sample_texts = False
                continue
            elif in_sample_texts and line.strip() == ']':
                in_sample_texts = False
                continue
                
            if in_sample_texts:
                match = text_pattern.match(line)
                if match:
                    text = match.group(1)
                    # To avoid massive memory leak, we update counts directly
                    ngrams = get_ngrams(text)
                    ngram_counts.update(ngrams)
                    
    print(f"Lectura completa. Tiempo: {time.time()-t0:.1f}s. Líneas totales: {total_lines:,}")
    return ngram_counts

def main():
    print("=" * 60)
    print(" PASO 5: CONSTRUCCIÓN DEL LÉXICO DE CONDICIONES EXPERIMENTALES ")
    print("=" * 60)
    
    if not INPUT_FILE.exists():
        print(f"❌ Error: No se encontró {INPUT_FILE}")
        return
        
    counts = process_large_json(INPUT_FILE)
    
    print(f"Total de n-gramos extraídos: {sum(counts.values()):,}")
    print(f"N-gramos únicos: {len(counts):,}")
    
    # Filter by min frequency to reduce noise and file size
    MIN_FREQ = 3
    filtered_counts = {k: v for k, v in counts.items() if v >= MIN_FREQ}
    print(f"N-gramos con frecuencia >= {MIN_FREQ}: {len(filtered_counts):,}")
    
    # Guardar índice
    print(f"Guardando diccionario léxico en {OUTPUT_FILE.name}...")
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        json.dump(filtered_counts, f, indent=2, ensure_ascii=False)
        
    print("✅ ¡Léxico construido con éxito!")

if __name__ == "__main__":
    main()
