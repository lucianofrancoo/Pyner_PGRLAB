# Query_generator — Frontend Integration (npm run dev)

**Propósito**: Instrucciones simples para que un frontend developer implemente una interfaz web que consume la API de generación de consultas booleanas NCBI.

---

## **SETUP RÁPIDO (5 minutos)**

### **Paso 1: Inicia el servidor backend**
```bash
cd /home/lahumada/disco1/Pyner_PGRLAB
python3 Query_generator/wrapper_api.py
```
Output esperado:
```
INFO:     Uvicorn running on http://127.0.0.1:8001
Press CTRL+C to quit
```

### **Paso 2: Crea tu proyecto frontend**
```bash
# Opción A: Vite + React
npm create vite@latest mi-buscador -- --template react
cd mi-buscador
npm install
npm run dev

# Opción B: Next.js (más fácil si ya conoces React)
npx create-next-app@latest
npm run dev
```

### **Paso 3: Consume la API desde tu frontend**
```javascript
// Ejemplo básico con fetch
async function generateQuery(userQuery) {
  const response = await fetch('http://127.0.0.1:8001/generate', {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify({ 
      query: userQuery,
      top_k: 5,
      expand: false 
    })
  });
  return await response.json();
}

// Uso
const result = await generateQuery("Arabidopsis root drought");
console.log(result.boolean_query);  // (("Arabidopsis"[Organism]) AND ("root" OR "drought"))
```

---

## **¿Qué es esta API?**

- **URL**: `http://127.0.0.1:8001/generate`
- **Método**: POST
- **Propósito**: Convierte consultas en lenguaje natural → consulta booleana NCBI lista para usar
- **Ejemplo**: 
  - Input: `"me gustaría buscar arabidopsis en sus raíces bajo sequía"`
  - Output: `(("Arabidopsis" OR "Arabidopsis thaliana")[Organism]) AND ("root" OR "drought")`

---

## **Backend: Cómo funciona**

- Input (natural language) → embeddings → FAISS vector search → KB enrichment (organisms, estrategias) → extracción de keywords → generador de boolean query
- Los keywords se extraen automáticamente del texto del usuario (soporta inglés y español)
---

## **API Reference**

### **POST /generate**

**Request:**
```json
{
  "query": "Mouse RNAseq liver",
  "top_k": 5,
  "expand": false
}
```

**Response:**
```json
{
  "input": "Mouse RNAseq liver",
  "boolean_query": "(\"Mus musculus\"[Organism]) AND (\"RNA-Seq\" OR \"RNA sequencing\")[Strategy]",
  "suggestions": [
    {
      "query_text": "Research studies on Mus musculus",
      "query_type": "organism",
      "similarity_score": 0.52,
      "kb_data": {"organism": "Mus musculus"}
    }
  ],
  "note": "Filtered 2 irrelevant results"
}
```

---

## **Ejemplo Frontend Completo (React)**

### **1. Component React Simple**
```jsx
// QueryBuilder.jsx
import { useState } from 'react';

export default function QueryBuilder() {
  const [query, setQuery] = useState('');
  const [result, setResult] = useState(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleGenerate = async (e) => {
    e.preventDefault();
    setLoading(true);
    setError(null);

    try {
      const response = await fetch('http://127.0.0.1:8001/generate', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ 
          query: query,
          top_k: 5,
          expand: false 
        })
      });

      if (!response.ok) throw new Error('API error');
      const data = await response.json();
      setResult(data);
    } catch (err) {
      setError(err.message);
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container">
      <h1>🔍 NCBI Query Generator</h1>
      
      <form onSubmit={handleGenerate}>
        <textarea
          value={query}
          onChange={(e) => setQuery(e.target.value)}
          placeholder="Ej: Arabidopsis en sus raíces bajo sequía..."
          rows={4}
        />
        <button type="submit" disabled={loading}>
          {loading ? 'Generando...' : 'Generar Query'}
        </button>
      </form>

      {error && <p className="error">Error: {error}</p>}

      {result && (
        <div className="results">
          <h2>Resultado:</h2>
          
          <div className="boolean-box">
            <h3>Boolean Query (copiar a NCBI):</h3>
            <code>{result.boolean_query}</code>
            <button onClick={() => navigator.clipboard.writeText(result.boolean_query)}>
              📋 Copiar
            </button>
          </div>

          {result.suggestions.length > 0 && (
            <div className="suggestions">
              <h3>Sugerencias:</h3>
              <ul>
                {result.suggestions.map((s, i) => (
                  <li key={i}>
                    <strong>{s.query_type.toUpperCase()}</strong> 
                    {s.query_text} 
                    <span className="score">{s.similarity_score.toFixed(3)}</span>
                  </li>
                ))}
              </ul>
            </div>
          )}
        </div>
      )}
    </div>
  );
}
```

### **2. Estilos (CSS)**
```css
.container {
  max-width: 900px;
  margin: 0 auto;
  padding: 20px;
  font-family: Arial, sans-serif;
}

h1 { color: #333; }

textarea {
  width: 100%;
  padding: 12px;
  border: 1px solid #ddd;
  border-radius: 4px;
  font-size: 16px;
}

button {
  margin-top: 10px;
  padding: 10px 20px;
  background: #007bff;
  color: white;
  border: none;
  border-radius: 4px;
  cursor: pointer;
  font-size: 16px;
}

button:hover { background: #0056b3; }
button:disabled { background: #ccc; cursor: not-allowed; }

.boolean-box {
  background: #f0f0f0;
  padding: 15px;
  border-radius: 4px;
  margin: 20px 0;
}

.boolean-box code {
  display: block;
  background: #fff;
  padding: 15px;
  border-radius: 4px;
  overflow-x: auto;
  font-family: 'Courier New', monospace;
  word-break: break-all;
  margin: 10px 0;
}

.suggestions {
  margin-top: 20px;
}

.suggestions li {
  padding: 8px;
  border-bottom: 1px solid #eee;
  display: flex;
  justify-content: space-between;
  align-items: center;
}

.score {
  background: #e7f3ff;
  padding: 4px 8px;
  border-radius: 3px;
  font-size: 12px;
  color: #0056b3;
}

.error {
  color: red;
  padding: 10px;
  background: #ffe7e7;
  border-radius: 4px;
}
```

---

## **Pasos por pasos: Setup + Deploy Local**

### **Para Linux/Mac:**
```bash
# 1. Terminal 1: Inicia el backend
cd /home/lahumada/disco1/Pyner_PGRLAB
python3 Query_generator/wrapper_api.py

# 2. Terminal 2: Crea tu frontend
npm create vite@latest mi-buscador -- --template react
cd mi-buscador
npm install
npm run dev

# 3. Abre http://localhost:5173 en el navegador
```

### **Para Windows (PowerShell):**
```powershell
# 1. PowerShell 1
cd C:\ruta\a\Pyner_PGRLAB
python Query_generator\wrapper_api.py

# 2. PowerShell 2
npm create vite@latest mi-buscador -- --template react
cd mi-buscador
npm install
npm run dev
```

---

## **Solución de Errores Comunes**

| Error | Solución |
|-------|----------|
| `CORS error` | Asegúrate de que `wrapper_api.py` esté corriendo + usa `http://127.0.0.1:8001` exacto |
| `Connection refused` | Backend no está iniciado. Verifica que Python esté corriendo en puerto 8001 |
| `No keywords extracted` | Normal si la consulta no tiene términos biológicos conocidos. La query se genera igual con organismos. |
| `Model not found (Ollama)` | Es OK - usa fallback local. No necesitas Ollama instalado |

---

## **Tus Responsabilidades + Backend**

### **Tú (Frontend developer):**
- ✅ HTML/CSS/JavaScript atractivo
- ✅ Inputs para el usuario 
- ✅ Llamadas fetch a `http://127.0.0.1:8001/generate`
- ✅ Mostrar resultados + botón para copiar

### **Backend (ya hecho):**
- ✅ Extrae keywords del texto natural (inglés + español)
- ✅ Busca en vector DB (FAISS)
- ✅ Valida organismos contra KB
- ✅ Genera boolean query listo para NCBI
- ✅ API en puerto 8001

---

## **Estructura de Archivos Backend (FYI)**

```
Query_generator/
  ├── generate_boolean_query.py    # CLI + función core
  ├── wrapper_api.py               # FastAPI server (puerto 8001)
  ├── FRONTEND_INTEGRATION.md      # Este archivo
  └── phases/
      ├── phase1/output/stage3_knowledge_base.json
      └── phase2/data/pyner_vectors.faiss
```

---

## **Requirements Backend**
```python
# Esto ya está instalado pero por si acaso:
fastapi
uvicorn
sentence-transformers
faiss-cpu  # o faiss-gpu
```

Para verificar:
```bash
python3 -c "import fastapi, uvicorn, faiss; print('✅ All OK')"
```

---

**¿Preguntas?** Revisa los logs de `wrapper_api.py` si algo falla. El backend logea todas las queries y errores.

12) Notes for frontend design team
- We only provide backend logic and boolean generation. The frontend should:
  - Request the boolean from the backend, display to user, and use it to call NCBI or show in the SRA Advanced Search form.
  - Optionally let users edit the boolean before submission.

Contact: The KB and FAISS index used are under `Query_generator/phases/phase1` and `Query_generator/phases/phase2` respectively.
