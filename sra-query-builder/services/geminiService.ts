
import { GoogleGenAI, Type } from "@google/genai";
import { SRAInputs, GeneratedQuery } from "../types";

export const generateSRAQuery = async (inputs: SRAInputs): Promise<GeneratedQuery> => {
  const ai = new GoogleGenAI({ apiKey: process.env.API_KEY });
  
  const prompt = `
    Act as a senior bioinformatics expert. 
    Analyze the following natural language request for SRA (Sequence Read Archive) data:
    REQUEST: "${inputs.naturalQuery}"

    TASK:
    1. Extract biological entities: Organism, Tissues, Conditions, Strategies (RNA-Seq, etc).
    2. Generate a professional NCBI E-Search boolean query.
    3. Create a mock dataset of 12 realistic entries that satisfy the request.

    JSON SCHEMA REQUIREMENTS:
    - natural_query: A technical summary of what was understood.
    - esearch_query: A valid NCBI search string (e.g. "Arabidopsis thaliana"[Organism] AND "RNA-Seq"[Strategy]...).
    - results: Array of objects with:
       - bioproject, title, organism, tissue, strategy.
       - libraryCount (10-2000), experimentCount (5-100).
       - timePoints (e.g. [0, 12, 24, 48] hours or days).
       - cultivar, isTimeSeries (boolean).
       - geneStats: Array of 2 objects { condition: string, up: number, down: number }.

    Ensure the mock data is diverse and reflects high-quality peer-reviewed studies.
  `;

  try {
    const response = await ai.models.generateContent({
      model: 'gemini-3-flash-preview',
      contents: prompt,
      config: {
        responseMimeType: 'application/json',
        responseSchema: {
          type: Type.OBJECT,
          properties: {
            natural_query: { type: Type.STRING },
            esearch_query: { type: Type.STRING },
            results: {
              type: Type.ARRAY,
              items: {
                type: Type.OBJECT,
                properties: {
                  bioproject: { type: Type.STRING },
                  title: { type: Type.STRING },
                  organism: { type: Type.STRING },
                  tissue: { type: Type.STRING },
                  strategy: { type: Type.STRING },
                  libraryCount: { type: Type.NUMBER },
                  experimentCount: { type: Type.NUMBER },
                  timePoints: { type: Type.ARRAY, items: { type: Type.NUMBER } },
                  cultivar: { type: Type.STRING },
                  isTimeSeries: { type: Type.BOOLEAN },
                  geneStats: {
                    type: Type.ARRAY,
                    items: {
                      type: Type.OBJECT,
                      properties: {
                        condition: { type: Type.STRING },
                        up: { type: Type.NUMBER },
                        down: { type: Type.NUMBER }
                      }
                    }
                  }
                }
              }
            }
          },
          required: ['natural_query', 'esearch_query', 'results']
        }
      }
    });

    return JSON.parse(response.text || '{}') as GeneratedQuery;
  } catch (error) {
    throw new Error("NCBI NLP Interface failed to interpret query. Check syntax.");
  }
};
