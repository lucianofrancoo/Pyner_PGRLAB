#!/bin/bash

# Test Pyner Query Generator desde terminal
# No necesita navegador

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║          PYNER QUERY GENERATOR - TERMINAL TEST                 ║"
echo "╚════════════════════════════════════════════════════════════════╝"
echo ""

# 1. Health Check
echo "1️⃣  Health Check..."
curl -s http://localhost:8000/ | jq .
echo ""

# 2. Stats
echo "2️⃣  System Stats..."
curl -s http://localhost:8000/stats | jq .
echo ""

# 3. Queries de prueba
echo "3️⃣  Test Queries..."
echo ""

# Query 1
echo "   📝 Query 1: Arabidopsis thaliana drought stress"
curl -s -X POST http://localhost:8000/generate \
  -H "Content-Type: application/json" \
  -d '{"query": "Arabidopsis thaliana drought stress RNA-Seq", "use_llm": false}' | jq '.ncbi_query'
echo ""

# Query 2
echo "   📝 Query 2: Homo sapiens cancer"
curl -s -X POST http://localhost:8000/generate \
  -H "Content-Type: application/json" \
  -d '{"query": "Homo sapiens cancer RNA-Seq", "use_llm": false}' | jq '.ncbi_query'
echo ""

# Query 3
echo "   📝 Query 3: Mus musculus liver transcriptome"
curl -s -X POST http://localhost:8000/generate \
  -H "Content-Type: application/json" \
  -d '{"query": "Mus musculus liver hepatocyte transcriptome RNA-Seq ChIP-Seq", "use_llm": false}' | jq '.ncbi_query'
echo ""

echo "╔════════════════════════════════════════════════════════════════╗"
echo "║                      ✅ TEST COMPLETE                          ║"
echo "╚════════════════════════════════════════════════════════════════╝"
