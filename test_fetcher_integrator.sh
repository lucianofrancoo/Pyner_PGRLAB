#!/bin/bash
#
# Query Generator → Fetcher_NCBI Integration (Interactive)
# ========================================================
# Starts Phase 3 query generator, confirms query, then fetches all results
# and saves them as a CSV table.
#

set -e

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
OUTPUT_FILE="$ROOT_DIR/test_fetcher_results.csv"

printf "\n==========================================\n"
printf "PYNER INTEGRATION TEST\n"
printf "==========================================\n\n"

read -r -p "Ingresa tu consulta en lenguaje natural: " USER_INPUT
if [ -z "$USER_INPUT" ]; then
    echo "ERROR: Input vacio. Abortando."
    exit 1
fi

printf "\n--- Generando query con Phase 3 ---\n\n"

cd "$ROOT_DIR/Query_generator/phases/phase3"
GEN_OUTPUT=$(python api/main.py --quick "$USER_INPUT")

echo "$GEN_OUTPUT"

QUERY=$(printf "%s\n" "$GEN_OUTPUT" | awk '/NCBI Query:/{getline; gsub(/^ +/, ""); print; exit}')

if [ -z "$QUERY" ]; then
    echo "ERROR: No se pudo extraer el query."
    exit 1
fi


read -r -p "Esta bien este query? (s/si): " CONFIRM
CONFIRM=$(echo "$CONFIRM" | tr '[:upper:]' '[:lower:]')
if [ -n "$CONFIRM" ] && [ "$CONFIRM" != "s" ] && [ "$CONFIRM" != "si" ] && [ "$CONFIRM" != "y" ] && [ "$CONFIRM" != "yes" ]; then
    echo "ERROR: Cancelado por el usuario."
    exit 1
fi

printf "\n--- Consultando NCBI (sin limite de resultados) ---\n\n"

cd "$ROOT_DIR/Fetcher_NCBI"
python main.py -q "$QUERY" -o "$OUTPUT_FILE" --min-bioprojects 30

printf "\nOK: Resultados guardados en: %s\n" "$OUTPUT_FILE"
