#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BACKEND_DIR="$ROOT_DIR/backend"
FRONTEND_DIR="$ROOT_DIR/frontend"
BACKEND_VENV="$BACKEND_DIR/.venv"

BACKEND_HOST="${BACKEND_HOST:-0.0.0.0}"
BACKEND_PORT="${BACKEND_PORT:-8010}"
OLLAMA_URL="${OLLAMA_URL:-http://localhost:11434}"

is_port_in_use() {
  local port="$1"
  ss -ltn 2>/dev/null | awk 'NR>1 {print $4}' | grep -Eq "(^|[:.])${port}$"
}

resolve_backend_port() {
  local desired_port="$1"
  local candidate="$desired_port"
  local attempts=0

  while (( attempts < 20 )); do
    if ! is_port_in_use "$candidate"; then
      echo "$candidate"
      return 0
    fi
    candidate=$((candidate + 1))
    attempts=$((attempts + 1))
  done

  return 1
}

# Si no se pasa OLLAMA_MODEL por entorno, pedirlo en modo interactivo.
if [[ -z "${OLLAMA_MODEL:-}" ]]; then
  if [[ -t 0 ]]; then
    DEFAULT_OLLAMA_MODEL="qwen2.5:14b"
    OLLAMA_MODELS=()

    if command -v ollama >/dev/null 2>&1; then
      while IFS= read -r MODEL_NAME; do
        [[ -n "$MODEL_NAME" ]] && OLLAMA_MODELS+=("$MODEL_NAME")
      done < <(ollama list 2>/dev/null | awk 'NR>1 {print $1}')
    fi

    echo "[config] Modelo Ollama por defecto: ${DEFAULT_OLLAMA_MODEL}"
    if (( ${#OLLAMA_MODELS[@]} > 0 )); then
      echo "[config] Modelos disponibles en tu Ollama:"
      for i in "${!OLLAMA_MODELS[@]}"; do
        printf "  [%d] %s\n" "$((i + 1))" "${OLLAMA_MODELS[$i]}"
      done
      echo "  [manual] Escribe otro nombre de modelo"
      read -r -p "Elige número o nombre [${DEFAULT_OLLAMA_MODEL}]: " INPUT_OLLAMA_MODEL

      if [[ -z "${INPUT_OLLAMA_MODEL}" ]]; then
        OLLAMA_MODEL="$DEFAULT_OLLAMA_MODEL"
      elif [[ "${INPUT_OLLAMA_MODEL}" =~ ^[0-9]+$ ]] && (( INPUT_OLLAMA_MODEL >= 1 && INPUT_OLLAMA_MODEL <= ${#OLLAMA_MODELS[@]} )); then
        OLLAMA_MODEL="${OLLAMA_MODELS[$((INPUT_OLLAMA_MODEL - 1))]}"
      else
        OLLAMA_MODEL="${INPUT_OLLAMA_MODEL}"
      fi
    else
      echo "[config] No pude listar modelos con 'ollama list'."
      read -r -p "Elige modelo Ollama (ej: qwen3.5:9b) [${DEFAULT_OLLAMA_MODEL}]: " INPUT_OLLAMA_MODEL
      OLLAMA_MODEL="${INPUT_OLLAMA_MODEL:-$DEFAULT_OLLAMA_MODEL}"
    fi
  else
    OLLAMA_MODEL="qwen2.5:14b"
  fi
fi

if [[ ! -d "$BACKEND_DIR" || ! -d "$FRONTEND_DIR" ]]; then
  echo "Error: no se encontraron backend/ y frontend/ dentro de minero_web."
  exit 1
fi

if [[ ! -d "$BACKEND_VENV" ]]; then
  echo "[setup] Creando entorno virtual del backend..."
  python3 -m venv "$BACKEND_VENV"
fi

ORIGINAL_BACKEND_PORT="$BACKEND_PORT"
if RESOLVED_BACKEND_PORT="$(resolve_backend_port "$BACKEND_PORT")"; then
  BACKEND_PORT="$RESOLVED_BACKEND_PORT"
  if [[ "$BACKEND_PORT" != "$ORIGINAL_BACKEND_PORT" ]]; then
    echo "[config] Puerto $ORIGINAL_BACKEND_PORT ocupado. Usando backend port $BACKEND_PORT."
  fi
else
  echo "Error: no encontré puerto libre para backend (probé desde $BACKEND_PORT)."
  exit 1
fi

echo "[setup] Instalando dependencias del backend..."
"$BACKEND_VENV/bin/pip" install -r "$BACKEND_DIR/requirements.txt"

if [[ ! -d "$FRONTEND_DIR/node_modules" ]]; then
  echo "[setup] Instalando dependencias del frontend..."
  (cd "$FRONTEND_DIR" && npm install)
fi

cleanup() {
  echo
  echo "[stop] Cerrando servicios..."
  [[ -n "${BACKEND_PID:-}" ]] && kill "$BACKEND_PID" 2>/dev/null || true
  [[ -n "${FRONTEND_PID:-}" ]] && kill "$FRONTEND_PID" 2>/dev/null || true
}
trap cleanup EXIT INT TERM

echo "[start] Backend en http://${BACKEND_HOST}:${BACKEND_PORT}"
echo "[start] Ollama URL: ${OLLAMA_URL}"
echo "[start] Ollama model: ${OLLAMA_MODEL}"
(
  cd "$BACKEND_DIR"
  exec env \
    OLLAMA_HOST="$OLLAMA_URL" \
    OLLAMA_BASE_URL="$OLLAMA_URL" \
    OLLAMA_MODEL="$OLLAMA_MODEL" \
    "$BACKEND_VENV/bin/uvicorn" main:app --host "$BACKEND_HOST" --port "$BACKEND_PORT" --reload
) &
BACKEND_PID=$!

BACKEND_HEALTH_URL="http://127.0.0.1:${BACKEND_PORT}/api/minero/health"
BACKEND_READY=0
for _ in {1..40}; do
  if curl -fsS -m 2 "$BACKEND_HEALTH_URL" >/dev/null 2>&1; then
    BACKEND_READY=1
    break
  fi
  if ! kill -0 "$BACKEND_PID" >/dev/null 2>&1; then
    break
  fi
  sleep 1
done

if [[ "$BACKEND_READY" -ne 1 ]]; then
  echo "[error] Backend no respondió en ${BACKEND_HEALTH_URL}."
  echo "[error] Revisa logs del backend arriba (o ejecuta uvicorn manualmente para ver el traceback)."
  exit 1
fi
echo "[start] Backend listo: ${BACKEND_HEALTH_URL}"

if [[ -n "${VITE_MINERO_API_URL:-}" ]]; then
  echo "[start] Frontend con API override: ${VITE_MINERO_API_URL}"
  (
    cd "$FRONTEND_DIR"
    exec env VITE_MINERO_API_URL="$VITE_MINERO_API_URL" npm run dev
  ) &
else
  FRONTEND_PROXY_TARGET="http://127.0.0.1:${BACKEND_PORT}"
  echo "[start] Frontend con proxy dinámico: /api -> ${FRONTEND_PROXY_TARGET}"
  (
    cd "$FRONTEND_DIR"
    exec env VITE_PROXY_TARGET="$FRONTEND_PROXY_TARGET" npm run dev
  ) &
fi
FRONTEND_PID=$!

wait "$BACKEND_PID" "$FRONTEND_PID"
