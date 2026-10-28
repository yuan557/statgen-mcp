#!/bin/bash
# Start llama.cpp inference server on a CPU instance (e.g. c7i.2xlarge)
# Uses Qwen3-4B GGUF Q4 - ~2.5GB download, ~5-10 tok/s on 8 vCPUs
# OpenAI-compatible API on port 8000, same as the GPU/vLLM setup
#
# Usage: bash deploy/start_llama_cpp.sh
# First run downloads the model from HuggingFace (~2.5GB)

set -e

MODEL_FILE="Qwen3-1.7B-Q4_K_M.gguf"
MODEL_REPO="unsloth/Qwen3-1.7B-GGUF"
N_THREADS=${N_THREADS:-8}   # match vCPU count; override with: N_THREADS=4 bash start_llama_cpp.sh

source llama_env/bin/activate

# download model if not already present
if [ ! -f "$MODEL_FILE" ]; then
    echo "Downloading $MODEL_FILE from HuggingFace..."
    huggingface-cli download "$MODEL_REPO" "$MODEL_FILE" --local-dir .
fi

echo "Starting llama.cpp server on port 8000 (threads: $N_THREADS)"
python -m llama_cpp.server \
    --model "$MODEL_FILE" \
    --n_ctx 8192 \
    --n_threads "$N_THREADS" \
    --host 0.0.0.0 \
    --port 8000
