#!/bin/bash
# Start vLLM + GPU monitor on EC2 (g4dn.xlarge, T4 16GB)
# Usage: bash deploy/start_vllm.sh
# Model is downloaded from HuggingFace on first run (~2.5GB)

source vllm_env/bin/activate

# start GPU monitor in background
uvicorn gpu_monitor:app --app-dir deploy --host 0.0.0.0 --port 8001 &
MONITOR_PID=$!
echo "GPU monitor started on port 8001 (PID: $MONITOR_PID)"

# start vLLM (foreground)
vllm serve Qwen/Qwen3-4B-AWQ \
    --quantization awq \
    --dtype half \
    --max-model-len 8192 \
    --gpu-memory-utilization 0.85 \
    --max-num-seqs 8 \
    --host 0.0.0.0 \
    --port 8000

# cleanup: kill monitor when vLLM exits
kill $MONITOR_PID 2>/dev/null
