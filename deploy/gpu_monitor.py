"""
GPU monitor API: serves NVIDIA GPU stats and vLLM metrics.
Runs alongside vLLM on EC2 (port 8001).
"""

from contextlib import asynccontextmanager
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pynvml import (
    nvmlInit,
    nvmlShutdown,
    nvmlDeviceGetHandleByIndex,
    nvmlDeviceGetMemoryInfo,
    nvmlDeviceGetPowerUsage,
    nvmlDeviceGetName,
    nvmlDeviceGetTemperature,
    nvmlDeviceGetUtilizationRates,
    nvmlDeviceGetClockInfo,
    NVML_TEMPERATURE_GPU,
    NVML_CLOCK_SM,
    NVML_CLOCK_MEM,
)
import httpx

VLLM_METRIC_URL = "http://localhost:8000/metrics"

handle = None


@asynccontextmanager
async def lifespan(app: FastAPI):
    global handle
    nvmlInit()
    handle = nvmlDeviceGetHandleByIndex(0)
    yield
    nvmlShutdown()


app = FastAPI(title="GPU Monitor", lifespan=lifespan)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["GET"],
    allow_headers=["*"],
)


@app.get("/stats")
def gpu_stats():
    mem = nvmlDeviceGetMemoryInfo(handle)
    util = nvmlDeviceGetUtilizationRates(handle)
    temp = nvmlDeviceGetTemperature(handle, NVML_TEMPERATURE_GPU)
    power = nvmlDeviceGetPowerUsage(handle)
    name = nvmlDeviceGetName(handle)
    clock_sm = nvmlDeviceGetClockInfo(handle, NVML_CLOCK_SM)
    clock_mem = nvmlDeviceGetClockInfo(handle, NVML_CLOCK_MEM)

    # pynvml returns bytes in older versions
    if isinstance(name, bytes):
        name = name.decode()

    return {
        "gpu_name": name,
        "gpu_util": util.gpu,
        "mem_util": util.memory,
        "vram_used_mb": mem.used // (1024 * 1024),
        "vram_total_mb": mem.total // (1024 * 1024),
        "vram_free_mb": mem.free // (1024 * 1024),
        "temp_c": temp,
        "power_w": round(power / 1000, 1),
        "clock_sm_mhz": clock_sm,
        "clock_mem_mhz": clock_mem,
    }


@app.get("/vllm_metrics")
async def vllm_metrics():
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(VLLM_METRIC_URL, timeout=5.0)
            resp.raise_for_status()
    except Exception as e:
        return {"error": str(e)}

    parsed = {}
    for line in resp.text.splitlines():
        if line.startswith("#"):
            continue
        if "vllm:kv_cache_usage_perc" in line:
            parsed["kv_cache_usage"] = _parse_value(line)
        elif "vllm:num_requests_running" in line:
            parsed["requests_running"] = int(_parse_value(line))
        elif "vllm:num_requests_waiting" in line:
            parsed["requests_waiting"] = int(_parse_value(line))
        elif "vllm:num_requests_swapped" in line:
            parsed["requests_swapped"] = int(_parse_value(line))
        elif "vllm:prompt_tokens_total" in line and "counter" not in line:
            parsed["prompt_tokens_total"] = int(_parse_value(line))
        elif "vllm:generation_tokens_total" in line and "counter" not in line:
            parsed["generation_tokens_total"] = int(_parse_value(line))

    return parsed


@app.get("/health")
def health():
    return {"status": "ok"}


def _parse_value(line: str) -> float:
    parts = line.strip().rsplit(" ", 1)
    if len(parts) == 2:
        try:
            return float(parts[1])
        except ValueError:
            return 0.0
    return 0.0
