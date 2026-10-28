#!/bin/bash
# Quick smoke test for the LLM endpoint.
# Usage: bash deploy/test_llm.sh <ec2-ip>
# Or set base_url in ~/.statgen-mcp/config.toml and run without args.

HOST=${1:-localhost}
PORT=${2:-8000}

curl -s http://$HOST:$PORT/v1/chat/completions \
  -H "Content-Type: application/json" \
  -d '{
    "messages": [{"role": "user", "content": "Reply with one word: working"}],
    "max_tokens": 200
  }' | python3 -c "import sys,json; r=json.load(sys.stdin); print(r['choices'][0]['message']['content'])"
