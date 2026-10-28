"""
Agentic pipeline client: drives statgen-mcp tools using a self-hosted LLM.
Connects to the local statgen-mcp server via stdio, then calls the configured
LLM endpoint (e.g. Qwen3-4B-AWQ on EC2 via vLLM) in an agentic loop.

Setup:
    pip install -e ".[llm]"

Configure ~/.statgen-mcp/config.toml:
    [llm]
    base_url = "http://<ec2-ip>:8000/v1"
    model    = "Qwen/Qwen3-4B-AWQ"
    api_key  = "none"

Usage:
    python deploy/run_pipeline.py "run data QC on /path/to/data with maf 0.01 and geno 0.05"
    python deploy/run_pipeline.py --dry-run "run PCA on /path/to/qc_output"
"""

import argparse
import asyncio
import json
import sys

from openai import AsyncOpenAI
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from statgen_mcp.config import get_config

SYSTEM_PROMPT = """\
You are a statistical genetics assistant with access to PLINK2 and bcftools tools.
Complete the user's analysis task by calling tools in the correct pipeline order.
Think step by step. After each tool result, decide what to do next.
If a tool returns an error, diagnose it and adjust your approach.
When the task is complete, summarize what was done and where the output files are.\
"""


async def run(task: str, dry_run: bool = False) -> None:
    cfg = get_config()

    if not cfg.llm.base_url:
        print("Error: [llm] base_url not configured in ~/.statgen-mcp/config.toml")
        print("See config.example.toml for the [llm] section.")
        sys.exit(1)

    llm = AsyncOpenAI(base_url=cfg.llm.base_url, api_key=cfg.llm.api_key, timeout=300.0)
    model = cfg.llm.model

    if dry_run:
        task = task + " (use dry_run=true for all tool calls)"

    server_params = StdioServerParameters(
        command="statgen-mcp",
        args=[],
    )

    async with stdio_client(server_params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()

            tools_result = await session.list_tools()
            openai_tools = [_to_openai_tool(t) for t in tools_result.tools]

            print(f"statgen-mcp connected  ({len(openai_tools)} tools)")
            print(f"model: {model}")
            print(f"task:  {task}\n")
            print("-" * 60)

            messages = [
                {"role": "system", "content": SYSTEM_PROMPT},
                {"role": "user",   "content": task},
            ]

            step = 0
            max_steps = 20
            while step < max_steps:
                response = await llm.chat.completions.create(
                    model=model,
                    messages=messages,
                    tools=openai_tools,
                    tool_choice="auto",
                    max_tokens=512,
                )

                msg = response.choices[0].message

                # Build assistant message dict explicitly (safer across SDK versions)
                assistant_entry: dict = {"role": "assistant", "content": msg.content or ""}
                if msg.tool_calls:
                    assistant_entry["tool_calls"] = [
                        {
                            "id": tc.id,
                            "type": "function",
                            "function": {
                                "name": tc.function.name,
                                "arguments": tc.function.arguments,
                            },
                        }
                        for tc in msg.tool_calls
                    ]
                messages.append(assistant_entry)

                if msg.content:
                    print(f"\n[model]: {msg.content[:300]}")

                # No tool calls — model is done
                if not msg.tool_calls:
                    print("\n--- Done ---")
                    print(msg.content)
                    break

            else:
                print(f"\n--- Stopped: reached {max_steps} steps ---")

                # Execute tool calls
                step += 1
                for tc in msg.tool_calls:
                    name = tc.function.name
                    try:
                        args = json.loads(tc.function.arguments)
                    except json.JSONDecodeError:
                        args = {}

                    print(f"\nStep {step}: {name}")
                    print(f"  args: {json.dumps(args, indent=4)}")

                    result = await session.call_tool(name, args)
                    result_text = (
                        result.content[0].text
                        if result.content
                        else "(no output)"
                    )

                    # Print a preview — full result goes into messages
                    preview = result_text[:400]
                    if len(result_text) > 400:
                        preview += f"\n  ... ({len(result_text)} chars total)"
                    print(f"  result: {preview}")

                    # Truncate long results before storing to keep context bounded
                    stored_text = result_text[:800] if len(result_text) > 800 else result_text
                    messages.append({
                        "role": "tool",
                        "tool_call_id": tc.id,
                        "content": stored_text,
                    })


def _strip_schema_descriptions(schema: dict) -> dict:
    """Remove description fields from JSON schema to reduce token count."""
    if not isinstance(schema, dict):
        return schema
    result = {k: v for k, v in schema.items() if k != "description"}
    if "properties" in result:
        result["properties"] = {
            k: _strip_schema_descriptions(v)
            for k, v in result["properties"].items()
        }
    return result


def _to_openai_tool(tool) -> dict:
    return {
        "type": "function",
        "function": {
            "name": tool.name,
            "description": (tool.description or "")[:120],
            "parameters": _strip_schema_descriptions(tool.inputSchema),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a statgen-mcp pipeline with a self-hosted LLM")
    parser.add_argument("task", nargs="+", help="Task description in natural language")
    parser.add_argument("--dry-run", action="store_true", help="Pass dry_run=true to all tool calls")
    args = parser.parse_args()

    asyncio.run(run(task=" ".join(args.task), dry_run=args.dry_run))


if __name__ == "__main__":
    main()
