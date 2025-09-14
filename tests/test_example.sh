#!/usr/bin/env bash
set -e

# Test 1: check pixi.toml exists
[ -f "pixi.toml" ] || { echo "pixi.toml missing!"; exit 1; }

# Test 2: simple string check
echo "pixi is great" | grep -q "pixi" || { echo "String test failed"; exit 1; }

echo "✅ Shell tests passed"
