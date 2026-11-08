SHELL := /bin/bash

PYTHON ?= python3
CC ?= cc
GM ?= gm
LAKE ?= lake
OUT := graphenium-out

.DEFAULT_GOAL := verify
.PHONY: verify generator c-search datalog lean clean

# Run all four verification layers in order.
verify:
	@$(MAKE) --no-print-directory generator
	@$(MAKE) --no-print-directory c-search
	@$(MAKE) --no-print-directory datalog
	@$(MAKE) --no-print-directory lean
	@echo
	@echo "All MUB-6 verification layers passed."

# 1. Exact/numerically certified fact and certificate generation.
generator:
	$(PYTHON) python/generate_facts.py

# 2. Exhaustive C search for the 48 Bjorck biunimodular sequences.
c-search: | $(OUT)
	$(CC) -O3 -o $(OUT)/f6_bjorck_search tools/f6_bjorck_search.c -lm
	$(OUT)/f6_bjorck_search > $(OUT)/f6_solutions.txt 2> $(OUT)/f6_search.log
	@count=$$(wc -l < "$(OUT)/f6_solutions.txt" | tr -d ' '); \
	if [ "$$count" -ne 48 ]; then \
		echo "error: expected 48 Bjorck sequences, got $$count" >&2; \
		exit 1; \
	fi; \
	echo "C exhaustive search: $$count Bjorck sequences"

# 3. Graphenium embedded Datalog clique and K_{6,6} verification.
datalog: | $(OUT)
	GM="$(GM)" PYTHON="$(PYTHON)" ./graphenium/verify_embedded.sh

# 4. Lean 4 kernel-checked core build (Mathlib-free).
lean:
	cd lean4 && $(LAKE) build

$(OUT):
	mkdir -p "$(OUT)"

clean:
	rm -rf "$(OUT)" lean4/.lake
