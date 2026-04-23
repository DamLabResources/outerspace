# Conda environment for this clone (e.g. /home/.../outerspace/venv)
VENV := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/venv)
PY := $(VENV)/bin/python
PIP := $(VENV)/bin/pip

DIR_RB_READS := /data/share/nonn-lab/rachel-test-crispr/reads
P1 := $(DIR_RB_READS)/409-4_S1_L001_R1_001.fastq.gz
P2 := $(DIR_RB_READS)/409-4_S1_L001_R2_001.fastq.gz

run:
	bin/main.py $(P1) $(P2)

# Default: unit tests only (see pyproject.toml / tox.ini)
test:
	$(PY) -m tox -e py

ruff:
	$(PY) -m tox -e ruff

format-check: ruff

venv:
	conda create -y -p $(VENV) python=3.12
	$(PIP) install -e ".[dev,pipeline]"
	$(VENV)/bin/pre-commit install

files:
	@ls $(P1)
	@ls $(P2)

py:
	find original tests outerspace bin -name '*.py' -type f | grep -v checkpoint

md:
	find . -type f -name '*.md'

vim:
	vim -p outerspace/extraction_attempt.py bin/main.py

run_script:
	cd outerspace; python extraction_attempt.py

reads:
	find $(DIR_RB_READS)/reads

hello:
	ls -l $(DIR_RB_READS)

clean:
	rm -f outerspace/.ipynb_checkpoints/extraction_attempt-checkpoint.py
	rm -f testing.cfg
	rm -rf outerspace.egg-info

clobber:
	make clean
	rm -rf outdir
	rm -rf venv

RB:
	findseq rb.cfg -1 reads_sample/409-4_S1_L002_R1_001.fastq.gz -2 reads_sample/409-4_S1_L002_R2_001.fastq.gz -o 409-4_S1_L002_R1_R2_output.csv
    
# Unit-test coverage (same selection as `make test` / tox -e py)
coverage:
	$(PY) -m coverage run -m pytest

report: 
	$(PY) -m coverage report -m
