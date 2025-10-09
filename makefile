DIR_RB_READS := /data/share/nonn-lab/rachel-test-crispr/reads
P1 := $(DIR_RB_READS)/409-4_S1_L001_R1_001.fastq.gz
P2 := $(DIR_RB_READS)/409-4_S1_L001_R2_001.fastq.gz

run:
	bin/main.py $(P1) $(P2)

test:
	pytest

black:
	black outerspace
	black tests

venv:
	# Create a new conda environment in the venv directory
	conda create -p ./venv python=3.12 pytest
	conda run -p ./venv pip install .

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
    
# Running coverage on pytest test scripts
# Do this first- then report below
coverage:
	coverage run -m pytest tests/test_*.py

report: 
	coverage report -m
