
import os
import os.path
import re
import glob
import yaml
import hashlib # md5

configfile: "config/config.yaml"

include: "workflow/rules/archr.smk" 

rule all:
	input:
		'out/archr/Atlas/ATAC/enrich_motifs.rds'