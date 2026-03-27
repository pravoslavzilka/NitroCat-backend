import json
from query import query_enzymes
from uniprot import enrich_results

output   = query_enzymes("CC(O)C", "CC(=O)C")
enriched = enrich_results(output)

# Print full data for the first result
print(json.dumps(enriched["result"], indent=2))
