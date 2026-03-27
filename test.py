import requests, json

resp = requests.get("https://rest.uniprot.org/uniprotkb/P00330.json", timeout=10)
data = resp.json()

for c in data.get("comments", []):
    if c["commentType"] == "BIOPHYSICOCHEMICAL PROPERTIES":
        print(json.dumps(c, indent=2))