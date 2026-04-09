import time
import hashlib
import os
import requests
from concurrent.futures import ThreadPoolExecutor
import json
import requests, re
import xml.etree.ElementTree as ET

# Set these after registering at brenda-enzymes.org
BRENDA_EMAIL    = os.getenv("BRENDA_EMAIL", "")
BRENDA_PASSWORD = os.getenv("BRENDA_PASSWORD", "")


# ── UniProt ───────────────────────────────────────────────────────────────────

def _fetch_uniprot(uniprot_id: str) -> dict:
    try:
        resp = requests.get(
            f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json",
            timeout=10,
        )
        if resp.status_code == 404:
            return {"_error": f"{uniprot_id} not found in UniProt"}
        if resp.status_code != 200:
            return {"_error": f"UniProt returned HTTP {resp.status_code}"}
        return resp.json()
    except requests.exceptions.Timeout:
        return {"_error": "UniProt request timed out"}
    except Exception as e:
        return {"_error": str(e)}


def _parse_uniprot(uniprot_id: str, data: dict) -> dict:
    if "_error" in data:
        return {"uniprot_id": uniprot_id, "error": data["_error"]}

    result = {"uniprot_id": uniprot_id}

    # Protein name — try recommendedName first, fall back to submittedName
    try:
        result["protein_name"] = (
            data["proteinDescription"]["recommendedName"]["fullName"]["value"]
        )
    except (KeyError, IndexError):
        try:
            result["protein_name"] = (
                data["proteinDescription"]["submittedNames"][0]["fullName"]["value"]
            )
        except (KeyError, IndexError):
            result["protein_name"] = None

    # Organism
    try:
        result["organism"] = data["organism"]["scientificName"]
    except KeyError:
        result["organism"] = None

    # EC number — try recommendedName, then submittedNames, then protein existence annotation
    ec = None
    try:
        ec = data["proteinDescription"]["recommendedName"]["ecNumbers"][0]["value"]
    except (KeyError, IndexError):
        pass
    if not ec:
        try:
            ec = data["proteinDescription"]["submittedNames"][0]["ecNumbers"][0]["value"]
        except (KeyError, IndexError):
            pass
    if not ec:
        try:
            name = result.get("protein_name", "") or ""
            import re
            match = re.search(r"EC[: ]+(\d+\.\d+\.\d+\.[\d-]+)", name)
            if match:
                ec = match.group(1)
        except Exception:
            pass
    result["ec_number"] = ec

    result["kcat"] = None
    result["km"]   = None

    for comment in data.get("comments", []):
        ctype = comment.get("commentType", "")

        if ctype == "KINETICS":
            kms, kcats = [], []
            for item in comment.get("michaelisConstants", []):
                kms.append({
                    "value":     item.get("constant"),
                    "unit":      item.get("unit"),
                    "substrate": item.get("substrate"),
                })
            for item in comment.get("maximumVelocities", []):
                kcats.append({
                    "value": item.get("velocity"),
                    "unit":  item.get("unit"),
                })
            if kms:
                result["km"] = kms
            if kcats:
                result["kcat"] = kcats

    # Gene name
    try:
        result["gene_name"] = data["genes"][0]["geneName"]["value"]
    except (KeyError, IndexError):
        result["gene_name"] = None

    # Sequence length
    try:
        result["sequence_length"] = data["sequence"]["length"]
    except KeyError:
        result["sequence_length"] = None

    # Comments — function, optimal pH, optimal temp
    result["function"]     = None
    result["optimal_ph"]   = None
    result["optimal_temp"] = None

    for comment in data.get("comments", []):
        ctype = comment.get("commentType", "")

        if ctype == "FUNCTION" and result["function"] is None:
            try:
                result["function"] = comment["texts"][0]["value"]
            except (KeyError, IndexError):
                pass

        if ctype == "BIOPHYSICOCHEMICAL PROPERTIES":

            # Km values
            kinetics = comment.get("kineticParameters", {})
            km_list  = []
            for mc in kinetics.get("michaelisConstants", []):
                km_list.append({
                    "value":     mc.get("constant"),
                    "unit":      mc.get("unit"),
                    "substrate": mc.get("substrate"),
                })
            if km_list:
                result["km"] = km_list

            # Vmax
            vmax_list = []
            for vmax in kinetics.get("maximumVelocities", []):
                vmax_list.append({
                    "value":  vmax.get("velocity"),
                    "unit":   vmax.get("unit"),
                    "enzyme": vmax.get("enzyme"),
                })
            if vmax_list:
                result["vmax"] = vmax_list

            # pH optimum
            if result["optimal_ph"] is None:
                try:
                    result["optimal_ph"] = comment["phDependence"]["texts"][0]["value"]
                except (KeyError, IndexError):
                    pass

            # Temperature optimum
            if result["optimal_temp"] is None:
                try:
                    result["optimal_temp"] = comment["temperatureDependence"]["texts"][0]["value"]
                except (KeyError, IndexError):
                    pass

    # Links
    result["alphafold_url"] = f"https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
    result["uniprot_url"]   = f"https://www.uniprot.org/uniprotkb/{uniprot_id}"

    return result


# ── BRENDA JSON ───────────────────────────────────────────────────────────────

def _fetch_brenda_json(ec: str, uniprot_id: str) -> dict:
    """
    Fetch enzyme properties from the local BRENDA JSON file via brenda.py.
    Returns a dict with km_q1_mM, kcat_q1_per_s, temperature_median, ph_median.
    Falls back gracefully if brenda.py or the JSON file is unavailable.
    """
    if not ec:
        return {"brenda_json_note": "No EC number — cannot query BRENDA JSON"}

    try:
        from brenda import get_enzyme_properties
    except ImportError:
        return {"brenda_json_note": "brenda.py not found — place it in the same directory"}

    try:
        props = get_enzyme_properties(ec, uniprot_id)
        return {
            "brenda_km_q1_mM":           props.get("km_q1_mM"),
            "brenda_kcat_q1_per_s":      props.get("kcat_q1_per_s"),
            "brenda_temperature_median": props.get("temperature_median"),
            "brenda_ph_median":          props.get("ph_median"),
            "brenda_enzyme_name":        props.get("enzyme_name"),
            "cofactors":                 props.get("cofactor_names", []),
        }
    except ValueError as e:
        # Protein not found in BRENDA JSON — not an error, just missing data
        return {"brenda_json_note": str(e)}
    except Exception as e:
        return {"brenda_json_note": f"BRENDA JSON query failed: {e}"}


# ── BRENDA SOAP (legacy) ──────────────────────────────────────────────────────

def _fetch_brenda(ec_number: str) -> dict:
    if not BRENDA_EMAIL or not BRENDA_PASSWORD:
        return {"brenda_note": "BRENDA credentials not set — register free at brenda-enzymes.org"}
    if not ec_number:
        return {"brenda_note": "No EC number — cannot query BRENDA"}
    try:
        from zeep import Client
    except ImportError:
        return {"brenda_note": "zeep not installed — run: pip install zeep"}

    try:
        pw_hash = hashlib.sha256(BRENDA_PASSWORD.encode()).hexdigest()
        client  = Client("https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl")

        params = f"{BRENDA_EMAIL},{pw_hash},ecNumber*{ec_number}#"

        result = {}

        try:
            rows = client.service.getTurnoverNumber(params)
            result["kcat"] = [
                {"value": r.get("turnoverNumber"), "unit": "1/s",
                 "substrate": r.get("substrate"), "organism": r.get("organism")}
                for r in (rows or []) if r.get("turnoverNumber")
            ][:5]
        except Exception:
            result["kcat"] = None

        try:
            rows = client.service.getKmValue(params)
            result["km"] = [
                {"value": r.get("kmValue"), "unit": "mM",
                 "substrate": r.get("substrate"), "organism": r.get("organism")}
                for r in (rows or []) if r.get("kmValue")
            ][:5]
        except Exception:
            result["km"] = None

        try:
            rows = client.service.getKcatKmValue(params)
            result["kcat_km"] = [
                {"value": r.get("kcatKmValue"), "unit": "1/(mM·s)",
                 "substrate": r.get("substrate"), "organism": r.get("organism")}
                for r in (rows or []) if r.get("kcatKmValue")
            ][:5]
        except Exception:
            result["kcat_km"] = None

        try:
            rows = client.service.getPhOptimum(params)
            result["ph_optimum"] = [
                {"value": r.get("phOptimum"), "organism": r.get("organism")}
                for r in (rows or []) if r.get("phOptimum")
            ][:3]
        except Exception:
            result["ph_optimum"] = None

        try:
            rows = client.service.getTemperatureOptimum(params)
            result["temp_optimum"] = [
                {"value": r.get("temperatureOptimum"), "unit": "°C",
                 "organism": r.get("organism")}
                for r in (rows or []) if r.get("temperatureOptimum")
            ][:3]
        except Exception:
            result["temp_optimum"] = None

        return result

    except Exception as e:
        return {"brenda_note": f"BRENDA query failed: {e}"}


# ── Single enzyme enrichment ──────────────────────────────────────────────────

def _enrich_single(enzyme: dict) -> dict:
    uniprot_id = enzyme["uniprot"]

    # UniProt — always fetch first
    raw    = _fetch_uniprot(uniprot_id)
    result = {**enzyme, **_parse_uniprot(uniprot_id, raw)}

    # SABIO-RK — only fetch if UniProt has no Km
    if not result.get("km"):
        time.sleep(1.0)
        sabio = _fetch_sabiork(uniprot_id)
        for key in ("kcat", "km", "ph_optimum", "temp_optimum", "sabio_entries"):
            if not result.get(key) and sabio.get(key):
                result[key] = sabio[key]

    # BRENDA JSON — fill in missing km/kcat/temp/pH using local flat file
    ec = result.get("ec_number")
    brenda_data = _fetch_brenda_json(ec, uniprot_id)
    result.update(brenda_data)

    # Back-fill km and kcat from BRENDA JSON if still missing
    if not result.get("km") and result.get("brenda_km_q1_mM") is not None:
        result["km"] = [{"value": result["brenda_km_q1_mM"], "unit": "mM", "source": "BRENDA_Q1"}]
    if not result.get("kcat") and result.get("brenda_kcat_q1_per_s") is not None:
        result["kcat"] = [{"value": result["brenda_kcat_q1_per_s"], "unit": "1/s", "source": "BRENDA_Q1"}]
    if not result.get("optimal_temp") and result.get("brenda_temperature_median") is not None:
        result["optimal_temp"] = result["brenda_temperature_median"]
    if not result.get("optimal_ph") and result.get("brenda_ph_median") is not None:
        result["optimal_ph"] = result["brenda_ph_median"]

    return result


# ── Public function ───────────────────────────────────────────────────────────

def enrich_results(query_output: dict, workers: int = 10) -> dict:
    if query_output["status"] == "error":
        return query_output

    enriched = []
    for enzyme in query_output["result"]:
        result = {**enzyme}

        # UniProt
        raw = _fetch_uniprot(enzyme["uniprot"])
        result.update(_parse_uniprot(enzyme["uniprot"], raw))

        # SABIO-RK
        sabio = _fetch_sabiork(enzyme["uniprot"])
        result.update(sabio)

        # BRENDA JSON — always run, fills missing fields
        ec = result.get("ec_number")
        brenda_data = _fetch_brenda_json(ec, enzyme["uniprot"])
        result.update(brenda_data)

        # Back-fill from BRENDA JSON where other sources returned nothing
        if not result.get("km") and result.get("brenda_km_q1_mM") is not None:
            result["km"] = [{"value": result["brenda_km_q1_mM"], "unit": "mM", "source": "BRENDA_Q1"}]
        if not result.get("kcat") and result.get("brenda_kcat_q1_per_s") is not None:
            result["kcat"] = [{"value": result["brenda_kcat_q1_per_s"], "unit": "1/s", "source": "BRENDA_Q1"}]
        if not result.get("optimal_temp") and result.get("brenda_temperature_median") is not None:
            result["optimal_temp"] = result["brenda_temperature_median"]
        if not result.get("optimal_ph") and result.get("brenda_ph_median") is not None:
            result["optimal_ph"] = result["brenda_ph_median"]

        enriched.append(result)

    return {
        "status":   query_output["status"],
        "result":   enriched,
        "comments": query_output["comments"],
    }


def _fetch_sabiork(uniprot_id: str) -> dict:
    """Fetch kinetic parameters from SABIO-RK."""

    try:
        id_resp = requests.get(
            "https://sabiork.h-its.org/sabioRestWebServices/searchKineticLaws/entryIDs",
            params={"format": "txt", "q": f"UniprotID:{uniprot_id}"},
            timeout=15,
        )
        text = id_resp.text.strip()

        if id_resp.status_code != 200 or not text or "no data found" in text.lower():
            return {
                "kcat": None, "km": None,
                "ph_optimum": None, "temp_optimum": None,
                "sabio_note": "Not in SABIO-RK",
            }

        entry_ids = text.split("\n")

    except Exception as e:
        return {"sabio_note": f"SABIO-RK lookup failed: {e}"}

    sbml_resp = None
    for attempt in range(3):
        try:
            time.sleep(0.5 * (attempt + 1))
            sbml_resp = requests.get(
                "https://sabiork.h-its.org/sabioRestWebServices/kineticLaws",
                params={"kinlawids": ",".join(entry_ids[:10])},
                timeout=15,
            )
            if sbml_resp.status_code == 200:
                break
        except Exception as e:
            if attempt == 2:
                return {"sabio_note": f"SABIO-RK SBML fetch failed: {e}"}

    if not sbml_resp or sbml_resp.status_code != 200:
        return {"sabio_note": f"SABIO-RK SBML fetch failed (status {sbml_resp.status_code if sbml_resp else 'no response'})"}

    ET.register_namespace("sbrk", "http://sabiork.h-its.org")
    ns = {
        "sbml": "http://www.sbml.org/sbml/level3/version1/core",
        "sbrk": "http://sabiork.h-its.org",
    }

    try:
        root  = ET.fromstring(sbml_resp.text)
        model = root.find("sbml:model", ns)
    except ET.ParseError as e:
        return {"sabio_note": f"SABIO-RK XML parse error: {e}"}

    kms, kcats, phs, temps = [], [], [], []

    for rxn in model.findall(".//sbml:reaction", ns):
        kl = rxn.find("sbml:kineticLaw", ns)
        if kl is not None:
            for param in kl.findall(".//sbml:localParameter", ns):
                pid  = param.get("id", "").lower()
                val  = param.get("value")
                unit = param.get("units", "")
                unit = unit.replace("swedgeone", "s⁻¹").replace("M", "M")
                if val is None:
                    continue
                try:
                    val = float(val)
                except ValueError:
                    continue
                if "km" in pid:
                    kms.append({"value": val, "unit": unit})
                elif "kcat" in pid:
                    kcats.append({"value": val, "unit": unit})

        temp_el = rxn.find(".//sbrk:startValueTemperature", ns)
        ph_el   = rxn.find(".//sbrk:startValuepH", ns)

        if temp_el is not None and temp_el.text:
            try:
                temps.append(float(temp_el.text))
            except ValueError:
                pass
        if ph_el is not None and ph_el.text:
            try:
                phs.append(float(ph_el.text))
            except ValueError:
                pass

    return {
        "kcat":          kcats            if kcats  else None,
        "km":            kms              if kms    else None,
        "ph_optimum":    list(set(phs))   if phs    else None,
        "temp_optimum":  list(set(temps)) if temps  else None,
        "sabio_entries": len(entry_ids),
    }


if __name__ == "__main__":
    query_output = {
    "status": "success",
    "result": [
        {"uniprot": "P00326"},
        {"uniprot": "P07327"},
        {"uniprot": "Q8KLT9"},
        ],
        "comments": []
    }

    enriched = enrich_results(query_output)

    for enzyme in enriched["result"]:
        print(enzyme)