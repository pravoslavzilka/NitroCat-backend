import numpy as np
from uniprot import _fetch_uniprot, _parse_uniprot, _fetch_sabiork, _fetch_brenda_json
import re

def _extract_float(value) -> float | None:
    """
    Extract first number from a value that may be a float, int, or string
    like 'Optimum pH is 8.' or '37 degrees' or '6.5 - 7.5'.
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    # Extract first number (including decimals) from string
    match = re.search(r"\d+\.?\d*", str(value))
    return float(match.group()) if match else None


def _unify_enzyme(raw_result: dict) -> dict:
    temp = (
        raw_result.get("brenda_temperature_median")
        or raw_result.get("temp_optimum")
        or raw_result.get("optimal_temp")
    )
    if isinstance(temp, list):
        temp = temp[0].get("value") if temp else None

    ph = (
        raw_result.get("brenda_ph_median")
        or raw_result.get("ph_optimum")
        or raw_result.get("optimal_ph")
    )
    if isinstance(ph, list):
        ph = ph[0].get("value") if ph else None

    km = raw_result.get("brenda_km_q1_mM")
    if km is None:
        km_list = raw_result.get("km")
        if isinstance(km_list, list) and km_list:
            km = km_list[0].get("value")

    kcat = raw_result.get("brenda_kcat_q1_per_s")
    if kcat is None:
        kcat_list = raw_result.get("kcat")
        if isinstance(kcat_list, list) and kcat_list:
            kcat = kcat_list[0].get("value")

    return {
        "uniprot":         raw_result.get("uniprot_id") or raw_result.get("uniprot"),
        "protein_name":    raw_result.get("protein_name"),
        "gene_name":       raw_result.get("gene_name"),
        "organism":        raw_result.get("organism"),
        "ec_number":       raw_result.get("ec_number"),
        "sequence_length": raw_result.get("sequence_length"),
        "function":        raw_result.get("function"),
        "alphafold_url":   raw_result.get("alphafold_url"),
        "uniprot_url":     raw_result.get("uniprot_url"),
        "temperature":     _extract_float(temp),
        "ph":              _extract_float(ph),
        "km_mM":           _extract_float(km),
        "kcat_per_s":      _extract_float(kcat),
    }


def _compute_group_stats(enzymes: list[dict]) -> dict:
    """
    Compute group-level statistics across all enzymes.
    Returns median temperature, median pH, Q1 of KM and Kcat.
    """
    def collect(field):
        return [e[field] for e in enzymes if e.get(field) is not None]

    temps  = collect("temperature")
    phs    = collect("ph")
    kms    = collect("km_mM")
    kcats  = collect("kcat_per_s")

    def stats(values: list) -> dict | None:
        if not values:
            return None
        arr = np.array(values)
        return {
            "median": float(np.median(arr)),
            "q1":     float(np.percentile(arr, 25)),
            "q3":     float(np.percentile(arr, 75)),
            "mean":   float(np.mean(arr)),
            "min":    float(np.min(arr)),
            "max":    float(np.max(arr)),
            "n":      len(values),
        }

    return {
        "temperature": stats(temps),
        "ph":          stats(phs),
        "km_mM":       stats(kms),
        "kcat_per_s":  stats(kcats),
    }


def enrich_results(query_output: dict) -> dict:
    if query_output["status"] == "error":
        return query_output

    raw_results = []

    for enzyme in query_output["result"]:
        result = {**enzyme}

        # UniProt
        raw = _fetch_uniprot(enzyme["uniprot"])
        result.update(_parse_uniprot(enzyme["uniprot"], raw))

        # SABIO-RK
        sabio = _fetch_sabiork(enzyme["uniprot"])
        result.update(sabio)

        # BRENDA JSON
        ec = result.get("ec_number")
        result.update(_fetch_brenda_json(ec, enzyme["uniprot"]))

        raw_results.append(result)

    # Unify each enzyme into clean structure
    enzymes = [_unify_enzyme(r) for r in raw_results]

    # Compute group-level stats across all enzymes
    group_stats = _compute_group_stats(enzymes)

    return {
        "status":   query_output["status"],
        "comments": query_output["comments"],
        "enzymes":  enzymes,
        "group_stats": {
            "n_enzymes":   len(enzymes),
            "temperature": group_stats["temperature"],
            "ph":          group_stats["ph"],
            "km_mM":       group_stats["km_mM"],
            "kcat_per_s":  group_stats["kcat_per_s"],
        }
    }