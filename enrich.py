import numpy as np
from uniprot import _fetch_uniprot, _parse_uniprot, _fetch_sabiork, _fetch_brenda_json
import re

def _extract_float(value) -> float | None:
   
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return float(value)
    # Extract first number (including decimals) from string
    match = re.search(r"\d+\.?\d*", str(value))
    return float(match.group()) if match else None


def _unify_enzyme(raw_result: dict) -> dict:

    def resolve(value):
        """
        Reduce any value to a scalar:
        - None          → None
        - float/int/str → as-is (let _extract_float handle strings)
        - list of dicts → first dict's 'value' key
        - list of scalars → first element
        """
        if value is None:
            return None
        if isinstance(value, (int, float, str)):
            return value
        if isinstance(value, list):
            if not value:
                return None
            first = value[0]
            if isinstance(first, dict):
                return first.get("value")
            return first
        return value

    def first_not_none(*values):
        """Return first value that is not None, regardless of truthiness."""
        for v in values:
            if v is not None:
                return v
        return None

    temp = resolve(first_not_none(
        raw_result.get("brenda_temperature_median"),
        raw_result.get("temp_optimum"),
        raw_result.get("optimal_temp"),
    ))

    ph = resolve(first_not_none(
        raw_result.get("brenda_ph_median"),
        raw_result.get("ph_optimum"),
        raw_result.get("optimal_ph"),
    ))

    km = resolve(first_not_none(
        raw_result.get("brenda_km_q1_mM"),
        raw_result.get("km"),
    ))

    kcat = resolve(first_not_none(
        raw_result.get("brenda_kcat_q1_per_s"),
        raw_result.get("kcat"),
    ))

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
        "cofactors":       raw_result.get("cofactors", []),
    }


def _compute_group_stats(enzymes: list[dict]) -> dict:
    def collect(field):
        return [e[field] for e in enzymes if e.get(field) is not None]

    temps  = collect("temperature")
    phs    = collect("ph")
    kms    = collect("km_mM")
    kcats  = collect("kcat_per_s")

    # Collect all cofactors across all enzymes with counts
    from collections import Counter
    cofactor_counts = Counter(
        cofactor
        for e in enzymes
        for cofactor in e.get("cofactors", [])
    )

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
        "cofactors":   [
            {"name": name, "count": count}
            for name, count in cofactor_counts.most_common()
        ],
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

        # SABIO-RK — only fill fields that are still None (never overwrite UniProt data)
        sabio = _fetch_sabiork(enzyme["uniprot"])
        for key, val in sabio.items():
            if val is not None and not result.get(key):
                result[key] = val

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
            "cofactors":   group_stats["cofactors"],
        }
    }