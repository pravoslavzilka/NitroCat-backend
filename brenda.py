import json
import re
import numpy as np

with open("brenda_2026_1.json") as f:
    brenda = json.load(f)


def find_protein_id(ec: str, uniprot_id: str) -> tuple[str, dict]:
    ec_entry = brenda["data"].get(ec)
    if not ec_entry:
        raise ValueError(f"EC {ec} not found in BRENDA")
    for pid, protein in ec_entry.get("protein", {}).items():
        if uniprot_id in protein.get("accessions", []):
            return pid, ec_entry
    raise ValueError(f"UniProt ID {uniprot_id} not found under EC {ec}")


def parse_value_substrate(raw: str) -> tuple[str, str]:
    """Parse '0.77 {ethanol}' → ('0.77', 'ethanol')"""
    match = re.match(r"([\d.eE+\-]+(?:\s*-\s*[\d.eE+\-]+)?)\s*(?:\{(.+)\})?", raw.strip())
    if match:
        return match.group(1), match.group(2) or ""
    return raw.strip(), ""


def get_entries(ec_entry: dict, field: str, protein_id: str) -> list:
    return [
        e for e in ec_entry.get(field, [])
        if protein_id in e.get("proteins", [])
    ]

def get_entries_with_fallback(ec_entry: dict, field: str, protein_id: str, organism: str) -> list:
    """
    First try protein-specific entries.
    If none found, fall back to all entries for this organism.
    """
    # Try protein-level first
    entries = [
        e for e in ec_entry.get(field, [])
        if protein_id in e.get("proteins", [])
    ]
    if entries:
        return entries

    # Fall back to organism-level — find all protein IDs for this organism
    organism_pids = {
        pid for pid, p in ec_entry.get("protein", {}).items()
        if p.get("organism", "").lower() == organism.lower()
    }

    entries = [
        e for e in ec_entry.get(field, [])
        if organism_pids & set(e.get("proteins", []))
    ]
    if entries:
        print(f"  [{field}] no protein-specific data, using organism-level ({len(entries)} entries)")
        return entries

    # Last resort — return everything for this EC
    all_entries = ec_entry.get(field, [])
    if all_entries:
        print(f"  [{field}] no organism-specific data, using EC-level ({len(all_entries)} entries)")
    return all_entries


def summarize(ec: str, uniprot_id: str) -> dict:
    protein_id, ec_entry = find_protein_id(ec, uniprot_id)
    protein_info = ec_entry["protein"][protein_id]
    organism = protein_info["organism"]

    result = {
        "ec":           ec,
        "uniprot":      uniprot_id,
        "protein_id":   protein_id,
        "organism":     organism,
        "enzyme_name":  ec_entry.get("recommended_name", ""),

        "temperature_optimum":   [],
        "temperature_range":     [],
        "temperature_stability": [],
        "ph_optimum":            [],
        "ph_range":              [],
        "ph_stability":          [],
        "km_value":              [],
        "turnover_number":       [],
        "kcat_km_value":         [],
        "ki_value":              [],
        "ic50_value":            [],
        "specific_activity":     [],
    }

    # Temperature and pH — use fallback
    for field, val_key in [
        ("temperature_optimum",   "value_celsius"),
        ("temperature_range",     "range_celsius"),
        ("temperature_stability", "value_celsius"),
        ("ph_optimum",            "value"),
        ("ph_range",              "range"),
        ("ph_stability",          "value"),
    ]:
        for e in get_entries_with_fallback(ec_entry, field, protein_id, organism):
            val, _ = parse_value_substrate(e["value"])
            result[field].append({
                val_key:   val,
                "comment": e.get("comment", ""),
            })

    # Kinetics — keep strict protein-level filter
    for e in get_entries(ec_entry, "km_value", protein_id):
        val, substrate = parse_value_substrate(e["value"])
        result["km_value"].append({
            "value_mM":  val,
            "substrate": substrate,
            "comment":   e.get("comment", ""),
        })

    for e in get_entries(ec_entry, "turnover_number", protein_id):
        val, substrate = parse_value_substrate(e["value"])
        result["turnover_number"].append({
            "value_per_s": val,
            "substrate":   substrate,
            "comment":     e.get("comment", ""),
        })

    for e in get_entries(ec_entry, "kcat_km_value", protein_id):
        val, substrate = parse_value_substrate(e["value"])
        result["kcat_km_value"].append({
            "value_mM_per_s": val,
            "substrate":      substrate,
            "comment":        e.get("comment", ""),
        })

    for e in get_entries(ec_entry, "ki_value", protein_id):
        val, substrate = parse_value_substrate(e["value"])
        result["ki_value"].append({
            "value_mM":  val,
            "inhibitor": substrate,
            "comment":   e.get("comment", ""),
        })

    for e in get_entries(ec_entry, "specific_activity", protein_id):
        val, _ = parse_value_substrate(e["value"])
        result["specific_activity"].append({
            "value_umol_per_min_per_mg": val,
            "comment": e.get("comment", ""),
        })

    return result



def print_summary(ec: str, uniprot_id: str):
    data = summarize(ec, uniprot_id)

    print(f"\n{'='*60}")
    print(f"Enzyme:    {data['enzyme_name']} (EC {data['ec']})")
    print(f"UniProt:   {data['uniprot']}")
    print(f"Organism:  {data['organism']}")
    print(f"Optimal temp:  {data['temperature_optimum']}")
    print(f"ph Optimal:  {data['ph_optimum']}")
    print(f"{'='*60}")

    sections = [
        ("TEMPERATURE OPTIMUM (°C)",    "temperature_optimum",  "value_celsius"),
        ("TEMPERATURE RANGE (°C)",      "temperature_range",    "range_celsius"),
        ("TEMPERATURE STABILITY (°C)",  "temperature_stability","value_celsius"),
        ("PH OPTIMUM",                  "ph_optimum",           "value"),
        ("PH RANGE",                    "ph_range",             "range"),
        ("PH STABILITY",                "ph_stability",         "value"),
        ("KM VALUE (mM)",               "km_value",             "value_mM"),
        ("TURNOVER NUMBER (1/s)",       "turnover_number",      "value_per_s"),
        ("KCAT/KM (mM/s)",              "kcat_km_value",        "value_mM_per_s"),
        ("KI VALUE (mM)",               "ki_value",             "value_mM"),
        ("IC50 VALUE (mM)",             "ic50_value",           "value_mM"),
        ("SPECIFIC ACTIVITY (µmol/min/mg)", "specific_activity","value_umol_per_min_per_mg"),
    ]

    for title, field, val_key in sections:
        entries = data[field]
        if not entries:
            continue
        print(f"\n--- {title} ({len(entries)} entries) ---")
        for e in entries:
            val     = e[val_key]
            sub     = e.get("substrate") or e.get("inhibitor") or ""
            comment = e.get("comment", "")
            sub_str     = f"  [{sub}]"     if sub     else ""
            comment_str = f"  — {comment}" if comment else ""
            print(f"  {val}{sub_str}{comment_str}")





def get_km_q1_summary(ec: str, uniprot_id: str) -> dict:
    data = summarize(ec, uniprot_id)
    km_entries = data["km_value"]

    if not km_entries:
        print("No KM values found.")
        return {}

    # Parse values to float, skip unparseable
    parsed = []
    for e in km_entries:
        try:
            parsed.append((float(e["value_mM"]), e["substrate"], e["comment"]))
        except ValueError:
            continue

    values = np.array([p[0] for p in parsed])

    # Stats
    q1      = np.percentile(values, 25)
    median  = np.percentile(values, 50)
    mean    = np.mean(values)
    std     = np.std(values)

    # Keep only entries in the first quartile (below Q1)
    q1_entries = [(v, s, c) for v, s, c in parsed if v <= q1]

    print(f"\n--- KM VALUE SUMMARY for {ec} / {uniprot_id} ---")
    print(f"  Total entries : {len(parsed)}")
    print(f"  Mean          : {mean:.4f} mM")
    print(f"  Std dev       : {std:.4f} mM")
    print(f"  Median (Q2)   : {median:.4f} mM")
    print(f"  Q1 threshold  : {q1:.4f} mM")
    print(f"  Entries in Q1 : {len(q1_entries)}")
    print(f"\n  --- First Quartile entries (≤ {q1:.4f} mM) ---")
    for val, substrate, comment in sorted(q1_entries):
        sub_str     = f"  [{substrate}]" if substrate else ""
        comment_str = f"  — {comment}"   if comment  else ""
        print(f"  {val:.5f} mM{sub_str}{comment_str}")

    return {
        "mean":       mean,
        "std":        std,
        "median":     median,
        "q1":         q1,
        "q1_entries": q1_entries,
    }


def get_optimal_conditions(ec: str, uniprot_id: str) -> dict:
    data = summarize(ec, uniprot_id)

    def extract_numeric(entries: list, val_key: str) -> list[float]:
        values = []
        for e in entries:
            raw = e[val_key]
            # handle ranges like "6.5 - 7.5" → take midpoint
            if "-" in str(raw):
                parts = raw.split("-")
                try:
                    values.append((float(parts[0]) + float(parts[1])) / 2)
                except ValueError:
                    continue
            else:
                try:
                    values.append(float(raw))
                except ValueError:
                    continue
        return values

    def summarize_values(values: list[float], label: str, unit: str) -> dict:
        if not values:
            print(f"  No {label} data found.")
            return {}
        arr = np.array(values)
        result = {
            "mean":   float(np.mean(arr)),
            "median": float(np.median(arr)),
            "std":    float(np.std(arr)),
            "min":    float(np.min(arr)),
            "max":    float(np.max(arr)),
            "q1":     float(np.percentile(arr, 25)),
            "q3":     float(np.percentile(arr, 75)),
            "n":      len(values),
        }
        print(f"\n  --- {label} ({unit}) ---")
        print(f"  N entries  : {result['n']}")
        print(f"  Mean       : {result['mean']:.2f}")
        print(f"  Median     : {result['median']:.2f}")
        print(f"  Std dev    : {result['std']:.2f}")
        print(f"  Range      : {result['min']:.2f} — {result['max']:.2f}")
        print(f"  IQR        : {result['q1']:.2f} — {result['q3']:.2f}")
        return result

    print(f"\n{'='*60}")
    print(f"OPTIMAL CONDITIONS: EC {ec} / UniProt {uniprot_id}")
    print(f"Enzyme   : {data['enzyme_name']}")
    print(f"Organism : {data['organism']}")
    print(f"{'='*60}")

    # --- Temperature ---
    temp_opt    = extract_numeric(data["temperature_optimum"],   "value_celsius")
    temp_range  = extract_numeric(data["temperature_range"],     "range_celsius")
    temp_stable = extract_numeric(data["temperature_stability"], "value_celsius")

    temp_all = temp_opt + temp_range + temp_stable

    temp_stats = summarize_values(temp_opt,    "Temperature Optimum",   "°C")
    summarize_values(temp_range,               "Temperature Range",     "°C")
    summarize_values(temp_stable,              "Temperature Stability", "°C")
    combined_temp = summarize_values(temp_all, "Temperature Combined",  "°C")

    # --- pH ---
    ph_opt    = extract_numeric(data["ph_optimum"],   "value")
    ph_range  = extract_numeric(data["ph_range"],     "range")
    ph_stable = extract_numeric(data["ph_stability"], "value")

    ph_all = ph_opt + ph_range + ph_stable

    ph_stats = summarize_values(ph_opt,        "pH Optimum",   "")
    summarize_values(ph_range,                 "pH Range",     "")
    summarize_values(ph_stable,                "pH Stability", "")
    combined_ph = summarize_values(ph_all,     "pH Combined",  "")

    # --- Best estimate ---
    print(f"\n{'='*60}")
    print(f"BEST ESTIMATE (median of all available data)")
    print(f"{'='*60}")

    if temp_all:
        print(f"  Temperature : {np.median(temp_all):.1f} °C  "
              f"(range {min(temp_all):.1f} — {max(temp_all):.1f} °C)")
    else:
        print("  Temperature : no data")

    if ph_all:
        print(f"  pH          : {np.median(ph_all):.2f}  "
              f"(range {min(ph_all):.2f} — {max(ph_all):.2f})")
    else:
        print("  pH          : no data")

    return {
        "temperature": {
            "optimum":   temp_stats,
            "combined":  combined_temp,
            "all_values": temp_all,
        },
        "ph": {
            "optimum":   ph_stats,
            "combined":  combined_ph,
            "all_values": ph_all,
        },
        "best_estimate": {
            "temperature_median": float(np.median(temp_all)) if temp_all else None,
            "ph_median":          float(np.median(ph_all))   if ph_all   else None,
        }
    }


def get_enzyme_properties(ec: str, uniprot_id: str) -> dict:
    data = summarize(ec, uniprot_id)

    def normalize_to_value(entries: list, possible_keys: list[str]) -> list[dict]:
        """Rename whichever key exists to 'value' for uniform processing."""
        result = []
        for e in entries:
            for key in possible_keys:
                if key in e:
                    result.append({"value": e[key], "comment": e.get("comment", "")})
                    break
        return result

    def q1_values(entries: list) -> float | None:
        values = []
        for e in entries:
            try:
                values.append(float(e["value"]))
            except (ValueError, TypeError):
                continue
        if not values:
            return None
        return float(np.percentile(values, 25))

    def median_values(entries: list) -> float | None:
        values = []
        for e in entries:
            raw = str(e["value"])
            if "-" in raw:
                parts = raw.split("-")
                try:
                    values.append((float(parts[0]) + float(parts[1])) / 2)
                except ValueError:
                    continue
            else:
                try:
                    values.append(float(raw))
                except ValueError:
                    continue
        if not values:
            return None
        return float(np.median(values))

    # Normalize all temp/pH entries to use "value" key
    all_temp = normalize_to_value(data["temperature_optimum"],   ["value_celsius"]) + \
               normalize_to_value(data["temperature_range"],     ["range_celsius"]) + \
               normalize_to_value(data["temperature_stability"], ["value_celsius"])

    all_ph   = normalize_to_value(data["ph_optimum"],   ["value"]) + \
               normalize_to_value(data["ph_range"],     ["range"]) + \
               normalize_to_value(data["ph_stability"], ["value"])

    all_km   = normalize_to_value(data["km_value"],        ["value_mM"])
    all_kcat = normalize_to_value(data["turnover_number"], ["value_per_s"])

    result = {
        "ec":                 ec,
        "uniprot":            uniprot_id,
        "organism":           data["organism"],
        "enzyme_name":        data["enzyme_name"],
        "km_q1_mM":           q1_values(all_km),
        "kcat_q1_per_s":      q1_values(all_kcat),
        "temperature_median": median_values(all_temp),
        "ph_median":          median_values(all_ph),
    }

    print(f"\n{'='*50}")
    print(f"Enzyme     : {result['enzyme_name']} (EC {ec})")
    print(f"UniProt    : {uniprot_id}")
    print(f"Organism   : {result['organism']}")
    print(f"{'='*50}")
    print(f"KM   Q1    : {result['km_q1_mM']:.5f} mM"          if result['km_q1_mM']           is not None else "KM   Q1    : no data")
    print(f"Kcat Q1    : {result['kcat_q1_per_s']:.3f} 1/s"    if result['kcat_q1_per_s']       is not None else "Kcat Q1    : no data")
    print(f"Temp median: {result['temperature_median']:.1f} °C" if result['temperature_median']  is not None else "Temp median: no data")
    print(f"pH   median: {result['ph_median']:.2f}"             if result['ph_median']           is not None else "pH   median: no data")

    return result



if __name__ ==  "__main__":
    
    # Usage
    props = get_enzyme_properties("1.1.1.1", "P00326")

    print(props)






