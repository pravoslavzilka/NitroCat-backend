import hashlib, time
from zeep import Client, Settings

# --- Credentials ---
EMAIL = "mafia@nitroduck.tech"
PASSWORD = hashlib.sha256("a3I5ps$IX29!".encode("utf-8")).hexdigest()

# --- Connect ---
WSDL = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
settings = Settings(strict=False)
client = Client(WSDL, settings=settings)


EC       = "1.1.1.1"
ORGANISM = "Homo sapiens"

EC       = "1.1.1.1"
ORGANISM = "Homo sapiens"

client = Client(
    "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl",
    settings=Settings(strict=False)
)

# ── KM Values ───────────────────────────────────────────────────────────────
time.sleep(1)
km_values = client.service.getKmValue(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    kmValue="*", kmValueMaximum="*", substrate="*",
    commentary="*", ligandStructureId="*", literature="*"
)
print(f"\n=== KM VALUES ({len(km_values)} entries) ===")
for e in km_values:
    print(dict(e))

# ── Turnover Number (kcat) ───────────────────────────────────────────────────
time.sleep(1)
kcat = client.service.getTurnoverNumber(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    turnoverNumber="*", turnoverNumberMaximum="*", substrate="*",
    commentary="*", ligandStructureId="*", literature="*"
)
print(f"\n=== KCAT VALUES ({len(kcat)} entries) ===")
for e in kcat:
    print(dict(e))

# ── pH Optimum ───────────────────────────────────────────────────────────────
time.sleep(1)
ph_opt = client.service.getPhOptimum(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    phOptimum="*", phOptimumMaximum="*",
    commentary="*", literature="*"
)
print(f"\n=== PH OPTIMUM ({len(ph_opt)} entries) ===")
for e in ph_opt:
    print(dict(e))

# ── pH Range ─────────────────────────────────────────────────────────────────
time.sleep(1)
ph_range = client.service.getPhRange(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    phRange="*", phRangeMaximum="*",
    commentary="*", literature="*"
)
print(f"\n=== PH RANGE ({len(ph_range)} entries) ===")
for e in ph_range:
    print(dict(e))

# ── Temperature Optimum ──────────────────────────────────────────────────────
time.sleep(1)
temp_opt = client.service.getTemperatureOptimum(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    temperatureOptimum="*", temperatureOptimumMaximum="*",
    commentary="*", literature="*"
)
print(f"\n=== TEMPERATURE OPTIMUM ({len(temp_opt)} entries) ===")
for e in temp_opt:
    print(dict(e))

# ── Temperature Range ────────────────────────────────────────────────────────
time.sleep(1)
temp_range = client.service.getTemperatureRange(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    temperatureRange="*", temperatureRangeMaximum="*",
    commentary="*", literature="*"
)
print(f"\n=== TEMPERATURE RANGE ({len(temp_range)} entries) ===")
for e in temp_range:
    print(dict(e))

# ── Inhibitors ───────────────────────────────────────────────────────────────
time.sleep(1)
inhibitors = client.service.getInhibitors(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    inhibitor="*", commentary="*",
    ligandStructureId="*", literature="*"
)
print(f"\n=== INHIBITORS ({len(inhibitors)} entries) ===")
for e in inhibitors:
    print(dict(e))

# ── Cofactors ────────────────────────────────────────────────────────────────
time.sleep(1)
cofactors = client.service.getCofactor(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    cofactor="*", commentary="*",
    ligandStructureId="*", literature="*"
)
print(f"\n=== COFACTORS ({len(cofactors)} entries) ===")
for e in cofactors:
    print(dict(e))

# ── Sequences ────────────────────────────────────────────────────────────────
time.sleep(1)
sequences = client.service.getSequence(
    email=EMAIL, password=PASSWORD,
    ecNumber=EC, organism=ORGANISM,
    firstAccessionCode="*", sequence="*",
    noOfAminoAcids="*", source="*", id=""
)
print(f"\n=== SEQUENCES ({len(sequences)} entries) ===")
for e in sequences:
    print(dict(e))