import requests

def fetch_genome_id(symbol, species="malus_domestica"):
    """
    Given a gene symbol (e.g., MdBPM2), return genome ID and CDS sequence.
    """
    server = "https://plants.ensembl.org"
    ext = f"/lookup/symbol/{species}/{symbol}?expand=1"

    headers = {"Content-Type": "application/json"}

    response = requests.get(server + ext, headers=headers)
    if not response.ok:
        print(f"Error fetching {symbol}: {response.status_code}")
        return None, None

    data = response.json()
    genome_id = data.get("id")
    # Fetch CDS sequence
    seq_ext = f"/sequence/id/{genome_id}?type=cds"
    seq_response = requests.get(server + seq_ext, headers=headers)
    if not seq_response.ok:
        print(f"Error fetching CDS for {genome_id}")
        return genome_id, None

    cds_seq = seq_response.text.strip()
    return genome_id, cds_seq

# Example usage
symbol = "MdBPM2"
genome_id, cds = fetch_genome_id(symbol)
print(f"Symbol: {symbol}, Genome ID: {genome_id}")
print(f"CDS length: {len(cds) if cds else 'not found'}")
