import requests
import time

sequence = "MENFQKVEKIGEGTYGVVYKARNKLTGEVVALKKIRLDTETEGVPSTAIREIS" # Replace this with homepage input

# 1. Submit to EBI BLAST API
job = requests.post(
    "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run",
    data={
        "email": "your@email.com",
        "sequence": sequence,
        "database": "uniprotkb_swissprot",  # or "uniprotkb"
        "stype": "protein",
        "program": "blastp"
    }
)

if job.status_code == 200:
    job_id = job.text
    print(f"BLAST job submitted: {job_id}")
else:
    raise Exception("Failed to submit BLAST job")

# 2. Poll for completion
status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
while True:
    status = requests.get(status_url).text
    if status == "FINISHED":
        print("BLAST job finished.")
        break
    elif status in ["RUNNING", "PENDING"]:
        print("Waiting for job to finish...")
        time.sleep(5)
    else:
        raise Exception(f"BLAST job failed with status: {status}")

# 3. Get results (e.g. UniProt accession)
result_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/tsv"
tsv_result = requests.get(result_url).text

# Each row corresponds to a hit; columns include accession
lines = tsv_result.strip().split('\n')
header = lines[0].split('\t')
accession_index = header.index("Accession")

top_hit = lines[1].split('\t')[accession_index]
print(f"Top matching UniProt accession: {top_hit}")

# 4. Query AlphaFold
af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{top_hit}"
af_response = requests.get(af_url)

if af_response.status_code == 200:
    pdb_url = af_response.json()[0]["pdbUrl"]
    pdb_response = requests.get(pdb_url)
    if pdb_response.status_code == 200:
        with open(f"pdb_outputs/{top_hit}.pdb", "w") as f:
            f.write(pdb_response.text)
    print(f"✅ AlphaFold structure downloaded for {top_hit}")
else:
    print(f"❌ No AlphaFold structure found for {top_hit}")
