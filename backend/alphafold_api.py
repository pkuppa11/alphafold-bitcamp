from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import requests, time, os
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # or specify your frontend URL like ["http://localhost:5173"]
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class SequenceInput(BaseModel):
    sequence: str

@app.post("/predict")
def predict_structure(data: SequenceInput):
    sequence = data.sequence

    # BLAST submission
    job = requests.post(
        "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run",
        data={
            "email": "your@email.com",
            "sequence": sequence,
            "database": "uniprotkb_swissprot",
            "stype": "protein",
            "program": "blastp"
        }
    )

    if job.status_code != 200:
        raise HTTPException(status_code=500, detail="BLAST job submission failed")

    job_id = job.text

    # Poll for job completion
    status_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/status/{job_id}"
    while True:
        status = requests.get(status_url).text
        if status == "FINISHED":
            break
        elif status in ["RUNNING", "PENDING"]:
            time.sleep(5)
        else:
            raise HTTPException(status_code=500, detail=f"BLAST job failed: {status}")

    # Get BLAST results
    result_url = f"https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/result/{job_id}/tsv"
    tsv_result = requests.get(result_url).text
    lines = tsv_result.strip().split('\n')
    header = lines[0].split('\t')
    accession_index = header.index("Accession")
    top_hit = lines[1].split('\t')[accession_index]

    # Query AlphaFold
    af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{top_hit}"
    af_response = requests.get(af_url)

    if af_response.status_code == 200:
        pdb_url = af_response.json()[0]["pdbUrl"]
        pdb_response = requests.get(pdb_url)
        if pdb_response.status_code == 200:
            os.makedirs("pdb_outputs", exist_ok=True)
            pdb_path = f"pdb_outputs/{top_hit}.pdb"
            with open(pdb_path, "w") as f:
                f.write(pdb_response.text)
            return {"message": "Structure downloaded", "accession": top_hit, "pdb_file": pdb_path}
        else:
            raise HTTPException(status_code=500, detail="PDB download failed")
    else:
        raise HTTPException(status_code=404, detail="No AlphaFold structure found")
    
    return 1
