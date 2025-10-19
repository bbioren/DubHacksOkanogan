# backend/main.py
from fastapi import FastAPI
from pydantic import BaseModel
from typing import List
from genebot import extract_traits, search_pubmed, summarize_genes

app = FastAPI()

class TraitRequest(BaseModel):
    query: str

@app.post("/run_genebot")
def run_genebot(request: TraitRequest):
    traits = extract_traits(request.query)
    abstracts = search_pubmed(traits)
    summary = summarize_genes(abstracts, traits)
    return {
        "traits": traits,
        "abstracts": abstracts,
        "summary": summary
    }
