from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from genebot_logic import extract_traits, search_pubmed, summarize_genes
from fastapi.staticfiles import StaticFiles

app = FastAPI()

# Allow CORS from the frontend dev server during development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173", "http://localhost:3000", "*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Serve the built frontend (if present) so a single host can serve API + app
# Place the built Vite files in `frontend/dist` (run `npm run build` in frontend)
try:
    app.mount("/", StaticFiles(directory="frontend/dist", html=True), name="static")
except Exception:
    # If the dist folder doesn't exist (dev environment), ignore and let Vite serve frontend
    pass

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
