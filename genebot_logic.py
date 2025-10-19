# genebot_logic.py
import ssl
import certifi
from Bio import Entrez
from google import genai

# -----------------------------
# CONFIGURATION
# -----------------------------
client = genai.Client()  # GEMINI_API_KEY must be set in your environment
Entrez.email = "ben.bioren@gmail.com"
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

# -----------------------------
# PRELOADED ABSTRACTS
# -----------------------------
# genebot_logic.py
import ssl
import certifi
from Bio import Entrez
from google import genai

# -----------------------------
# CONFIGURATION
# -----------------------------
client = genai.Client()  # GEMINI_API_KEY must be set in your environment
Entrez.email = "ben.bioren@gmail.com"
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

# -----------------------------
# PRELOADED ABSTRACTS
# -----------------------------
PRELOADED_ABSTRACTS = {
    "Wood flexibility": "...",
    "Drooping branches": "...",
    "Weeping habit": "...",
    "Annual fruiting": "...",
    "Fruit size": "...",
    "Sweetness": "...",
    "Disease resistance": "..."
}

# -----------------------------
# GEMINI UTILITIES
# -----------------------------
def extract_traits(user_query):
    prompt = f"""
    You are a plant genetics assistant. Extract 5-10 apple tree traits from this request:
    "{user_query}"
    Output only short keywords (1-3 words each) separated by commas. Do NOT include explanations.
    """
    response = client.models.generate_content(model="gemini-2.5-flash", contents=prompt)
    if not response.text:
        return []
    return [kw.strip() for kw in response.text.split(",") if kw.strip()]

def summarize_genes(abstracts, traits):
    if abstracts.strip() in ["", "No abstracts found."]:
        prompt = f"Suggest genes and CRISPR edits for traits: {traits}"
    else:
        prompt = f"Given abstracts:\n{abstracts}\nSuggest genes and CRISPR edits for traits: {traits}"
    response = client.models.generate_content(model="gemini-2.5-flash", contents=prompt)
    return response.text if response.text else "No response from Gemini."

# -----------------------------
# PUBMED SEARCH
# -----------------------------
def search_pubmed(keywords, max_results=5):
    collected_abstracts = []
    for kw in keywords:
        query = f"{kw} AND apple malus domestica"
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            ids = Entrez.read(handle)["IdList"]
            if ids:
                handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
                collected_abstracts.append(handle.read())
            elif kw in PRELOADED_ABSTRACTS:
                collected_abstracts.append(PRELOADED_ABSTRACTS[kw])
        except Exception:
            if kw in PRELOADED_ABSTRACTS:
                collected_abstracts.append(PRELOADED_ABSTRACTS[kw])
    return "\n---\n".join(collected_abstracts) if collected_abstracts else "No abstracts found."
