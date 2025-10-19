from google import genai

from Bio import Entrez
Entrez.email = "your_email@example.com"


# The client gets the API key from the environment variable `GEMINI_API_KEY`.
client = genai.Client()

response = client.models.generate_content(
    model="gemini-2.5-flash", contents="Explain how AI works in a few words"
)

def extract_traits(query):
    prompt = f"""
    Extract plant trait keywords and related biological pathways from this request:
    "{query}"
    Output as a comma-separated list of keywords for database search.
    """
    response = client.models.generate_content(
        model="gemini-2.5-flash", contents=prompt
    )

    return [kw.strip() for kw in response.text.split(",")]

print(extract_traits("I want apple trees that stay short and flower earlier."))



print(response.text)


def search_pubmed(keywords, max_results=5):
    query = " OR ".join(keywords) + " AND apple"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    ids = Entrez.read(handle)["IdList"]
    if not ids: return []
    handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
    return handle.read()
