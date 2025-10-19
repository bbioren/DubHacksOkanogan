# genebot_streamlit.py
import ssl
import certifi
from Bio import Entrez
from google import genai
import streamlit as st

# -----------------------------
# CONFIGURATION
# -----------------------------
client = genai.Client()  # GEMINI_API_KEY must be set in your environment
Entrez.email = "ben.bioren@gmail.com"

# SSL fix for PubMed
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

# -----------------------------
# PRELOADED APPLE GENOME SNIPPETS
# -----------------------------
PRELOADED_ABSTRACTS = {
    "Wood flexibility": """
Genes involved in secondary cell wall biosynthesis, lignin and cellulose metabolism, and regulatory transcription factors (MYB, NAC) influence branch strength and flexibility in apple trees. Modifications in lignin biosynthesis can alter rigidity and drooping behavior.
""",
    "Drooping branches": """
Weeping habit in apple trees is linked to auxin transport and cell wall-modifying enzymes. Genes like PIN-FORMED (PIN) auxin transporters and expansins contribute to branch angle and drooping tendencies.
""",
    "Weeping habit": """
Regulation of branch orientation in apple trees is influenced by gravitropism-related genes and cytoskeletal regulators. Alterations in these genes can produce drooping or weeping forms.
"""
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
    response = client.models.generate_content(
        model="gemini-2.5-flash",
        contents=prompt
    )
    if not response.text:
        return []
    return [kw.strip() for kw in response.text.split(",") if kw.strip()]

def summarize_genes(abstracts, traits):
    prompt = f"""
    You are a plant bioinformatics assistant. Given the following abstracts, 
    list genes or gene families relevant to these traits: {traits}.
    For each gene, summarize its known function and suggest a potential CRISPR edit 
    (knockout, upregulation, promoter edit, etc.).

    Abstracts:
    {abstracts}
    """
    response = client.models.generate_content(
        model="gemini-2.5-flash",
        contents=prompt
    )
    return response.text if response.text else "No response from Gemini."

# -----------------------------
# PUBMED SEARCH
# -----------------------------
def search_pubmed(keywords, max_results=5, truncate_chars=500):
    collected_abstracts = []

    for kw in keywords:
        query = f"{kw} AND apple"
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        ids = Entrez.read(handle)["IdList"]
        if ids:
            handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
            text = handle.read()
            collected_abstracts.append(text[:truncate_chars])
        elif kw in PRELOADED_ABSTRACTS:
            collected_abstracts.append(PRELOADED_ABSTRACTS[kw][:truncate_chars])

    if not collected_abstracts:
        return "No abstracts found."
    return "\n---\n".join(collected_abstracts)

# -----------------------------
# STREAMLIT UI
# -----------------------------
st.set_page_config(page_title="Apple GeneBot", layout="wide")
st.title("🍏 Apple GeneBot")
st.write(
    """
    Enter desired apple tree traits in natural language, and GeneBot will:
    1. Extract traits.
    2. Search PubMed and preloaded literature.
    3. Suggest candidate genes and potential CRISPR edits.
    """
)

user_input = st.text_area("Enter apple tree traits:", placeholder="e.g., short trees that fruit on new wood")

if st.button("Run GeneBot"):
    if not user_input.strip():
        st.warning("Please enter some traits!")
    else:
        with st.spinner("Extracting traits..."):
            traits = extract_traits(user_input)
        if not traits:
            st.error("No traits found. Check your input or API key.")
        else:
            st.success(f"Traits extracted: {', '.join(traits)}")

            with st.spinner("Searching PubMed and preloaded literature..."):
                abstracts = search_pubmed(traits)
            st.text_area("Abstracts (truncated)", abstracts, height=300)

            with st.spinner("Summarizing genes with Gemini..."):
                summary = summarize_genes(abstracts, traits)
            st.text_area("Suggested Genes & CRISPR Edits", summary, height=400)
