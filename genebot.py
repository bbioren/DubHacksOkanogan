# genebot_streamlit_full.py
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
# PRELOADED APPLE GENOME SNIPPETS (richer)
# -----------------------------
PRELOADED_ABSTRACTS = {
    "Wood flexibility": """
Genes involved in secondary cell wall biosynthesis: Lignin biosynthesis genes (CCR, CAD), Cellulose Synthase genes (CESA), transcription factors MYB46, NAC. Alterations in these genes modulate wood rigidity and flexibility. CRISPR upregulation increases rigidity, knockout/downregulation increases flexibility.
""",
    "Drooping branches": """
Genes affecting branch angle and weeping habit: PIN-FORMED (PIN) auxin transporters, EXPANSINs, TIR1 auxin receptor, LAZY1. Modifying these genes changes branch orientation and drooping.
""",
    "Weeping habit": """
Genes regulating branch gravitropism and orientation: LAZY1, TAC1, cytoskeletal regulators. CRISPR edits can induce or reduce weeping forms.
""",
    "Annual fruiting": """
Genes controlling flowering and fruiting cycles: FT (Flowering Locus T), TFL1, SOC1, AP1, LFY. CRISPR upregulation promotes yearly fruiting; knockouts can delay flowering.
""",
    "Fruit size": """
Genes controlling fruit growth: Cell cycle regulators (CYCD3, CYCB1), Auxin response factors (ARF), Gibberellin pathway genes (GA20ox, GA3ox). Modulating these can increase or reduce fruit size.
""",
    "Sweetness": """
Sugar metabolism genes: Sucrose Synthase (SUS), Invertases (INV), Hexose Transporters (HT). Upregulation increases sugar accumulation; knockouts decrease sweetness.
""",
    "Disease resistance": """
Resistance genes: R-genes (NBS-LRR family), Pathogenesis-Related proteins (PR1, PR5), WRKY transcription factors. CRISPR can enhance or modify resistance.
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
    """
    Summarize candidate genes and suggest CRISPR edits.
    Fallback if abstracts are sparse.
    """
    if abstracts.strip() in ["", "No abstracts found."]:
        # Fallback: ask Gemini directly to suggest genes for the traits
        prompt = f"""
        You are a plant genetics expert. For each apple tree trait below, suggest 2-5 relevant genes or gene families,
        explain their function, and suggest potential CRISPR edits (knockout, upregulation, promoter edit, etc.).

        Traits:
        {traits}
        """
    else:
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
def search_pubmed(keywords, max_results=5):
    collected_abstracts = []

    for kw in keywords:
        query = f"{kw} AND apple malus domestica"
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
            ids = Entrez.read(handle)["IdList"]
            if ids:
                handle = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
                text = handle.read()
                collected_abstracts.append(text)
            elif kw in PRELOADED_ABSTRACTS:
                collected_abstracts.append(PRELOADED_ABSTRACTS[kw])
        except Exception as e:
            # fallback to preloaded abstracts if PubMed fails
            if kw in PRELOADED_ABSTRACTS:
                collected_abstracts.append(PRELOADED_ABSTRACTS[kw])

    if not collected_abstracts:
        return "No abstracts found."
    return "\n---\n".join(collected_abstracts)

# -----------------------------
# STREAMLIT UI
# -----------------------------
st.set_page_config(page_title="Apple GeneBot", layout="wide")
st.title("🍏 Apple GeneBot (Hackathon Version)")

st.write(
    """
    Enter desired apple tree traits in natural language. GeneBot will:
    1. Extract traits.
    2. Search PubMed + preloaded literature.
    3. Suggest candidate genes and potential CRISPR edits.
    If no literature is found, it will generate genes directly from the traits.
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
            st.text_area("Abstracts (full)", abstracts, height=300)

            with st.spinner("Summarizing genes with Gemini..."):
                summary = summarize_genes(abstracts, traits)
            st.text_area("Suggested Genes & CRISPR Edits", summary, height=400)
