# genebot_streamlit_enhanced.py
import ssl
import certifi
from Bio import Entrez
from google import genai
import streamlit as st

# -----------------------------
# CONFIGURATION
# -----------------------------
client = genai.Client()
Entrez.email = "ben.bioren@gmail.com"
ssl._create_default_https_context = lambda: ssl.create_default_context(cafile=certifi.where())

# -----------------------------
# PRELOADED APPLE GENOME SNIPPETS
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
    if abstracts.strip() in ["", "No abstracts found."]:
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
            if kw in PRELOADED_ABSTRACTS:
                collected_abstracts.append(PRELOADED_ABSTRACTS[kw])

    if not collected_abstracts:
        return "No abstracts found."
    return "\n---\n".join(collected_abstracts)

# -----------------------------
# CUSTOM CSS
# -----------------------------
st.set_page_config(page_title="CRISPR Apples", page_icon="🍏", layout="wide")

st.markdown("""
    <style>
    /* Main container styling */
    .main {
        background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
    }
    
    /* Header styling */
    .header-container {
        background: linear-gradient(120deg, #84fab0 0%, #8fd3f4 100%);
        padding: 2rem;
        border-radius: 15px;
        margin-bottom: 2rem;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
    }
    
    .main-title {
        font-size: 3rem;
        font-weight: 700;
        color: #2d3748;
        text-align: center;
        margin: 0;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.1);
    }
    
    .subtitle {
        font-size: 1.2rem;
        color: #4a5568;
        text-align: center;
        margin-top: 0.5rem;
    }
    
    /* Card styling */
    .info-card {
        background: white;
        padding: 1.5rem;
        border-radius: 12px;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        margin-bottom: 1.5rem;
        border-left: 4px solid #84fab0;
    }
    
    .result-card {
        background: white;
        padding: 1.5rem;
        border-radius: 12px;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        margin-top: 1rem;
    }
    
    /* Trait badges */
    .trait-badge {
        display: inline-block;
        background: linear-gradient(120deg, #84fab0 0%, #8fd3f4 100%);
        color: #2d3748;
        padding: 0.5rem 1rem;
        border-radius: 20px;
        margin: 0.3rem;
        font-weight: 600;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    }
    
    /* Button styling */
    .stButton>button {
        background: linear-gradient(120deg, #84fab0 0%, #8fd3f4 100%);
        color: #2d3748;
        font-weight: 600;
        border: none;
        padding: 0.75rem 2rem;
        border-radius: 25px;
        font-size: 1.1rem;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        transition: transform 0.2s;
    }
    
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(0, 0, 0, 0.15);
    }
    
    /* Section headers */
    .section-header {
        font-size: 1.5rem;
        font-weight: 600;
        color: #2d3748;
        margin-bottom: 1rem;
        padding-bottom: 0.5rem;
        border-bottom: 3px solid #84fab0;
    }
    
    /* Examples styling */
    .example-box {
        background: #f7fafc;
        padding: 1rem;
        border-radius: 8px;
        border-left: 3px solid #8fd3f4;
        margin-top: 1rem;
        font-style: italic;
        color: #4a5568;
    }
    
    /* Step indicators */
    .step-indicator {
        display: inline-block;
        background: #8fd3f4;
        color: white;
        width: 30px;
        height: 30px;
        border-radius: 50%;
        text-align: center;
        line-height: 30px;
        font-weight: bold;
        margin-right: 0.5rem;
    }
    </style>
""", unsafe_allow_html=True)

# -----------------------------
# STREAMLIT UI
# -----------------------------

# Header
st.markdown("""
    <div class="header-container">
        <h1 class="main-title">🍏 CRISPR Apples</h1>
        <p class="subtitle">AI-Powered Genetic Fruit Designer</p>
    </div>
""", unsafe_allow_html=True)

# Introduction card
st.markdown("""
    <div class="info-card">
        <h3 style="color: #2d3748; margin-top: 0;">How It Works</h3>
        <p><span class="step-indicator">1</span> Describe your desired apple tree traits in natural language</p>
        <p><span class="step-indicator">2</span> GeneBot extracts and analyzes key traits using AI</p>
        <p><span class="step-indicator">3</span> Searches scientific literature and genomic databases</p>
        <p><span class="step-indicator">4</span> Suggests candidate genes and CRISPR editing strategies</p>
    </div>
""", unsafe_allow_html=True)

# Main input section
col1, col2 = st.columns([2, 1])

with col1:
    st.markdown('<p class="section-header">🧬 Describe Your Ideal Apple Tree</p>', unsafe_allow_html=True)
    user_input = st.text_area(
        "",
        placeholder="Example: I want compact trees with large, sweet fruit that resist fire blight and produce annually...",
        height=150,
        label_visibility="collapsed"
    )

with col2:
    st.markdown('<p class="section-header">💡 Example Traits</p>', unsafe_allow_html=True)
    st.markdown("""
        <div class="example-box">
        • Dwarf growth habit<br>
        • Disease resistance<br>
        • Large fruit size<br>
        • High sugar content<br>
        • Cold hardiness<br>
        • Early fruiting<br>
        • Weeping branches
        </div>
    """, unsafe_allow_html=True)

# Center the button
col_center = st.columns([1, 1, 1])
with col_center[1]:
    run_button = st.button("🚀 Analyze Traits", use_container_width=True)

# Results section
if run_button:
    if not user_input.strip():
        st.warning("⚠️ Please enter some traits to analyze!")
    else:
        # Trait extraction
        with st.spinner("🔍 Extracting traits from your description..."):
            traits = extract_traits(user_input)
        
        if not traits:
            st.error("❌ No traits found. Please check your input or API configuration.")
        else:
            # Display extracted traits
            st.markdown('<div class="result-card">', unsafe_allow_html=True)
            st.markdown('<p class="section-header">✅ Extracted Traits</p>', unsafe_allow_html=True)
            traits_html = "".join([f'<span class="trait-badge">{trait}</span>' for trait in traits])
            st.markdown(traits_html, unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)

            # Literature search
            with st.spinner("📚 Searching PubMed and genomic databases..."):
                abstracts = search_pubmed(traits)
            
            # Display abstracts in expandable section
            with st.expander("📄 View Scientific Literature (Click to expand)", expanded=False):
                st.text_area("", abstracts, height=300, label_visibility="collapsed")

            # Gene summary
            with st.spinner("🧬 Analyzing genes and generating CRISPR recommendations..."):
                summary = summarize_genes(abstracts, traits)
            
            st.markdown('<div class="result-card">', unsafe_allow_html=True)
            st.markdown('<p class="section-header">🎯 Candidate Genes & CRISPR Editing Strategies</p>', unsafe_allow_html=True)
            st.markdown(summary)
            st.markdown('</div>', unsafe_allow_html=True)

# Footer
st.markdown("""
    <div style="text-align: center; margin-top: 3rem; padding: 2rem; color: #718096;">
        <p>Powered by Gemini AI & NCBI PubMed | For Research Purposes</p>
    </div>
""", unsafe_allow_html=True)