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
# CUSTOM CSS
# -----------------------------
st.markdown("""
<style>
    /* Main container styling */
    .main {
        background: linear-gradient(135deg, #f5f7fa 0%, #e8f5e9 100%);
    }
    
    /* Header styling */
    .header-container {
        background: linear-gradient(135deg, #2e7d32 0%, #66bb6a 100%);
        padding: 2rem;
        border-radius: 15px;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        margin-bottom: 2rem;
        text-align: center;
    }
    
    .header-title {
        color: white;
        font-size: 3rem;
        font-weight: bold;
        margin: 0;
        text-shadow: 2px 2px 4px rgba(0, 0, 0, 0.2);
    }
    
    .header-subtitle {
        color: #e8f5e9;
        font-size: 1.2rem;
        margin-top: 0.5rem;
    }
    
    /* Info box styling */
    .info-box {
        background: white;
        padding: 1.5rem;
        border-radius: 10px;
        border-left: 5px solid #66bb6a;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
        margin-bottom: 2rem;
    }
    
    /* Result card styling */
    .result-card {
        background: white;
        padding: 1.5rem;
        border-radius: 10px;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
        margin-bottom: 1.5rem;
    }
    
    .result-header {
        color: #2e7d32;
        font-size: 1.5rem;
        font-weight: bold;
        margin-bottom: 1rem;
        padding-bottom: 0.5rem;
        border-bottom: 2px solid #e8f5e9;
    }
    
    /* Trait badge styling */
    .trait-badge {
        display: inline-block;
        background: linear-gradient(135deg, #66bb6a 0%, #81c784 100%);
        color: white;
        padding: 0.5rem 1rem;
        border-radius: 20px;
        margin: 0.25rem;
        font-weight: 500;
        box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    }
    
    /* Button styling */
    .stButton>button {
        background: linear-gradient(135deg, #2e7d32 0%, #66bb6a 100%);
        color: white;
        font-size: 1.1rem;
        font-weight: bold;
        padding: 0.75rem 2rem;
        border-radius: 25px;
        border: none;
        box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        transition: all 0.3s ease;
    }
    
    .stButton>button:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(0, 0, 0, 0.15);
    }
    
    /* Text area styling */
    .stTextArea textarea {
        border-radius: 10px;
        border: 2px solid #e0e0e0;
        font-size: 1rem;
    }
    
    .stTextArea textarea:focus {
        border-color: #66bb6a;
        box-shadow: 0 0 0 2px rgba(102, 187, 106, 0.2);
    }
    
    /* Example section */
    .example-section {
        background: #f1f8e9;
        padding: 1rem;
        border-radius: 8px;
        margin-top: 1rem;
        font-size: 0.9rem;
    }
    
    .example-title {
        color: #2e7d32;
        font-weight: bold;
        margin-bottom: 0.5rem;
    }
</style>
""", unsafe_allow_html=True)

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
        query = f"{kw} AND apple"
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
# STREAMLIT UI
# -----------------------------
st.set_page_config(
    page_title="Apple GeneBot",
    page_icon="🍏",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# Header
st.markdown("""
<div class="header-container">
    <h1 class="header-title">🍏 Apple GeneBot</h1>
    <p class="header-subtitle">AI-Powered Gene Discovery for Apple Tree Breeding</p>
</div>
""", unsafe_allow_html=True)

# Info section
st.markdown("""
<div class="info-box">
    <h3 style="color: #2e7d32; margin-top: 0;">How It Works</h3>
    <p style="margin-bottom: 0;">
    Describe your ideal apple tree in natural language, and GeneBot will:
    </p>
    <ol style="margin-top: 0.5rem;">
        <li><strong>Extract key traits</strong> from your description</li>
        <li><strong>Search scientific literature</strong> (PubMed + curated database)</li>
        <li><strong>Identify candidate genes</strong> and suggest CRISPR edits</li>
    </ol>
</div>
""", unsafe_allow_html=True)

# Main input area
col1, col2 = st.columns([2, 1])

with col1:
    user_input = st.text_area(
        "Describe your ideal apple tree:",
        placeholder="e.g., I want a compact tree with large, sweet fruit that resists disease and produces fruit every year",
        height=150,
        key="main_input"
    )
    
with col2:
    st.markdown("""
    <div class="example-section">
        <div class="example-title">💡 Example Traits:</div>
        • Short/dwarf trees<br>
        • Large fruit size<br>
        • High sweetness<br>
        • Disease resistance<br>
        • Annual fruiting<br>
        • Flexible branches<br>
        • Weeping habit
    </div>
    """, unsafe_allow_html=True)

# Run button (centered)
col1, col2, col3 = st.columns([1, 1, 1])
with col2:
    run_button = st.button("🧬 Analyze Traits", use_container_width=True)

# Results section
if run_button:
    if not user_input.strip():
        st.warning("⚠️ Please enter some traits to analyze!")
    else:
        # Trait extraction
        with st.spinner("🔍 Extracting traits..."):
            traits = extract_traits(user_input)
        
        if not traits:
            st.error("❌ No traits found. Please check your input or API configuration.")
        else:
            # Display extracted traits
            st.markdown('<div class="result-card">', unsafe_allow_html=True)
            st.markdown('<div class="result-header">📋 Extracted Traits</div>', unsafe_allow_html=True)
            traits_html = "".join([f'<span class="trait-badge">{trait}</span>' for trait in traits])
            st.markdown(traits_html, unsafe_allow_html=True)
            st.markdown('</div>', unsafe_allow_html=True)
            
            # Literature search
            with st.spinner("📚 Searching scientific literature..."):
                abstracts = search_pubmed(traits)
            
            # Display abstracts in expandable section
            with st.expander("📄 View Literature Abstracts", expanded=False):
                st.text_area("Full Abstracts", abstracts, height=300, key="abstracts")
            
            # Gene analysis
            with st.spinner("🧬 Analyzing genes with AI..."):
                summary = summarize_genes(abstracts, traits)
            
            # Display gene recommendations
            st.markdown('<div class="result-card">', unsafe_allow_html=True)
            st.markdown('<div class="result-header">🎯 Gene Recommendations & CRISPR Strategies</div>', unsafe_allow_html=True)
            st.markdown(summary)
            st.markdown('</div>', unsafe_allow_html=True)
            
            # Success message
            st.success("✅ Analysis complete! Review the gene recommendations above.")

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #666; padding: 1rem;">
    <p style="margin: 0;">Built with Streamlit, Gemini AI, and PubMed</p>
    <p style="margin: 0; font-size: 0.9rem;">For research and educational purposes</p>
</div>
""", unsafe_allow_html=True)