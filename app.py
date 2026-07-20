import streamlit as st
import plotly.express as px
import pandas as pd
import random
import re

# --- DESIGN & VISUAL INTERFACE ARCHITECTURE --- 
st.set_page_config( 
    page_title="CRISPR Computational Workstation", 
    layout="wide", 
    initial_sidebar_state="expanded" 
)

# Inject production-grade clinical workspace styling layers
st.markdown(""" 
    <style> 
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@500;700&display=swap'); 
    
    /* Primary Layout Foundations */
    html, body, [data-testid="stAppViewContainer"] { 
        background-color: #0e1626 !important; 
        font-family: 'Inter', sans-serif !important; 
        color: #f1f5f9 !important; 
    } 
    p, span, label { color: #f1f5f9 !important; } 
    
    /* Technical Display Card Containers */
    .dashboard-card { 
        background-color: #152238; 
        border: 1px solid #22314d; 
        border-radius: 8px; 
        padding: 24px; 
        margin-bottom: 20px; 
    } 
    
    /* Clinical Summary Metric Blocks */
    div[data-testid="stMetric"] { 
        background-color: #152238 !important; 
        border-radius: 6px !important; 
        padding: 16px !important; 
        border: 1px solid #22314d !important; 
        border-left: 4px solid #38bdf8 !important; 
    } 
    div[data-testid="stMetricLabel"] div, div[data-testid="stMetricLabel"] p { 
        color: #94a3b8 !important; 
        font-size: 0.8rem !important; 
        text-transform: uppercase !important; 
        letter-spacing: 0.05em !important; 
        font-weight: 600 !important; 
    } 
    div[data-testid="stMetricValue"] div, div[data-testid="stMetricValue"] span { 
        color: #38bdf8 !important; 
        font-family: 'JetBrains Mono', monospace !important; 
        font-weight: 700 !important; 
        font-size: 1.6rem !important; 
    } 
    
    /* Strict Typography Rules */
    h1 {
        font-family: 'Inter', sans-serif !important;
        font-weight: 700 !important;
        color: #ffffff !important;
        letter-spacing: -0.02em;
    }
    h2, h3, h4 { 
        font-family: 'Inter', sans-serif !important; 
        font-weight: 500 !important; 
        color: #f1f5f9 !important; 
    } 
    
    /* Text Input Elements */
    textarea { 
        background-color: #1c2d4a !important; 
        border: 1px solid #2d3f66 !important; 
        color: #ffffff !important; 
        font-family: 'JetBrains Mono', monospace !important; 
        border-radius: 6px !important; 
    } 
    
    /* Refined Primary Operational Button */
    .stButton>button { 
        background: #38bdf8 !important; 
        color: #0e1626 !important; 
        font-weight: 600 !important; 
        letter-spacing: 0.01em; 
        border-radius: 6px !important; 
        border: none !important; 
        padding: 12px 24px !important; 
        transition: background-color 0.15s ease-in-out !important; 
    } 
    .stButton>button:hover { 
        background: #0ea5e9 !important; 
        color: #0e1626 !important; 
    } 
    
    /* Professional Sidebar Framing */
    [data-testid="stSidebar"] { 
        background-color: #090e1a !important; 
        border-right: 1px solid #22314d; 
    } 
    </style> 
""", unsafe_allow_html=True)


# --- COMPUTATIONAL BIOINFORMATICS ENGINE --- 
class CRISPRAnalyzerEngine:
    """
    Object-Oriented Analysis Engine isolating biological data computation 
    from presentation layer components.
    """
    def __init__(self, seed: int = 42, num_loci: int = 50):
        self.target_base = "ATCGATCGATCGATCGATCG"
        self.bases = ['A', 'T', 'C', 'G']
        self.seed = seed
        self.num_loci = num_loci
        self.off_target_db = self._generate_variant_matrix()

    def _generate_variant_matrix(self) -> list:
        random.seed(self.seed)
        matrix_bank = []
        for _ in range(self.num_loci):
            sequence_list = list(self.target_base)
            num_mutations = random.choice([0, 1, 1, 2, 2, 3, 4, 5])
            mutation_positions = random.sample(range(20), min(num_mutations, 20))
            
            for pos in mutation_positions:
                current_base = sequence_list[pos]
                alternatives = [b for b in self.bases if b != current_base]
                sequence_list[pos] = random.choice(alternatives)
                
            matrix_bank.append("".join(sequence_list))
        return matrix_bank

    @staticmethod
    def clean_sequence(input_text: str) -> str:
        cleaned = input_text.upper().replace("\n", "").replace("\r", "").strip()
        return re.sub(r'[^ATCG]', '', cleaned)

    def calculate_safety_score(self, grna: str) -> float:
        total_penalty = 0
        for off_target in self.off_target_db:
            mismatches = 0
            seed_penalty = 0
            
            for pos in range(20):
                if grna[pos] != off_target[pos]:
                    mismatches += 1
                    if pos >= 10:
                        seed_penalty += 2
                    else:
                        seed_penalty += 10
                        
            if mismatches == 0:
                total_penalty += 25
            elif mismatches == 1:
                total_penalty += max(5, 15 - seed_penalty)
            elif mismatches == 2:
                total_penalty += max(2, 8 - seed_penalty)
                
        return round(max(0, 100 - total_penalty), 2)

    def process_strand(self, raw_dna: str) -> list:
        sequence = self.clean_sequence(raw_dna)
        pipeline_results = []
        
        for i in range(23, len(sequence)):
            if sequence[i:i+2] == "GG":
                grna_target = sequence[i-23:i-3]
                if len(grna_target) == 20:
                    gc_content = ((grna_target.count('G') + grna_target.count('C')) / 20) * 100
                    efficiency = max(0, 100 - (abs(50 - gc_content) * 2))
                    safety = self.calculate_safety_score(grna_target)
                    composite = (efficiency * 0.4) + (safety * 0.6)
                    
                    pipeline_results.append({
                        "Position": i - 23,
                        "gRNA Sequence": grna_target,
                        "PAM": "N" + sequence[i:i+2],
                        "Efficiency Score": round(efficiency, 2),
                        "Safety Score": safety,
                        "Global Composite Score": round(composite, 2)
                    })
        return sorted(pipeline_results, key=lambda x: x['Global Composite Score'], reverse=True)


# --- INTERACTION TERMINAL CONTROL LAYER --- 
with st.sidebar: 
    st.markdown("### Data Control Panel") 
    input_method = st.radio("Genomic Data Source Type:", ["Manual Text Entry", "Upload FASTA File"]) 
    st.markdown("---") 
    st.markdown("<small style='color: #64748b;'>System Engine v4.0.0<br>Architecture: Object-Oriented</small>", unsafe_allow_html=True) 

# Instantiate analytical core
analyzer = CRISPRAnalyzerEngine()

st.title("CRISPR Computational Alignment & Target Analysis Workstation") 
st.markdown("A professional bioinformatics terminal designed for calculating multi-vector guide RNA cutting efficiency profiles and off-target sequence safety indices.") 

with st.container(): 
    st.markdown('<div class="dashboard-card">', unsafe_allow_html=True) 
    sequence_to_process = "" 
    
    if input_method == "Manual Text Entry": 
        default_sequence = "ATCGATCGATCGATCGATCGATCGGGATCGATCGATCGATCGATCGGGATCGATCG" 
        sequence_to_process = st.text_area("Input Raw Sequence Base Pairs:", value=default_sequence, height=120) 
    else: 
        uploaded_file = st.file_uploader("Upload Target (.fasta, .txt):", type=["fasta", "fa", "txt"]) 
        if uploaded_file is not None: 
            file_contents = uploaded_file.read().decode("utf-8") 
            sequence_lines = [line for line in file_contents.splitlines() if not line.startswith(">")] 
            sequence_to_process = "".join(sequence_lines) 
    st.markdown('</div>', unsafe_allow_html=True) 

# Process logic
if sequence_to_process: 
    if st.button("Run Comprehensive Diagnostics", type="primary"): 
        with st.spinner("Processing structural alignment data..."): 
            analysis_results = analyzer.process_strand(sequence_to_process) 
            
            if len(analysis_results) == 0: 
                st.error("Processing Fault: No valid PAM anchors identified. Check input requirements.") 
                st.session_state.analysis_df = None 
            else: 
                st.session_state.analysis_df = pd.DataFrame(analysis_results) 

# Render summary layouts if data persists in state
if "analysis_df" in st.session_state and st.session_state.analysis_df is not None: 
    df = st.session_state.analysis_df 
    
    st.write("### Diagnostic Target Summary Analytics") 
    m1, m2, m3 = st.columns(3) 
    m1.metric("Loci Isolated", f"{len(df)} sequences") 
    m2.metric("Optimal Efficiency Found", f"{df['Efficiency Score'].max()}%") 
    m3.metric("Peak Safety Confidence", f"{df['Safety Score'].max()}%") 

    col1, col2 = st.columns(2) 
    
    with col1: 
        st.write("#### Core Variance Scatter Map") 
        fig = px.scatter( 
            df, 
            x="Efficiency Score", 
            y="Safety Score", 
            size="Global Composite Score", 
            color="Global Composite Score", 
            hover_data=["Position", "gRNA Sequence"], 
            labels={ 
                "Efficiency Score": "Cutting Efficiency (%)", 
                "Safety Score": "Safety Threshold (%)", 
                "Global Composite Score": "Composite Score" 
            }, 
            color_continuous_scale=px.colors.sequential.ice 
        ) 
        fig.update_layout( 
            plot_bgcolor="rgba(0,0,0,0)", 
            paper_bgcolor="rgba(0,0,0,0)", 
            font_color="#f1f5f9", 
            margin=dict(l=20, r=20, t=20, b=20) 
        ) 
        
        fig.update_xaxes(showgrid=True, gridcolor="rgba(243, 244, 246, 0.05)", zeroline=False, range=[-5, 110])
        fig.update_yaxes(showgrid=True, gridcolor="rgba(243, 244, 246, 0.05)", zeroline=False, range=[-5, 110])
        
        fig.update_traces(
            hovertemplate="<b>Position:</b> %{customdata[0]}<br>" +
                          "<b>Efficiency:</b> %{x}%<br>" +
                          "<b>Safety:</b> %{y}%<br>" +
                          "<b>gRNA:</b> %{customdata[1]}<extra></extra>"
        )
        st.plotly_chart(fig, use_container_width=True)
        
    with col2:
        st.write("#### Ranked Optimization Log")
        st.dataframe(
            df.set_index("Position"), 
            use_container_width=True,
            column_config={
                "gRNA Sequence": st.column_config.TextColumn("gRNA Sequence"),
                "PAM": st.column_config.TextColumn("PAM"),
                "Efficiency Score": st.column_config.NumberColumn("Efficiency", format="%.1f%%"),
                "Safety Score": st.column_config.NumberColumn("Safety", format="%.1f%%"),
                "Global Composite Score": st.column_config.ProgressColumn(
                    "Global Composite", 
                    format="%.1f", 
                    min_value=0.0, 
                    max_value=100.0,
                )
            }
        )
        
        # Operational Data Exporter
        csv_data = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Export Analytics Log to CSV",
            data=csv_data,
            file_name="crispr_target_export.csv",
            mime="text/csv"
        )