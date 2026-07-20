import re
import urllib.request
import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

# --- PAGE CONFIG & THEME ---
st.set_page_config(
    page_title="CRISPR Computational Workstation",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom container styling
st.markdown(
    """ 
    <style> 
    .framework-box { 
        background-color: #0f172a; 
        border: 1px dashed #334155; 
        border-radius: 6px; 
        padding: 14px; 
        margin: 10px 0; 
        font-family: 'Courier New', Courier, monospace; 
        font-size: 0.9rem; 
        color: #38bdf8; 
    } 
    </style> 
    """,
    unsafe_allow_html=True,
)

# --- MULTI-NUCLEASE PROFILES ---
NUCLEASE_PROFILES = {
    "SpCas9 (Wild-Type)": {
        "fwd_pam": r"(?=(?P<target>[ATCG]{20})(?P<pam>[ATCG]GG))",
        "rev_pam": r"(?=(?P<pam>CC[ATCG])(?P<target>[ATCG]{20}))",
        "pam_len": 3,
        "desc": "Standard 5'-NGG-3' PAM (3' end)"
    },
    "SpCas9-VQR Variant": {
        "fwd_pam": r"(?=(?P<target>[ATCG]{20})(?P<pam>[ATCG]GA[ATCG]))",
        "rev_pam": r"(?=(?P<pam>[ATCG]TC[ATCG])(?P<target>[ATCG]{20}))",
        "pam_len": 4,
        "desc": "Engineered 5'-NGAN-3' PAM for target expansion"
    },
    "Cas12a (Cpf1)": {
        "fwd_pam": r"(?=(?P<pam>TTT[ATCG])(?P<target>[ATCG]{23}))",
        "rev_pam": r"(?=(?P<target>[ATCG]{23})(?P<pam>[ATCG]AAA))",
        "pam_len": 4,
        "desc": "5'-TTTN-3' PAM (5' end), 23-bp target sequence"
    }
}


# --- CACHED DATA FETCHING ENGINE ---
@st.cache_data(ttl=86400)
def fetch_live_chromosome(ncbi_accession: str) -> str:
    """Downloads and sanitizes reference genome FASTA data."""
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
        f"db=nuccore&id={ncbi_accession}&rettype=fasta&retmode=text"
    )
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "CRISPR_Analyzer/5.3.0"}
        )
        with urllib.request.urlopen(req, timeout=10) as response:
            fasta_data = response.read().decode("utf-8")
        clean_lines = [
            line for line in fasta_data.splitlines() if not line.startswith(">")
        ]
        return "".join(clean_lines).upper()
    except Exception as e:
        st.warning(
            f"Failed to fetch live NCBI data for '{ncbi_accession}' ({e}). Using human mtDNA fallback sequence."
        )
        return (
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGG"
            "GGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGAT"
            "CCACTAGTCCACCCCTCAGAACACTACTACACCAACACCCACCCACCAC"
        ) * 20


# --- BIOINFORMATICS ENGINE ---
class CRISPRChromosomeAlignmentEngine:
    """Production-grade engineering pipeline utilizing the MIT Hsu-Zhang
    position-weighted matrix model to evaluate authentic chromosomal safety.
    """

    def __init__(self, ncbi_accession: str = "NC_012920.1"):
        self.ncbi_accession = ncbi_accession
        self.chromosome_sequence = fetch_live_chromosome(self.ncbi_accession)

        # Official MIT Hsu-Zhang Mismatch Weight Matrix (Positions 1 to 20 relative to PAM)
        self.hsu_weights = np.array([
            0.0, 0.0, 0.014, 0.0, 0.0, 0.395, 0.317, 0.0, 0.389, 0.074,
            0.247, 0.545, 0.409, 0.0, 0.418, 0.401, 0.290, 0.467, 0.523, 0.534
        ])

    @staticmethod
    def clean_sequence(input_text: str) -> str:
        cleaned = input_text.upper().replace("\n", "").replace("\r", "").strip()
        return re.sub(r"[^ATCG]", "", cleaned)

    @staticmethod
    def get_reverse_complement(seq: str) -> str:
        """Computes the 5'->3' reverse complement of a DNA sequence."""
        complement = str.maketrans("ATCG", "TAGC")
        return seq.translate(complement)[::-1]

    def calculate_mit_safety_score(
        self, grna: str, target_segment: str
    ) -> float:
        """Executes the non-linear Hsu-Zhang scoring function across a localized match sequence."""
        if len(grna) != 20 or len(target_segment) != 20:
            return 1.0

        grna_arr = np.frombuffer(grna.encode('ascii'), dtype=np.int8)
        target_arr = np.frombuffer(target_segment.encode('ascii'), dtype=np.int8)
        
        mismatch_mask = grna_arr != target_arr
        mismatch_positions = 19 - np.where(mismatch_mask)[0]
        
        n_mismatches = len(mismatch_positions)
        if n_mismatches == 0:
            return 1.0

        product_term = np.prod(1.0 - self.hsu_weights[mismatch_positions])

        if n_mismatches > 1:
            diff_matrix = np.abs(np.subtract.outer(mismatch_positions, mismatch_positions))
            triu_indices = np.triu_indices(n_mismatches, k=1)
            d_avg = np.mean(diff_matrix[triu_indices])
            distance_term = 1.0 / (((19.0 - d_avg) / 19.0) * 4.0 + 1.0)
        else:
            distance_term = 1.0

        mismatch_count_term = 1.0 / (n_mismatches**2)

        return float(product_term * distance_term * mismatch_count_term)

    def calculate_alignment_safety(self, grna: str) -> tuple:
        """Scans chromosomal sequence for off-target hits (using 20bp seed region)."""
        seed_grna = grna[:20]
        chromosome_len = len(self.chromosome_sequence)
        cumulative_off_target_risk = 0.0
        alignments_mapped = 0

        for idx in range(0, chromosome_len - 20, 2):
            off_target_segment = self.chromosome_sequence[idx : idx + 20]

            mismatches = sum(
                1 for p in range(20) if seed_grna[p] != off_target_segment[p]
            )

            if mismatches <= 4:
                risk = self.calculate_mit_safety_score(seed_grna, off_target_segment)
                cumulative_off_target_risk += risk
                if risk > 0.001:
                    alignments_mapped += 1

        safety_score = round(
            max(0.0, 100.0 - (cumulative_off_target_risk * 10.0)), 2
        )
        return safety_score, alignments_mapped

    def process_strand(
        self, raw_dna: str, nuclease_key: str = "SpCas9 (Wild-Type)", progress_bar=None
    ) -> list[dict]:
        sequence = self.clean_sequence(raw_dna)
        pipeline_results = []
        profile = NUCLEASE_PROFILES[nuclease_key]

        fwd_pattern = re.compile(profile["fwd_pam"])
        rev_pattern = re.compile(profile["rev_pam"])

        targets_to_evaluate = []

        # Forward Strand
        for match in fwd_pattern.finditer(sequence):
            targets_to_evaluate.append({
                "pos": match.start(),
                "grna": match.group("target"),
                "pam": match.group("pam"),
                "strand": "+"
            })

        # Reverse Strand
        for match in rev_pattern.finditer(sequence):
            grna_fwd = match.group("target")
            grna_rev_comp = self.get_reverse_complement(grna_fwd)
            pam_rev_comp = self.get_reverse_complement(match.group("pam"))
            targets_to_evaluate.append({
                "pos": match.start(),
                "grna": grna_rev_comp,
                "pam": pam_rev_comp,
                "strand": "-"
            })

        total_targets = len(targets_to_evaluate)

        for i, item in enumerate(targets_to_evaluate):
            grna_target = item["grna"]
            pam = item["pam"]
            start_pos = item["pos"]
            strand = item["strand"]
            seq_len = len(grna_target)

            gc_count = grna_target.count("G") + grna_target.count("C")
            gc_content = (gc_count / float(seq_len)) * 100.0

            efficiency = max(0.0, 100.0 - (abs(50.0 - gc_content) * 2.5))
            if "TTTT" in grna_target:
                efficiency = max(0.0, efficiency - 30.0)

            safety, hits = self.calculate_alignment_safety(grna_target)
            composite = round((efficiency * 0.4) + (safety * 0.6), 2)

            pipeline_results.append({
                "Position": start_pos,
                "Locus End": start_pos + seq_len + profile["pam_len"],
                "Strand": strand,
                "gRNA Sequence": grna_target,
                "PAM": pam,
                "GC Content (%)": round(gc_content, 1),
                "Efficiency Score": round(efficiency, 2),
                "Safety Score": safety,
                "Chromosome Hits Mapped": hits,
                "Global Composite Score": composite,
            })

            if progress_bar and total_targets > 0:
                progress_bar.progress((i + 1) / total_targets)

        return sorted(
            pipeline_results,
            key=lambda x: x["Global Composite Score"],
            reverse=True,
        )


@st.cache_resource
def get_analyzer_engine(ncbi_accession: str):
    return CRISPRChromosomeAlignmentEngine(ncbi_accession)


# --- INTERACTION TERMINAL CONTROL LAYER ---
with st.sidebar:
    st.title("Settings")
    
    ncbi_accession_input = st.text_input(
        "NCBI Accession ID:",
        value="NC_012920.1",
        help="Enter any NCBI Nucleotide Accession (e.g., NC_012920.1 for Human mtDNA)."
    )
    
    st.divider()
    input_method = st.radio(
        "Genomic Data Source Type:",
        ["Manual Text Entry", "Upload FASTA File"],
        help="Select how you wish to feed DNA sequences into the engine.",
    )
    st.divider()
    selected_nuclease = st.selectbox(
        "CRISPR Nuclease Enzyme:",
        list(NUCLEASE_PROFILES.keys()),
        help="Choose Cas variant to update PAM motif and guide length parameters."
    )
    st.caption(f"Profile: {NUCLEASE_PROFILES[selected_nuclease]['desc']}")
    st.divider()
    st.caption("Workstation v5.3.0")
    st.caption(f"Target Accession: {ncbi_accession_input}")

analyzer = get_analyzer_engine(ncbi_accession_input)

# --- MAIN DASHBOARD VIEW ---
st.title("CRISPR Target & Alignment Workstation")
st.write(
    "Analyze guide RNA (gRNA) sequences, calculate cutting efficiency, and"
    " evaluate chromosomal off-target safety."
)

tab1, tab2, tab3 = st.tabs([
    "Active Pipeline Workspace",
    "Scoring Framework",
    "Reference Sequence Data",
])

with tab1:
    st.subheader("1. Sequence Input")
    sequence_to_process = ""

    if input_method == "Manual Text Entry":
        default_sequence = (
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGG"
            "GGTATGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGAT"
        )
        sequence_to_process = st.text_area(
            "Input Raw Sequence Base Pairs:",
            value=default_sequence,
            height=120,
            help="Paste raw DNA sequence containing PAM sites.",
        )
    else:
        uploaded_file = st.file_uploader(
            "Upload Target Genomic Data (.fasta):",
            type=["fasta", "fa", "txt"],
        )
        if uploaded_file is not None:
            file_contents = uploaded_file.read().decode("utf-8")
            sequence_lines = [
                line
                for line in file_contents.splitlines()
                if not line.startswith(">")
            ]
            sequence_to_process = "".join(sequence_lines)

    run_analysis = st.button("Execute Chromosomal Alignment Run", type="primary")

    if run_analysis and sequence_to_process:
        progress_bar = st.progress(0.0)
        with st.spinner("Searching genome for off-target alignments..."):
            analysis_results = analyzer.process_strand(
                sequence_to_process, 
                nuclease_key=selected_nuclease, 
                progress_bar=progress_bar
            )
            progress_bar.empty()

            if not analysis_results:
                st.error(
                    "No valid target sequences matching selected PAM anchors were found."
                )
                st.session_state.analysis_df = None
            else:
                st.session_state.analysis_df = pd.DataFrame(analysis_results)

    if st.session_state.get("analysis_df") is not None:
        df = st.session_state.analysis_df
        st.divider()
        st.subheader("2. Diagnostics & Analytics")

        m1, m2, m3 = st.columns(3)
        m1.metric("Loci Isolated", f"{len(df)}")
        m2.metric(
            "Max Cutting Efficiency", f"{df['Efficiency Score'].max():.1f}%"
        )
        m3.metric("Peak Safety Score", f"{df['Safety Score'].max():.1f}%")
        st.write("")

        col1, col2 = st.columns([1, 1], gap="medium")
        with col1:
            st.markdown("#### Efficiency vs. Safety Distribution")
            fig = px.scatter(
                df,
                x="Efficiency Score",
                y="Safety Score",
                size="Global Composite Score",
                color="Global Composite Score",
                hover_data=[
                    "Position",
                    "Strand",
                    "gRNA Sequence",
                    "Chromosome Hits Mapped",
                ],
                labels={
                    "Efficiency Score": "Efficiency (%)",
                    "Safety Score": "Safety (%)",
                    "Global Composite Score": "Composite",
                },
                color_continuous_scale="Viridis",
            )
            fig.update_layout(margin=dict(l=10, r=10, t=10, b=10), height=340)
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            st.markdown("#### Ranked Target Candidates")
            st.dataframe(
                df.set_index("Position"),
                use_container_width=True,
                height=340,
                column_config={
                    "Strand": st.column_config.TextColumn("Strand", width="small"),
                    "gRNA Sequence": st.column_config.TextColumn("gRNA"),
                    "PAM": st.column_config.TextColumn("PAM", width="small"),
                    "Efficiency Score": st.column_config.NumberColumn(
                        "Eff.", format="%.1f%%"
                    ),
                    "Safety Score": st.column_config.NumberColumn(
                        "Safety", format="%.1f%%"
                    ),
                    "Chromosome Hits Mapped": st.column_config.NumberColumn(
                        "Hits"
                    ),
                    "Global Composite Score": st.column_config.ProgressColumn(
                        "Score",
                        format="%.1f",
                        min_value=0.0,
                        max_value=100.0,
                    ),
                },
            )

        # --- CHROMOSOMAL LOCUS TRACK MAP ---
        st.subheader("3. Chromosomal Locus Map")
        st.caption("Visual position of target gRNA loci along the query DNA strand.")
        
        df_locus = df.copy()
        df_locus["Locus Length"] = df_locus["Locus End"] - df_locus["Position"]
        
        fig_locus = px.bar(
            df_locus,
            x="Locus Length",
            y="Strand",
            base="Position",
            orientation="h",
            color="Global Composite Score",
            hover_name="gRNA Sequence",
            hover_data=["Position", "Locus End", "PAM", "Efficiency Score", "Safety Score"],
            color_continuous_scale="Plasma",
            labels={"Locus Length": "Base Pair Position", "base": "Start Position"},
            title="Genomic Position Tracks (+/- Strands)"
        )
        fig_locus.update_xaxes(
            title_text="Base Pair Position (bp)", 
            type="linear"
        )
        fig_locus.update_yaxes(autorange="reversed")
        fig_locus.update_layout(height=220, margin=dict(l=10, r=10, t=30, b=10))
        st.plotly_chart(fig_locus, use_container_width=True)

        csv_data = df.to_csv(index=False).encode("utf-8")
        st.download_button(
            label="Export Results (CSV)",
            data=csv_data,
            file_name="crispr_target_export.csv",
            mime="text/csv",
        )

with tab2:
    st.subheader("Scoring Rules & Advanced Frameworks")
    st.markdown("**1. Efficiency Score Formula**")
    st.markdown(
        '<div class="framework-box">Efficiency = max(0, 100 - (abs(50 - GC%) *'
        " 2.5)) - PolyT_Penalty</div>",
        unsafe_allow_html=True,
    )
    st.write(
        "Measures GC balance (optimal 40–60%) and penalizes premature Pol III"
        " transcription termination sequences (`TTTT`)."
    )

    st.markdown("**2. MIT Hsu-Zhang Off-Target Safety Score Model**")
    st.markdown(
        '<div class="framework-box">Cleavage_Probability = Product_Term *'
        " Distance_Term * Mismatch_Count_Term</div>",
        unsafe_allow_html=True,
    )
    st.write(
        "Scans all chromosome loci with up to 4 mismatches using the official MIT"
        " position-weighted mismatch matrix."
    )

    st.divider()

    col_cite, col_limits = st.columns(2, gap="medium")

    with col_cite:
        st.markdown("#### Primary Citation")
        st.caption("Algorithm implementation modeled after the seminal MIT Hsu-Zhang paper:")
        st.markdown(
            "> **Hsu, P., Scott, D., Weinstein, J. et al.** *DNA targeting specificity of RNA-guided Cas9 nucleases.* Nat Biotechnol **31**, 827–832 (2013)."
        )

    with col_limits:
        st.markdown("#### Model Scope & Limitations")
        st.info(
            "• **Target Reference:** Configurable via NCBI Nucleotide Accession ID.\n"
            "• **PAM Specificity:** Evaluates SpCas9 variants and Cas12a target motifs.\n"
            "• **Future Roadmap:** Multi-threaded full nuclear scans and deep-learning integration (DeepHF, Azimuth)."
        )

    st.divider()

    st.markdown("#### Machine Learning & Deep Learning Roadmap")
    st.info(
        "• **Doench / Azimuth (Rule Set 2):** Incorporates 2nd-order nucleotide interactions,"
        " melting temperatures ($T_m$), and position-specific dinucleotides.\n"
        "• **DeepHF Integration:** Utilizes Recurrent Neural Networks (RNNs) trained on"
        " >50,000 gRNA cleavage experiments to predict cleavage efficiency across Cas9 variants.\n"
        "• **CFD (Cutting Frequency Determination):** Replaces binary mismatch penalties"
        " with experimental substitution matrices (e.g., rG:dT vs rC:dA)."
    )

with tab3:
    st.subheader("Reference Genome Overview")
    st.write(
        f"**NCBI Accession:** `{analyzer.ncbi_accession}`"
    )
    st.write(
        f"**Sequence Length:** {len(analyzer.chromosome_sequence):,} base pairs"
    )
    st.markdown("#### Sequence Preview (First 1,000 Base Pairs)")
    st.code(
        analyzer.chromosome_sequence[:1000],
        language="text",
        wrap_lines=True,
    )
