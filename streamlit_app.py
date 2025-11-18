import streamlit as st
import pandas as pd
import os
from pathlib import Path
import requests
import stmol
from stmol import showmol
import py3Dmol

st.set_page_config(
    page_title="PhIP-seq 3D Visualization",
    layout="wide"
)

st.title("PhIP-seq 3D Visualization")
st.markdown("Interactive 3D PhIP-Seq enriched epitope visualization")

@st.cache_data
def load_epitope_data():
    try:
        df = pd.read_csv('data/epitope_data.csv')
        return df
    except FileNotFoundError:
        return pd.DataFrame()

@st.cache_data
def fetch_pdb_content(pdb_id):
    local_path = f'data/pdb_files/{pdb_id}.pdb'
    if os.path.exists(local_path):
        with open(local_path, 'r') as f:
            return f.read()
    
    try:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.text
        else:
            return None
    except Exception as e:
        st.error(f"Error fetching PDB {pdb_id}: {str(e)}")
        return None

st.sidebar.header("Controls")

epitope_df = load_epitope_data()

if not epitope_df.empty:
    available_pdbs = sorted(epitope_df['pdb_id'].unique())
    selected_pdb = st.sidebar.selectbox("Select Structure:", available_pdbs)
    
    min_zscore = float(epitope_df['zscore'].min())
    max_zscore = float(epitope_df['zscore'].max())
    zscore_range = st.sidebar.slider(
        "Z-score Range:", 
        min_value=min_zscore, 
        max_value=max_zscore, 
        value=(min_zscore, max_zscore),
        step=0.1
    )
    
    bg_color = st.sidebar.selectbox(
        "Background:",
        ["white", "black"]
    )
    
    filtered_df = epitope_df[
        (epitope_df['pdb_id'] == selected_pdb) &
        (epitope_df['zscore'] >= zscore_range[0]) &
        (epitope_df['zscore'] <= zscore_range[1])
    ]
    
    col1, col2 = st.columns([1, 1])
    with col1:
        st.metric("Total PhIP-Seq enriched epitopes", len(filtered_df))
    with col2:
        st.metric("PDB ID", selected_pdb)
    
    st.subheader(f"3D Molecular Viewer - {selected_pdb}")
    
    if not filtered_df.empty:
        with st.spinner(f'Loading PDB structure for {selected_pdb}...'):
            pdb_content = fetch_pdb_content(selected_pdb)
        
        if pdb_content:
            st.success(f"Structure loaded: {selected_pdb}")
            
            col1, col2 = st.columns([3, 1])
            
            with col1:
                xyzview = py3Dmol.view(width=800, height=600)
                xyzview.addModel(pdb_content, 'pdb')
                xyzview.setStyle({'cartoon': {'color': 'lightgray'}})
                
                for _, row in filtered_df.iterrows():
                    xyzview.addStyle(
                        {'resi': list(range(int(row['start']), int(row['end'])+1))},
                        {'cartoon': {'color': 'red'}}
                    )
                
                xyzview.setBackgroundColor(bg_color)
                xyzview.zoomTo()
                
                showmol(xyzview, height=600, width=800)
            
            with col2:
                st.markdown("### Epitope Legend")
                st.markdown("---")
                
                epitope_options = []
                epitope_details = {}
                
                for _, row in filtered_df.iterrows():
                    peptide_id = str(row['peptide_id'])
                    # Remove .0 if present (for integer IDs)
                    if peptide_id.endswith('.0'):
                        peptide_id = peptide_id[:-2]
                    
                    option_label = f"Peptide ID {peptide_id} (Z={row['zscore']:.2f})"
                    epitope_options.append(option_label)
                    epitope_details[option_label] = {
                        'peptide_id': peptide_id,
                        'zscore': row['zscore'],
                        'start': int(row['start']),
                        'end': int(row['end']),
                        'length': int(row['length'])
                    }
                
                epitope_options.insert(0, "View Summary")
                
                selected_epitope = st.selectbox(
                    "Enriched epitopes list:",
                    epitope_options,
                    index=0
                )
                
                if selected_epitope == "View Summary":
                    st.markdown(f"**Total epitopes:** {len(filtered_df)}")
                else:
                    details = epitope_details[selected_epitope]
                    st.markdown(f"**Peptide ID:** {details['peptide_id']}")
                    st.markdown(f"**Z-score:** {details['zscore']:.2f}")
                
        else:
            st.error(f"Could not load structure for {selected_pdb}")
    else:
        st.warning("No epitopes found with current filters")

else:
    st.warning("No epitope data found.")
    st.info("Check the pipeline logs.")

st.markdown("---")
st.markdown("**Generated by PhIP-Seq pipeline (SPHERES Lab Team) | Powered by Streamlit | Development Phase (Prototype)**")
