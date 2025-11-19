import streamlit as st
import pandas as pd
import os
from pathlib import Path
import requests
import stmol
from stmol import showmol
import py3Dmol;

st.set_page_config(
    page_title="PhIP-seq 3D Visualization",
    layout="wide"
)

AA_MAP = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def get_pdb_residue_mapping(pdb_content, fmt='pdb'):
    chains = {}
    
    lines = pdb_content.splitlines()
    
    if fmt == 'pdb':
        for line in lines:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                res_name = line[17:20].strip()
                chain_id = line[21]
                try:
                    res_seq = int(line[22:26])
                except ValueError:
                    continue 
                
                if chain_id not in chains:
                    chains[chain_id] = {'seq': '', 'ids': []}
                
                if res_name in AA_MAP:
                    chains[chain_id]['seq'] += AA_MAP[res_name]
                    chains[chain_id]['ids'].append(res_seq)

    elif fmt == 'cif':
        loop_indices = {}
        in_loop = False
        
        for line in lines:
            line = line.strip()
            if line.startswith('loop_'):
                in_loop = True
                loop_indices = {}
                current_col_idx = 0
                continue
            
            if in_loop and line.startswith('_atom_site.'):
                key = line.split('.')[1].strip()
                loop_indices[key] = current_col_idx
                current_col_idx += 1
                continue
            
            if in_loop and (line.startswith('ATOM') or line.startswith('HETATM')):
                if 'label_atom_id' in loop_indices and 'label_comp_id' in loop_indices and 'label_asym_id' in loop_indices and 'label_seq_id' in loop_indices:
                    parts = line.split()
                    try:
                        atom_idx = loop_indices['label_atom_id']
                        if len(parts) > atom_idx and parts[atom_idx] == 'CA':
                            res_name = parts[loop_indices['label_comp_id']]
                            chain_id = parts[loop_indices['label_asym_id']]
                            res_seq = int(parts[loop_indices['label_seq_id']])
                            
                            if chain_id not in chains:
                                chains[chain_id] = {'seq': '', 'ids': []}
                            
                            if res_name in AA_MAP:
                                chains[chain_id]['seq'] += AA_MAP[res_name]
                                chains[chain_id]['ids'].append(res_seq)
                    except (ValueError, IndexError):
                        continue
    
    return chains


def find_sequence_locations(peptide_seq, pdb_mapping, first_match_only=False):
    matches = []
    if not peptide_seq or peptide_seq == 'Unknown' or pd.isna(peptide_seq):
        return matches
        
    peptide_seq = str(peptide_seq).strip()
    
    found_exact = False
    for chain_id, data in pdb_mapping.items():
        seq = data['seq']
        start = 0
        while True:
            idx = seq.find(peptide_seq, start)
            if idx == -1:
                break
            
            end_idx = idx + len(peptide_seq)
            residues = data['ids'][idx:end_idx]
            matches.append({'chain': chain_id, 'resi': residues, 'method': 'exact'})
            found_exact = True
            
            if first_match_only:
                return matches
            start = idx + 1
            
    if found_exact:
        return matches


    WINDOW_SIZE = 7
    if len(peptide_seq) < WINDOW_SIZE:
        WINDOW_SIZE = len(peptide_seq)
        
    fragment_matches = {} 
    
    for i in range(len(peptide_seq) - WINDOW_SIZE + 1):
        fragment = peptide_seq[i : i + WINDOW_SIZE]
        
        for chain_id, data in pdb_mapping.items():
            seq = data['seq']
            idx = seq.find(fragment)
            if idx != -1:
                if chain_id not in fragment_matches:
                    fragment_matches[chain_id] = set()
                
                res_ids = data['ids'][idx : idx + WINDOW_SIZE]
                fragment_matches[chain_id].update(res_ids)

    if fragment_matches:
        for chain_id, res_set in fragment_matches.items():
            matches.append({
                'chain': chain_id, 
                'resi': sorted(list(res_set)),
                'method': 'partial'
            })
            if first_match_only:
                break
                
    return matches

st.title("PhIP-seq 3D Visualization")
st.markdown("Interactive 3D PhIP-Seq enriched epitope visualization")

@st.cache_data
def load_epitope_data():
    try:
        df = pd.read_csv('data/epitope_data.csv')
        if 'sample' in df.columns:
            df['sample'] = df['sample'].astype(str)
        return df
    except FileNotFoundError:
        return pd.DataFrame()

@st.cache_data
def fetch_pdb_content(pdb_id):
    local_path_pdb = f'data/pdb_files/{pdb_id}.pdb'
    if os.path.exists(local_path_pdb):
        with open(local_path_pdb, 'r') as f:
            return f.read(), 'pdb'
            
    local_path_cif = f'data/pdb_files/{pdb_id}.cif'
    if os.path.exists(local_path_cif):
        with open(local_path_cif, 'r') as f:
            return f.read(), 'cif'
    
    try:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.text, 'pdb'
        else:
            return None, None
    except Exception as e:
        st.error(f"Error fetching PDB {pdb_id}: {str(e)}")
        return None, None

st.sidebar.header("Controls")

epitope_df = load_epitope_data()

if not epitope_df.empty:
    available_pdbs = sorted(epitope_df['pdb_id'].unique())
    selected_pdb = st.sidebar.selectbox("Select Structure:", available_pdbs)
    
    pdb_subset_df = epitope_df[epitope_df['pdb_id'] == selected_pdb]
    
    if 'sample' in pdb_subset_df.columns:
        available_samples = ["View Structure Only"] + sorted(pdb_subset_df['sample'].unique())
        selected_sample = st.sidebar.selectbox("Select Sample:", available_samples)
    else:
        selected_sample = "View Structure Only"
        st.sidebar.warning("No sample information found in data.")
    
    bg_color = st.sidebar.selectbox("Background:", ["white", "black"])
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("**Visualization Options**")
    highlight_all_chains = st.sidebar.checkbox(
        "Highlight all chains (Multimer)", 
        value=True,
        help="Highlights the sequence in all protein chains. Uncheck to see only in one chain."
    )
    
    if selected_sample and selected_sample != "View Structure Only":
        filtered_df = pdb_subset_df[
            (pdb_subset_df['sample'] == selected_sample)
        ]
    else:
        filtered_df = pd.DataFrame(columns=pdb_subset_df.columns)
    
    col1, col2 = st.columns(2)
    with col1:
        st.metric("Total PhIP-Seq enriched epitopes", len(filtered_df))
    with col2:
        st.metric("PDB/CIF ID", selected_pdb)
    
    st.subheader(f"3D Molecular Viewer - {selected_pdb}")
    
    with st.spinner(f'Loading PDB structure for {selected_pdb}...'):
        pdb_content, pdb_format = fetch_pdb_content(selected_pdb)
    
    if pdb_content:
        view_col, details_col = st.columns([3, 1])
        
        pdb_mapping = get_pdb_residue_mapping(pdb_content, pdb_format)
        
        selected_epitope_id = None
        
        with details_col:
            st.markdown("### Epitope Selection")
            st.markdown("Select an epitope below to isolate it in the viewer.")
            
            if not filtered_df.empty:
                epitope_options = []
                epitope_map = {} 
                
                for _, row in filtered_df.iterrows():
                    peptide_id = str(row['peptide_id'])
                    if peptide_id.endswith('.0'):
                        peptide_id = peptide_id[:-2]
                    
                    option_label = f"Peptide ID {peptide_id} (Z={row['zscore']:.2f})"
                    epitope_options.append(option_label)
                    epitope_map[option_label] = row
                
                epitope_options.insert(0, "View All Enriched")
                
                selected_option = st.selectbox(
                    "Select Epitope:",
                    epitope_options,
                    index=0
                )
                
                if selected_option != "View All Enriched":
                    selected_row = epitope_map[selected_option]
                    selected_epitope_id = selected_row['peptide_id']
                    
                    st.markdown("---")
                    st.markdown(f"**Peptide ID:** {selected_epitope_id}")
                    st.markdown(f"**Z-score:** {selected_row['zscore']:.2f}")
                    st.markdown(f"**Sequence:**")
                    st.code(selected_row.get('sequence', 'N/A'))

                else:
                    st.markdown("---")
                    st.info(f"Showing all {len(filtered_df)} enriched epitopes.")

            else:
                st.markdown("*No epitopes available.*")
            
            st.markdown("---")
            with st.expander("PDB/CIF Structure Info", expanded=False):
                for chain, data in pdb_mapping.items():
                    if data['ids']:
                        st.text(f"Chain {chain}: Residues {min(data['ids'])} - {max(data['ids'])}")
                    else:
                        st.text(f"Chain {chain}: No residues found")

        with view_col:
            xyzview = py3Dmol.view(width=800, height=600)
            xyzview.addModel(pdb_content, pdb_format)
            xyzview.setStyle({'cartoon': {'color': 'lightgray'}})
            
            if not filtered_df.empty:
                
                if selected_epitope_id is not None:
                    rows_to_highlight = filtered_df[filtered_df['peptide_id'] == selected_epitope_id]
                    highlight_color = 'green'
                else:
                    rows_to_highlight = filtered_df
                    highlight_color = 'red'
                
                found_any = False
                for _, row in rows_to_highlight.iterrows():
                    peptide_seq = row.get('sequence', '')
                    
                    matches = find_sequence_locations(
                        peptide_seq, 
                        pdb_mapping, 
                        first_match_only=(not highlight_all_chains)
                    )
                    
                    if matches:
                        found_any = True
                        for match in matches:
                            xyzview.addStyle(
                                {'chain': match['chain'], 'resi': match['resi']},
                                {'cartoon': {'color': highlight_color}}
                            )

                if selected_epitope_id is not None:
                    if found_any:
                        st.success(f"Highlighted peptide {selected_epitope_id}")
                    else:
                        st.warning(f"Peptide {selected_epitope_id} not found. Check 'PDB Structure Coverage Info' on the right to see valid residue ranges.")

            elif selected_sample != "View Structure Only":
                st.info(f"No epitopes found for {selected_sample}. Showing base structure.")
            
            xyzview.setBackgroundColor(bg_color)
            xyzview.zoomTo()
            showmol(xyzview, height=600, width=800)
            
    else:
        st.error(f"Could not load structure for {selected_pdb}")

else:
    st.warning("No epitope data found.")
    st.info("Check the pipeline logs for more details.")

st.markdown("---")
st.markdown("**Generated by PhIP-Seq pipeline (SPHERES Lab Team) | Powered by Streamlit | Development Phase (Prototype)**")
