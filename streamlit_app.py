import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import os
import base64
import json
import requests

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

    elif fmt == 'cif' or fmt == 'mmcif':
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
    bg_color_hex = "#ffffff" if bg_color == "white" else "#000000"
    
    st.sidebar.markdown("---")
    st.sidebar.markdown("**Visualization Options**")
    
    spin_on_load = st.sidebar.checkbox("Auto-Spin", value=False)
    
    highlight_all_chains = st.sidebar.checkbox(
        "Highlight all chains (Multimer)", 
        value=True,
        help="Highlights the sequence in all protein chains."
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
                
                selected_option = st.selectbox("Select Epitope:", epitope_options, index=0)
                
                if selected_option != "View All Enriched":
                    selected_row = epitope_map[selected_option]
                    selected_epitope_id = selected_row['peptide_id']
                    st.markdown("---")
                    st.markdown(f"**Peptide ID:** {selected_epitope_id}")
                    st.markdown(f"**Z-score:** {selected_row['zscore']:.2f}")
                    st.code(selected_row.get('sequence', 'N/A'))
                else:
                    st.markdown("---")
                    st.info(f"Showing all {len(filtered_df)} enriched epitopes.")
            else:
                st.markdown("*No epitopes available.*")
            
            st.markdown("---")
            with st.expander("Structure Info", expanded=False):
                for chain, data in pdb_mapping.items():
                    if data['ids']:
                        st.text(f"Chain {chain}: {min(data['ids'])} - {max(data['ids'])}")

        with view_col:
            molstar_selections = []
            

            for chain_id in pdb_mapping.keys():
                molstar_selections.append({
                    'struct_asym_id': chain_id,
                    'color': {'r': 220, 'g': 220, 'b': 220}, 
                    'focus': False
                })
            
            if not filtered_df.empty:
                if selected_epitope_id is not None:
                    rows_to_highlight = filtered_df[filtered_df['peptide_id'] == selected_epitope_id]
                    color_rgb = {'r': 0, 'g': 255, 'b': 0} 
                else:
                    rows_to_highlight = filtered_df
                    color_rgb = {'r': 255, 'g': 0, 'b': 0}
                
                for _, row in rows_to_highlight.iterrows():
                    peptide_seq = row.get('sequence', '')
                    matches = find_sequence_locations(
                        peptide_seq, 
                        pdb_mapping, 
                        first_match_only=(not highlight_all_chains)
                    )
                    
                    for match in matches:
                        residues = sorted(match['resi'])
                        if not residues: continue
                        
                        ranges = []
                        start = residues[0]
                        prev = residues[0]
                        
                        for r in residues[1:]:
                            if r != prev + 1:
                                ranges.append((start, prev))
                                start = r
                            prev = r
                        ranges.append((start, prev))
                        
                        for r_start, r_end in ranges:
                            molstar_selections.append({
                                'struct_asym_id': match['chain'],
                                'start_residue_number': r_start,
                                'end_residue_number': r_end,
                                'color': color_rgb,
                                'focus': False
                            })

            b64_pdb = base64.b64encode(pdb_content.encode()).decode()
            

            html_code = (
                f'<!DOCTYPE html>'
                f'<html lang="en">'
                f'<head>'
                f'    <meta charset="utf-8" />'
                f'    <meta name="viewport" content="width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0">'
                f'    <link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.2.0/build/pdbe-molstar.css">'
                f'    <script type="text/javascript" src="https://cdn.jsdelivr.net/npm/pdbe-molstar@3.2.0/build/pdbe-molstar-plugin.js"></script>'
                f'    <style>'
                f'        * {{ margin: 0; padding: 0; box-sizing: border-box; }}'
                f'        .viewer-section {{ position: relative; width: 100%; height: 600px; }}'
                f'        #myViewer {{ position: absolute; top: 0; left: 0; right: 0; bottom: 0; }}'
                f'    </style>'
                f'</head>'
                f'<body>'
                f'    <div class="viewer-section">'
                f'        <div id="myViewer"></div>'
                f'    </div>'
                f'    <script>'
                f'        var pdbBase64 = "{b64_pdb}";'
                f'        var pdbFormat = "{pdb_format}";'
                f'        var bgColorHex = "{bg_color_hex}";'
                f'        function hexToRgb(hex) {{'
                f'            var result = /^#?([a-f0-9]{{2}})([a-f0-9]{{2}})([a-f0-9]{{2}})$/i.exec(hex);'
                f'            return result ? {{ r: parseInt(result[1], 16), g: parseInt(result[2], 16), b: parseInt(result[3], 16) }} : {{r: 255, g: 255, b: 255}};'
                f'        }}'
                f'        function b64toBlob(b64Data, contentType="", sliceSize=512) {{'
                f'            var byteCharacters = atob(b64Data);'
                f'            var byteArrays = [];'
                f'            for (var offset = 0; offset < byteCharacters.length; offset += sliceSize) {{'
                f'                var slice = byteCharacters.slice(offset, offset + sliceSize);'
                f'                var byteNumbers = new Array(slice.length);'
                f'                for (var i = 0; i < slice.length; i++) {{ byteNumbers[i] = slice.charCodeAt(i); }}'
                f'                var byteArray = new Uint8Array(byteNumbers);'
                f'                byteArrays.push(byteArray);'
                f'            }}'
                f'            return new Blob(byteArrays, {{type: contentType}});'
                f'        }}'
                f'        var blob = b64toBlob(pdbBase64, "text/plain");'
                f'        var blobUrl = URL.createObjectURL(blob);'
                f'        var viewerInstance = new PDBeMolstarPlugin();'
                f'        var options = {{'
                f'            customData: {{ url: blobUrl, format: pdbFormat, binary: false }},'
                f'            bgColor: hexToRgb(bgColorHex),'
                f'            hideControls: true,'
                f'            hideCanvasControls: ["selection", "controlInfo"]'
                f'        }};'
                f'        var viewerContainer = document.getElementById("myViewer");'
                f'        viewerInstance.render(viewerContainer, options);'
                f'        viewerInstance.events.loadComplete.subscribe(() => {{'
                f'            var selections = {json.dumps(molstar_selections)};'
                f'            if(selections.length > 0) {{'
                f'                try {{ viewerInstance.visual.select({{ data: selections }}); }} catch(err) {{ console.error("Selection error", err); }}'
                f'            }}'
                f'            var spin = {str(spin_on_load).lower()};'
                f'            if(spin) {{ viewerInstance.visual.toggleSpin(true); }}'
                f'        }});'
                f'    </script>'
                f'</body>'
                f'</html>'
            )
            
            components.html(html_code, height=600)
            
    else:
        st.error(f"Could not load structure for {selected_pdb}")

else:
    st.warning("No epitope data found.")

st.markdown("---")
st.markdown("*Generated by PhIP-Seq pipeline (SPHERES Lab Team) | Powered by Streamlit & MolStar | Development Phase (Prototype)*")
