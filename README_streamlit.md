# PhIP-seq Streamlit Dashboard

Interactive dashboard for PhIP-seq epitope analysis with stmol 3D visualization.

## Features

- stmol Integration: Native Streamlit molecular viewer
- Interactive 3D PhIP-Seq enriched epitope visualization
- Red highlighting for all enriched epitopes
- Informational epitope dropdown for details

## Quick Start

pip install -r requirements.txt
chmod +x deploy_streamlit.sh
./deploy_streamlit.sh

Then open: http://localhost:8501

## Usage

1. Select Structure: Choose PDB from dropdown
2. Adjust Z-score range filter
3. Choose background color (white/black)
4. Use epitope dropdown to view individual epitope details
5. All epitopes are highlighted in red on 3D structure
