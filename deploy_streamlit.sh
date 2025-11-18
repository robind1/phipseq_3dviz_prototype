#!/bin/bash
echo "Starting PhIP-seq Streamlit Dashboard with stmol"

if ! command -v streamlit &> /dev/null; then
    echo "Installing dependencies..."
    pip install -r requirements.txt
fi

mkdir -p .streamlit
cp config.toml .streamlit/

echo "Starting app on port 8501..."
echo "Access: http://localhost:8501"

streamlit run streamlit_app.py --server.port=8501
