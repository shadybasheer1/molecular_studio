# Molecular Studio — free starter project

A Streamlit cheminformatics app using RDKit for molecular descriptors and structures, py3Dmol for a browser-based 3D viewer, and ReportLab for PDF export.

## Local setup on Windows (recommended: Miniconda)

1. Install Miniconda for Windows: https://www.anaconda.com/docs/getting-started/miniconda/install#windows
2. Open **Anaconda Prompt**.
3. Run:

```bash
conda create -n molecular-studio -c conda-forge python=3.12 rdkit -y
conda activate molecular-studio
pip install streamlit py3Dmol pandas reportlab pillow
```

4. Extract this ZIP to a folder, e.g. `C:\MolecularStudio`.
5. In Anaconda Prompt, run:

```bash
cd C:\MolecularStudio
streamlit run app.py
```

6. Open the local URL printed in the terminal (usually http://localhost:8501).

## Deploy for free

1. Create a GitHub repository and upload `app.py`, `requirements.txt`, and `README.md`.
2. Go to https://share.streamlit.io/ and sign in with GitHub.
3. Choose **Create app**, select the repository, branch `main`, and file `app.py`.
4. Click **Deploy**. Streamlit Community Cloud provides a shareable `streamlit.app` URL, subject to its platform limits.

## Notes

- The 3D viewer loads the 3Dmol.js library in the browser, so the viewer needs internet access.
- RDKit calculates descriptors from the supplied SMILES. Descriptor values are model-based and may differ from experimental measurements.
- Sigma/pi bond-component counts use a documented Kekulized-structure convention and should be checked against the convention expected in your course.
- This is an educational starter project. Test all examples and generated PDFs before submitting it.
