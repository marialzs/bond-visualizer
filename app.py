import streamlit as st
import py3Dmol

st.set_page_config(page_title="BondViz:Protein–Ligand Interaction Explorer", layout="wide")
st.title("BondViz:Protein–Ligand Interaction Explorer")

uploaded_pdb = st.file_uploader("Upload a PDB file", type=["pdb", "cif"])
lig_text = st.text_input("Specify ligand ID (optional):")
if not uploaded_pdb:
    st.info("Drop a .pdb to begin.")
    st.stop()

pdb_file = uploaded_pdb.read().decode("utf-8")

view = py3Dmol.view(width=800, height=600)
view.addModel(pdb_file, "pdb")
view.setStyle({'cartoon': {'color': '#607D8B'}})
view.addStyle({'resn': lig_text}, {'stick': {}})
view.zoomTo()

st.components.v1.html(view._make_html(), height=600, width=800)