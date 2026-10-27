import io
import math
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Draw
from rdkit.Chem import Crippen, Lipinski
from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import mm
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
)

APP_TITLE = "Molecular Studio"

EXAMPLES = {
    "Water": ("O", "A small inorganic molecule used as a reference for polarity and hydrogen bonding."),
    "Methane": ("C", "The simplest alkane, with a tetrahedral carbon centre."),
    "Ethane": ("CC", "A saturated hydrocarbon containing a carbon-carbon single bond."),
    "Ethanol": ("CCO", "An alcohol containing a hydroxyl group; commonly used as a baseline example."),
    "Acetone": ("CC(=O)C", "A ketone with a polar carbonyl group."),
    "Acetic acid": ("CC(=O)O", "A carboxylic acid with hydrogen-bonding functionality."),
    "Benzene": ("c1ccccc1", "An aromatic six-membered carbon ring."),
    "Toluene": ("Cc1ccccc1", "An aromatic ring bearing a methyl substituent."),
    "Phenol": ("Oc1ccccc1", "An aromatic compound containing a hydroxyl group."),
    "Aniline": ("Nc1ccccc1", "An aromatic amine."),
    "Glycine": ("NCC(=O)O", "The simplest proteinogenic amino acid in its neutral SMILES representation."),
    "Caffeine": ("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "A nitrogen-containing heterocyclic compound found in coffee and tea."),
    "Aspirin": ("CC(=O)Oc1ccccc1C(=O)O", "Acetylsalicylic acid, containing ester and carboxylic-acid groups."),
    "Glucose": ("C(C1C(C(C(C(O1)O)O)O)O)O", "A six-carbon sugar represented here as a cyclic structure."),
    "Lactic acid": ("CC(O)C(=O)O", "An alpha-hydroxy carboxylic acid.")
}

st.set_page_config(page_title=APP_TITLE, page_icon="🧪", layout="wide")

st.markdown("""
<style>
.stApp {background: linear-gradient(135deg, #0b1220 0%, #101b2d 55%, #101827 100%); color: #edf4ff;}
[data-testid="stSidebar"] {background: #0c1524;}
h1, h2, h3 {color: #e9f8ff !important;}
p, label, .stMarkdown {color: #d7e2f0;}
div[data-testid="stMetric"] {background: #14243a; border: 1px solid #24415d; padding: 14px; border-radius: 14px;}
div[data-testid="stMetricLabel"] {color: #9fd9e8;}
div[data-testid="stMetricValue"] {color: #ffffff;}
.stButton button, .stDownloadButton button {border-radius: 10px; border: 1px solid #32c6bd; background: #123b4b; color: white;}
.stButton button:hover, .stDownloadButton button:hover {border-color: #78f0df; color: #ffffff;}
</style>
""", unsafe_allow_html=True)

st.title("🧪 Molecular Studio")
st.caption("Interactive cheminformatics laboratory · Molecular structures · Descriptors · Educational reports")
st.info("Educational tool: calculated descriptors depend on the input structure and software conventions. They are not experimental measurements or a substitute for expert chemical review.")

def count_bond_components(mol):
    """Count sigma/pi bond components using a Kekulized copy.
    Each bonded atom pair contributes one sigma component; each double bond adds one pi,
    each triple bond adds two pi. Aromatic bonds are kekulized first where possible.
    Hydrogens are added so X-H bonds are included.
    """
    kek = Chem.Mol(mol)
    try:
        Chem.Kekulize(kek, clearAromaticFlags=True)
    except Exception:
        # Fallback: count sigma for each bond and infer pi from non-aromatic bond order.
        sigma = mol.GetNumBonds()
        pi = 0
        for b in mol.GetBonds():
            if b.GetIsAromatic():
                continue
            pi += max(0, int(round(b.GetBondTypeAsDouble())) - 1)
        return sigma, pi, "Fallback convention: aromatic pi components may be undercounted."
    pi = 0
    for b in kek.GetBonds():
        order = int(round(b.GetBondTypeAsDouble()))
        pi += max(0, order - 1)
    with_h = Chem.AddHs(kek)
    sigma = with_h.GetNumBonds()
    return sigma, pi, "Kekulized representation: each atom-pair bond contributes one σ; each additional bond order contributes a π component. Explicit/added X–H bonds are included."

def safe_descriptor(label, fn, unit="", digits=3):
    try:
        value = fn()
        if isinstance(value, float):
            value = round(value, digits)
        return value, unit
    except Exception:
        return "Not available", unit

def analyze_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse this SMILES. Check the spelling, brackets, ring numbers, and bond symbols.")
    Chem.SanitizeMol(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    sigma, pi, bond_note = count_bond_components(mol)
    n_atoms = mol.GetNumAtoms()
    heavy = mol.GetNumHeavyAtoms()
    heavy_bonds = sum(1 for b in mol.GetBonds()
                      if mol.GetAtomWithIdx(b.GetBeginAtomIdx()).GetAtomicNum() > 1
                      and mol.GetAtomWithIdx(b.GetEndAtomIdx()).GetAtomicNum() > 1)
    carbon_count = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6)
    sp3_c = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetHybridization() == Chem.HybridizationType.SP3)
    frac_csp3 = (sp3_c / carbon_count) if carbon_count else 0.0
    props = [
        ("Molecular formula", formula, ""),
        ("Molecular weight", round(Descriptors.MolWt(mol), 3), "g/mol"),
        ("Exact molecular weight", round(Descriptors.ExactMolWt(mol), 4), "Da"),
        ("LogP (Crippen)", round(Crippen.MolLogP(mol), 3), ""),
        ("TPSA", round(rdMolDescriptors.CalcTPSA(mol), 2), "Å²"),
        ("Hydrogen bond donors", rdMolDescriptors.CalcNumHBD(mol), ""),
        ("Hydrogen bond acceptors", rdMolDescriptors.CalcNumHBA(mol), ""),
        ("Rotatable bonds", rdMolDescriptors.CalcNumRotatableBonds(mol), ""),
        ("Ring count", rdMolDescriptors.CalcNumRings(mol), ""),
        ("Aromatic rings", rdMolDescriptors.CalcNumAromaticRings(mol), ""),
        ("Aliphatic rings", rdMolDescriptors.CalcNumAliphaticRings(mol), ""),
        ("Number of atoms", n_atoms, "heavy + explicit atoms in input graph; implicit H excluded"),
        ("Heavy atoms", heavy, ""),
        ("Heavy atom bonds", heavy_bonds, ""),
        ("Sigma (σ) bond components", sigma, ""),
        ("Pi (π) bond components", pi, ""),
        ("Total bond components (σ + π)", sigma + pi, ""),
        ("Fraction Csp3", round(frac_csp3, 3), ""),
        ("Molar refractivity", round(Crippen.MolMR(mol), 3), "")
    ]
    return mol, formula, props, sigma, pi, bond_note

def interpretation(name, value, unit):
    key = name.lower()
    if key == "molecular formula":
        return "Shows the elemental composition represented by this molecular structure."
    if key == "molecular weight":
        return "Average molecular mass calculated using standard atomic weights."
    if key == "exact molecular weight":
        return "Monoisotopic mass calculated from exact isotopic masses."
    if "logp" in key:
        return "Estimated partitioning descriptor; higher values generally indicate greater lipophilicity. It is a model descriptor, not a direct solubility measurement."
    if key == "tpsa":
        return "Topological polar surface area; a descriptor related to polar functionality and molecular interactions."
    if "donors" in key:
        return "Counts groups able to donate hydrogen bonds under the descriptor's structural rules."
    if "acceptors" in key:
        return "Counts atoms/groups able to accept hydrogen bonds under the descriptor's structural rules."
    if "rotatable" in key:
        return "Estimates conformational flexibility from rotatable bonds under RDKit's definition."
    if "aromatic rings" in key:
        return "Counts aromatic rings identified by RDKit."
    if "aliphatic rings" in key:
        return "Counts non-aromatic rings identified by RDKit."
    if key == "ring count":
        return "Counts rings using RDKit's ring-perception conventions."
    if "heavy atom bonds" in key:
        return "Counts bonds connecting non-hydrogen atoms."
    if "sigma" in key:
        return "Bond-component count under the convention documented below; hydrogens are included when added for counting."
    if "pi" in key:
        return "Additional bond components associated with multiple bonding in the Kekulized representation."
    if "total bond" in key:
        return "Sum of the calculated sigma and pi bond components; not the number of unique atom-pair connections."
    if "fraction csp3" in key:
        return "Fraction of carbon atoms assigned sp³ hybridization by RDKit; zero when the molecule has no carbon atoms."
    if "molar refractivity" in key:
        return "Calculated descriptor related to molecular volume and electronic polarizability."
    if "number of atoms" in key:
        return "Counts atoms explicitly present in the input molecular graph; implicit hydrogens are excluded."
    if "heavy atoms" in key:
        return "Counts all non-hydrogen atoms in the molecular graph."
    return "A calculated structural descriptor that helps characterize the molecule."

def draw_2d(mol):
    m = Chem.Mol(mol)
    try:
        Chem.rdDepictor.Compute2DCoords(m)
    except Exception:
        pass
    return Draw.MolToImage(m, size=(800, 500), kekulize=True)

def generate_3d_html(mol):
    """Create HTML for py3Dmol viewer. Returns None if conformer generation fails."""
    m = Chem.AddHs(Chem.Mol(mol))
    params = AllChem.ETKDGv3()
    params.randomSeed = 2026
    status = AllChem.EmbedMolecule(m, params)
    if status != 0:
        return None, "Could not generate a 3D conformer for this structure."
    try:
        if AllChem.MMFFHasAllMoleculeParams(m):
            AllChem.MMFFOptimizeMolecule(m, maxIters=200)
        else:
            AllChem.UFFOptimizeMolecule(m, maxIters=200)
    except Exception:
        pass
    conf = m.GetConformer()
    lines = []
    for atom in m.GetAtoms():
        p = conf.GetAtomPosition(atom.GetIdx())
        symbol = atom.GetSymbol()
        lines.append(f"{symbol} {p.x:.5f} {p.y:.5f} {p.z:.5f}")
    xyz = f"{m.GetNumAtoms()}\nMolecular Studio conformer\n" + "\n".join(lines)
    html = f"""
    <html><head><script src="https://3Dmol.org/build/3Dmol-min.js"></script></head>
    <body style="margin:0;background:#101827;">
    <div id="viewer" style="width:100%;height:430px;position:relative;"></div>
    <script>
    const element = document.getElementById('viewer');
    const viewer = $3Dmol.createViewer(element, {{backgroundColor:'#101827'}});
    const xyz = {xyz!r};
    viewer.addModel(xyz, 'xyz');
    viewer.setStyle({{}}, {{stick:{{radius:0.16}}, sphere:{{scale:0.23}}}});
    viewer.zoomTo();
    viewer.render();
    viewer.zoom(0.9);
    </script></body></html>
    """
    return html, None

def make_pdf(mol, smiles, canonical, formula, props, bond_note, student_name, reg_no, image_bytes):
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=A4, rightMargin=16*mm, leftMargin=16*mm, topMargin=15*mm, bottomMargin=15*mm)
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="SmallText", parent=styles["BodyText"], fontSize=8, leading=10))
    styles.add(ParagraphStyle(name="SectionTitle", parent=styles["Heading2"], textColor=colors.HexColor("#087e8b"), spaceBefore=10, spaceAfter=6))
    story = []
    story.append(Paragraph("CHEMINFORMATICS VIRTUAL LABORATORY", styles["Title"]))
    story.append(Paragraph("Molecular Analysis Report", styles["Heading2"]))
    story.append(Spacer(1, 6))
    story.append(Paragraph(f"<b>Student name:</b> {student_name or 'Not provided'}", styles["BodyText"]))
    story.append(Paragraph(f"<b>Registration number:</b> {reg_no or 'Not provided'}", styles["BodyText"]))
    story.append(Paragraph(f"<b>Input SMILES:</b> {smiles}", styles["BodyText"]))
    story.append(Paragraph(f"<b>Canonical SMILES:</b> {canonical}", styles["BodyText"]))
    story.append(Paragraph(f"<b>Molecular formula:</b> {formula}", styles["BodyText"]))
    story.append(Spacer(1, 10))
    story.append(Paragraph("1. Molecular Structure", styles["SectionTitle"]))
    story.append(Image(io.BytesIO(image_bytes), width=160*mm, height=100*mm))
    story.append(Paragraph("2. Molecular Properties and Interpretation", styles["SectionTitle"]))
    table_data = [["Property", "Value", "Unit", "Interpretation"]]
    for name, value, unit in props:
        interp = interpretation(name, value, unit)
        table_data.append([
            Paragraph(str(name), styles["SmallText"]),
            Paragraph(str(value), styles["SmallText"]),
            Paragraph(str(unit), styles["SmallText"]),
            Paragraph(interp, styles["SmallText"])
        ])
    table = Table(table_data, colWidths=[36*mm, 25*mm, 25*mm, 88*mm], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0,0), (-1,0), colors.HexColor("#123b4b")),
        ("TEXTCOLOR", (0,0), (-1,0), colors.white),
        ("GRID", (0,0), (-1,-1), 0.35, colors.HexColor("#b7c8d6")),
        ("VALIGN", (0,0), (-1,-1), "TOP"),
        ("LEFTPADDING", (0,0), (-1,-1), 4),
        ("RIGHTPADDING", (0,0), (-1,-1), 4),
        ("TOPPADDING", (0,0), (-1,-1), 4),
        ("BOTTOMPADDING", (0,0), (-1,-1), 4),
        ("ROWBACKGROUNDS", (0,1), (-1,-1), [colors.white, colors.HexColor("#f2f7fa")]),
    ]))
    story.append(table)
    story.append(Spacer(1, 8))
    story.append(Paragraph("Bond-counting convention", styles["SectionTitle"]))
    story.append(Paragraph(bond_note, styles["SmallText"]))
    story.append(Paragraph("3. Point-wise Interpretation", styles["SectionTitle"]))
    for name, value, unit in props:
        story.append(Paragraph(f"<b>{name} ({value}{(' ' + unit) if unit else ''}):</b> {interpretation(name, value, unit)}", styles["SmallText"]))
        story.append(Spacer(1, 3))
    story.append(Paragraph("4. Conclusion and Learning Note", styles["SectionTitle"]))
    story.append(Paragraph(
        "The calculated descriptors should be interpreted together with the 2D connectivity and generated 3D conformer to understand composition, bonding, polarity, flexibility, and structural character. These values are computational estimates based on the submitted molecular representation; experimental measurements may differ with conditions and methodology.",
        styles["BodyText"]))
    doc.build(story)
    buffer.seek(0)
    return buffer.getvalue()

with st.sidebar:
    st.header("Molecule input")
    selected_example = st.selectbox("Choose an example molecule", ["Custom SMILES"] + list(EXAMPLES.keys()))
    default_smiles = EXAMPLES[selected_example][0] if selected_example in EXAMPLES else "CCO"
    smiles = st.text_input("SMILES notation", value=default_smiles, help="Example: CCO for ethanol; c1ccccc1 for benzene.")
    student_name = st.text_input("Student name (optional)")
    reg_no = st.text_input("Registration number (optional)")
    analyze = st.button("Analyze molecule", type="primary", use_container_width=True)
    st.caption("Choose an example or enter a valid SMILES string.")

if "analysis_result" not in st.session_state:
    st.session_state.analysis_result = None
if analyze:
    try:
        mol, formula, props, sigma, pi, bond_note = analyze_smiles(smiles.strip())
        canonical = Chem.MolToSmiles(mol, canonical=True)
        image = draw_2d(mol)
        img_buf = io.BytesIO()
        image.save(img_buf, format="PNG")
        img_bytes = img_buf.getvalue()
        st.session_state.analysis_result = {
            "mol": mol, "formula": formula, "props": props, "canonical": canonical,
            "smiles": smiles.strip(), "sigma": sigma, "pi": pi, "bond_note": bond_note,
            "image_bytes": img_bytes, "student_name": student_name, "reg_no": reg_no
        }
    except Exception as exc:
        st.session_state.analysis_result = None
        st.error(f"Could not analyze this SMILES: {exc}")

result = st.session_state.analysis_result
if result is None:
    st.subheader("Explore example molecules")
    st.write("Select a molecule from the sidebar, then click **Analyze molecule**.")
    cols = st.columns(3)
    for i, (name, (smi, desc)) in enumerate(EXAMPLES.items()):
        with cols[i % 3]:
            st.markdown(f"### {name}")
            st.code(smi, language=None)
            st.caption(desc)
    st.stop()

mol = result["mol"]
props = result["props"]
prop_map = {name: value for name, value, unit in props}
st.success(f"Analysis complete · {result['formula']} · Canonical SMILES: {result['canonical']}")

m1, m2, m3, m4, m5 = st.columns(5)
m1.metric("Formula", result["formula"])
m2.metric("Molecular weight", f"{prop_map['Molecular weight']} g/mol")
m3.metric("LogP", str(prop_map["LogP (Crippen)"]))
m4.metric("TPSA", f"{prop_map['TPSA']} Å²")
m5.metric("Heavy atoms", str(prop_map["Heavy atoms"]))

tab1, tab2, tab3, tab4 = st.tabs(["Structures", "Properties", "Interpretation", "Export report"])

with tab1:
    left, right = st.columns(2)
    with left:
        st.subheader("2D molecular structure")
        st.image(result["image_bytes"], use_container_width=True)
        st.download_button("Download 2D structure (PNG)", data=result["image_bytes"],
                           file_name="molecular_structure.png", mime="image/png")
    with right:
        st.subheader("Interactive 3D ball-and-stick model")
        html, error = generate_3d_html(mol)
        if error:
            st.warning(error)
        else:
            components.html(html, height=450, scrolling=False)
        st.caption("The 3D model is a generated conformer, not an experimentally determined structure.")

with tab2:
    st.subheader("Complete molecular properties")
    rows = []
    for name, value, unit in props:
        rows.append({"Property": name, "Calculated value": value, "Unit / convention": unit,
                     "Interpretation": interpretation(name, value, unit)})
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
    st.caption(result["bond_note"])
    st.download_button("Download properties as CSV", data=pd.DataFrame(rows).to_csv(index=False).encode("utf-8"),
                       file_name="molecular_properties.csv", mime="text/csv")

with tab3:
    st.subheader("Point-wise molecular interpretation")
    for name, value, unit in props:
        with st.expander(f"{name}: {value}{(' ' + unit) if unit else ''}", expanded=name in ["Molecular formula", "LogP (Crippen)", "TPSA"]):
            st.write(interpretation(name, value, unit))
    st.subheader("Overall conclusion")
    st.write(f"The submitted structure has formula **{result['formula']}** and average molecular weight **{prop_map['Molecular weight']} g/mol**. Interpret its descriptors together with the 2D connectivity and 3D conformer. Descriptor values are computational estimates, not experimental measurements.")
    st.warning("Sigma/pi counts use the displayed convention and may depend on aromatic kekulization. They should be checked against the convention required by your course or report.")

with tab4:
    st.subheader("Generate a molecular analysis PDF")
    st.write("The report is generated from the current analysis. Student fields are optional.")
    if st.button("Prepare PDF report"):
        try:
            pdf_bytes = make_pdf(mol, result["smiles"], result["canonical"], result["formula"], props,
                                 result["bond_note"], result["student_name"], result["reg_no"], result["image_bytes"])
            st.session_state.pdf_bytes = pdf_bytes
        except Exception as exc:
            st.error(f"PDF generation failed: {exc}")
    if st.session_state.get("pdf_bytes"):
        st.download_button("Download Molecular Analysis Report (PDF)",
                           data=st.session_state.pdf_bytes,
                           file_name="molecular_analysis_report.pdf",
                           mime="application/pdf")
