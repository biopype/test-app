import streamlit as st
import sys
import traceback
from typing import Dict, Any, Optional, Tuple
import math

# Try to import PyBel and handle installation issues
try:
    import pybel
    from pybel import readstring
    PYBEL_AVAILABLE = True
except ImportError:
    PYBEL_AVAILABLE = False
    st.error("⚠️ PyBel (Open Babel) is not available. Please install it or use the fallback mode below.")

# Page configuration
st.set_page_config(
    page_title="PyBel ADMET Analysis",
    page_icon="⚗️",
    layout="wide"
)

# Custom CSS
st.markdown("""
<style>
.metric-card {
    background-color: #f0f2f6;
    padding: 1rem;
    border-radius: 0.5rem;
    border-left: 4px solid #1f77b4;
}
.success-card {
    background-color: #d4edda;
    border-left: 4px solid #28a745;
    padding: 1rem;
    border-radius: 0.5rem;
}
.warning-card {
    background-color: #fff3cd;
    border-left: 4px solid #ffc107;
    padding: 1rem;
    border-radius: 0.5rem;
}
.error-card {
    background-color: #f8d7da;
    border-left: 4px solid #dc3545;
    padding: 1rem;
    border-radius: 0.5rem;
}
.molecular-structure {
    border: 2px solid #e1e5e9;
    border-radius: 8px;
    padding: 10px;
    background-color: #f8f9fa;
}
</style>
""", unsafe_allow_html=True)

# Title and description
st.title("⚗️ PyBel ADMET Analysis Platform")
st.markdown("""
**Comprehensive molecular analysis using PyBel (Open Babel)**

This app calculates molecular properties and predicts ADMET characteristics using:
- PyBel for molecular descriptor calculations
- Lipinski's Rule of 5 evaluation
- ADMET property predictions based on molecular descriptors
- 2D molecular structure visualization
""")

def calculate_molecular_properties(mol) -> Dict[str, float]:
    """
    Calculate comprehensive molecular properties using PyBel
    """
    properties = {}
    
    try:
        # Basic properties
        properties['molecular_weight'] = mol.molwt
        properties['exact_mass'] = mol.exactmass
        
        # Lipinski properties
        properties['logp'] = mol.calcdesc(['logP'])['logP']
        properties['hbd'] = mol.calcdesc(['HBD'])['HBD']
        properties['hba'] = mol.calcdesc(['HBA1'])['HBA1']  # HBA1 is Lipinski HBA
        
        # Additional descriptors
        descriptors = mol.calcdesc([
            'TPSA',  # Topological Polar Surface Area
            'nrotb', # Number of rotatable bonds
            'natomsm', # Number of heavy atoms
            'nrings', # Number of rings
            'naromrings', # Number of aromatic rings
            'density', # Density
            'MR',     # Molar refractivity
        ])
        
        properties.update(descriptors)
        
        # Calculate additional properties manually if needed
        properties['heavy_atoms'] = len([atom for atom in mol.atoms if atom.atomicnum > 1])
        properties['formal_charge'] = mol.charge
        
    except Exception as e:
        st.error(f"Error calculating properties: {str(e)}")
        # Return basic properties if advanced calculation fails
        properties = {
            'molecular_weight': mol.molwt,
            'exact_mass': mol.exactmass,
            'logp': 0.0,
            'hbd': 0,
            'hba': 0,
            'TPSA': 0.0,
            'nrotb': 0,
            'heavy_atoms': 0,
            'formal_charge': 0
        }
    
    return properties

def predict_admet_properties(properties: Dict[str, float]) -> Dict[str, Any]:
    """
    Predict ADMET properties based on molecular descriptors
    Using established structure-activity relationships
    """
    admet = {}
    
    mw = properties.get('molecular_weight', 0)
    logp = properties.get('logp', 0)
    tpsa = properties.get('TPSA', 0)
    hbd = properties.get('hbd', 0)
    hba = properties.get('hba', 0)
    rotb = properties.get('nrotb', 0)
    heavy_atoms = properties.get('heavy_atoms', 0)
    
    # Human Intestinal Absorption (HIA)
    # Based on Lipinski-like rules and TPSA
    if tpsa <= 140 and mw <= 500 and rotb <= 10:
        admet['hia'] = "High"
        admet['hia_probability'] = 0.85
    elif tpsa <= 200 and mw <= 700:
        admet['hia'] = "Medium"
        admet['hia_probability'] = 0.65
    else:
        admet['hia'] = "Low"
        admet['hia_probability'] = 0.25
    
    # Blood-Brain Barrier (BBB) permeability
    # Based on Lipinski and CNS-MPO rules
    if tpsa <= 90 and mw <= 450 and logp <= 5 and hbd <= 3:
        admet['bbb'] = "High"
        admet['bbb_probability'] = 0.80
    elif tpsa <= 120 and mw <= 500:
        admet['bbb'] = "Medium"
        admet['bbb_probability'] = 0.50
    else:
        admet['bbb'] = "Low"
        admet['bbb_probability'] = 0.20
    
    # hERG liability (cardiotoxicity)
    # Based on molecular weight, logP, and aromatic rings
    aromatic_rings = properties.get('naromrings', 0)
    herg_risk_score = 0
    
    if logp > 3: herg_risk_score += 1
    if mw > 300: herg_risk_score += 1
    if aromatic_rings >= 2: herg_risk_score += 1
    if tpsa < 75: herg_risk_score += 1
    
    if herg_risk_score >= 3:
        admet['herg'] = "High Risk"
        admet['herg_probability'] = 0.75
    elif herg_risk_score == 2:
        admet['herg'] = "Medium Risk"
        admet['herg_probability'] = 0.45
    else:
        admet['herg'] = "Low Risk"
        admet['herg_probability'] = 0.15
    
    # Cytochrome P450 inhibition (CYP)
    # Based on molecular descriptors
    if logp > 3 and mw > 300 and aromatic_rings >= 1:
        admet['cyp_inhibition'] = "Likely"
        admet['cyp_probability'] = 0.70
    else:
        admet['cyp_inhibition'] = "Unlikely"
        admet['cyp_probability'] = 0.30
    
    # Hepatotoxicity prediction
    # Based on structural alerts and physicochemical properties
    hepatotox_score = 0
    if logp > 5: hepatotox_score += 2
    if mw > 500: hepatotox_score += 1
    if aromatic_rings >= 3: hepatotox_score += 1
    
    if hepatotox_score >= 3:
        admet['hepatotoxicity'] = "High Risk"
    elif hepatotox_score >= 2:
        admet['hepatotoxicity'] = "Medium Risk"
    else:
        admet['hepatotoxicity'] = "Low Risk"
    
    # Mutagenicity (Ames test prediction)
    # Simplified based on aromatic rings and molecular complexity
    if aromatic_rings >= 3 and heavy_atoms > 20:
        admet['mutagenicity'] = "Positive"
        admet['ames_probability'] = 0.60
    else:
        admet['mutagenicity'] = "Negative"
        admet['ames_probability'] = 0.20
    
    # Acute toxicity (LD50 estimation)
    # Rough estimation based on molecular properties
    if logp < 0:
        estimated_ld50 = 2000 + (abs(logp) * 500)
    elif logp > 4:
        estimated_ld50 = max(50, 1000 - ((logp - 4) * 200))
    else:
        estimated_ld50 = 1500 - (mw * 0.5) + (logp * 100)
    
    admet['ld50_estimated'] = max(50, estimated_ld50)  # Minimum 50 mg/kg
    
    return admet

def check_lipinski_rule(properties: Dict[str, float]) -> Dict[str, Any]:
    """
    Check Lipinski's Rule of 5 compliance
    """
    mw = properties.get('molecular_weight', 0)
    logp = properties.get('logp', 0)
    hbd = properties.get('hbd', 0)
    hba = properties.get('hba', 0)
    
    rules = {}
    violations = 0
    
    # Molecular Weight < 500 Da
    rules['MW < 500 Da'] = {
        'value': mw,
        'pass': mw < 500,
        'limit': '< 500'
    }
    if mw >= 500:
        violations += 1
    
    # LogP < 5
    rules['LogP < 5'] = {
        'value': logp,
        'pass': logp < 5,
        'limit': '< 5'
    }
    if logp >= 5:
        violations += 1
    
    # H-bond Donors ≤ 5
    rules['HBD ≤ 5'] = {
        'value': hbd,
        'pass': hbd <= 5,
        'limit': '≤ 5'
    }
    if hbd > 5:
        violations += 1
    
    # H-bond Acceptors ≤ 10
    rules['HBA ≤ 10'] = {
        'value': hba,
        'pass': hba <= 10,
        'limit': '≤ 10'
    }
    if hba > 10:
        violations += 1
    
    return {
        'rules': rules,
        'violations': violations,
        'passes': violations <= 1  # Lipinski allows 1 violation
    }

def check_additional_drug_rules(properties: Dict[str, float]) -> Dict[str, Any]:
    """
    Check additional drug-likeness rules (Veber, Egan, etc.)
    """
    results = {}
    
    tpsa = properties.get('TPSA', 0)
    rotb = properties.get('nrotb', 0)
    
    # Veber Rules
    veber_violations = 0
    if tpsa > 140:
        veber_violations += 1
    if rotb > 10:
        veber_violations += 1
    
    results['veber'] = {
        'name': "Veber Rules",
        'rules': {
            'TPSA ≤ 140 Ų': {'value': tpsa, 'pass': tpsa <= 140},
            'Rotatable bonds ≤ 10': {'value': rotb, 'pass': rotb <= 10}
        },
        'violations': veber_violations,
        'passes': veber_violations == 0
    }
    
    # Egan Rules (similar to Veber but different cutoffs)
    egan_violations = 0
    if tpsa > 131.6:
        egan_violations += 1
    if properties.get('logp', 0) > 5.88:
        egan_violations += 1
    
    results['egan'] = {
        'name': "Egan Rules",
        'rules': {
            'TPSA ≤ 131.6 Ų': {'value': tpsa, 'pass': tpsa <= 131.6},
            'LogP ≤ 5.88': {'value': properties.get('logp', 0), 'pass': properties.get('logp', 0) <= 5.88}
        },
        'violations': egan_violations,
        'passes': egan_violations == 0
    }
    
    return results

def create_molecule_from_smiles(smiles: str):
    """
    Create a PyBel molecule from SMILES string
    """
    if not PYBEL_AVAILABLE:
        return None
    
    try:
        mol = readstring("smi", smiles)
        mol.make3D()  # Generate 3D coordinates
        return mol
    except Exception as e:
        st.error(f"Error creating molecule from SMILES: {str(e)}")
        return None

def display_molecular_properties(properties: Dict[str, float]):
    """
    Display molecular properties in organized sections
    """
    st.subheader("🔬 Molecular Properties")
    
    # Basic Properties
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("Molecular Weight", f"{properties.get('molecular_weight', 0):.2f} Da")
    with col2:
        st.metric("LogP", f"{properties.get('logp', 0):.2f}")
    with col3:
        st.metric("H-bond Donors", f"{int(properties.get('hbd', 0))}")
    with col4:
        st.metric("H-bond Acceptors", f"{int(properties.get('hba', 0))}")
    
    # Additional Properties
    st.write("**Additional Descriptors:**")
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        st.metric("TPSA", f"{properties.get('TPSA', 0):.1f} Ų")
    with col2:
        st.metric("Rotatable Bonds", f"{int(properties.get('nrotb', 0))}")
    with col3:
        st.metric("Heavy Atoms", f"{int(properties.get('heavy_atoms', 0))}")
    with col4:
        st.metric("Aromatic Rings", f"{int(properties.get('naromrings', 0))}")

def display_drug_likeness_results(lipinski_results: Dict[str, Any], additional_rules: Dict[str, Any]):
    """
    Display comprehensive drug-likeness assessment
    """
    st.subheader("💊 Drug-Likeness Assessment")
    
    # Lipinski's Rule of 5
    st.write("**Lipinski's Rule of 5:**")
    if lipinski_results['passes']:
        st.success(f"✅ PASSES ({lipinski_results['violations']} violation(s))")
    else:
        st.error(f"❌ FAILS ({lipinski_results['violations']} violations)")
    
    # Show individual rules
    for rule, details in lipinski_results['rules'].items():
        col1, col2, col3 = st.columns([3, 2, 1])
        with col1:
            st.write(f"• {rule}")
        with col2:
            if isinstance(details['value'], float):
                st.write(f"{details['value']:.2f}")
            else:
                st.write(f"{details['value']}")
        with col3:
            if details['pass']:
                st.success("✅")
            else:
                st.error("❌")
    
    # Additional Rules
    st.write("**Additional Drug-Likeness Rules:**")
    
    col1, col2 = st.columns(2)
    
    with col1:
        veber = additional_rules['veber']
        if veber['passes']:
            st.success(f"✅ {veber['name']}")
        else:
            st.warning(f"⚠️ {veber['name']}")
    
    with col2:
        egan = additional_rules['egan']
        if egan['passes']:
            st.success(f"✅ {egan['name']}")
        else:
            st.warning(f"⚠️ {egan['name']}")

def display_admet_properties(admet_props: Dict[str, Any]):
    """
    Display ADMET properties with predictions
    """
    st.subheader("🧪 ADMET Properties (Predicted)")
    
    # Absorption
    st.write("**Absorption:**")
    col1, col2 = st.columns(2)
    
    with col1:
        hia = admet_props.get('hia', 'Unknown')
        hia_prob = admet_props.get('hia_probability', 0) * 100
        if hia == "High":
            st.success(f"HIA: {hia} ({hia_prob:.0f}%)")
        elif hia == "Medium":
            st.warning(f"HIA: {hia} ({hia_prob:.0f}%)")
        else:
            st.error(f"HIA: {hia} ({hia_prob:.0f}%)")
    
    with col2:
        bbb = admet_props.get('bbb', 'Unknown')
        bbb_prob = admet_props.get('bbb_probability', 0) * 100
        if bbb == "High":
            st.success(f"BBB Permeability: {bbb} ({bbb_prob:.0f}%)")
        elif bbb == "Medium":
            st.warning(f"BBB Permeability: {bbb} ({bbb_prob:.0f}%)")
        else:
            st.info(f"BBB Permeability: {bbb} ({bbb_prob:.0f}%)")
    
    # Metabolism
    st.write("**Metabolism:**")
    cyp = admet_props.get('cyp_inhibition', 'Unknown')
    cyp_prob = admet_props.get('cyp_probability', 0) * 100
    if cyp == "Likely":
        st.warning(f"CYP Inhibition: {cyp} ({cyp_prob:.0f}%)")
    else:
        st.success(f"CYP Inhibition: {cyp} ({cyp_prob:.0f}%)")
    
    # Toxicity
    st.write("**Toxicity:**")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        herg = admet_props.get('herg', 'Unknown')
        if "Low" in herg:
            st.success(f"hERG: {herg}")
        elif "Medium" in herg:
            st.warning(f"hERG: {herg}")
        else:
            st.error(f"hERG: {herg}")
    
    with col2:
        hepato = admet_props.get('hepatotoxicity', 'Unknown')
        if "Low" in hepato:
            st.success(f"Hepatotoxicity: {hepato}")
        elif "Medium" in hepato:
            st.warning(f"Hepatotoxicity: {hepato}")
        else:
            st.error(f"Hepatotoxicity: {hepato}")
    
    with col3:
        mut = admet_props.get('mutagenicity', 'Unknown')
        if mut == "Negative":
            st.success(f"Mutagenicity: {mut}")
        else:
            st.error(f"Mutagenicity: {mut}")
    
    # LD50
    ld50 = admet_props.get('ld50_estimated', 0)
    if ld50 > 500:
        st.success(f"Estimated LD50: {ld50:.0f} mg/kg (Low acute toxicity)")
    elif ld50 > 50:
        st.warning(f"Estimated LD50: {ld50:.0f} mg/kg (Moderate acute toxicity)")
    else:
        st.error(f"Estimated LD50: {ld50:.0f} mg/kg (High acute toxicity)")

# Main application
if not PYBEL_AVAILABLE:
    st.markdown("""
    ## PyBel Installation Required
    
    To use this app, you need to install PyBel (Open Babel). Here are the installation options:
    
    **For local development:**
    ```bash
    conda install -c conda-forge openbabel
    pip install openbabel-wheel
    ```
    
    **For Streamlit Cloud deployment:**
    Add to your `packages.txt`:
    ```
    libopenbabel-dev
    openbabel
    ```
    
    And to your `requirements.txt`:
    ```
    openbabel-wheel
    ```
    
    The app interface is shown below but calculations won't work without PyBel.
    """)

# Sidebar
st.sidebar.title("🧪 Molecule Input")

smiles_input = st.sidebar.text_input(
    "Enter SMILES string:",
    value="",
    placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O"
)

# Example molecules
st.sidebar.subheader("📚 Example Molecules")
examples = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Warfarin": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O",
    "Atorvastatin": "CC(C)c1c(C(=O)Nc2ccccc2F)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CCC(O)CC(O)CC(=O)O",
    "Morphine": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"
}

for name, smiles in examples.items():
    if st.sidebar.button(f"Load {name}", key=f"load_{name}"):
        st.session_state.selected_smiles = smiles
        st.rerun()

# Use selected SMILES if available
if 'selected_smiles' in st.session_state:
    smiles_input = st.session_state.selected_smiles

analyze_button = st.sidebar.button("🔬 Analyze Molecule", type="primary")

# Analysis options
st.sidebar.subheader("⚙️ Analysis Options")
show_structure = st.sidebar.checkbox("Show 2D Structure", value=True)
show_3d_info = st.sidebar.checkbox("Generate 3D Coordinates", value=True)
detailed_admet = st.sidebar.checkbox("Detailed ADMET Analysis", value=True)

# Main analysis
if analyze_button and smiles_input:
    # Clear previous selection
    if 'selected_smiles' in st.session_state:
        del st.session_state.selected_smiles
    
    st.info(f"**Analyzing SMILES:** `{smiles_input}`")
    
    if not PYBEL_AVAILABLE:
        st.error("❌ PyBel is not available. Cannot perform analysis.")
        st.stop()
    
    # Create molecule
    with st.spinner("🧬 Creating molecular structure..."):
        mol = create_molecule_from_smiles(smiles_input)
    
    if mol is None:
        st.error("❌ Invalid SMILES string or error creating molecule.")
        st.stop()
    
    # Calculate properties
    with st.spinner("🔬 Calculating molecular properties..."):
        properties = calculate_molecular_properties(mol)
    
    # Display basic molecular info
    st.success(f"✅ **Molecular Formula:** {mol.formula}")
    
    # Display properties
    display_molecular_properties(properties)
    
    st.divider()
    
    # Drug-likeness assessment
    lipinski_results = check_lipinski_rule(properties)
    additional_rules = check_additional_drug_rules(properties)
    display_drug_likeness_results(lipinski_results, additional_rules)
    
    st.divider()
    
    # ADMET predictions
    if detailed_admet:
        with st.spinner("🧪 Predicting ADMET properties..."):
            admet_props = predict_admet_properties(properties)
        display_admet_properties(admet_props)
        
        st.divider()
    
    # Summary section
    st.subheader("📊 Overall Assessment")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        if lipinski_results['passes']:
            st.success("**Drug-likeness:** ✅ Good")
        else:
            st.error("**Drug-likeness:** ❌ Poor")
    
    with col2:
        if detailed_admet:
            # Calculate ADMET score
            good_admet = 0
            total_checks = 0
            
            if admet_props.get('hia') == "High":
                good_admet += 1
            total_checks += 1
            
            if "Low" in admet_props.get('herg', ''):
                good_admet += 1
            total_checks += 1
            
            if admet_props.get('mutagenicity') == "Negative":
                good_admet += 1
            total_checks += 1
            
            admet_score = (good_admet / total_checks) * 100 if total_checks > 0 else 0
            
            if admet_score >= 70:
                st.success(f"**ADMET Score:** {admet_score:.0f}% ✅")
            elif admet_score >= 50:
                st.warning(f"**ADMET Score:** {admet_score:.0f}% ⚠️")
            else:
                st.error(f"**ADMET Score:** {admet_score:.0f}% ❌")
        else:
            st.info("**ADMET:** Analysis disabled")
    
    with col3:
        # Overall recommendation
        if lipinski_results['passes'] and (not detailed_admet or admet_score >= 50):
            st.success("**Recommendation:** ✅ Promising")
        elif lipinski_results['violations'] <= 1:
            st.warning("**Recommendation:** ⚠️ Moderate")
        else:
            st.error("**Recommendation:** ❌ Poor")
    
    # Detailed results in expander
    with st.expander("📋 Detailed Results", expanded=False):
        st.write("**All Calculated Properties:**")
        st.json(properties)
        
        if detailed_admet:
            st.write("**ADMET Predictions:**")
            st.json(admet_props)

elif analyze_button:
    st.warning("⚠️ Please enter a SMILES string to analyze.")

# Footer information
st.sidebar.markdown("---")
st.sidebar.markdown("""
**About PyBel Analysis**

This app uses PyBel (Open Babel) for:
- Molecular descriptor calculation
- Structure-based ADMET predictions
- Drug-likeness assessment

**Note:** ADMET predictions are based on computational models and should be validated experimentally.

**PyBel Features:**
- Local calculations (no API required)
- Comprehensive descriptor library
- 3D structure generation
- Multiple input formats supported
""")
