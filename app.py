import streamlit as st
import requests
import json
import time
from typing import Dict, Any, Optional

# Page configuration
st.set_page_config(
    page_title="ADMETlab 2.0 Drug Analysis",
    page_icon="🧬",
    layout="wide"
)

# Custom CSS for better styling
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
</style>
""", unsafe_allow_html=True)

# Title and description
st.title("🧬 ADMETlab 2.0 Drug Analysis Platform")
st.markdown("""
**Analyze molecular properties and ADMET characteristics using ADMETlab 2.0 API**

Enter a SMILES string to get comprehensive drug-likeness and ADMET predictions including:
- Molecular properties (MW, LogP, H-bond donors/acceptors)
- Lipinski's Rule of 5 compliance
- ADMET properties (Absorption, Distribution, Metabolism, Excretion, Toxicity)
""")

def make_api_request(smiles: str) -> Optional[Dict[str, Any]]:
    """
    Send SMILES to ADMETlab 2.0 API and return the response
    """
    url = "https://admetmesh.scbdd.com/service/predict"
    
    headers = {
        'Content-Type': 'application/json',
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    payload = {
        "smiles": smiles
    }
    
    try:
        with st.spinner("🔬 Analyzing molecular properties..."):
            response = requests.post(url, json=payload, headers=headers, timeout=30)
            response.raise_for_status()
            return response.json()
    except requests.exceptions.RequestException as e:
        st.error(f"API request failed: {str(e)}")
        return None
    except json.JSONDecodeError as e:
        st.error(f"Failed to parse API response: {str(e)}")
        return None

def extract_molecular_properties(data: Dict[str, Any]) -> Dict[str, float]:
    """
    Extract molecular weight, logP, H-bond donors and acceptors from API response
    """
    properties = {}
    
    # Common property mappings (adjust based on actual API response structure)
    property_mappings = {
        'molecular_weight': ['MW', 'molecular_weight', 'mol_weight', 'mw'],
        'logp': ['LogP', 'logp', 'clogp', 'log_p'],
        'hbd': ['HBD', 'hbd', 'h_bond_donors', 'num_hbd'],
        'hba': ['HBA', 'hba', 'h_bond_acceptors', 'num_hba']
    }
    
    # Search through the data structure
    def find_property(data, possible_keys):
        if isinstance(data, dict):
            for key in possible_keys:
                if key in data:
                    return data[key]
            # Search recursively
            for value in data.values():
                result = find_property(value, possible_keys)
                if result is not None:
                    return result
        return None
    
    for prop, keys in property_mappings.items():
        value = find_property(data, keys)
        if value is not None:
            try:
                properties[prop] = float(value)
            except (ValueError, TypeError):
                properties[prop] = None
    
    return properties

def extract_admet_properties(data: Dict[str, Any]) -> Dict[str, Any]:
    """
    Extract ADMET properties from API response
    """
    admet_props = {}
    
    # ADMET property mappings
    admet_mappings = {
        'hia': ['HIA', 'hia', 'human_intestinal_absorption'],
        'herg': ['hERG', 'herg', 'herg_blocker', 'cardiotoxicity'],
        'ames': ['Ames', 'ames', 'ames_test', 'mutagenicity'],
        'carcinogenicity': ['Carcinogenicity', 'carcinogenicity', 'carc'],
        'ld50': ['LD50', 'ld50', 'acute_toxicity']
    }
    
    def find_admet_property(data, possible_keys):
        if isinstance(data, dict):
            for key in possible_keys:
                if key in data:
                    return data[key]
            for value in data.values():
                result = find_admet_property(value, possible_keys)
                if result is not None:
                    return result
        return None
    
    for prop, keys in admet_mappings.items():
        value = find_admet_property(data, keys)
        admet_props[prop] = value
    
    return admet_props

def check_lipinski_rule(properties: Dict[str, float]) -> Dict[str, Any]:
    """
    Check Lipinski's Rule of 5 compliance
    """
    mw = properties.get('molecular_weight')
    logp = properties.get('logp')
    hbd = properties.get('hbd')
    hba = properties.get('hba')
    
    rules = {}
    violations = 0
    
    if mw is not None:
        rules['MW < 500 Da'] = {
            'value': mw,
            'pass': mw < 500,
            'limit': '< 500'
        }
        if mw >= 500:
            violations += 1
    
    if logp is not None:
        rules['LogP < 5'] = {
            'value': logp,
            'pass': logp < 5,
            'limit': '< 5'
        }
        if logp >= 5:
            violations += 1
    
    if hbd is not None:
        rules['HBD ≤ 5'] = {
            'value': hbd,
            'pass': hbd <= 5,
            'limit': '≤ 5'
        }
        if hbd > 5:
            violations += 1
    
    if hba is not None:
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

def display_molecular_properties(properties: Dict[str, float]):
    """
    Display molecular properties in a formatted section
    """
    st.subheader("🔬 Molecular Properties")
    
    col1, col2, col3, col4 = st.columns(4)
    
    with col1:
        mw = properties.get('molecular_weight')
        if mw is not None:
            st.metric("Molecular Weight", f"{mw:.2f} Da")
        else:
            st.metric("Molecular Weight", "N/A")
    
    with col2:
        logp = properties.get('logp')
        if logp is not None:
            st.metric("LogP", f"{logp:.2f}")
        else:
            st.metric("LogP", "N/A")
    
    with col3:
        hbd = properties.get('hbd')
        if hbd is not None:
            st.metric("H-bond Donors", f"{int(hbd)}")
        else:
            st.metric("H-bond Donors", "N/A")
    
    with col4:
        hba = properties.get('hba')
        if hba is not None:
            st.metric("H-bond Acceptors", f"{int(hba)}")
        else:
            st.metric("H-bond Acceptors", "N/A")

def display_lipinski_results(lipinski_results: Dict[str, Any]):
    """
    Display Lipinski's Rule of 5 results
    """
    st.subheader("⚖️ Lipinski's Rule of 5")
    
    if not lipinski_results['rules']:
        st.warning("Insufficient data to evaluate Lipinski's Rule of 5")
        return
    
    # Overall result
    if lipinski_results['passes']:
        st.success(f"✅ **PASSES** Lipinski's Rule of 5 ({lipinski_results['violations']} violation(s))")
    else:
        st.error(f"❌ **FAILS** Lipinski's Rule of 5 ({lipinski_results['violations']} violations)")
    
    # Individual rules
    st.write("**Individual Rule Assessment:**")
    
    for rule, details in lipinski_results['rules'].items():
        col1, col2, col3 = st.columns([3, 2, 1])
        
        with col1:
            st.write(f"**{rule}**")
        
        with col2:
            if isinstance(details['value'], float):
                st.write(f"{details['value']:.2f} ({details['limit']})")
            else:
                st.write(f"{details['value']} ({details['limit']})")
        
        with col3:
            if details['pass']:
                st.success("✅ Pass")
            else:
                st.error("❌ Fail")

def display_admet_properties(admet_props: Dict[str, Any]):
    """
    Display ADMET properties in formatted sections
    """
    st.subheader("🧪 ADMET Properties")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**Absorption & Distribution**")
        
        hia = admet_props.get('hia')
        if hia is not None:
            if str(hia).lower() in ['yes', 'true', '1', 'high', 'good']:
                st.success(f"HIA (Human Intestinal Absorption): {hia}")
            else:
                st.warning(f"HIA (Human Intestinal Absorption): {hia}")
        else:
            st.info("HIA: Data not available")
        
        st.write("**Toxicity**")
        
        herg = admet_props.get('herg')
        if herg is not None:
            if str(herg).lower() in ['no', 'false', '0', 'low', 'non-blocker']:
                st.success(f"hERG Blocker: {herg}")
            else:
                st.error(f"hERG Blocker: {herg}")
        else:
            st.info("hERG: Data not available")
    
    with col2:
        st.write("**Mutagenicity & Carcinogenicity**")
        
        ames = admet_props.get('ames')
        if ames is not None:
            if str(ames).lower() in ['negative', 'no', 'false', '0']:
                st.success(f"Ames Test: {ames}")
            else:
                st.error(f"Ames Test: {ames}")
        else:
            st.info("Ames Test: Data not available")
        
        carc = admet_props.get('carcinogenicity')
        if carc is not None:
            if str(carc).lower() in ['no', 'false', '0', 'non-carcinogenic']:
                st.success(f"Carcinogenicity: {carc}")
            else:
                st.error(f"Carcinogenicity: {carc}")
        else:
            st.info("Carcinogenicity: Data not available")
        
        st.write("**Acute Toxicity**")
        ld50 = admet_props.get('ld50')
        if ld50 is not None:
            try:
                ld50_val = float(ld50)
                if ld50_val > 500:
                    st.success(f"LD50: {ld50_val:.2f} mg/kg (Low toxicity)")
                elif ld50_val > 50:
                    st.warning(f"LD50: {ld50_val:.2f} mg/kg (Moderate toxicity)")
                else:
                    st.error(f"LD50: {ld50_val:.2f} mg/kg (High toxicity)")
            except (ValueError, TypeError):
                st.info(f"LD50: {ld50}")
        else:
            st.info("LD50: Data not available")

# Main application interface
st.sidebar.title("Input")
smiles_input = st.sidebar.text_input(
    "Enter SMILES string:",
    value="",
    placeholder="e.g., CC(=O)Oc1ccccc1C(=O)O"
)

# Example molecules
st.sidebar.subheader("Example Molecules")
examples = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Ibuprofen": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Warfarin": "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"
}

for name, smiles in examples.items():
    if st.sidebar.button(f"Load {name}", key=f"load_{name}"):
        st.session_state.selected_smiles = smiles
        st.rerun()

# Use selected SMILES if available
if 'selected_smiles' in st.session_state:
    smiles_input = st.session_state.selected_smiles

analyze_button = st.sidebar.button("🔬 Analyze Molecule", type="primary")

# Main analysis
if analyze_button and smiles_input:
    # Clear previous selection
    if 'selected_smiles' in st.session_state:
        del st.session_state.selected_smiles
    
    st.info(f"**Analyzing SMILES:** `{smiles_input}`")
    
    # Make API request
    api_response = make_api_request(smiles_input)
    
    if api_response:
        # Show raw API response in expandable section
        with st.expander("📋 Raw API Response", expanded=False):
            st.json(api_response)
        
        # Extract and display molecular properties
        molecular_props = extract_molecular_properties(api_response)
        display_molecular_properties(molecular_props)
        
        st.divider()
        
        # Check and display Lipinski's Rule
        lipinski_results = check_lipinski_rule(molecular_props)
        display_lipinski_results(lipinski_results)
        
        st.divider()
        
        # Extract and display ADMET properties
        admet_props = extract_admet_properties(api_response)
        display_admet_properties(admet_props)
        
        # Summary section
        st.subheader("📊 Summary")
        
        col1, col2 = st.columns(2)
        
        with col1:
            if lipinski_results['passes']:
                st.success("**Drug-likeness:** Passes Lipinski's Rule ✅")
            else:
                st.error("**Drug-likeness:** Fails Lipinski's Rule ❌")
        
        with col2:
            # Simple ADMET score based on available data
            good_admet = 0
            total_admet = 0
            
            if admet_props.get('hia') and str(admet_props['hia']).lower() in ['yes', 'true', '1', 'high', 'good']:
                good_admet += 1
            if admet_props.get('hia') is not None:
                total_admet += 1
            
            if admet_props.get('herg') and str(admet_props['herg']).lower() in ['no', 'false', '0', 'low', 'non-blocker']:
                good_admet += 1
            if admet_props.get('herg') is not None:
                total_admet += 1
            
            if admet_props.get('ames') and str(admet_props['ames']).lower() in ['negative', 'no', 'false', '0']:
                good_admet += 1
            if admet_props.get('ames') is not None:
                total_admet += 1
            
            if total_admet > 0:
                admet_score = (good_admet / total_admet) * 100
                if admet_score >= 70:
                    st.success(f"**ADMET Score:** {admet_score:.0f}% ✅")
                elif admet_score >= 50:
                    st.warning(f"**ADMET Score:** {admet_score:.0f}% ⚠️")
                else:
                    st.error(f"**ADMET Score:** {admet_score:.0f}% ❌")
            else:
                st.info("**ADMET Score:** Insufficient data")

elif analyze_button:
    st.warning("Please enter a SMILES string to analyze.")

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown("""
**About ADMETlab 2.0**

This app uses the ADMETlab 2.0 API to predict molecular properties and ADMET characteristics. Results are for research purposes only.

[ADMETlab 2.0 Website](https://admetmesh.scbdd.com/)
""")
