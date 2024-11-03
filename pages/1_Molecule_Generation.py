import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import io
import base64
import json
from models.llm_chain import generate_smiles

st.set_page_config(layout="wide")
# Inject custom CSS to style the sidebar
custom_style = """
    <style>
    /* Sidebar background color */
    [data-testid="stSidebar"] {
        background-color: #1a1f2d;
        color:white;
    }

    /* Main content background color (container) */
    .main,{
        background-color: #eff3f8;
        color: black;  /* Change text color for better contrast */
    }

    /* Container background color */
    [data-testid="stAppViewBlockContainer"]{
        background-color: #fefefe;
        padding: 20px;
        border-radius: 10px;
        color: black; /* Set text color to white for better contrast */
    }
    div[data-testid="stSidebarNav"] li div a span{
        color:white;
        # margin-left: 1rem;
        padding: 0.5rem;
        width: 300px;
        border-radius: 0.5rem;
    }
    div[data-testid="stSidebarNav"] li div::focus-visible {
        background-color: rgba(151, 166, 195, 0.15);
    }
    [data-testid="stSidebarNav"]::before {
                content: "💊 MedGuide AI 🧑🏻‍⚕️";
                margin-left: 20px;
                margin-top: 20px;
                font-size: 30px;
                position: relative;
                top: -50px;
            }
    
    </style>
"""

# Apply the CSS
st.markdown(custom_style, unsafe_allow_html=True)


def call_nvidia_api(payload, test_mode=False):
    invoke_url = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate"
    headers = {
        # "Authorization": "Bearer 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",
        "Authorization": "Bearer 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",
        "Accept": "application/json",
    }
    session = requests.Session()
    
    try:
        if test_mode:
            response = session.get(invoke_url, headers=headers)
        else:
            response = session.post(invoke_url, headers=headers, json=payload)
        
        response.raise_for_status()
        return response.json() if not test_mode else {"status": "API is responsive"}
    except requests.exceptions.HTTPError as e:
        if e.response.status_code == 500:
            raise Exception("NVIDIA API server error. The server might be temporarily unavailable or overloaded.")
        else:
            raise
    except requests.exceptions.RequestException as e:
        raise Exception(f"Error connecting to NVIDIA API: {str(e)}")

def mol_to_img(mol):
    img = Draw.MolToImage(mol)
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    return base64.b64encode(buffered.getvalue()).decode()

# Main app title
st.title("Molecule Generation and Visualization")

# Create a container for inputs
    
# 80% width container
cols = st.columns([4, 2])  # To center the content (left, main, right)
with cols[0]:
    container = st.container(border=True)
    with container:
        st.markdown("<h3 style='text-align: center;'>Input Parameters</h3>", unsafe_allow_html=True)
        # Input parameters inside the container
        smi = st.text_input("SMILES", value="[H][C@@]12Cc3c[nH]c4cccc(C1=C[C@H](NC(=O)N(CC)CC)CN2C)c34")
        col = st.columns([1,1])
        with col[0]:
            algorithm = st.selectbox("Algorithm", ["CMA-ES"])
        with col[1]:
            num_molecules = st.number_input("Number of Molecules", min_value=1, max_value=100, value=20)
        property_name = st.selectbox("Property Name", ["QED","plogP"])
        minimize = st.checkbox("Minimize", value=False)
        min_similarity = st.slider("Minimum Similarity", min_value=0.0, max_value=1.0, value=0.3)
        with col[0]:
            particles = st.number_input("Particles", min_value=1, max_value=100, value=30)
        with col[1]:
            iterations = st.number_input("Iterations", min_value=1, max_value=100, value=10)

with cols[1]:
    chat = st.container(height=560)
    with chat:
        st.markdown("<h3 style='text-align: center;'>MoleculeMirth</h3>", unsafe_allow_html=True)
        # Initialize chat history
        # if "messages" not in st.session_state:
        #     st.session_state.messages = []

        # Display chat messages from history on app rerun
        # for message in st.session_state.messages:
        #     with st.chat_message(message["role"]):
        #         st.markdown(message["content"])

        # React to user input
        st.markdown("------")
        if prompt := st.chat_input("How can I support your research today?"):
            # Display user message in chat message container
            st.chat_message("user").markdown(prompt)
            # Add user message to chat history
            # st.session_state.messages.append({"role": "user", "content": prompt})

            # Generate and add AI response
            with st.spinner("🧠 Thinking..."):
                response = generate_smiles(prompt)
                # print(response)
                # st.write(response)

                if response:
                    assistant_message = {
                        "role":"assistant",
                        "content":response
                    }
                st.chat_message("assistant").markdown(response)
                # st.session_state.messages.append(assistant_message)


# Button to generate molecules
if st.button("Generate Molecules"):
    payload = {
        "algorithm": algorithm,
        "num_molecules": num_molecules,
        "property_name": property_name,
        "minimize": minimize,
        "min_similarity": min_similarity,
        "particles": particles,
        "iterations": iterations,
        "smi": smi
    }

    try:
        response_body = call_nvidia_api(payload)
        
        if 'molecules' in response_body:
            json_string = response_body['molecules']
            data = json.loads(json_string)
            
            st.subheader(f"Generated Molecules ({len(data)})")
            
            img_cols = st.columns(4)  # Create 4 columns for images
            for idx, mol_data in enumerate(data):
                if 'sample' in mol_data:
                    mol = Chem.MolFromSmiles(mol_data['sample'])
                    if mol:
                        AllChem.Compute2DCoords(mol)
                        img = mol_to_img(mol)
                        img_cols[idx % 4].image(f"data:image/png;base64,{img}", caption=f"Score: {round(mol_data['score'], 2)}", use_column_width=True)
                        if 'properties' in mol_data:
                            for prop, value in mol_data['properties'].items():
                                img_cols[idx % 4].write(f"{prop}: {value}")
                    else:
                        img_cols[idx % 4].write(f"Molecule {idx+1}: Invalid SMILES")
                else:
                    img_cols[idx % 4].write(f"Molecule {idx+1}: No SMILES data")
        else:
            st.warning("No molecules generated. Check the API response.")
    
    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        st.error("Troubleshooting steps:")
        st.error("1. Check your internet connection")
        st.error("2. Verify that the API key is correct")
        st.error("3. Try again later as the server might be temporarily unavailable")
        st.error("4. Contact NVIDIA support if the issue persists")

# Footer or additional information
# st.info("This app generates molecules using the NVIDIA MolMIM API and visualizes them using RDKit.")
