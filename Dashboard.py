import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import io
import base64
import json

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

# import streamlit as st

# st.set_page_config(layout="wide")

# Add custom CSS for container style and hover effect
container_style = """
    <style>
    .container-box {
        background-color: #f0f2f6;
        border-radius: 10px;
        padding: 30px;  /* Adjust padding for container height */
        margin: 10px;
        cursor: pointer;
        text-align: center;
        font-size: 20px;
        color: black;
        transition: transform 0.3s ease;
        height: 200px;
        # display: flex;
        justify-content: center;
        align-items: center;
        position: relative;
    }

    .container-box:hover {
        background-color: #d9e2f3;
        transform: scale(1.05);
    }

    .container-image {
        max-height: 100px;
        max-width: 100px;
        border-radius: 10px;
    }

    /* App name styling */
    .app-name {
        font-size: 36px;
        color: #1a1f2d;
        text-align: center;
        margin-bottom: 30px;
        font-weight: bold;
    }
    </style>
"""

st.markdown(container_style, unsafe_allow_html=True)

# Add app name at the top
st.markdown('<div class="app-name">MedGuide AI Dashboard</div>', unsafe_allow_html=True)


# Define callback functions for each container
def operation_1():
    st.switch_page("./pages/0_MedGuide_Chat.py")

def operation_2():
    st.switch_page("./pages/1_Molecule_Generation.py")

def operation_3():
    st.switch_page("./pages/2_Reasearch.py")

def operation_4():
    st.switch_page("./pages/3_Physio_Planner.py")

def operation_5():
    st.switch_page("./pages/4_MediRegistry.py")

# Function to convert image to base64
def image_to_base64(image_path):
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode()

# Function to handle container clicks
def container_click(container_name, operation, image_path):
    # Convert the local image to a base64 string
    base64_image = image_to_base64(image_path)
    image_data_url = f"data:image/png;base64,{base64_image}"

    st.markdown(f"""
        <div class="container-box">
            <img src="{image_data_url}" class="container-image"/>
            <div>{container_name}</div>
        </div>
    """, unsafe_allow_html=True)

    if st.button(f"Click {container_name}", key=f"btn_{container_name}"):
        operation()

# Image URLs or paths for containers
image_1 = "img//icons//consultation.png"  # Replace with your image URL or local path
image_2 = "img//icons//thymine.png"
image_3 = "img//icons//formula.png"
image_4 = "img//icons//program.png"
image_5 = "img//icons//capsules.png"

# Create columns for the first row with 3 containers
col1, col2, col3 = st.columns(3)

with col1:
    container_click("MedGuide Chat", operation_1,image_1)

with col2:
    container_click("Molecule Generation", operation_2,image_2)

with col3:
    container_click("Research", operation_3,image_3)

# Create columns for the second row with 2 containers
col4, col5 = st.columns(2)

with col4:
    container_click("Physio Planner", operation_4,image_4)

with col5:
    container_click("MediScan", operation_5,image_5)


