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


