import streamlit as st
from utils.ayurvedic_utils import analyze_dosha, get_remedies
from models.llm_chain import create_ayurvedic_llm_chain
import json

# Custom CSS with Ayurvedic theme
custom_style = """
    <style>
    /* Sidebar styling with earthy colors */
    [data-testid="stSidebar"] {
        background-color: #1a1f2d;
        color:white;
    }
    
    /* Main content styling */
    .main {
        background-color: #eff3f8;
        color: black;
    }
    
    /* Container styling */
    [data-testid="stAppViewBlockContainer"] {
        background-color: #fefefe;
        padding: 20px;
        border-radius: 10px;
        box-shadow: 0 2px 6px rgba(139, 69, 19, 0.1);
        color:black;
    }

    /* Sidebar navigation styling */
    div[data-testid="stSidebarNav"] li div a span {
        color: white;
        padding: 0.5rem;
        width: 300px;
        border-radius: 0.5rem;
    }

    div[data-testid="stSidebarNav"] li div::focus-visible {
        background-color: rgba(151, 166, 195, 0.15);
    }

    /* Custom header for sidebar */
    [data-testid="stSidebarNav"]::before {
        content: "💊 MedGuide AI 🧑🏻‍⚕️";
        margin-left: 20px;
        margin-top: 20px;
        font-size: 30px;
        position: relative;
        top: -50px;
        # color: #F5DEB3;
    }

    /* Custom styling for headers */
    h1, h2, h3 {
        font-family: 'Helvetica Neue', sans-serif;
    }

    /* Custom card styling */
    .dosha-card {
        # background-color: #FFF;
        padding: 20px;
        border-radius: 10px;
        border-left: 5px solid #8B4513;
        margin: 10px 0;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    </style>
"""

st.markdown(custom_style, unsafe_allow_html=True)
import streamlit as st
import pandas as pd
import os
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate
from langchain_groq import ChatGroq
from dotenv import load_dotenv
from datetime import datetime, timedelta
import json

# Load environment variables
load_dotenv()



# Enhanced Custom CSS
st.markdown("""
<style>
    .stSelectbox, .stMultiSelect {
        margin-bottom: 1rem;
    }
    .exercise-card {
        background-color: #f8f9fa;
        padding: 1.5rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }
    .main-header {
        color: #1f77b4;
        margin-bottom: 2rem;
    }
    .ai-response {
        background-color: #e3f2fd;
        padding: 1.5rem;
        border-radius: 0.5rem;
        margin: 1rem 0;
        border-left: 4px solid #1976d2;
    }
    .day-card {
        background-color: white;
        padding: 1rem;
        border-radius: 0.5rem;
        margin: 0.5rem 0;
        border: 1px solid #ddd;
    }
    .muscle-tag {
        background-color: #e0e0e0;
        padding: 0.2rem 0.5rem;
        border-radius: 0.25rem;
        margin-right: 0.5rem;
        font-size: 0.9em;
    }
    .weekly-schedule {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
        gap: 1rem;
        padding: 1rem;
    }
</style>
""", unsafe_allow_html=True)
