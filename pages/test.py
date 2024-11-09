import streamlit as st


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
from langchain_google_genai import ChatGoogleGenerativeAI
from dotenv import load_dotenv
from datetime import datetime, timedelta
import json

# Load environment variables
load_dotenv()

api_key = os.getenv("GROQ_API_KEY")


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



import streamlit as st
from groq import Groq
from langchain_core.messages import HumanMessage, SystemMessage
from langchain_groq import ChatGroq
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate
import base64
from io import BytesIO
from PIL import Image

# Initialize Streamlit app with custom styling
st.title("Medicine Image Analyzer 💊")

# Add custom CSS
st.markdown("""
    <style>
    .stButton>button {
        width: 100%;
        background-color: #ff4b4b;
        color: white;
    }
    .stTextInput>div>div>input {
        background-color: #f0f2f6;
    }
    </style>
    """, unsafe_allow_html=True)

def resize_image(image, max_size):
    """Resize image while maintaining aspect ratio"""
    img = Image.open(image)
    img.thumbnail(max_size, Image.Resampling.LANCZOS)
    return img

def compress_image(image, quality=85):
    """Compress image to reduce file size"""
    buffer = BytesIO()
    image.save(buffer, format="JPEG", quality=quality, optimize=True)
    return buffer

def encode_image(uploaded_file, max_size=(800, 800), quality=85):
    """Convert uploaded image to base64 string with automatic size optimization"""
    if uploaded_file is not None:
        try:
            # First attempt with original size
            resized_img = resize_image(uploaded_file, max_size)
            compressed_buffer = compress_image(resized_img, quality)
            base64_str = base64.b64encode(compressed_buffer.getvalue()).decode()
            
            # If base64 string is too large, reduce size gradually
            current_size = max_size
            current_quality = quality
            
            while len(base64_str) > 250000 and (current_size[0] > 200 or current_quality > 30):
                if current_size[0] > 200:
                    current_size = (current_size[0] // 1.5, current_size[1] // 1.5)
                    uploaded_file.seek(0)  # Reset file pointer
                    resized_img = resize_image(uploaded_file, current_size)
                else:
                    current_quality -= 10
                
                compressed_buffer = compress_image(resized_img, current_quality)
                base64_str = base64.b64encode(compressed_buffer.getvalue()).decode()
                
                if current_size[0] <= 200 and current_quality <= 30:
                    st.warning("Image has been significantly compressed to meet size requirements. Quality may be reduced.")
            
            return f"data:image/jpeg;base64,{base64_str}"
        except Exception as e:
            raise ValueError(f"Error processing image: {str(e)}")
    return None

def create_langchain_agent():
    """Create LangChain agent with Groq"""
    try:
        api_key = os.getenv("GROQ_API_KEY")
        llm = ChatGroq(
            groq_api_key=api_key,
            model_name="llama-3.2-11b-vision-preview",
            max_tokens=4000
        )
        
        template = """just give me a correct medicine name in english
        
        Image: {image_data}
        """

        prompt = PromptTemplate(
            input_variables=["image_data"],
            template=template
        )
        
        return LLMChain(llm=llm, prompt=prompt)
    except Exception as e:
        raise Exception(f"Error creating LangChain agent: {str(e)}")

def analyze_medicine(chain, image_data, uploaded_file):
    """Analyze medicine using the LangChain agent with automatic retry"""
    try:
        response = chain.run(image_data=image_data)
        return response
    except Exception as e:
        if "rate_limit_exceeded" in str(e):
            try:
                # Retry with more aggressive compression
                st.info("Reducing image size for analysis...")
                uploaded_file.seek(0)  # Reset file pointer
                smaller_image_data = encode_image(
                    uploaded_file,
                    max_size=(400, 400),  # Smaller max size
                    quality=60  # Lower quality
                )
                return chain.run(image_data=smaller_image_data)
            except Exception as retry_e:
                if "rate_limit_exceeded" in str(retry_e):
                    try:
                        # Final attempt with minimum size
                        st.info("Making final attempt with minimum image size...")
                        uploaded_file.seek(0)
                        smallest_image_data = encode_image(
                            uploaded_file,
                            max_size=(200, 200),  # Minimum size
                            quality=30  # Minimum quality
                        )
                        return chain.run(image_data=smallest_image_data)
                    except Exception as final_e:
                        return "Error: Unable to analyze image even after maximum compression. Please try a different image."
                return f"Error during retry: {str(retry_e)}"
        return f"Error analyzing image: {str(e)}"

# Sidebar for API key
api_key = st.sidebar.text_input("Enter your Groq API key", type="password")

# Main app layout
st.header("Upload Medicine Image")
uploaded_file = st.file_uploader("Choose an image...", type=['jpg', 'jpeg', 'png'])

if uploaded_file is not None:
    try:
        # Display uploaded image
        col1, col2 = st.columns(2)
        with col1:
            st.image(uploaded_file, caption="Uploaded Medicine Image", use_column_width=True)
        
        # Analyze button
        with col2:
            if st.button("Analyze Medicine"):
                with st.spinner("Analyzing image..."):
                    try:
                        # Initial conversion with standard settings
                        image_data = encode_image(uploaded_file)
                        
                        # Create LangChain agent
                        chain = create_langchain_agent()
                        
                        # Get analysis with automatic retry
                        analysis = analyze_medicine(chain, image_data, uploaded_file)
                        
                        # Display results
                        st.markdown("### Analysis Results")
                        st.markdown(analysis)
                        
                        # Add download button for results
                        st.download_button(
                            label="Download Analysis",
                            data=analysis,
                            file_name="medicine_analysis.txt",
                            mime="text/plain"
                        )
                    except Exception as e:
                        st.error(f"Error during analysis: {str(e)}")
    except Exception as e:
        st.error(f"Error processing upload: {str(e)}")

# Footer
st.markdown("---")
st.markdown("""
    <div style='text-align: center'>
        <p>⚠️ This tool is for informational purposes only. Always consult with a healthcare professional before taking any medication.</p>
    </div>
    """, unsafe_allow_html=True)