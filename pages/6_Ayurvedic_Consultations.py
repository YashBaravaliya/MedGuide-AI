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

# Initialize session state
if "messages" not in st.session_state:
    st.session_state.messages = []
if "user_profile" not in st.session_state:
    st.session_state.user_profile = None

# Create LLM chain
llm_chain = create_ayurvedic_llm_chain()

def generate_consultation(query, user_profile):
    try:
        # Analyze dosha based on query and user profile
        dosha_analysis = analyze_dosha(query, user_profile)
        remedies = get_remedies(dosha_analysis)
        
        consultation_response = llm_chain.run(
            user_input=query,
            user_profile=json.dumps(user_profile) if user_profile else "Not provided",
            consultation_history=str(st.session_state.messages[-5:] if len(st.session_state.messages) > 0 else "No previous consultation")
        )
        
        return consultation_response, dosha_analysis, remedies
    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        return None, None, None

# Streamlit UI
st.title("🌿 Ayurvedic Consultation Guide")

# Create two columns
col1, col2 = st.columns([3, 2])

# Column 1: Chatbot
with col1:
    st.subheader("Consultation Chat")
    
    # Chat interface
    user_input = st.chat_input("Describe your concerns...")
    
    # Display chat messages
    for message in st.session_state.messages:
        with st.chat_message(message["role"]):
            st.markdown(message["content"])
            if "dosha_analysis" in message:
                with st.expander("View Dosha Analysis"):
                    st.json(message["dosha_analysis"])
    
    if user_input:
        # Add user message
        st.session_state.messages.append({"role": "user", "content": user_input})
        
        # Generate response
        with st.spinner("🧘‍♀️ Analyzing..."):
            response, dosha_analysis, remedies = generate_consultation(user_input, st.session_state.user_profile)
            
            if response:
                assistant_message = {
                    "role": "assistant",
                    "content": response,
                    "dosha_analysis": dosha_analysis,
                    "remedies": remedies
                }
                
                st.session_state.messages.append(assistant_message)
                
                with st.chat_message("assistant"):
                    st.markdown(response)
                    if dosha_analysis:
                        with st.expander("🔍 Detailed Analysis"):
                            st.json(dosha_analysis)
            else:
                st.error("I apologize, I couldn't generate a proper consultation. Please try rephrasing your concern.")

# Column 2: Ayurvedic Information
with col2:
    st.subheader("Ayurvedic Profile & Recommendations")
    
    # User Profile Section
    with st.expander("👤 Your Profile", expanded=True):
        if not st.session_state.user_profile:
            st.write("Please complete your profile for better consultation:")
            age = st.number_input("Age", 1, 120, 25)
            gender = st.selectbox("Gender", ["Male", "Female", "Other"])
            prakriti = st.multiselect("Known Prakriti (Body Constitution)", 
                                    ["Vata", "Pitta", "Kapha"],
                                    max_selections=2)
            
            if st.button("Save Profile"):
                st.session_state.user_profile = {
                    "age": age,
                    "gender": gender,
                    "prakriti": prakriti
                }
        else:
            st.write(f"Age: {st.session_state.user_profile['age']}")
            st.write(f"Gender: {st.session_state.user_profile['gender']}")
            st.write(f"Prakriti: {', '.join(st.session_state.user_profile['prakriti'])}")
            if st.button("Update Profile"):
                st.session_state.user_profile = None
    
    # Current Recommendations
    if st.session_state.messages and "remedies" in st.session_state.messages[-1]:
        remedies = st.session_state.messages[-1]["remedies"]
        st.markdown("### 📝 Current Recommendations")
        
        if remedies.get("herbs"):
            with st.expander("🌿 Recommended Herbs", expanded=True):
                for herb in remedies["herbs"]:
                    st.markdown(f"**{herb['name']}**")
                    st.write(herb['benefits'])
        
        if remedies.get("lifestyle"):
            with st.expander("🌅 Lifestyle Suggestions", expanded=True):
                for suggestion in remedies["lifestyle"]:
                    st.markdown(f"- {suggestion}")
        
        if remedies.get("diet"):
            with st.expander("🍽️ Dietary Guidelines", expanded=True):
                st.write(remedies["diet"])
    
    # Disclaimer
    st.markdown("""
    ---
    *Note: This is an AI-powered Ayurvedic consultation guide. 
    Always consult with a qualified Ayurvedic practitioner for 
    personalized medical advice.*
    """)