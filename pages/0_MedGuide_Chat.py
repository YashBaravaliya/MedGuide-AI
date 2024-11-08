import streamlit as st
from utils.pubmed_utils import fetch_pubmed_articles
from utils.medicine_utils import search_medicine
from utils.pdf_utils import generate_pdf
from utils.agents import HealthcareSearchAgent
from models.llm_chain import create_llm_chain
from ui.custom_css import apply_custom_css
import requests
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import io
import base64
import json

# Apply custom CSS
# apply_custom_css()
custom_style = """
    <style>
    /* Sidebar background color */
    [data-testid="stSidebar"] {
        background-color: #1a1f2d;
        color:white;
    }
    
    /* Main content background color (container) */
    .main {
        # background-color: #eff3f8;
        color: black;  /* Change text color for better contrast */
    }
    
    /* Container background color */
    [data-testid="stAppViewBlockContainer"] {
        # background-color: #fefefe;
        padding: 20px;
        border-radius: 10px;
        color: black; /* Set text color to white for better contrast */
    }
    div[data-testid="stSidebarNav"] li div a span {
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
st.markdown(custom_style, unsafe_allow_html=True)

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = []

# Create LLM chain
llm_chain = create_llm_chain()

def generate_answer(question):
    try:
        pubmed_results = fetch_pubmed_articles(question)
        medicine_results = search_medicine(question)
        
        if medicine_results:
            primary_medicine = medicine_results["primary"]
            alternative_medicines = medicine_results.get("alternatives", [])
            
            pubmed_formatted = "\n\n".join([f"Title: {article['meta']['title']}\nAbstract: {article['content']}" for article in pubmed_results])
            medicine_formatted = json.dumps(primary_medicine, indent=2)
            alternative_formatted = json.dumps(alternative_medicines, indent=2)

            response = llm_chain.run(
                question=question,
                pubmed_results=pubmed_formatted,
                medicine_info=medicine_formatted,
                alternative_medicines=alternative_formatted
            )

            image_data = []
            if primary_medicine.get('img_url'):
                image_data.append(("Primary Medicine", primary_medicine['img_url']))
            for alt_med in alternative_medicines:
                if alt_med.get('img_url'):
                    image_data.append((f"Alternative: {alt_med['Name']}", alt_med['img_url']))

            pdf_buffer = generate_pdf(primary_medicine, alternative_medicines, response)

            return response, image_data, True, pdf_buffer, medicine_results
        else:
            general_response = llm_chain.run(
                question=question,
                pubmed_results="\n\n".join([f"Title: {article['meta']['title']}\nAbstract: {article['content']}" for article in pubmed_results]),
                medicine_info="No specific medicine information found in the database.",
                alternative_medicines="No alternative medicines found."
            )
            return general_response, None, False, None, None
    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        return None, None, False, None, None
    
def seach_query(query):
    healthcare_search = HealthcareSearchAgent()
    result = healthcare_search.search(query)
    return result
    
# Streamlit UI
st.title("MedGuide Chat")

# Create two columns
col1, col2 = st.columns([4, 2])

# Column 1: Chatbot
with col1:
    st.subheader("Chatbot")
    
    # Chat input
    user_input = st.chat_input("Type your message here...")

    # Display chat messages
    for message in reversed(st.session_state.messages):  # try to revresed the massages
        with st.chat_message(message["role"]):
            st.markdown(message["content"])
            if "image_data" in message:
                for caption, url in message["image_data"]:
                    st.image(url, caption=caption, use_column_width=True)
    
    if user_input:
        # Add user message to chat history
        st.session_state.messages.append({"content": user_input,"role": "user"})
        
        # Generate and add AI response
        with st.spinner("🧠 Thinking..."):
            answer, image_data, is_medicine_query, pdf_buffer, medicine_results = generate_answer(user_input)
            
            if answer:
                # Create assistant message with text and image data
                assistant_message = {
                    "role": "assistant",
                    "content": answer,
                    "user_input": user_input,
                    # "image_data": image_data if image_data else []
                    "medicine_results": medicine_results if is_medicine_query else None
                }
                
                # Add assistant message to chat history
                st.session_state.messages.append(assistant_message)
                
                # Display assistant response in chat message container
                with st.chat_message("assistant"):
                    st.markdown(answer)
                    if is_medicine_query and image_data:
                        for caption, url in image_data:
                            st.image(url, caption=caption, use_column_width=True)
                    if is_medicine_query and pdf_buffer:
                        st.download_button(
                            label="Download Full Medicine Information",
                            data=pdf_buffer,
                            file_name=f"{medicine_results['primary']['Name']}.pdf",
                            mime="application/pdf"
                        )
                    elif not is_medicine_query:
                        st.info("No specific medicine information found in the database. This response is based on general AI knowledge.")
            else:
                st.error("Sorry, I couldn't generate an answer. Please try rephrasing your question.")
        
        # Force a rerun to update the chat display
        st.rerun()

# Column 2: Medication Information
with col2:
    st.subheader("Medication Information")
    
    if st.session_state.messages and "medicine_results" in st.session_state.messages[-1]:
        medicine_results = st.session_state.messages[-1]["medicine_results"]
        if medicine_results:
            primary_medicine = medicine_results["primary"]
            st.write(f"### {primary_medicine['Name']}")
            st.image(primary_medicine['img_url'], use_column_width=True)
            
            st.write("#### Details:")
            # st.write(primary_medicine)
            st.write(f"**Uses:** {primary_medicine['Uses']}")
            st.write(f"**Composition:** {primary_medicine['Composition']}")
            st.write(f"**Manufacturer:** {primary_medicine['Manufacturer']}")
            st.write(f"**Side_effects:** {primary_medicine['Side_effects']}")
            # st.write(f"**Dosage:** {primary_medicine.get('Dosage', 'N/A')}")
            # st.write(f"**Side Effects:** {primary_medicine.get('Side Effects', 'N/A')}")
            # st.write(f"**Precautions:** {primary_medicine.get('Precautions', 'N/A')}")
            
            if medicine_results.get("alternatives"):
                st.write("### Alternative Medicines")
                for alt_med in medicine_results["alternatives"]:
                    st.write(f"**{alt_med['Name']}**")
                    st.image(alt_med['img_url'], width=150)
                    # st.write(f"Category: {alt_med.get('Category', 'N/A')}")
    else:
        st.write("No specific medication information is available for the current query. Ask about a medicine to see details here.")

    try:
        search = st.session_state.messages[-1]["user_input"]
        print(f"this is seach: {search}")
        st.write(search)
        st.subheader("Search for medical sources related to the current query:")
        # Initialize the search system
        healthcare_search = HealthcareSearchAgent()
        result = healthcare_search.search(search)
        
        # Print results
        if result["status"] == "success":
            # st.subheader("Top 3 Medical Sources:")
            for i, res in enumerate(result["results"], 1):
                st.markdown(f"{i}. URL: {res['url'].split(' ')[1]}")
                st.markdown(f"Summary: {res['summary']}")
        else:
            st.error(f"Error: {result['message']}")
    except Exception as e:
        pass