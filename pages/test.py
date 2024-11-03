import streamlit as st

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = []

# Function to generate AI response (placeholder)
def get_ai_response(user_input):
    return f"AI: I received your message: '{user_input}'"

# Layout
st.title("Medical Chatbot with Split Layout")

# Create two columns
col1, col2 = st.columns([4, 2])

# Column 1: Chatbot
with col1:
    st.subheader("Chatbot")
    # Chat input
    user_input = st.chat_input("Type your message here...")

    # Display chat messages
    for message in reversed(st.session_state.messages):
        with st.chat_message(message["role"]):
            st.write(message["content"])

    if user_input:
        # Add user message to chat history
        st.session_state.messages.append({"role": "user", "content": user_input})

        # Generate and add AI response
        ai_response = get_ai_response(user_input)
        st.session_state.messages.append({"role": "assistant", "content": ai_response})

        # Force a rerun to update the chat display
        st.rerun()

# Column 2: Medication Information
with col2:
    st.subheader("Medication Information")
    
    # Sample medication data
    medications = {
        "Aspirin": "Used for pain relief and reducing inflammation.",
        "Ibuprofen": "Used for pain, fever, and inflammation.",
        "Paracetamol": "Used for pain relief and fever reduction.",
    }
    
    # Display medication information
    selected_medication = st.selectbox("Select a medication:", list(medications.keys()))
    st.write(f"**{selected_medication}**: {medications[selected_medication]}")
    
    st.write("This space includes information on various medications. You can select a medication from the list to learn more about it.")
