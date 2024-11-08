import streamlit as st
from models.llm_chain import gym_llm_chain

CHAT_TEMPLATE = {
    "variables": ["user_question"],
    "content": """
    You are a knowledgeable fitness assistant. Please provide helpful and accurate information about:
    {user_question}
    
    Base your response on scientific principles and best practices in fitness and exercise science.
    """
}

def render_fitness_assistant():
    st.title("Fitness Assistant")
    
    # Initialize chat history
    if 'chat_history' not in st.session_state:
        st.session_state.chat_history = []
    
    # Chat input
    handle_user_input()
    # Display chat history
    display_chat_history()

def display_chat_history():
    for message in st.session_state.chat_history:
        role = message["role"]
        content = message["content"]
        with st.container():
            if role == "user":
                st.markdown(
                    f'<div class="chat-message user-message">👤 You: {content}</div>', 
                    unsafe_allow_html=True
                )
            else:
                st.markdown(
                    f'<div class="chat-message bot-message">🤖 Assistant: \n{content}</div>', 
                    unsafe_allow_html=True
                )

def handle_user_input():
    user_input = st.text_input(
        "Ask anything about fitness, exercises, or your workout plan:",
        key="user_input"
    )
    
    if st.button("Send", key="send_button"):
        if user_input:
            process_user_message(user_input)

def process_user_message(user_input):
    # Add user message to chat history
    st.session_state.chat_history.append({"role": "user", "content": user_input})
    
    try:
        chain = gym_llm_chain(CHAT_TEMPLATE)
        response = chain.run(user_question=user_input)
        
        # Add assistant response to chat history
        st.session_state.chat_history.append({"role": "assistant", "content": response})
        
        # Rerun to update chat display
        # st.experimental_rerun()
    except Exception as e:
        st.error(f"Error: {str(e)}")