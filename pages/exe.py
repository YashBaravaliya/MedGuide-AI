import pandas as pd
import streamlit as st
from langchain.chains import LLMChain
from langchain.prompts import PromptTemplate
from langchain_groq import ChatGroq
from dotenv import load_dotenv
import os
# from config import GROQ_API_KEY

# Load environment variables
load_dotenv()


@st.cache_data
def load_exercise_data(path):
    return pd.read_csv(path)

def create_llm_chain(template):
    GROQ_API_KEY = os.getenv('GROQ_API_KEY')
    llm = ChatGroq(
        groq_api_key=GROQ_API_KEY,
        model_name="gemma-7b-it",
        temperature=0.5,
        max_tokens=500
    )
    
    prompt = PromptTemplate(
        input_variables=template["variables"],
        template=template["content"]
    )
    
    return LLMChain(llm=llm, prompt=prompt)

# src/components/exercise_browser.py
import streamlit as st
import pandas as pd

def render_exercise_browser(df):
    with st.sidebar:
        st.title("Exercise Filters")
        st.markdown("---")
        
        force_options = ["All"] + list(df['force'].dropna().unique())
        level_options = ["All"] + list(df['level'].unique())
        mechanic_options = ["All"] + list(df['mechanic'].dropna().unique())
        equipment_options = ["All"] + list(df['equipment'].dropna().unique())
        
        filters = {
            "force": st.selectbox("Force Type", force_options, key='browser_force'),
            "level": st.selectbox("Experience Level", level_options, key='browser_level'),
            "mechanic": st.selectbox("Movement Type", mechanic_options, key='browser_mechanic'),
            "equipment": st.selectbox("Equipment Required", equipment_options, key='browser_equipment'),
            "muscles": st.multiselect("Target Muscles", df['primaryMuscles'].unique(), key='browser_muscles')
        }
    
    st.title("Exercise Browser")
    
    filtered_data = filter_exercises(df, filters)
    display_exercises(filtered_data)

def filter_exercises(df, filters):
    filtered_data = df.copy()
    
    if filters["force"] != "All":
        filtered_data = filtered_data[filtered_data['force'] == filters["force"]]
    if filters["level"] != "All":
        filtered_data = filtered_data[filtered_data['level'] == filters["level"]]
    if filters["mechanic"] != "All":
        filtered_data = filtered_data[filtered_data['mechanic'] == filters["mechanic"]]
    if filters["equipment"] != "All":
        filtered_data = filtered_data[filtered_data['equipment'] == filters["equipment"]]
    if filters["muscles"]:
        filtered_data = filtered_data[filtered_data['primaryMuscles'].isin(filters["muscles"])]
    
    return filtered_data

def display_exercises(filtered_data):
    if not filtered_data.empty:
        col1, col2 = st.columns([2, 1])
        with col1:
            st.subheader(f"Found {len(filtered_data)} exercises matching your criteria")
        with col2:
            if st.button("🎲 Random Exercise"):
                filtered_data = filtered_data.sample(n=1)
        
        for idx, row in filtered_data.iterrows():
            with st.container():
                st.markdown("---")
                cols = st.columns([3, 1])
                with cols[0]:
                    st.subheader(row['name'])
                    st.markdown(f"**Level:** {row['level']} | **Equipment:** {row['equipment']}")
                    
                    instructions = row['instructions']
                    if isinstance(instructions, str):
                        try:
                            instructions = eval(instructions)
                        except:
                            instructions = [instructions]
                    
                    st.markdown("**Instructions:**")
                    # for i, instruction in enumerate(instructions, 1):
                    #     st.markdown(f"{i}. {instruction}")
                
                with cols[1]:
                    images = row['images']
                    if isinstance(images, str):
                        image_list = images.split(", ")
                        for img_path in image_list:
                            st.image(
                                f"Assets/exercises/{img_path.strip()}", 
                                use_column_width=True,
                                caption="Exercise demonstration"
                            )
    else:
        st.warning("No exercises found matching your criteria. Try adjusting the filters.")

# src/components/workout_planner.py
# import streamlit as st
from datetime import datetime
# from utils import create_llm_chain

WORKOUT_TEMPLATE = {
    "variables": ["muscles", "level", "equipment", "days_per_week", "special_requirements"],
    "content": """
    Create a personalized weekly workout plan based on the following information:
    Target Muscles: {muscles}
    Experience Level: {level}
    Available Equipment: {equipment}
    Days per Week: {days_per_week}
    Special Requirements: {special_requirements}

    Please provide a detailed plan with:
    1. Specific exercises for each day
    2. Sets and reps recommendations
    3. Rest periods
    4. Tips for progression
    """
}

def render_workout_planner(df):
    st.title("Personal Workout Planner")
    
    col1, col2 = st.columns([1, 1])
    
    with col1:
        plan_inputs = {
            "target_muscles": st.multiselect(
                "Select Target Muscles",
                df['primaryMuscles'].unique(),
                key='planner_muscles'
            ),
            "experience_level": st.select_slider(
                "Experience Level",
                options=['Beginner', 'Intermediate', 'Advanced'],
                value='Intermediate'
            ),
            "available_equipment": st.multiselect(
                "Available Equipment",
                df['equipment'].unique(),
                # default=['Dumbbell', 'Bodyweight']
            )
        }
        
    with col2:
        plan_inputs.update({
            "days_per_week": st.slider(
                "Days per Week",
                min_value=1,
                max_value=7,
                value=4
            ),
            "special_requirements": st.text_area(
                "Special Requirements/Notes",
                placeholder="E.g., injuries, time constraints, specific goals..."
            )
        })
    
    if st.button("Generate Workout Plan"):
        generate_workout_plan(plan_inputs)

def generate_workout_plan(inputs):
    try:
        chain = create_llm_chain(WORKOUT_TEMPLATE)
        response = chain.run({
            "muscles": ", ".join(inputs["target_muscles"]),
            "level": inputs["experience_level"],
            "equipment": ", ".join(inputs["available_equipment"]),
            "days_per_week": inputs["days_per_week"],
            "special_requirements": inputs["special_requirements"]
        })
        
        st.markdown("### Your Personalized Workout Plan")
        st.markdown(response)
        
        st.download_button(
            "Download Plan",
            response,
            file_name=f"workout_plan_{datetime.now().strftime('%Y%m%d')}.txt",
            mime="text/plain"
        )
    except Exception as e:
        st.error(f"Error generating workout plan: {str(e)}")

# src/components/fitness_assistant.py
import streamlit as st
# from utils import create_llm_chain

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
    
    # Display chat history
    display_chat_history()
    
    # Chat input
    handle_user_input()

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
                    f'<div class="chat-message bot-message">🤖 Assistant: {content}</div>', 
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
        chain = create_llm_chain(CHAT_TEMPLATE)
        response = chain.run(user_question=user_input)
        
        # Add assistant response to chat history
        st.session_state.chat_history.append({"role": "assistant", "content": response})
        
        # Rerun to update chat display
        # st.experimental_rerun()
    except Exception as e:
        st.error(f"Error: {str(e)}")

# src/app.py
import streamlit as st
import pandas as pd
# from config import PAGE_CONFIG, CUSTOM_CSS, DATA_PATH
# from components.exercise_browser import render_exercise_browser
# from components.workout_planner import render_workout_planner
# from components.fitness_assistant import render_fitness_assistant
# from utils import load_exercise_data

def main():
    # Configure page
    # st.set_page_config(**PAGE_CONFIG)
    # st.markdown(CUSTOM_CSS, unsafe_allow_html=True)
    
    # Load data
    df = load_exercise_data('Assets\dist\exercises.csv')
    
    # Create tabs
    tab1, tab2, tab3 = st.tabs(["Exercise Browser", "Workout Planner", "Fitness Assistant"])
    
    with tab1:
        render_exercise_browser(df)
    
    with tab2:
        render_workout_planner(df)
    
    with tab3:
        render_fitness_assistant()
    
    # Footer
    st.markdown("---")
    st.markdown(
        """
        <div style='text-align: center; color: #666;'>
            Remember to consult with a fitness professional before starting a new exercise routine.
        </div>
        """, 
        unsafe_allow_html=True
    )

if __name__ == "__main__":
    main()