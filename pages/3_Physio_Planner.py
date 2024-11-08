import pandas as pd
import streamlit as st
from utils.exercise_browser import render_exercise_browser
from utils.workout_planner import render_workout_planner
from utils.fitness_assistant import render_fitness_assistant

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

@st.cache_data
def load_exercise_data(path):
    return pd.read_csv(path)

def main():

    st.write("---")
    st.header("🏋️‍♂️ Physio Planner")
    
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