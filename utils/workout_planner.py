import streamlit as st
from datetime import datetime
from models.llm_chain import gym_llm_chain

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
        chain = gym_llm_chain(WORKOUT_TEMPLATE)
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
