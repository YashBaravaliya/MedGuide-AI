import streamlit as st
import pandas as pd

def render_exercise_browser(df):
    # with st.expander("## Exercise Filters", expanded=True):
    st.title("Exercise Filters")
    # st.markdown("---")
    
    force_options = ["All"] + list(df['force'].dropna().unique())
    level_options = ["All"] + list(df['level'].unique())
    mechanic_options = ["All"] + list(df['mechanic'].dropna().unique())
    equipment_options = ["All"] + list(df['equipment'].dropna().unique())

    col1,col2 = st.columns(2)

    with col1:
        force = st.selectbox("Force Type",options=force_options,index=1)
        level = st.select_slider("Experience Level",options=level_options,value="beginner")
        equipment = st.selectbox("Equipment Required", equipment_options, key='browser_equipment',index=1)

    with col2:
        mechanic = st.radio("Movement Type",options=mechanic_options)
        muscles = st.multiselect("Target Muscles", df['primaryMuscles'].unique(), key='browser_muscles')

    filters = [force,level,mechanic,equipment,muscles]
        # filters = {
            # force,level,equipment,mechanic,muscles
            # "force": st.selectbox("Force Type", force_options, key='browser_force'),
            # "level": st.selectbox("Experience Level", level_options, key='browser_level'),
            # "mechanic": st.selectbox("Movement Type", mechanic_options, key='browser_mechanic'),
            # "equipment": st.selectbox("Equipment Required", equipment_options, key='browser_equipment'),
            # "muscles": st.multiselect("Target Muscles", df['primaryMuscles'].unique(), key='browser_muscles')
        # }
    
    # st.title("Exercise Browser")

    if st.button("set filters"):
        filtered_data = filter_exercises(df, filters)
        display_exercises(filtered_data)
        return filtered_data

def filter_exercises(df, filters):
    filtered_data = df.copy()
    
    if filters[0] != "All":
        filtered_data = filtered_data[filtered_data['force'] == filters[0]]
    if filters[1] != "All":
        filtered_data = filtered_data[filtered_data['level'] == filters[1]]
    if filters[2] != "All":
        filtered_data = filtered_data[filtered_data['mechanic'] == filters[2]]
    if filters[3] != "All":
        filtered_data = filtered_data[filtered_data['equipment'] == filters[3]]
    if filters[4]:
        filtered_data = filtered_data[filtered_data['primaryMuscles'].isin(filters[4])]
    
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
                    try:
                        for i, instruction in enumerate(instructions, 1):
                            st.markdown(f"{i}. {instruction}")
                    except Exception as e:
                        pass
                
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
