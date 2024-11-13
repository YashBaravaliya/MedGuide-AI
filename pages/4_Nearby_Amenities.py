import streamlit as st
import requests
import folium
from streamlit_folium import folium_static
import requests
from geopy.distance import geodesic
import pandas as pd
from folium.plugins import MarkerCluster
from streamlit_current_location import current_position

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


def get_medical_facilities(lat, lon, radius_km=2):
    """Fetch medical facilities from the OSM Overpass API within the specified radius."""
    try:
        radius_m = radius_km * 1000
        overpass_url = "https://overpass-api.de/api/interpreter"
        query = f"""
        [out:json][timeout:25];
        (
          node["amenity"="hospital"](around:{radius_m},{lat},{lon});
          node["amenity"="clinic"](around:{radius_m},{lat},{lon});
          node["amenity"="doctors"](around:{radius_m},{lat},{lon});
          node["healthcare"="hospital"](around:{radius_m},{lat},{lon});
          node["healthcare"="clinic"](around:{radius_m},{lat},{lon});
          node["healthcare"="doctor"](around:{radius_m},{lat},{lon});
          node["healthcare"="ayurvedic"](around:{radius_m},{lat},{lon});
          node["healthcare"="homeopathic"](around:{radius_m},{lat},{lon});
          node["healthcare"="nursing_home"](around:{radius_m},{lat},{lon});
        );
        out body;
        >;
        out skel qt;
        """
        response = requests.get(overpass_url, params={'data': query})
        data = response.json()

        facilities = []
        for element in data.get('elements', []):
            if 'tags' in element:
                name = element['tags'].get('name', 'Unnamed Medical Facility')
                facility_type = (
                    element['tags'].get('amenity', '') or
                    element['tags'].get('healthcare', 'Medical Facility')
                ).title()

                addr_parts = [
                    element['tags'].get('addr:street', ''),
                    element['tags'].get('addr:suburb', ''),
                    element['tags'].get('addr:city', '')
                ]
                address = ', '.join(part for part in addr_parts if part) or 'Address not available'

                facility = {
                    'id': f"{element['id']}",
                    'name': name,
                    'type': facility_type,
                    'lat': element['lat'],
                    'lon': element['lon'],
                    'phone': element['tags'].get('phone', 'N/A'),
                    'address': address,
                    'emergency': element['tags'].get('emergency', 'N/A'),
                    'distance': round(geodesic((lat, lon), (element['lat'], element['lon'])).kilometers, 2)
                }
                facilities.append(facility)

        return sorted(facilities, key=lambda x: x['distance'])

    except Exception as e:
        st.error(f"Error fetching medical facilities: {str(e)}")
        return []

def create_map(lat, lon, facilities, accuracy=None):
    """Create a Folium map with clustered markers for medical facilities."""
    m = folium.Map(location=[lat, lon], zoom_start=17)

    marker_cluster = MarkerCluster().add_to(m)

    # Your location marker
    folium.Marker(
        [lat, lon],
        popup='Your Location',
        icon=folium.Icon(color='red', icon='info-sign')
    ).add_to(m)

    # Accuracy circle
    if accuracy:
        folium.Circle(
            [lat, lon],
            radius=accuracy,
            color='red',
            fill=True,
            popup='Accuracy: location (low accuracy)'
        ).add_to(m)

    for facility in facilities:
        popup_html = f"""
        <div style="width: 200px">
            <b>{facility['name']}</b><br>
            Type: {facility['type']}<br>
            Distance: {facility['distance']} km<br>
            Phone: {facility['phone']}<br>
            Emergency: {facility['emergency']}<br>
            Address: {facility['address']}<br>
            <a href="https://www.google.com/maps/dir/?api=1&destination={facility['lat']},{facility['lon']}&travelmode=driving"
               target="_blank">Get Directions</a>
        </div>
        """

        folium.Marker(
            [facility['lat'], facility['lon']],
            popup=folium.Popup(popup_html, max_width=200),
            icon=folium.Icon(color='blue', icon='plus', prefix='fa')
        ).add_to(marker_cluster)

    return m

st.title("Nearby Medical Facilities")

with st.sidebar:
    st.write("Search Radius:")
    radius = st.slider("", 1, 30, 15,1)

# Get the current location
position = current_position()

# Process and display results if location is available
if position:
    lat, lon = position['latitude'], position['longitude']

    with st.spinner('Searching for medical facilities...'):
        facilities = get_medical_facilities(lat, lon, radius)

        if facilities:
            tab1, tab2 = st.tabs(["🗺️ Map", "📋 List"])

            with tab1:
                col1, col2 = st.columns([4, 2])
                with col1:
                    st.write(f"Found {len(facilities)} medical facilities within {radius} km")
                    m = create_map(lat, lon, facilities, position.get('accuracy'))
                    folium_static(m)

                with col2:
                    st.write("#### 🏥 Nearest Facilities")
                    # Divider with a customized style
                    st.markdown("---")
                    for facility in facilities[:5]:
                        # Title with bold and colored text for better visibility
                        st.markdown(f"**{facility['name']}**", unsafe_allow_html=True)
                        
                        # Add icons and color to the distance
                        st.markdown(f"🛣️ **Distance:** {facility['distance']} km", unsafe_allow_html=True)
                        
                        # Use icons for better readability
                        st.markdown(f"🏷️ **Type:** {facility['type']}", unsafe_allow_html=True)
                        
                        # Address with some space and added map preview link
                        st.markdown(f"📍 **Address:** {facility['address']}", unsafe_allow_html=True)
                        
                        # Get directions with enhanced link styling
                        st.markdown(f"[🗺️ Get Directions](https://www.google.com/maps/dir/?api=1&destination={facility['lat']},{facility['lon']}&travelmode=driving)", unsafe_allow_html=True)
                        
                        # Phone number with an icon, if available
                        if facility['phone'] != 'N/A':
                            st.markdown(f"📞 **Phone:** {facility['phone']}", unsafe_allow_html=True)
                        
                        # Divider with a customized style
                        st.markdown("---", unsafe_allow_html=True)


            with tab2:
                df = pd.DataFrame(facilities)
                st.dataframe(
                    df[['name', 'type', 'distance', 'phone', 'emergency', 'address']],
                    hide_index=True
                )
        else:
            st.warning('No medical facilities found in the specified radius.')
