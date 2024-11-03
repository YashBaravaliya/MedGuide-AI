import streamlit as st
from utils.pubmed_utils import fetch_pubmed_articles_for_research
from models.llm_chain import research_llm_chain
import json

# Apply custom CSS
custom_style = """
    <style>
    [data-testid="stSidebar"] {
        background-color: #1a1f2d;
        color: white;
    }
    
    .main {
        color: black;
    }
    
    [data-testid="stAppViewBlockContainer"] {
        padding: 20px;
        border-radius: 10px;
        color: black;
    }
    
    div[data-testid="stSidebarNav"] li div a span {
        color: white;
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

# Create LLM chain
llm_chain = research_llm_chain()

def generate_research_analysis(query):
    try:
        pubmed_results = fetch_pubmed_articles_for_research(query)
        
        pubmed_formatted = "\n\n".join([
            f"Title: {article['meta']['title']}\n"
            f"Abstract: {article['content']}"
            + (f"\nDOI: {article['meta']['doi']}" if article['meta'].get('doi') else "")
            + (f"\nPubMed ID: {article['meta']['pubmed_id']}" if article['meta'].get('pubmed_id') else "")
            for article in pubmed_results
        ])

        response = llm_chain.run(
            question=query,
            pubmed_results=pubmed_formatted
        )
        
        return response, pubmed_results
    except Exception as e:
        st.error(f"An error occurred: {str(e)}")
        return None, None

# Streamlit UI
st.title("Research Assistant")

# Create two columns
col1, col2 = st.columns([3, 2])

# Column 1: Research Analysis
with col1:
    st.subheader("Research Analysis")
    
    user_input = st.text_input("Enter your research query...",placeholder="West Nile virus")
    
    if st.button("Analyze"):
        if user_input:
            with st.spinner("🔬 Analyzing research..."):
                analysis, articles = generate_research_analysis(user_input)
                
                if analysis:
                    st.markdown(analysis)
                else:
                    st.error("Sorry, I couldn't generate an analysis. Please try another query.")

# Column 2: Article Details
with col2:
    st.subheader("Source Articles")
    
    if user_input and 'articles' in locals():
        for idx, article in enumerate(articles):
            with st.expander(f"Article {idx + 1}", expanded=True):
                title = article['meta']['title']
                clean_title = ' '.join([word for word in title.split() if not word.strip().isdigit()])
                
                # Display the main title with the first PubMed ID as a link if available
                if article['meta'].get('pubmed_id'):
                    pubmed_ids = article['meta']['pubmed_id'].split()
                    st.markdown(f"**[{clean_title}](https://pubmed.ncbi.nlm.nih.gov/{pubmed_ids[0]}/)**")
                else:
                    st.markdown(f"**{clean_title}**")
                
                if article['meta'].get('publication_date'):
                    st.markdown("#### Publication Date")
                    st.markdown(article['meta']['publication_date'])
                
                st.markdown("#### Abstract")
                st.markdown(article['content'])
                
                if article['meta'].get('keywords') and article['meta']['keywords']:
                    st.markdown("#### Keywords")
                    st.markdown(", ".join(article['meta']['keywords']))
                
                if article['meta'].get('pubmed_id'):
                    st.markdown("#### PubMed Links")
                    pubmed_ids = article['meta']['pubmed_id'].split()
                    for pubmed_id in pubmed_ids:
                        st.markdown(f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/")
    else:
        st.info("Enter a research query and click 'Analyze' to see article details.")