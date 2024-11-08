from langchain_community.tools import DuckDuckGoSearchResults
from langchain.chains import LLMChain
# from langchain_openai import ChatOpenAI
from langchain.prompts import PromptTemplate
from langchain_groq import ChatGroq
import os

class HealthcareSearchAgent:
    def __init__(self):
        GROQ_API_KEY = os.getenv('GROQ_API_KEY')
        # Initialize the search tool
        self.search_tool = DuckDuckGoSearchResults()
        
        # Initialize the LLM
        self.llm = ChatGroq(
            groq_api_key=GROQ_API_KEY,
            model_name="gemma-7b-it",
            temperature=0.5,
            max_tokens=500
        )
        
        # Create prompt template
        self.prompt = PromptTemplate(
            input_variables=["search_results"],
            template="""You are a healthcare information specialist. 
            Given the following search results, provide the 3 most relevant medical information sources 
            and their corresponding website links. Focus only on reputable medical sources like medical institutions,
            government health websites, and recognized medical organizations.
            
            Search Results: {search_results}
            
            Please provide the top 3 most relevant sources in this exact format:
            
            1. URL: [insert URL here]
            Summary: [insert brief summary here]
            
            2. URL: [insert URL here]
            Summary: [insert brief summary here]
            
            3. URL: [insert URL here]
            Summary: [insert brief summary here]"""
        )
        
        # Create the chain
        self.chain = LLMChain(llm=self.llm, prompt=self.prompt)
        
    
    def search(self, query: str) -> dict:
        try:
            # Perform the search
            search_results = self.search_tool.run(query)
            
            # Process results through the chain
            response = self.chain.run(search_results=search_results)
            
            # Parse the response into three separate results
            results = []
            sections = response.split("\n\n")
            
            for section in sections:
                if section.strip() and "URL:" in section:
                    url = section.split("URL:")[1].split("\n")[0].strip()
                    summary = section.split("Summary:")[1].strip()
                    results.append({
                        "url": url,
                        "summary": summary
                    })
            
            return {
                "results": results[:3],  # Ensure we only return top 3
                "status": "success"
            }
        except Exception as e:
            return {
                "status": "error",
                "message": str(e)
            }

def main():
    # Initialize the search system
    healthcare_search = HealthcareSearchAgent()
    
    # Example query
    query = "What are the symptoms of diabetes?"
    
    # Get results
    result = healthcare_search.search(query)
    
    # Print results
    if result["status"] == "success":
        print("\nTop 3 Medical Sources:\n")
        for i, res in enumerate(result["results"], 1):
            print(f"{i}. URL: {res['url']}")
            print(f"   Summary: {res['summary']}\n")
    else:
        print(f"Error: {result['message']}")

if __name__ == "__main__":
    main()