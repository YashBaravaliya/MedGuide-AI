from langchain_groq import ChatGroq
from langchain.prompts import PromptTemplate
from langchain.chains import LLMChain
import os

def create_llm_chain():
    GROQ_API_KEY = os.getenv('GROQ_API_KEY')
    
    llm = ChatGroq(
        groq_api_key=GROQ_API_KEY,
        model_name="gemma-7b-it",
        temperature=0.5,
        max_tokens=500
    )
    prompt_template = PromptTemplate(
        input_variables=["question", "pubmed_results", "medicine_info"],
        template="""
        **Question:** {question}
        
        **Medicine Information:**
        {medicine_info}

        **Related Research Articles:**
        {pubmed_results}

        **Guidelines for Answer:**
        - **Accuracy:** Ensure that all information provided is medically accurate and up-to-date.
        - **Clarity:** Explain complex medical terms in simple language that patients can easily understand.
        - **Patient Focus:** Tailor your response to prioritize patient safety and understanding. Provide clear instructions or advice where applicable.
        - **Detailed Explanation:** Offer a detailed explanation of the medicine, including its uses, potential side effects, and any interactions it may have.
        - **Simplified Summary:** Provide a brief, simplified summary of the information for patients who may need a quick overview.
        - **Professional Tone:** Maintain a professional, empathetic, and reassuring tone throughout your response.
        - **Alternative Ayurvedic Medicine:** Suggest an alternative medicine or treatment option if applicable, including its Ayurvedic plant name and their usecase which is most relative to patients problems.
        - **Ayurvedic Plant Composition:** Include the composition of the Ayurvedic plants in the alternative medicine, highlighting their key properties and benefits and how they can be used in the treatment.
        
        **Final Answer:**
        Please generate a response that combines the detailed medical information, relevant research insights, and patient-friendly advice based on the information provided above. Include information about the Ayurvedic plants in the composition and suggest an alternative medicine or treatment option at the end of your response. Ensure that the response is accurate, clear, and patient-focused. 
        """
    )

    return LLMChain(llm=llm, prompt=prompt_template)

def research_llm_chain():
    GROQ_API_KEY = os.getenv('GROQ_API_KEY')
    
    llm = ChatGroq(
        groq_api_key=GROQ_API_KEY,
        model_name="gemma-7b-it",
        temperature=0.5,
        max_tokens=1000  # Increased for more detailed research summaries
    )
    
    prompt_template = PromptTemplate(
        input_variables=["question", "pubmed_results"],
        template="""
        **Research Query:** {question}
        
        **Available Research Articles:**
        {pubmed_results}

        **Guidelines for Analysis:**
        1. **Comprehensive Summary:**
           - Synthesize the key findings from all provided research articles
           - Highlight any conflicting results or perspectives
           - Identify gaps in current research if apparent

        2. **Methodology Analysis:**
           - Comment on the research methods used in the studies
           - Highlight any limitations or strengths in the methodologies

        3. **Impact Assessment:**
           - Discuss the potential implications of the research findings
           - Consider both theoretical and practical applications

        4. **Future Directions:**
           - Suggest potential areas for future research based on current findings
           - Identify any unanswered questions or areas needing more investigation

        5. **Technical Details:**
           - Explain any complex technical terms or concepts in accessible language
           - Provide context for specialized terminology

        **Response Structure:**
        1. Brief overview of the research question
        2. Summary of key findings across all articles
        3. Analysis of methodologies and their implications
        4. Discussion of potential applications and impact
        5. Identification of future research directions
        6. Simplified explanation of technical concepts

        **Final Analysis:**
        Please provide a comprehensive analysis that synthesizes the research findings, evaluates methodologies, discusses implications, and suggests future directions. Ensure the response is both scholarly and accessible, suitable for researchers while remaining comprehensible to an informed general audience.
        """
    )

    return LLMChain(llm=llm, prompt=prompt_template)

def generate_smiles(chemical_description):
    MoleculeMirth = os.getenv("MoleculeMirth")
    llm = ChatGroq(
        groq_api_key=MoleculeMirth,
        model_name="llama3-groq-70b-8192-tool-use-preview",
        temperature=0.8,
        max_tokens=500
    )
    
    prompt_template = PromptTemplate(
        input_variables=["chemical_description"],
        template="""
        You are a highly skilled chemist and computational expert. Your task is to generate a SMILES (Simplified Molecular-Input Line-Entry System) representation for a given chemical bond or molecule.

        Input: {chemical_description}

        Please provide only the SMILES representation of the chemical bond or molecule.

        Example Input: "single bond between carbon and oxygen"
        Example Output: ```C-O```

        Now, please generate the SMILES representation for: {chemical_description}
        """
    )

    chain = LLMChain(llm=llm, prompt=prompt_template)
    return chain.run(chemical_description)


def create_ayurvedic_llm_chain():
    GROQ_API_KEY = os.getenv('GROQ_API_KEY')
    
    llm = ChatGroq(
        groq_api_key=GROQ_API_KEY,
        model_name="gemma-7b-it",
        temperature=0.7,  # Slightly higher for more natural conversation
        max_tokens=1000
    )
    
    prompt_template = PromptTemplate(
        input_variables=["user_input", "user_profile", "consultation_history"],
        template="""
        **Ayurvedic Consultation Context:**
        User Profile: {user_profile}
        Previous Consultation History: {consultation_history}
        Current Query: {user_input}

        **Consultation Guidelines:**
        1. **Initial Assessment:**
           - Analyze the user's prakriti (constitution) and current vikriti (imbalance)
           - Consider seasonal factors (ritucharya) and time of day
           - Evaluate lifestyle factors mentioned in the query

        2. **Dosha Analysis:**
           - Identify predominant doshas involved
           - Assess potential imbalances
           - Consider inter-dosha relationships

        3. **Recommendation Framework:**
           - Dietary suggestions based on dosha balance
           - Lifestyle modifications (dinacharya)
           - Appropriate herbs and spices
           - Exercise and yoga recommendations
           - Daily routine adjustments

        4. **Safety Considerations:**
           - Include relevant health disclaimers
           - Highlight any red flags requiring medical attention
           - Specify limitations of Ayurvedic guidance
           - Recommend professional consultation when necessary

        5. **Educational Component:**
           - Explain Ayurvedic concepts in simple terms
           - Provide rationale for recommendations
           - Include traditional wisdom with modern context
           - Offer practical implementation tips

        **Response Structure:**
        1. Acknowledge the user's concern
        2. Provide dosha-specific analysis
        3. Offer holistic recommendations
        4. Include necessary safety disclaimers
        5. Add educational insights
        6. Suggest follow-up questions if needed

        **Important Guidelines:**
        - Maintain a compassionate and supportive tone
        - Use clear, accessible language
        - Balance traditional wisdom with practical application
        - Always prioritize user safety
        - Encourage mindful self-observation
        - Respect the limits of virtual consultation

        **Medical Disclaimer:**
        This guidance is for educational purposes only and not a substitute for professional medical advice. Always consult qualified healthcare providers for medical conditions.

        **Response:**
        Based on the above guidelines, please provide appropriate Ayurvedic guidance that is both authentic and safe.
        """
    )

    return LLMChain(llm=llm, prompt=prompt_template)