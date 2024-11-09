# Medicine Recommendation & Research Chatbot

This project utilizes LangChain and Groq to build a chatbot that provides medically accurate answers to patient questions by integrating medicine information, related research articles, and Ayurvedic treatment suggestions.

## Project Overview

The model is designed to answer medical queries with a patient-focused approach, combining detailed information about medicines with relevant research articles and alternative Ayurvedic treatments. The model ensures clarity, accuracy, and a simplified explanation for patients to easily understand. It also suggests Ayurvedic plant-based alternatives and provides detailed compositions of these plants.

## Technologies Used

- **LangChain**: A framework to build language model applications. Used for chaining prompts and integrating LLMs.
- **Groq**: A hardware-accelerated model inference platform. It is used to run the language model, `gemma-7b-it`.
- **Python**: The main programming language used to build and run the application.
- **Environment Variables**: API keys are securely stored and accessed through environment variables.

## Steps to Create the Model

### 1. Set Up the Environment
Ensure you have the following dependencies installed:

```bash
pip install langchain groq
```

### 2. Create Environment Variables

You will need a valid **Groq API Key** to access the `gemma-7b-it` model. Set this key in your environment variables: https://console.groq.com/keys

```bash
export GROQ_API_KEY="your-groq-api-key"
```

### 3. Define the Prompt Template

The prompt template is defined to structure the information to be provided to the model. It includes the following components:

- **Question**: The user's medical question.
- **Medicine Information**: Detailed information about the medicine being discussed.
- **Related Research Articles**: Relevant PubMed articles.
- **Guidelines for Answer**: Instructions to ensure the response is accurate, clear, and patient-focused.

```python
prompt_template = PromptTemplate(
    input_variables=["question", "pubmed_results", "medicine_info"],
    template=""" ... """
)
```

### 4. Initialize the ChatGroq Model

The `ChatGroq` class is initialized with the following parameters:

- **groq_api_key**: Your Groq API key.
- **model_name**: The name of the model being used (`gemma-7b-it`).
- **temperature**: Controls the randomness of the response. Set to `0.5` for balanced responses.
- **max_tokens**: Limits the length of the generated response (500 tokens in this case).

```python
llm = ChatGroq(
    groq_api_key=GROQ_API_KEY,
    model_name="gemma-7b-it",
    temperature=0.5,
    max_tokens=500
)
```

### 5. Define the LLM Chain

The `LLMChain` is used to link the model and the prompt together, so that the inputs are passed through the prompt template to generate the final response.

```python
llm_chain = LLMChain(llm=llm, prompt=prompt_template)
```

### 6. Generate the Response

Once the LLM chain is created, you can use it to generate responses based on inputs such as the user's question, PubMed research results, and medicine information. Here's an example:

```python
response = llm_chain.run({
    "question": "What are the side effects of paracetamol?",
    "pubmed_results": "<Related PubMed articles>",
    "medicine_info": "<Detailed medicine information>"
})
```

This will generate a medically accurate and patient-friendly answer, including Ayurvedic alternatives where applicable.

## Example Output

### **Question**: What are the side effects of paracetamol?
### **Medicine Information**:
Paracetamol is a commonly used pain reliever. Some potential side effects include liver damage with high doses and allergic reactions such as skin rash.

### **Related Research Articles**:
- PubMed Article 1
- PubMed Article 2

### **Final Answer**:
Paracetamol is effective for relieving pain but can cause liver damage if taken in excess. It's essential to follow the prescribed dosage. If you're seeking an alternative treatment, **Ashwagandha** is an Ayurvedic plant with similar pain-relieving properties. It is commonly used in traditional medicine for managing stress and inflammation.
