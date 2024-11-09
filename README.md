# MedGuide AI Dashboard

MedGuide AI is a Streamlit-based application that integrates multiple healthcare modules, providing a versatile tool for medical information retrieval, molecular analysis, research support, and physiotherapy planning. This AI-driven tool leverages language models, PubMed research integration, and a custom FAISS index for optimized medicine suggestions, including Ayurvedic alternatives.

---

## Overview

![Dashboard](https://github.com/YashBaravaliya/MedGuide-AI/tree/main/img/ss/Dashboard.png)

The MedGuide AI dashboard offers the following main features:

- **MedGuide Chat**: Provides answers for medical questions using PubMed research articles and suggests Ayurvedic alternatives when available.
- **Molecule Generation**: Allows users to visualize molecular structures.
- **Research Support**: Assists in retrieving top scientific articles from PubMed related to healthcare topics.
- **Physio Planner**: Helps users create custom physiotherapy plans.
- **MediScan**: user can upload there medicine and get details about 


## 1. MedGuide Chat

MedGuide Chat is an interactive healthcare assistance application built with Streamlit. It helps users find reliable information on medications, symptoms, and related medical topics. The chatbot leverages advanced language models, medical databases, and PubMed resources to provide comprehensive and accurate responses.

### Features

#### 1.1. Interactive Chatbot

- **Conversational Interface**: Type questions related to medicines, symptoms, or general healthcare inquiries, and receive detailed, contextual answers.
- **Medicine-Specific Responses**: If a medication is found, the chatbot displays information on its uses, side effects, composition, and manufacturer.
- **Alternative Options**: If available, alternative medicines are suggested along with images.
- **Downloadable Reports**: Users can download a PDF report containing detailed information about the medication.

#### 1.2. Medication Information Panel

![img](https://github.com/YashBaravaliya/MedGuide-AI/tree/main/img/ss/medGuideChat1.png)
![img](https://github.com/YashBaravaliya/MedGuide-AI/tree/main/img/ssmedGuideChat2.png)
- **Primary Medicine Details**: View comprehensive information on the primary medicine found for the query, including:
  - Uses
  - Composition
  - Manufacturer
  - Side effects
- **Alternative Medicines**: Suggestions with images, giving users options beyond the primary medicine.
  
#### 1.3. Healthcare Source Search

![img](https://github.com/YashBaravaliya/MedGuide-AI/tree/main/img/ss/medGuideChat3.png)
- **HealthcareSearchAgent**: Searches medical sources relevant to the query and provides URLs and summaries of useful articles, making it easier to verify information and read further.

  
### Technology Stack

- **Streamlit**: For the web interface and chat interaction.
- **Groq API**: Used to run the language model with Groq's hardware acceleration. [more info](knowledge_base\medGuideChat.md)
- **Python**: Backend logic and data processing.
- **RDKit**: For handling chemical structures.
- **PubMed API**: To fetch reliable, research-backed medical articles.
- **Custom Libraries and Utils**:
  - `pubmed_utils` for article fetching
  - `medicine_utils` for medication search
  - `pdf_utils` for generating PDF reports
  - `llm_chain` for generating responses with language models
