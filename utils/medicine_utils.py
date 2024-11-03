import pandas as pd
import faiss
import pickle

def load_medicine_data():
    with open('model/sentence_transformer_model.pkl', 'rb') as f:
        model = pickle.load(f)
    index = faiss.read_index('model/medicine_faiss_index.bin')
    df = pd.read_pickle('model/medicine_df.pkl')
    return model, index, df

model, faiss_index, medicine_df = load_medicine_data()

def search_medicine(query):
    query_vector = model.encode([query])
    faiss.normalize_L2(query_vector)
    
    k = 3
    D, I = faiss_index.search(query_vector, k)
    
    results = []
    for i, dist in zip(I[0], D[0]):
        if dist > 0.5:
            result = medicine_df.iloc[i].to_dict()
            result['similarity_score'] = float(dist)
            
            if 'Composition' in result:
                ayurvedic_plants = [plant.strip() for plant in result['Composition'].split(',')]
                result['Ayurvedic_Plants'] = ayurvedic_plants
            
            results.append(result)
    
    if results:
        primary_medicine = results[0]
        alternative_medicines = results[1:4]
        return {
            "primary": primary_medicine,
            "alternatives": alternative_medicines
        }
    else:
        return None