# utils/ayurvedic_utils.py

from typing import Dict, List, Union, Optional
import re

class DoshaAnalyzer:
    def __init__(self):
        # Keywords associated with each dosha
        self.dosha_keywords = {
            'vata': [
                'anxiety', 'insomnia', 'dry', 'cold', 'light', 'irregular', 'constipation',
                'nervousness', 'unstable', 'restless', 'worry', 'joint pain', 'crackling',
                'stiffness', 'gas', 'bloating', 'variable appetite', 'weight loss'
            ],
            'pitta': [
                'anger', 'irritation', 'hot', 'sharp', 'acidic', 'burning', 'inflammation',
                'fever', 'rash', 'acidity', 'ulcer', 'heartburn', 'excessive hunger',
                'perfectionist', 'critical', 'competitive', 'skin problems', 'diarrhea'
            ],
            'kapha': [
                'lethargy', 'depression', 'heavy', 'slow', 'cold', 'congestion', 'mucus',
                'attachment', 'possessiveness', 'weight gain', 'water retention', 'diabetes',
                'sleeping excess', 'oversleeping', 'cough', 'allergies', 'asthma', 'slow digestion'
            ]
        }
        
        # Seasonal influences (Ritucharya)
        self.seasonal_effects = {
            'summer': {'vata': 0, 'pitta': 1.2, 'kapha': 0.8},
            'winter': {'vata': 1.2, 'pitta': 0.8, 'kapha': 1.1},
            'spring': {'vata': 0.8, 'pitta': 0.9, 'kapha': 1.2},
            'autumn': {'vata': 1.1, 'pitta': 0.9, 'kapha': 0.9},
            'monsoon': {'vata': 1.0, 'pitta': 1.1, 'kapha': 1.0}
        }

def analyze_dosha(query: str, user_profile: Optional[Dict] = None) -> Dict:
    """
    Analyze the dosha imbalances based on user query and profile.
    
    Args:
        query (str): User's health concern or query
        user_profile (dict): User's profile information including age, gender, known prakriti
        
    Returns:
        dict: Detailed dosha analysis
    """
    analyzer = DoshaAnalyzer()
    
    # Simple text preprocessing
    query = query.lower()
    # Remove special characters
    query = re.sub(r'[^\w\s]', '', query)
    # Split into words
    words = query.split()
    
    # Initialize dosha scores
    dosha_scores = {
        'vata': 0,
        'pitta': 0,
        'kapha': 0
    }
    
    # Analyze keywords in query
    for word in words:
        for dosha, keywords in analyzer.dosha_keywords.items():
            if word in keywords:
                dosha_scores[dosha] += 1
            # Check for partial matches (e.g., "anxious" matches "anxiety")
            for keyword in keywords:
                if keyword in word or word in keyword:
                    dosha_scores[dosha] += 0.5
    
    # Consider user profile if available
    if user_profile:
        # Age considerations
        age = user_profile.get('age', 30)
        if age < 30:  # Kapha predominant age
            dosha_scores['kapha'] += 0.5
        elif 30 <= age < 60:  # Pitta predominant age
            dosha_scores['pitta'] += 0.5
        else:  # Vata predominant age
            dosha_scores['vata'] += 0.5
        
        # Consider known prakriti
        if 'prakriti' in user_profile:
            for dosha in user_profile['prakriti']:
                dosha_scores[dosha.lower()] += 0.3
    
    # Determine primary and secondary doshas
    sorted_doshas = sorted(dosha_scores.items(), key=lambda x: x[1], reverse=True)
    
    # Prepare symptom lists without NLTK
    symptoms = {
        'vata_symptoms': [kw for kw in analyzer.dosha_keywords['vata'] if kw in query],
        'pitta_symptoms': [kw for kw in analyzer.dosha_keywords['pitta'] if kw in query],
        'kapha_symptoms': [kw for kw in analyzer.dosha_keywords['kapha'] if kw in query]
    }
    
    return {
        'primary_dosha': sorted_doshas[0][0],
        'secondary_dosha': sorted_doshas[1][0],
        'dosha_scores': dosha_scores,
        'detailed_analysis': symptoms
    }

def get_remedies(dosha_analysis: Dict) -> Dict:
    """
    Generate Ayurvedic remedies based on dosha analysis.
    
    Args:
        dosha_analysis (dict): Result from analyze_dosha function
        
    Returns:
        dict: Recommended remedies including herbs, lifestyle changes, and dietary advice
    """
    # Basic remedy templates
    remedy_templates = {
        'vata': {
            'herbs': [
                {'name': 'Ashwagandha', 'benefits': 'Calming, grounding, and strengthening for the nervous system'},
                {'name': 'Brahmi', 'benefits': 'Supports mental clarity and reduces anxiety'},
                {'name': 'Trikatu', 'benefits': 'Improves digestion and reduces gas'}
            ],
            'lifestyle': [
                'Maintain regular daily routine',
                'Practice gentle yoga',
                'Get adequate rest',
                'Avoid excessive travel or movement',
                'Use warm oil massage (abhyanga)'
            ],
            'diet': {
                'recommended': [
                    'Warm, cooked foods',
                    'Sweet fruits',
                    'Cooked vegetables',
                    'Warm milk',
                    'Ghee',
                    'Rice'
                ],
                'avoid': [
                    'Cold foods',
                    'Raw vegetables',
                    'Dry fruits',
                    'Carbonated drinks',
                    'Coffee'
                ],
                'general_guidelines': 'Favor warm, cooked, slightly oily foods. Include sweet, sour, and salty tastes. Avoid raw foods and irregular eating patterns.'
            }
        },
        'pitta': {
            'herbs': [
                {'name': 'Amalaki', 'benefits': 'Cooling and rejuvenating, balances stomach acid'},
                {'name': 'Brahmi', 'benefits': 'Cooling for mind and body, reduces inflammation'},
                {'name': 'Coriander', 'benefits': 'Cooling spice that aids digestion'}
            ],
            'lifestyle': [
                'Avoid excessive heat and sun exposure',
                'Practice cooling exercises like swimming',
                'Meditation for emotional balance',
                'Take breaks during work',
                'Moonlight walks'
            ],
            'diet': {
                'recommended': [
                    'Sweet fruits',
                    'Cooling vegetables',
                    'Coconut water',
                    'Milk',
                    'Rice',
                    'Mint tea'
                ],
                'avoid': [
                    'Spicy foods',
                    'Fermented foods',
                    'Sour fruits',
                    'Hot drinks',
                    'Alcohol'
                ],
                'general_guidelines': 'Favor cooling, sweet, bitter, and astringent foods. Avoid spicy, sour, and fermented foods. Eat at regular times.'
            }
        },
        'kapha': {
            'herbs': [
                {'name': 'Trikatu', 'benefits': 'Stimulates metabolism and reduces congestion'},
                {'name': 'Guggulu', 'benefits': 'Supports healthy weight management'},
                {'name': 'Ginger', 'benefits': 'Improves digestion and reduces mucus'}
            ],
            'lifestyle': [
                'Regular vigorous exercise',
                'Wake up early',
                'Dry massage with powder',
                'Keep active throughout the day',
                'Avoid daytime naps'
            ],
            'diet': {
                'recommended': [
                    'Light soups',
                    'Steamed vegetables',
                    'Spicy foods',
                    'Honey',
                    'Green tea',
                    'Bitter greens'
                ],
                'avoid': [
                    'Dairy products',
                    'Sweet foods',
                    'Cold drinks',
                    'Oily foods',
                    'Heavy meats'
                ],
                'general_guidelines': 'Favor light, dry, warm foods. Include pungent, bitter, and astringent tastes. Avoid heavy, cold, and sweet foods.'
            }
        }
    }
    
    # Get primary dosha remedies
    primary_dosha = dosha_analysis['primary_dosha']
    remedies = remedy_templates[primary_dosha].copy()
    
    # Add secondary dosha considerations
    secondary_dosha = dosha_analysis['secondary_dosha']
    if secondary_dosha:
        # Add one herb from secondary dosha
        remedies['herbs'].append(remedy_templates[secondary_dosha]['herbs'][0])
        # Add two lifestyle recommendations from secondary dosha
        remedies['lifestyle'].extend(remedy_templates[secondary_dosha]['lifestyle'][:2])
        
    # Add general recommendations
    remedies['general'] = [
        'Practice pranayama (breathing exercises)',
        'Maintain proper sleep schedule',
        'Stay hydrated with appropriate temperature water',
        'Practice mindful eating'
    ]
    
    return remedies