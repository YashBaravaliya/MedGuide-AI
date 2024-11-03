from pymed import PubMed

pubmed = PubMed(tool="StreamlitApp", email="dummyemail1@gmail.com")

def documentize(article):
    return {
        'content': article.abstract,
        'meta': {
            'title': article.title,
            'keywords': article.keywords
        }
    }

def fetch_pubmed_articles(queries: str):
    cleaned_queries = queries.strip().split('\n')
    articles = []
    for query in cleaned_queries:
        response = pubmed.query(query, max_results=1)
        documents = [documentize(article) for article in response]
        articles.extend(documents)
    return articles

def documentize_for_research(article):
    # Helper function to safely get attributes
    def safe_get(obj, attr):
        return getattr(obj, attr, None)
    
    # Helper function to safely get author names
    def safe_get_authors(article):
        authors = safe_get(article, 'authors')
        if authors:
            return [getattr(author, 'name', 'Unknown Author') for author in authors]
        return []

    return {
        'content': safe_get(article, 'abstract') or "Abstract not available",
        'meta': {
            'title': safe_get(article, 'title') or "Title not available",
            'keywords': safe_get(article, 'keywords') or [],
            'doi': safe_get(article, 'doi'),
            'pubmed_id': safe_get(article, 'pubmed_id'),
            'publication_date': str(safe_get(article, 'publication_date')) if safe_get(article, 'publication_date') else None,
            'authors': safe_get_authors(article)
        }
    }

def fetch_pubmed_articles_for_research(query: str, max_results: int = 5):
    try:
        response = pubmed.query(query, max_results=max_results)
        documents = [documentize_for_research(article) for article in response]
        return documents
    except Exception as e:
        print(f"Error fetching PubMed articles: {str(e)}")
        return []
