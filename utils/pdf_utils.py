import io
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from reportlab.lib.styles import getSampleStyleSheet
from PIL import Image as PILImage
import requests

def generate_pdf(primary_medicine, alternative_medicines, answer):
    pdf_file_name = primary_medicine['Name']
    buffer = io.BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter)
    styles = getSampleStyleSheet()
    story = []

    # Primary Medicine
    story.append(Paragraph(f"Primary Medicine: {primary_medicine['Name']}", styles['Title']))
    story.append(Spacer(1, 12))

    for key, value in primary_medicine.items():
        if key != 'img_url' and key != 'combined_text':
            if key == "Unnamed: 0":
                continue
            story.append(Paragraph(f"<b>{key}:</b> {value}", styles['Normal']))
            story.append(Spacer(1, 6))

    if primary_medicine.get('img_url'):
        try:
            img_response = requests.get(primary_medicine['img_url'])
            img = PILImage.open(io.BytesIO(img_response.content))
            img_width, img_height = img.size
            aspect = img_height / float(img_width)
            img_width = 300
            img_height = int(img_width * aspect)
            story.append(Image(io.BytesIO(img_response.content), width=img_width, height=img_height))
        except Exception as e:
            story.append(Paragraph(f"Error loading primary medicine image: {str(e)}", styles['Normal']))

    # Alternative Medicines
    story.append(Paragraph("Alternative Medicines", styles['Heading2']))
    story.append(Spacer(1, 12))

    for i, alt_med in enumerate(alternative_medicines, 1):
        story.append(Paragraph(f"Alternative {i}: {alt_med['Name']}", styles['Heading3']))
        for key, value in alt_med.items():
            if key != 'img_url' and key != 'combined_text':
                if key == "Unnamed: 0":
                    continue
                story.append(Paragraph(f"<b>{key}:</b> {value}", styles['Normal']))
                story.append(Spacer(1, 6))
        
        if alt_med.get('img_url'):
            try:
                img_response = requests.get(alt_med['img_url'])
                img = PILImage.open(io.BytesIO(img_response.content))
                img_width, img_height = img.size
                aspect = img_height / float(img_width)
                img_width = 200
                img_height = int(img_width * aspect)
                story.append(Image(io.BytesIO(img_response.content), width=img_width, height=img_height))
            except Exception as e:
                story.append(Paragraph(f"Error loading alternative medicine image: {str(e)}", styles['Normal']))
        
        story.append(Spacer(1, 12))

    doc.build(story)
    buffer.seek(0)
    return buffer