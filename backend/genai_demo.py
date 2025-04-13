from google import genai
from google.genai import types
from api_key import key

client = genai.Client(api_key=key)

response = client.models.generate_content(
    model="gemini-2.0-flash", contents="Explain how AI works in a few words",
    config=types.GenerateContentConfig(
            max_output_tokens=50,
            temperature=0.1
        )
)
print(response.text)