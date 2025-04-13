from google import genai
from google.genai import types
from api_key import key
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import requests, time, os
from fastapi.middleware.cors import CORSMiddleware

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # or specify your frontend URL like ["http://localhost:5173"]
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class UserInput(BaseModel):
    data: str

@app.post("/chat")
def generate_response(data : UserInput):
    client = genai.Client(api_key=key)
    response = client.models.generate_content(
        model="gemini-2.0-flash", contents=data,
        config=types.GenerateContentConfig(
            max_output_tokens=500,
            temperature=0.1
        )
    )

    return {"response" : response.text}
