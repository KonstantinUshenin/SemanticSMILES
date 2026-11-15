import os
from dotenv import load_dotenv
from openai import OpenAI, DefaultHttpxClient

if __name__ == "__main__":
    load_dotenv()
    OPENAI_API_KEY = os.getenv("OPENAI_API_KEY")
    OPENAI_API_BASE = os.getenv("OPENAI_API_BASE")

    client = OpenAI(
        base_url=OPENAI_API_BASE,
        api_key=OPENAI_API_KEY,
        http_client=DefaultHttpxClient(verify=False)
    )

    models = client.models.list()
    print([model.id for model in models.data])