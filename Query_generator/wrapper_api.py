from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
import uvicorn

from .generate_boolean_query import generate_boolean

app = FastAPI(title="Query Generator Wrapper")


class QueryRequest(BaseModel):
    query: str
    top_k: int = 5
    expand: bool = False


@app.post('/generate')
def generate(req: QueryRequest):
    try:
        out = generate_boolean(req.query, top_k=req.top_k, expand=req.expand)
        return out
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == '__main__':
    uvicorn.run('Query_generator.wrapper_api:app', host='0.0.0.0', port=8001, reload=False)
